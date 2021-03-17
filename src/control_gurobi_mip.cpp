#include <ros/ros.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/Twist.h>
#include <geometry_msgs/Pose2D.h>
#include <tf/transform_datatypes.h>
#include <nav_msgs/Path.h>
#include <visualization_msgs/Marker.h>

#include <yaml-cpp/yaml.h>
#include <fstream>
#include <math.h>
#include <chrono>
#include <algorithm>
#include <valarray>
#include <ifopt/problem.h>
#include <ifopt/snopt_solver.h>
#include "gurobi_c++.h"
#include "pusher/problem.h"
#include "pusher/CalculateControl.h"


class Control
{
public:
    Control(ros::NodeHandle n, GRBEnv envIn) : nh(n), env(envIn)
    {
        pusherVelPub = nh.advertise<geometry_msgs::Twist>("/cmd_vel", 1);
        pusherMovePub = nh.advertise<geometry_msgs::Twist>("/move_pusher", 1);

        pusherMoveSub = nh.subscribe("/move_pusher", 1, &Control::pusherMoveCallback, this);
        pusherOdomSub = nh.subscribe("/odom/pusher", 1, &Control::pusherOdomCallback, this);
        sliderOdomSub = nh.subscribe("/odom/slider", 1, &Control::sliderOdomCallback, this);
    
        pathNomiPub = nh.advertise<nav_msgs::Path>("/path/nominal", 10);
        pathMPCPub = nh.advertise<nav_msgs::Path>("/path/mpc", 10);
        markerPub = nh.advertise<visualization_msgs::Marker>( "/slider_marker", 0 );
        calculateControlClient = nh.serviceClient<pusher::CalculateControl>("/calculate_control");

        YAML::Node params = YAML::LoadFile("/home/mzwang/qsp_ws/src/pusher/Config/params.yaml");
        mu = params["mu"].as<double>();
        muGround = params["mu_g"].as<double>();
        m = params["m"].as<double>();
        xC = -params["length"].as<double>() / 2.0;
        width = params["width"].as<double>();
        mMax = params["mMax"].as<double>();
        fMax = muGround * m * 9.81;
        MPCSteps = params["n_step"].as<int>();
        tStep = params["t_step"].as<double>();
        QFWeight = params["QFinal"].as<double>();
        QWeight = params["Q"].as<double>();
        RWeight = params["R"].as<double>();
        debugInfo = params["debug_info"].as<bool>();

        QFinal[0] = 3*QFWeight; QFinal[1] = 3*QFWeight; QFinal[2] = 0.1*QFWeight; QFinal[3] = 0*QFWeight;
        Q[0] = 3*QWeight; Q[1] = 3*QWeight; Q[2] = 0.1*QWeight; Q[3] = 0*QWeight;
        R[0] = 1*RWeight; R[1] = 1*RWeight; R[2] = 0.01*RWeight;

        J.resize(2,3); 
        B.resize(3,2); 
        L.resize(3,3);
        L << 2/(fMax*fMax), 0, 0,
             0, 2/(fMax*fMax), 0,
             0, 0, 2/(mMax*mMax);

        state.resize(4*MPCSteps);
        control.resize(3*MPCSteps-3);
        for (int i = 0; i < MPCSteps; ++i) {
            state.segment(i*4, 4) << 0.05*i*tStep, 0,0,0;
        }
        for (int i = 0; i < MPCSteps-1; ++i) {
            control.segment(i*3, 3) << 0.3, 0.0, 0.0;
        }

        zeroVel.linear.x = 0;
        zeroVel.linear.y = 0;
        zeroVel.linear.z = 0;
        zeroVel.angular.x = 0;
        zeroVel.angular.y = 0;
        zeroVel.angular.z = 0;

        stepCounter = 0;

        marker.header.frame_id = "world";
        marker.ns = "my_namespace";
        marker.id = 0;
        marker.type = visualization_msgs::Marker::CUBE;
        marker.action = visualization_msgs::Marker::ADD;
        marker.scale.x = 0.01;
        marker.scale.y = 0.01;
        marker.scale.z = 0.01;
        marker.color.a = 0.5; // Don't forget to set the alpha!
        marker.color.r = 0.0;
        marker.color.g = 0.0;
        marker.color.b = 1.0;

        // Gurobi
        stateNominal.resize(MPCSteps*4);
        controlNominal.resize((MPCSteps-1)*3);

        // define the nominal trajectory
        float targetLV = 0.05;
        float targetLength = 0.5;
        totalSteps = int(targetLength/(targetLV*tStep));
        lineNomi.resize(4*totalSteps+4);
        for (int i = 0; i < totalSteps; ++i){
            lineNomi.segment(i*4, 4) << targetLV*i*tStep, 0, 0, 0;
        }
        lineNomi.tail(4) << targetLength, 0, 0, 0;
        lineControlNomi.resize(3*totalSteps);
        for (int i = 0; i < totalSteps-1; ++i){
            lineControlNomi.segment(i*3, 3) << 0.305, 0, 0, 0;
        }

        // circular tracking
        radius = 0.15;
        float targetAV = 0.2;
        int nPointsPi = 3.14/(targetAV*tStep);        
        eightNomi.resize(nPointsPi*4*4);
        for (int i = 0; i < nPointsPi; ++i){
            eightNomi.segment(i*4, 4) << radius * sin(targetAV*i*tStep), radius - radius*cos(targetAV*i*tStep), targetAV*i*tStep, 0;
            eightNomi.segment(4*nPointsPi + i*4, 4) << -radius * sin(targetAV*i*tStep), 3*radius - radius*cos(targetAV*i*tStep), 3.14 - targetAV*i*tStep, 0;
            eightNomi.segment(2*4*nPointsPi + i*4, 4) << radius * sin(targetAV*i*tStep), 3*radius + radius*cos(targetAV*i*tStep), - targetAV*i*tStep, 0;
            eightNomi.segment(3*4*nPointsPi + i*4, 4) << -radius * sin(targetAV*i*tStep), radius + radius*cos(targetAV*i*tStep), -3.14 + targetAV*i*tStep, 0;
        }


        // for recording
        actualState.resize(4*(totalSteps));
        actualState.setZero();
        solvingTime.resize(totalSteps);
        solvingTime.setZero();
        errors.resize(totalSteps);
        errors.setZero();
        timePath = "/home/mzwang/qsp_ws/src/pusher/logs/time_mip.txt";
        errorPath = "/home/mzwang/qsp_ws/src/pusher/logs/error_mip.txt";

        timer1 = nh.createTimer(ros::Duration(tStep), &Control::findSolution, this);
        // timerVel = nh.createTimer(ros::Duration(0.01), &Control::pubPusherVel, this);
    }

    void pusherOdomCallback(const nav_msgs::Odometry &msg) 
    {
        double xSliderToPusher = msg.pose.pose.position.x - sliderPose.x;
        double ySliderToPusher = msg.pose.pose.position.y - sliderPose.y;
        yC = -sin(sliderPose.theta) * xSliderToPusher + cos(sliderPose.theta) * ySliderToPusher;
        phi = atan(-yC/xC);
    }

    void sliderOdomCallback(const nav_msgs::Odometry &msg) 
    {
        sliderPose.x = msg.pose.pose.position.x;
        sliderPose.y = msg.pose.pose.position.y;
        tf::Pose pose;
        tf::poseMsgToTF(msg.pose.pose, pose);
        float newTheta = tf::getYaw(pose.getRotation());
        // To prevent sudden change from pi to -pi or vice versa
        if (abs(newTheta - sliderPose.theta) < 3){
            sliderPose.theta = newTheta;
        }

        marker.header.stamp = ros::Time();
        marker.pose = msg.pose.pose;
        marker.pose.position.z = 0;
    }

    void pusherMoveCallback(const geometry_msgs::Twist &msg){
        pusherVelPub.publish(msg);
        ros::Duration(tStep).sleep();
        pusherVelPub.publish(zeroVel);
    }
    
    void findSolution(const ros::TimerEvent&)
    {   
        // straight line
        if (stepCounter > totalSteps){
            timeFile.open(timePath, std::ios::app);
            if (timeFile.is_open()){
                timeFile << solvingTime.transpose() <<"\n";
            }
            else{
                std::cout << " WARNING: Unable to open the recording file.\n";
            }
            timeFile.close();
            errorFile.open(errorPath, std::ios::app);
            if (errorFile.is_open()){
                errorFile << errors.transpose() <<"\n";
            }
            else{
                std::cout << " WARNING: Unable to open the recording file.\n";
            }
            errorFile.close();
            ros::Duration(10).sleep();
        } else if (stepCounter*4+MPCSteps*4 > lineNomi.size()){
            stateNominal.head(MPCSteps*4 - 4) = stateNominal.tail(MPCSteps*4 - 4);
        } else {
            stateNominal = lineNomi.segment(stepCounter*4, MPCSteps*4);
            controlNominal = lineControlNomi.segment(stepCounter*3, (MPCSteps-1)*3);
        }

        // 8 shape
        // if (stepCounter*4+MPCSteps*4 > eightNomi.size()){
        //     stepCounter = 0;
        // }
        // stateNominal = eightNomi.segment(stepCounter*4, MPCSteps*4);
                
        std::cout << sliderPose.x - stateNominal(0) << ", " << sliderPose.y - stateNominal(1) << "," << sliderPose.theta - stateNominal(2) << ", " <<phi - stateNominal(3) << "\n";
        double error = sqrt((stateNominal(0) - sliderPose.x) * (stateNominal(0) - sliderPose.x) + 
                       (stateNominal(1) - sliderPose.y) * (stateNominal(1) - sliderPose.y));
        // Gurobi solver
        GRBModel model = GRBModel(env);
        GRBVar xBar[MPCSteps][4], uBar[MPCSteps-1][3], z[MPCSteps-1][3];

        // initial error
        xBar[0][0] = model.addVar(sliderPose.x - stateNominal(0), sliderPose.x - stateNominal(0), 0, GRB_CONTINUOUS, "state[0]");
        xBar[0][1] = model.addVar(sliderPose.y - stateNominal(1), sliderPose.y - stateNominal(1), 0, GRB_CONTINUOUS, "state[20]");
        xBar[0][2] = model.addVar(sliderPose.theta - stateNominal(2), sliderPose.theta - stateNominal(2), 0, GRB_CONTINUOUS, "state[40]");
        xBar[0][3] = model.addVar(phi - stateNominal(3), phi - stateNominal(3), 0, GRB_CONTINUOUS, "state[60]");

        for (int i = 1; i < MPCSteps; ++i){
            xBar[i][0] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "state["+std::to_string(i)+"]");
            xBar[i][1] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "state["+std::to_string(20+i)+"]");
            xBar[i][2] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "state["+std::to_string(40+i)+"]");
            xBar[i][3] = model.addVar(-0.70 - stateNominal(4*i+3), 0.70 - stateNominal(4*i+3), 0, GRB_CONTINUOUS, "state["+std::to_string(60+i)+"]");
        }

        for (int i = 0; i < MPCSteps-1; ++i){
            uBar[i][0] = model.addVar(-controlNominal(3*i), GRB_INFINITY, 0, GRB_CONTINUOUS, "control["+std::to_string(i)+"]");
            uBar[i][1] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "control["+std::to_string(19+i)+"]");
            uBar[i][2] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "control["+std::to_string(38+i)+"]");
            z[i][0] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "mode["+std::to_string(i)+"]");
            z[i][1] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "mode["+std::to_string(19+i)+"]");
            z[i][2] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "mode["+std::to_string(38+i)+"]");
        }

        // Set cost
        GRBQuadExpr cost = 0;
        cost.addTerms(QFinal, xBar[MPCSteps-1], xBar[MPCSteps-1], 4);
        for (int i = 0; i < MPCSteps-1; ++i){
            cost.addTerms(Q, xBar[i], xBar[i], 4);
            cost.addTerms(R, uBar[i], uBar[i], 3);
            if (i > 0){
                cost += 0.1*(z[i][0]-1)*(z[i][0]-1) + 0.1*z[i][1]*z[i][1]+0.1*z[i][2]*z[i][2];
            }
        }
        model.setObjective(cost, GRB_MINIMIZE);

        for (int i = 0; i < MPCSteps-1; ++i){
            double fn = controlNominal(i*3);
            double ft = controlNominal(i*3+1);
            double pdot = controlNominal(i*3+2);
            double theta = stateNominal(i*4+2);
            double phi = stateNominal(i*4+3);
            // dynamics A*xBar + B*uBar = x_k+1 - x
            double A[] = {0, 0, -2*fn*sin(theta)/pow(fMax,2) - 2*ft*cos(theta)/pow(fMax,2), 0, 
                          0, 0, 2*fn*cos(theta)/pow(fMax,2) - 2*ft*sin(theta)/pow(fMax,2), 0, 
                          0, 0, 0, 2*fn*xC*(pow(tan(phi),2) + 1)/pow(mMax,2), 
                          0, 0, 0, 0};
            double B[] = {2*cos(theta)/pow(fMax,2), -2*sin(theta)/pow(fMax,2), 0, 
                          2*sin(theta)/pow(fMax,2), 2*cos(theta)/pow(fMax,2), 0, 
                          2*xC*tan(phi)/pow(mMax,2), 2*xC/pow(mMax,2), 0, 
                          0, 0, 1};
            for (int j = 0; j < 4; ++j){
                GRBLinExpr dynamics = 0.0;
                dynamics.addTerms(A+j*4, xBar[i], 4);
                dynamics.addTerms(B+j*3, uBar[i], 3);
                model.addConstr(dynamics == (xBar[i+1][j] - xBar[i][j]));
            }
            // velocity mapping, bounds the ouptut velocity; pusherV - nominal velicity at step i
            Eigen::Vector2d pusherV;
            Eigen::MatrixXd dVdx(2,4), dVdu(2,3);
            pusherV << fn*(2*pow(xC,2)*pow(tan(phi),2)/pow(mMax,2) + 2/pow(fMax,2)) + 2*ft*pow(xC,2)*tan(phi)/pow(mMax,2), 
                       2*fn*pow(xC,2)*tan(phi)/pow(mMax,2) + ft*(2*pow(xC,2)/pow(mMax,2) + 2/pow(fMax,2)) - pdot*xC/pow(cos(phi),2); 
            dVdx << 0, 0, 0, 2*fn*pow(xC,2)*(2*pow(tan(phi),2) + 2)*tan(phi)/pow(mMax,2) + 2*ft*pow(xC,2)*(pow(tan(phi),2) + 1)/pow(mMax,2), 
                    0, 0, 0, 2*fn*pow(xC,2)*(pow(tan(phi),2) + 1)/pow(mMax,2) - 2*pdot*xC*sin(phi)/pow(cos(phi),3);
            dVdu << 2*pow(xC,2)*pow(tan(phi),2)/pow(mMax,2) + 2/pow(fMax,2), 2*pow(xC,2)*tan(phi)/pow(mMax,2), 0, 
                    2*pow(xC,2)*tan(phi)/pow(mMax,2), 2*pow(xC,2)/pow(mMax,2) + 2/pow(fMax,2), -xC/pow(cos(phi),2);
            GRBLinExpr vMappedX = pusherV(0), vMappedY = pusherV(1);
            for (int j = 0; j < 4; ++j){
                vMappedX += dVdx(0,j)*xBar[i][j];
                vMappedY += dVdx(1,j)*xBar[i][j];
            }
            for (int j = 0; j < 3; ++j){
                vMappedX += dVdu(0,j)*uBar[i][j];
                vMappedY += dVdu(1,j)*uBar[i][j];
            }
            model.addConstr(vMappedX >= 0);
            model.addConstr(vMappedX <= 0.3);
            model.addConstr(vMappedY >= -0.3);
            model.addConstr(vMappedY <= 0.3);
            // friction cone
            model.addConstr(ft+uBar[i][1] <= mu*(fn+uBar[i][0]));
            model.addConstr(ft+uBar[i][1] >= -mu*(fn+uBar[i][0]));
            // mode selection
            double M = 1.0;
            // stick
            model.addConstr(pdot+uBar[i][2] >= M*(z[i][0] - 1));
            model.addConstr(pdot+uBar[i][2] <= M*(-z[i][0] + 1));
            // up
            model.addConstr(pdot+uBar[i][2] >= 5*M*(z[i][1] - 1));
            model.addConstr(mu*(fn+uBar[i][0]) - (ft+uBar[i][1]) >= M*(z[i][1]-1));
            model.addConstr(mu*(fn+uBar[i][0]) - (ft+uBar[i][1]) <= M*(-z[i][1]+1));
            // down
            model.addConstr(pdot+uBar[i][2] <= 5*M*(-z[i][2] + 1));
            model.addConstr(mu*(fn+uBar[i][0]) + (ft+uBar[i][1]) >= M*(z[i][2]-1));
            model.addConstr(mu*(fn+uBar[i][0]) + (ft+uBar[i][1]) <= M*(-z[i][2]+1));
            model.addConstr(z[i][0]+z[i][1]+z[i][2] == 1);
        }
        // agglomerated
        // for (int i = 1; i < 5; ++i) {
        //     model.addConstr(z[1][0]==z[1+i][0]);
        //     model.addConstr(z[1][1]==z[1+i][1]);
        //     model.addConstr(z[1][2]==z[1+i][2]);
        //     model.addConstr(z[6][0]==z[6+i][0]);
        //     model.addConstr(z[6][1]==z[6+i][1]);
        //     model.addConstr(z[6][2]==z[6+i][2]);
        //     model.addConstr(z[11][0]==z[6+i][0]);
        //     model.addConstr(z[11][1]==z[6+i][1]);
        //     model.addConstr(z[11][2]==z[6+i][2]);
        // }
        // for (int i = 1; i < 4; ++i) {
        //     model.addConstr(z[16][0]==z[6+i][0]);
        //     model.addConstr(z[16][1]==z[6+i][1]);
        //     model.addConstr(z[16][2]==z[6+i][2]);
        // }

        // model.write("/home/mzwang/qsp_ws/src/pusher/model" + std::to_string(stepCounter) + ".lp");
        // getchar();
        auto tStart = std::chrono::system_clock::now();
        model.optimize();
        auto tEnd = std::chrono::system_clock::now();

        errors[stepCounter] = error;
        solvingTime[stepCounter] = std::chrono::duration<double>(tEnd - tStart).count();

        Eigen::Vector3d controlNow;
        controlNow << uBar[0][0].get(GRB_DoubleAttr_X) + controlNominal(0), 
                      uBar[0][1].get(GRB_DoubleAttr_X) + controlNominal(1), 
                      uBar[0][2].get(GRB_DoubleAttr_X) + controlNominal(2);
        double mode = z[0][1].get(GRB_DoubleAttr_X)*1 + z[0][2].get(GRB_DoubleAttr_X)*2;
        
        // convert pusher velocity from slider frame to world frame
        J << 1, 0, -yC,
            0, 1, xC;
        B.col(0) = J.transpose() * Eigen::Vector2d(1, 0);
        B.col(1) = J.transpose() * Eigen::Vector2d(0, 1);
        Eigen::MatrixXd Gc(2,3);
        Gc.setZero();
        Gc.leftCols(2) = J * L * B;
        Gc.col(2) << 0.0, -xC/(cos(phi)*cos(phi));
        Eigen::Vector2d vPusher = Gc * controlNow;
        if (vPusher(0) > 0.3){
            vPusher(0) = 0.3;
        }
        if (vPusher(1) > 0.3){
            vPusher(1) = 0.3;
        } else if (vPusher(1) < -0.3){
            vPusher(1) = -0.3;
        }
        geometry_msgs::Twist pusherVel;
        pusherVel.linear.x = cos(sliderPose.theta) * vPusher(0) - sin(sliderPose.theta) * vPusher(1);
        pusherVel.linear.y = sin(sliderPose.theta) * vPusher(0) + cos(sliderPose.theta) * vPusher(1);
        pusherVel.linear.z = 0;
        pusherVel.angular.x = 0;
        pusherVel.angular.y = 0;
        pusherVel.angular.z = 0;

        std::cout << "---------\n";

        // DEBUG
        if (debugInfo){
            nav_msgs::Path pathNomi, pathMPC;
            pathNomi.header.frame_id = "world";
            pathMPC.header.frame_id = "world";

            geometry_msgs::PoseStamped pose;
            pose.header.frame_id = "world";
            for (int i = 0; i < MPCSteps; ++i){
                pose.pose.position.x = stateNominal(i*4);
                pose.pose.position.y = stateNominal(i*4+1);
                pathNomi.poses.push_back(pose);
            }

            for (int i = 0; i < MPCSteps; ++i){
                pose.pose.position.x = state(i*4);
                pose.pose.position.y = state(i*4+1);
                pathMPC.poses.push_back(pose);
            }

            markerPub.publish( marker );
            pathNomiPub.publish(pathNomi);
            pathMPCPub.publish(pathMPC);
        }

        std::cout << error << "\n";
        std::cout << vPusher(0) << ", " << vPusher(1) << "\n";
        std::cout << mode << "\n";
        std::cout << controlNow.transpose() << "\n";
        pusherMovePub.publish(pusherVel);
        ++stepCounter;
    }
    
    // void pubPusherVel(const ros::TimerEvent&){
    //     pusherVelPub.publish(pusherVel);
    // }

private:
    ros::NodeHandle nh;
    ros::Subscriber sliderOdomSub, pusherOdomSub, pusherMoveSub;
    ros::Publisher pusherVelPub, pusherMovePub;
    ros::Publisher pathNomiPub, pathMPCPub, markerPub;
    ros::ServiceClient calculateControlClient;
    ros::Timer timer1, timerVel;

    visualization_msgs::Marker marker;

    geometry_msgs::Pose2D sliderPose;
    geometry_msgs::Twist zeroVel;
    double xC, yC, phi, fMax, mMax, c, mu, m, muGround, width;
    double QFWeight, QWeight, RWeight;
    int MPCSteps, stepCounter, totalSteps;
    double tStep;
    double radius;
    double alphaDot;
    bool debugInfo;

    double Q[4], QFinal[4], R[3];

    Eigen::VectorXd state, control;
    
    // Used to map force to velocity
    Eigen::MatrixXd J, L, B;

    // reference trajectory
    Eigen::VectorXd lineNomi, eightNomi, snakeNomi;
    Eigen::VectorXd lineControlNomi, eightControlNomi, snakeControlNomi;
    Eigen::VectorXd stateNominal, controlNominal;

    // for recording
    Eigen::VectorXd actualState, solvingTime, errors;
    std::ofstream timeFile, errorFile;
    std::string timePath, errorPath;
    GRBEnv env;
};

main(int argc, char *argv[])
{
    ros::init(argc, argv, "controller");
    ros::NodeHandle n("~");
    GRBEnv env = GRBEnv(true);
    env.start();
    Control control(n, env);
    ros::AsyncSpinner spinner(6);
    spinner.start();
    ros::waitForShutdown();
    return 0;
}
#include <ros/ros.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/Twist.h>
#include <geometry_msgs/Pose2D.h>
#include <tf/transform_datatypes.h>
#include <nav_msgs/Path.h>
#include <visualization_msgs/Marker.h>

#include <kortex_driver/SendTwistCommand.h>

#include <yaml-cpp/yaml.h>
#include <math.h>
#include <algorithm>
#include <ifopt/problem.h>
#include <ifopt/snopt_solver.h>
#include "pusher/problem.h"

class Control
{
public:
    Control(ros::NodeHandle n) : nh(n)
    {
        pusherOdomSub = nh.subscribe("/pose/pusher", 1, &Control::pusherOdomCallback, this);
        sliderOdomSub = nh.subscribe("/pose/slider", 1, &Control::sliderOdomCallback, this);
        gen3VelClient = nh.serviceClient<kortex_driver::SendTwistCommand>("/my_gen3/base/send_twist_command");
        pathNomiPub = nh.advertise<nav_msgs::Path>("/path/nominal", 10);
        pathMPCPub = nh.advertise<nav_msgs::Path>("/path/mpc", 10);
        markerPub = nh.advertise<visualization_msgs::Marker>( "/slider_marker", 0 );
        YAML::Node params = YAML::LoadFile("/home/pengchang/build_ws/src/planar_pushing/Config/params.yaml");
        mu = params["mu"].as<double>();
        muGround = params["mu_g"].as<double>();
        m = params["m"].as<double>();
        xC = -params["length"].as<double>() / 2.0;
        mMax = params["mMax"].as<double>();
        fMax = muGround * m * 9.81;
        MPCSteps = params["n_step"].as<int>();
        tStep = params["t_step"].as<double>();
        debugInfo = params["debug_info"].as<bool>();

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

        pusherVel.linear.x = 0;
        pusherVel.linear.y = 0;
        pusherVel.linear.z = 0;
        pusherVel.angular.x = 0;
        pusherVel.angular.y = 0;
        pusherVel.angular.z = 0;

        stepCounter = 0;

        // marker.header.frame_id = "base_link";
        // marker.ns = "my_namespace";
        // marker.id = 0;
        // marker.type = visualization_msgs::Marker::CUBE;
        // marker.action = visualization_msgs::Marker::ADD;
        // marker.scale.x = 0.01;
        // marker.scale.y = 0.01;
        // marker.scale.z = 0.01;
        // marker.color.a = 0.5; // Don't forget to set the alpha!
        // marker.color.r = 0.0;
        // marker.color.g = 0.0;
        // marker.color.b = 1.0;

        // define the nominal trajectory
        float targetLV = 0.02;
        lineNomi.resize(int(4*0.3/(targetLV*tStep))+4);
        for (int i = 0; i < 0.3/(targetLV*tStep); ++i){
            lineNomi.segment(i*4, 4) << 0.496+targetLV*i*tStep, 0.144, 0, 0;
        }
        lineNomi.tail(4) << 0.782, 0, 0, 0;

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

        ros::Duration(1).sleep();
        timer1 = nh.createTimer(ros::Duration(0.05), &Control::findSolution, this);
        // timerVel = nh.createTimer(ros::Duration(0.01), &Control::pubPusherVel, this);
    }

    void pusherOdomCallback(const geometry_msgs::Pose2D &msg) 
    {
        double xSliderToPusher = msg.x - sliderPose.x;
        double ySliderToPusher = msg.y - sliderPose.y;
        yC = -sin(sliderPose.theta) * xSliderToPusher + cos(sliderPose.theta) * ySliderToPusher;
        phi = atan(-yC/xC);
    }

    void sliderOdomCallback(const geometry_msgs::Pose2D &msg) 
    {
        sliderPose.x = msg.x;
        sliderPose.y = msg.y;
        float newTheta = msg.theta;
        // To prevent sudden change from pi to -pi or vice versa
        if (abs(newTheta - sliderPose.theta) < 3){
            sliderPose.theta = newTheta;
        }

        // marker.header.stamp = ros::Time();
        // marker.pose = msg.pose.pose;
        // marker.pose.position.z = 0;
    }
    
    void findSolution(const ros::TimerEvent&)
    {

        Eigen::VectorXd stateNominal(MPCSteps*4), controlNominal((MPCSteps-1)*3);

        for (int i = 0; i < MPCSteps-1; ++i){
            controlNominal.segment(i*3, 3) << 0,0,0;
        }

        // straight line
        if (stepCounter*4+MPCSteps*4 > lineNomi.size()){
            stateNominal.head(MPCSteps*4 - 4) = stateNominal.tail(MPCSteps*4 - 4);
        } else {
            stateNominal = lineNomi.segment(stepCounter*4, MPCSteps*4);
        }
        // 8 shape
        // if (stepCounter*4+MPCSteps*4 > eightNomi.size()){
        //     stepCounter = 0;
        // }
        // stateNominal = eightNomi.segment(stepCounter*4, MPCSteps*4);
        
        // stateNominal.head(4) << sliderPose.x, sliderPose.y, sliderPose.theta, phi;
        std::cout << "---------\n";
        std::cout << "current pose: " << sliderPose.x << ", " << sliderPose.y << ", " << sliderPose.theta << "\n";

        state.head(4) << sliderPose.x, sliderPose.y, sliderPose.theta, phi;

        ifopt::Problem nlp;
        nlp.AddVariableSet(std::make_shared<ifopt::ExVariables>(4*MPCSteps, "state", state));
        nlp.AddVariableSet(std::make_shared<ifopt::ExVariables>(3*(MPCSteps-1), "control", control));
        nlp.AddConstraintSet(std::make_shared<ifopt::ExConstraint>(10*(MPCSteps-1)));
        nlp.AddCostSet(std::make_shared<ifopt::ExCost>("cost", stateNominal, control));

        solver.Solve(nlp);

        Eigen::VectorXd variables = nlp.GetOptVariables()->GetValues();
        
        state.head(4 * MPCSteps - 4) = variables.segment(4, 4 * MPCSteps-4);
        state.tail(4) = variables.segment(4 * MPCSteps-4, 4);
        control.head(3 * MPCSteps - 6) = variables.segment(4 * MPCSteps+3, 3 * MPCSteps-6);
        Eigen::Vector3d controlNow = variables.segment(4 * MPCSteps, 3);
        // control.tail(3) << 0.5, 0, 0;
        
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
        if (vPusher(0) > 0.1){
            vPusher(0) = 0.1;
        }
        if (vPusher(1) > 0.05){
            vPusher(1) = 0.05;
        } else if (vPusher(1) < -0.05){
            vPusher(1) = -0.05;
        }
        float velX = cos(sliderPose.theta) * vPusher(0) - sin(sliderPose.theta) * vPusher(1);
        float velY = sin(sliderPose.theta) * vPusher(0) + cos(sliderPose.theta) * vPusher(1);

        // if not solved, reset control and state
        int status = solver.GetReturnStatus();
        if (status == 13){
            control.setZero();
        }

        std::cout << "reference goal: " << stateNominal.segment(4,4).transpose() << "\n";
        std::cout << "solved goal: " << state.head(4).transpose() << "\n";

        // DEBUG
        if (debugInfo){
            nlp.PrintCurrent();
            nav_msgs::Path pathNomi, pathMPC;
            pathNomi.header.frame_id = "base_link";
            pathMPC.header.frame_id = "base_link";

            geometry_msgs::PoseStamped pose;
            pose.header.frame_id = "base_link";
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

            // markerPub.publish(marker);
            pathNomiPub.publish(pathNomi);
            pathMPCPub.publish(pathMPC);
        }

        double error = sqrt((stateNominal(4) - sliderPose.x) * (stateNominal(4) - sliderPose.x) + 
                       (stateNominal(5) - sliderPose.y) * (stateNominal(5) - sliderPose.y));
        std::cout << error << "\n";
        std::cout << vPusher(0) << ", " << vPusher(1) << "\n";
        std::cout << controlNow.transpose() << "\n";

        std::cout << "xC, yC: " << xC << ", " << yC << "\n";

        // getchar();
        sendVelRequest(velX, velY);
        // ros::Duration(0.05).sleep();
        // sendVelRequest(0, 0);
        ++stepCounter;
    }
    
    void sendVelRequest(float x, float y)
    {
      kortex_driver::SendTwistCommand srv;
      srv.request.input.twist.linear_x = x;
      srv.request.input.twist.linear_y = y;
      srv.request.input.twist.linear_z = 0;
      srv.request.input.twist.angular_x = 0;
      srv.request.input.twist.angular_y = 0;
      srv.request.input.twist.angular_z = 0;
      gen3VelClient.call(srv);
    }

    // void pubPusherVel(const ros::TimerEvent&){
    //     pusherVelPub.publish(pusherVel);
    // }

private:
    ros::NodeHandle nh;
    ros::Subscriber sliderOdomSub, pusherOdomSub;
    ros::Publisher pusherVelPub;
    ros::Publisher pathNomiPub, pathMPCPub, markerPub;
    ros::ServiceClient gen3VelClient;
    ros::Timer timer1, timerVel;

    visualization_msgs::Marker marker;

    geometry_msgs::Pose2D sliderPose;
    geometry_msgs::Twist pusherVel;
    double xC, yC, phi, fMax, mMax, c, mu, m, muGround;
    int MPCSteps, stepCounter;
    double tStep;
    double radius;
    double alphaDot;
    bool debugInfo;

    Eigen::VectorXd state, control;
    
    // Used to map force to velocity
    Eigen::MatrixXd J, L, B;

    // reference trajectory
    Eigen::VectorXd lineNomi, eightNomi, snakeNomi;

    ifopt::SnoptSolver solver;
};

main(int argc, char *argv[])
{
    ros::init(argc, argv, "controller");
    ros::NodeHandle n("~");
    Control control(n);

    ros::AsyncSpinner spinner(6);
    spinner.start();
    ros::waitForShutdown();
    return 0;
}

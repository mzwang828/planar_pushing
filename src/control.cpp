#include <ros/ros.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/Twist.h>
#include <geometry_msgs/Pose2D.h>
#include <tf/transform_datatypes.h>
#include <nav_msgs/Path.h>

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
        pusherOdomSub = nh.subscribe("/odom/pusher", 1, &Control::pusherOdomCallback, this);
        sliderOdomSub = nh.subscribe("/odom/slider", 1, &Control::sliderOdomCallback, this);
        pusherVelPub = nh.advertise<geometry_msgs::Twist>("/cmd_vel", 1);
        pathNomiPub = nh.advertise<nav_msgs::Path>("/path/nominal", 10);
        pathMPCPub = nh.advertise<nav_msgs::Path>("/path/mpc", 10);

        YAML::Node params = YAML::LoadFile("/home/mzwang/qsp_ws/src/pusher/Config/params.yaml");
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

        // define the nominal trajectory
        float targetLV = 0.05;
        lineNomi.resize(int(4*10/(targetLV*tStep))+4);
        for (int i = 0; i < 10/(targetLV*tStep); ++i){
            lineNomi.segment(i*4, 4) << targetLV*i*tStep, 0, 0, 0;
        }
        lineNomi.tail(4) << 10, 0, 0, 0;

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
        timer1 = nh.createTimer(ros::Duration(0.05), &Control::findSolution, this);
        timerVel = nh.createTimer(ros::Duration(0.01), &Control::pubPusherVel, this);
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
    }
    
    void findSolution(const ros::TimerEvent&)
    {
        
        Eigen::VectorXd stateNominal(MPCSteps*4), controlNominal((MPCSteps-1)*3);

        for (int i = 0; i < MPCSteps-1; ++i){
            controlNominal.segment(i*3, 3) << 0,0,0;
        }

        // straight line
        // if (stepCounter*4+MPCSteps*4 > lineNomi.size()){
        //     stateNominal.head(MPCSteps*4 - 4) = stateNominal.tail(MPCSteps*4 - 4);
        // } else {
        //     stateNominal = lineNomi.segment(stepCounter*4, MPCSteps*4);
        // }
        // 8 shape
        if (stepCounter*4+MPCSteps*4 > eightNomi.size()){
            stepCounter = 0;
        }
        stateNominal = eightNomi.segment(stepCounter*4, MPCSteps*4);
        
        // stateNominal.head(4) << sliderPose.x, sliderPose.y, sliderPose.theta, phi;
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
        Eigen::Vector2d vPusher = Gc * control.segment(0, 3);
        if (vPusher(0) > 0.3){
            vPusher(0) = 0.3;
        }
        if (vPusher(1) > 0.3){
            vPusher(1) = 0.3;
        } else if (vPusher(1) < -0.3){
            vPusher(1) = -0.3;
        }
        pusherVel.linear.x = cos(sliderPose.theta) * vPusher(0) - sin(sliderPose.theta) * vPusher(1);
        pusherVel.linear.y = sin(sliderPose.theta) * vPusher(0) + cos(sliderPose.theta) * vPusher(1);

        // if not solved, reset control and state
        // int status = solver.GetReturnStatus();
        // if (status == 13){
        //     control.setZero();
        // }

        std::cout << "---------\n";

        // DEBUG
        if (debugInfo){
            nlp.PrintCurrent();
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

            pathNomiPub.publish(pathNomi);
            pathMPCPub.publish(pathMPC);
        }

        double error = sqrt((stateNominal(4) - sliderPose.x) * (stateNominal(4) - sliderPose.x) + 
                       (stateNominal(5) - sliderPose.y) * (stateNominal(5) - sliderPose.y));
        std::cout << error << "\n";
        std::cout << vPusher(0) << ", " << vPusher(1) << "\n";

        ++stepCounter;
    }
    
    void pubPusherVel(const ros::TimerEvent&){
        pusherVelPub.publish(pusherVel);
    }

private:
    ros::NodeHandle nh;
    ros::Subscriber sliderOdomSub, pusherOdomSub;
    ros::Publisher pusherVelPub, pathNomiPub, pathMPCPub;
    ros::Timer timer1, timerVel;

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

#include <ros/ros.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/Twist.h>
#include <geometry_msgs/Pose2D.h>
#include <tf/transform_datatypes.h>

#include <yaml-cpp/yaml.h>
#include <math.h>
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
        pusherVelPub = nh.advertise<geometry_msgs::Twist>("/cmd_vel", 10);

        YAML::Node params = YAML::LoadFile("/home/mzwang/qsp_ws/src/pusher/Config/params.yaml");
        mu = params["mu"].as<double>();
        muGround = params["mu_g"].as<double>();
        m = params["m"].as<double>();
        xC = -params["length"].as<double>() / 2.0;
        mMax = params["mMax"].as<double>();
        fMax = muGround * m * 9.81;
        c = fMax/mMax;
        MPCSteps = params["n_step"].as<int>();
        tStep = params["t_step"].as<double>();

        J.resize(2,3); 
        B.resize(3,2); 
        L.resize(3,3);
        L << 2/(fMax*fMax), 0, 0,
             0, 2/(fMax*fMax), 0,
             0, 0, 2/(mMax*mMax);
        
        state.resize(4*MPCSteps);
        control.resize(3*MPCSteps);
        for (int i = 0; i < MPCSteps; ++i) {
            state.segment(i*4, 4) << 0.05*i*tStep, 0,0,0;
            control.segment(i*3, 3) << 0.3, 0.0, 0.0;
        }
        
        // circular tracking
        radius = 0.20;
        alphaDot = 0.05/(2*3.14*radius);

        timer1 = nh.createTimer(ros::Duration(0.05), &Control::findSolution, this);
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
        sliderPose.theta = tf::getYaw(pose.getRotation());
    }
    
    /**
     * Pass in initial guess state and control, the nominal state and control and solve the problem
     **/
    void findSolution(const ros::TimerEvent&)
    {
        Eigen::VectorXd stateNominal(MPCSteps*4), controlNominal(MPCSteps*3);

        // straight line
        for (int i = 0; i < MPCSteps; ++i){
            stateNominal.segment(i*4, 4) << sliderPose.x + 0.05*i*tStep, 0,0,0;
            controlNominal.segment(i*3, 3) << 0.0, 0.0, 0.0;
        }

        // circular
        // double alpha = atan2(sliderPose.x, radius - sliderPose.y);
        // for (int i = 0; i < MPCSteps; ++i){
        //     stateNominal.segment(i*4, 4) << radius*sin(alpha + alphaDot*i*tStep), radius - radius*cos(alpha + alphaDot*i*tStep), alpha, 0;
        //     controlNominal.segment(i*3, 3) << 0.3, 0.0, 0.0;
        // }
        
        Eigen::VectorXd stateInit = stateNominal;
        stateInit.segment(0,4) << sliderPose.x, sliderPose.y, sliderPose.theta, phi;
        state.segment(0,4) << sliderPose.x, sliderPose.y, sliderPose.theta, phi;

        std::cout << sliderPose.x << ", " << sliderPose.y << ", " << sliderPose.theta << ", " << phi << "\n";
        std::cout << stateNominal.tail(4).transpose() << "\n";

        ifopt::Problem nlp;
        nlp.AddVariableSet(std::make_shared<ifopt::ExVariables>(4*MPCSteps, "state", state));
        nlp.AddVariableSet(std::make_shared<ifopt::ExVariables>(3*MPCSteps, "control", control));
        nlp.AddConstraintSet(std::make_shared<ifopt::ExConstraint>(10*(MPCSteps-1)));
        nlp.AddCostSet(std::make_shared<ifopt::ExCost>("cost", stateNominal, controlNominal));

        solver.Solve(nlp);
        Eigen::VectorXd variables = nlp.GetOptVariables()->GetValues();
        
        state = variables.segment(0, 4 * MPCSteps);
        control = variables.segment(4 * MPCSteps, 3 * MPCSteps);

        // convert pusher velocity from slider frame to world frame
        geometry_msgs::Twist pusherVel;

        J << 1, 0, -yC,
            0, 1, xC;
        B.col(0) = J.transpose() * Eigen::Vector2d(1, 0);
        B.col(1) = J.transpose() * Eigen::Vector2d(0, 1);
        Eigen::MatrixXd Gc(2,3);
        Gc.setZero();
        Gc.leftCols(2) = J * L * B;
        Gc.col(2) << 0.0, -xC/(cos(phi)*cos(phi));
        Eigen::Vector2d vPusher = Gc * control.segment(3, 3);
        if (vPusher(0) > 0.2){
            vPusher(0) = 0.2;
        }
        if (vPusher(1) > 0.2){
            vPusher(1) = 0.2;
        } else if (vPusher(1) < -0.2){
            vPusher(1) = -0.2;
        }
        pusherVel.linear.x = cos(sliderPose.theta) * vPusher(0) - sin(sliderPose.theta) * vPusher(1);
        pusherVel.linear.y = sin(sliderPose.theta) * vPusher(0) + cos(sliderPose.theta) * vPusher(1);

        // std::cout << control(3) << ", " << control(4) << ", " << control(5) << "\n";
        // std::cout << vPusher(0) << ", " << vPusher(1) << "\n";
        // std::cout << alpha << ", " << alpha + alphaDot*19*tStep << ", " << sliderPose.theta << "\n";
        std::cout << "---------\n";

        pusherVelPub.publish(pusherVel);
    }

    // void straightLine(const ros::TimerEvent&)
    // {
    //     Eigen::VectorXd stateNominal(MPCSteps*4), controlNominal(MPCSteps*3);
    //     for (int i = 0; i < MPCSteps; ++i){
    //         stateNominal.segment(i*4, 4) << sliderPose.x + 0.05*i*tStep, 0,0,0;
    //         controlNominal.segment(i*3, 3) << 0.3, 0.0, 0.0;
    //     }

    //     this->findSolution(state, control, stateNominal, controlNominal);

    //     // convert pusher velocity from slider frame to world frame
    //     geometry_msgs::Twist pusherVel, zeroVel;
    //     zeroVel.linear.x = 0.0;
    //     zeroVel.linear.y = 0.0;

    //     J << 1, 0, -yC,
    //         0, 1, xC;
    //     B.col(0) = J.transpose() * Eigen::Vector2d(1, 0);
    //     B.col(1) = J.transpose() * Eigen::Vector2d(0, 1);
    //     Eigen::MatrixXd Gc(2,3);
    //     Gc.leftCols(2) = J * L * B;
    //     Gc.col(2) << 0.0, -xC/(cos(phi)*cos(phi));
    //     Eigen::Vector2d vPusher = Gc * control.segment(3, 3);
    //     pusherVel.linear.x = cos(sliderPose.theta) * vPusher(0) - sin(sliderPose.theta) * vPusher(1);
    //     pusherVel.linear.y = sin(sliderPose.theta) * vPusher(0) + cos(sliderPose.theta) * vPusher(1);

    //     std::cout << control(3) << ", " << control(4) << ", " << control(5) << "\n";
    //     std::cout << vPusher(0) << ", " << vPusher(1) << "\n";
    //     std::cout << "---------\n";
    //     // getchar();

    //     pusherVelPub.publish(pusherVel);
    //     // ros::Duration(tStep).sleep();
    //     // pusherVelPub.publish(zeroVel);
    // }

private:
    ros::NodeHandle nh;
    ros::Subscriber sliderOdomSub, pusherOdomSub;
    ros::Publisher pusherVelPub;
    ros::Timer timer1;

    geometry_msgs::Pose2D sliderPose;
    double xC, yC, phi, fMax, mMax, c, mu, m, muGround;
    int MPCSteps;
    double tStep;
    double radius;
    double alphaDot;

    Eigen::VectorXd state, control;
    
    // Used to map force to velocity
    Eigen::MatrixXd J, L, B;

    ifopt::SnoptSolver solver;
};

main(int argc, char *argv[])
{
    ros::init(argc, argv, "controller");
    ros::NodeHandle n("~");
    Control control(n);

    ros::AsyncSpinner spinner(4);
    spinner.start();
    ros::waitForShutdown();
    return 0;
}

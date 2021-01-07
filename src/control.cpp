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
        px = params["length"].as<double>() / 2.0;
        mMax = params["mMax"].as<double>();
        fMax = muGround * m * 9.81;
        c = fMax/mMax;
    }

    void pusherOdomCallback(const nav_msgs::Odometry &msg) 
    {
        double xSliderToPusher = msg.pose.pose.position.x - sliderPose.x;
        double ySliderToPusher = msg.pose.pose.position.y - sliderPose.y;
        py = -sin(sliderPose.theta) * xSliderToPusher + cos(sliderPose.theta) * ySliderToPusher;
        yt = (mu*c*c - px*py + mu*px*px)/(c*c + py*py - mu*px*py);
        yb = (-mu*c*c - px*py - mu*px*px)/(c*c + py*py + mu*px*py);
    }

    void sliderOdomCallback(const nav_msgs::Odometry &msg) 
    {
        sliderPose.x = msg.pose.pose.position.x;
        sliderPose.y = msg.pose.pose.position.y;
        tf::Pose pose;
        tf::poseMsgToTF(msg.pose.pose, pose);
        sliderPose.theta = tf::getYaw(pose.getRotation());
    }

    void findSolution(const Eigen::VectorXd &stateIn, const Eigen::VectorXd &controlIn, int MPCSteps)
    {
        ifopt::Problem nlp;
        Eigen::VectorXd initState = stateIn;
        initState.segment(0,4) << sliderPose.x, sliderPose.y, sliderPose.theta, py;
        nlp.AddVariableSet(std::make_shared<ifopt::ExVariables>(4*MPCSteps, "state", initState));
        nlp.AddVariableSet(std::make_shared<ifopt::ExVariables>(2*MPCSteps, "control", controlIn));
        nlp.AddConstraintSet(std::make_shared<ifopt::ExConstraint>(4*3*(MPCSteps-1)));
        nlp.AddCostSet(std::make_shared<ifopt::ExCost>("cost", initState, controlIn));

        ifopt::SnoptSolver solver;
        solver.Solve(nlp);
        Eigen::VectorXd variables = nlp.GetOptVariables()->GetValues();
        Eigen::Map<Eigen::MatrixXd> state(variables.segment(0, 4 * MPCSteps).data(), 4, MPCSteps);
        Eigen::Map<Eigen::MatrixXd> control(variables.segment(4 * MPCSteps, 2 * MPCSteps).data(), 2, MPCSteps);
        
        // convert pusher velocity from slider frame to world frame
        geometry_msgs::Twist pusherVel;
        pusherVel.linear.x = cos(sliderPose.theta) * control(0,0) - sin(sliderPose.theta) * control(1,0);
        pusherVel.linear.y = sin(sliderPose.theta) * control(0,0) + cos(sliderPose.theta) * control(1,0);
        pusherVelPub.publish(pusherVel);
    }

    // for test use
    void straightLine()
    {
        int MPCSteps = 35;
        Eigen::VectorXd stateNominal(MPCSteps*4), controlNominal(MPCSteps*2);
        for (int i = 0; i < MPCSteps; ++i){
            stateNominal.segment(i*4, 4) << sliderPose.x + 0.05*i*0.03, 0,0,0;
            controlNominal.segment(i*2, 2) << 0.05, 0;
        }

        this->findSolution(stateNominal, controlNominal, MPCSteps);
    }

private:
    ros::NodeHandle nh;
    ros::Subscriber sliderOdomSub, pusherOdomSub;
    ros::Publisher pusherVelPub;

    geometry_msgs::Pose2D sliderPose;
    double px, py, fMax, mMax, c, mu, m, muGround;
    double yt, yb;
};

main(int argc, char *argv[])
{
    ros::init(argc, argv, "controller");
    ros::NodeHandle n("~");
    Control control(n);

    ros::spin();
    return 0;
}

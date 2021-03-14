#!/usr/bin/env python3
import gurobipy as gp 
from gurobipy import GRB 
import numpy as np

import rospy
from geometry_msgs.msg import Twist
from nav_msgs.msg import Odometry
from geometry_msgs.msg import Pose2D
import yaml

# parameters
with open ('/home/mzwang/qsp_ws/src/pusher/Config/params.yaml') as pfile:
    params = yaml.full_load(pfile)
mu = params.get('mu')
muG = params.get('mu_g')
m = params.get('m')
length = params.get('length') # along x direction 
width = params.get('width') # along y direction
mMax = params.get('mMax') # numerical integrated
fMax = mu*m*9.81
xC = - length/2
# params for MPC
nStep = params.get('n_step')
tStep = params.get('t_step')
QFinal = params.get('QFinal') * np.diag(np.array([3,3,0.1,0]))
Q = params.get('Q') * np.diag(np.array([3,3,0.1,0]))
R = params.get('R') * np.diag(np.array([1,1,0.01]))

L = np.array([[2/(fMax*fMax), 0, 0], [0, 2/(fMax*fMax), 0], [0, 0, 2/(mMax*mMax)]])
# define the nominal trajectory
targetLV = 0.05
totalStep = int(10/(targetLV*tStep))
lineNomi = np.zeros((4, totalStep))
lineControlNomi = np.zeros((3,totalStep))
for i in range(totalStep):
    lineNomi[0, i] = targetLV*i*tStep
lineControlNomi[0, :] = 0.305

sliderPose = Pose2D()
sliderPose.theta = 0
yC = 0
currentPhi = 0
stepCounter = 0

class Control:
    def __init__(self):
        rospy.Subscriber("/odom/pusher", Odometry, self.pusherOdomCallback)
        rospy.Subscriber("/odom/slider", Odometry, self.sliderOdomCallback)
        rospy.Subscriber("/move_pusher", Twist, self.pusherMoveCallback)
        self.pusherMovePub = rospy.Publisher("/move_pusher", Twist, queue_size=1)
        self.pusherVelPub = rospy.Publisher("/cmd_vel", Twist, queue_size=1)
        rospy.Timer(rospy.Duration(0.01), self.find_solution)
        self.zeroVel = Twist()
        self.zeroVel.linear.x = 0
        self.zeroVel.linear.y = 0
        self.zeroVel.linear.y = 0
        self.zeroVel.angular.x = 0
        self.zeroVel.angular.y = 0
        self.zeroVel.angular.z = 0

    def pusherMoveCallback(self, data):
        self.pusherVelPub.publish(data)
        rospy.sleep(0.05)
        self.pusherVelPub.publish(self.zeroVel)

    def pusherOdomCallback(self, data):
        global yC, currentPhi
        xSliderToPusher = data.pose.pose.position.x - sliderPose.x
        ySliderToPusher = data.pose.pose.position.y - sliderPose.y
        yC = -np.sin(sliderPose.theta) * xSliderToPusher + np.cos(sliderPose.theta) * ySliderToPusher
        currentPhi = np.arctan(-yC/xC)

    def sliderOdomCallback(self, data):
        global sliderPose
        sliderPose.x = data.pose.pose.position.x
        sliderPose.y = data.pose.pose.position.y
        qx = data.pose.pose.orientation.x
        qy = data.pose.pose.orientation.y
        qz = data.pose.pose.orientation.z
        qw = data.pose.pose.orientation.w
        t3 = +2.0 * (qw * qz + qx * qy)
        t4 = +1.0 - 2.0 * (qy * qy + qz * qz)
        newTheta = np.arctan2(t3,t4)
        if (abs(newTheta - sliderPose.theta) < 3):
            sliderPose.theta = newTheta

    def find_solution(self, event):
        global stepCounter
        # Get nominal trajectory
        stateNomi = lineNomi[:, stepCounter:stepCounter+nStep]
        controlNomi = lineControlNomi[:, stepCounter:stepCounter+nStep-1]
        # Set bounds
        xBarLB = np.empty((4, nStep))
        xBarLB[0, :] = -GRB.INFINITY
        xBarLB[1, :] = -GRB.INFINITY
        xBarLB[2, :] = -GRB.INFINITY
        xBarLB[3, :] = -0.35 - stateNomi[3, :] # so the pusher not exceeds the slider's edge
        xBarUB = np.empty((4, nStep))
        xBarUB[0, :] = GRB.INFINITY
        xBarUB[1, :] = GRB.INFINITY
        xBarUB[2, :] = GRB.INFINITY
        xBarUB[3, :] = 0.35 - stateNomi[3, :]
        uBarLB = np.empty((3, nStep-1))
        uBarLB[0, :] = -controlNomi[0, :]
        uBarLB[1, :] = -GRB.INFINITY
        uBarLB[2, :] = -GRB.INFINITY
        uBarUB = np.empty((3, nStep-1))
        uBarUB[0, :] = GRB.INFINITY
        uBarUB[1, :] = GRB.INFINITY
        uBarUB[2, :] = GRB.INFINITY
        # set current error state x_bar_0
        xBarLB[:, 0] = np.array([sliderPose.x-stateNomi[0,0], sliderPose.y-stateNomi[1,0], sliderPose.theta-stateNomi[2,0], currentPhi - stateNomi[3,0]])
        xBarUB[:, 0] = np.array([sliderPose.x-stateNomi[0,0], sliderPose.y-stateNomi[1,0], sliderPose.theta-stateNomi[2,0], currentPhi - stateNomi[3,0]])
        # pusher velocity bounds
        pusherVLB = np.array([0, -0.3])
        pusherVUB = np.array([0.3, 0.3])
        # define the QP
        gm = gp.Model("qp")
        xBar = gm.addMVar((4, nStep), lb = xBarLB, ub = xBarUB)
        uBar = gm.addMVar((3, nStep-1), lb = uBarLB, ub = uBarUB)
        cost = xBar[:, -1]@QFinal@xBar[:, -1]
        for i in range(nStep - 1):
            cost = cost + xBar[:, i]@Q@xBar[:, i] + uBar[:, i]@R@uBar[:, i]
        gm.setObjective(cost)

        for i in range(nStep-1):
            fn = controlNomi[0,i]
            ft = controlNomi[1,i]
            pdot = controlNomi[2,i]
            theta = stateNomi[2,i]
            phi = stateNomi[3,i]

            # dynamics
            A = np.array([[0, 0, -2*fn*np.sin(theta)/fMax**2 - 2*ft*np.cos(theta)/fMax**2, 0], 
                        [0, 0, 2*fn*np.cos(theta)/fMax**2 - 2*ft*np.sin(theta)/fMax**2, 0], 
                        [0, 0, 0, 2*fn*xC*(np.tan(phi)**2 + 1)/mMax**2], 
                        [0, 0, 0, 0]])
            B = np.array([[2*np.cos(theta)/fMax**2, -2*np.sin(theta)/fMax**2, 0], 
                        [2*np.sin(theta)/fMax**2, 2*np.cos(theta)/fMax**2, 0], 
                        [2*xC*np.tan(phi)/mMax**2, 2*xC/mMax**2, 0], 
                        [0, 0, 1]])
            gm.addConstr(xBar[:,i+1] - xBar[:,i] == tStep*(A@xBar[:,i] + B@uBar[:,i]))
            # velocity mapping
            pusherV = np.array([fn*(2*xC**2*np.tan(phi)**2/mMax**2 + 2/fMax**2) + 2*ft*xC**2*np.tan(phi)/mMax**2, 
                                2*fn*xC**2*np.tan(phi)/mMax**2 + ft*(2*xC**2/mMax**2 + 2/fMax**2) - pdot*xC/np.cos(phi)**2])
            dVdx = np.array([[0, 0, 0, 2*fn*xC**2*(2*np.tan(phi)**2 + 2)*np.tan(phi)/mMax**2 + 2*ft*xC**2*(np.tan(phi)**2 + 1)/mMax**2], 
                            [0, 0, 0, 2*fn*xC**2*(np.tan(phi)**2 + 1)/mMax**2 - 2*pdot*xC*np.sin(phi)/np.cos(phi)**3]])
            dVdu = np.array([[2*xC**2*np.tan(phi)**2/mMax**2 + 2/fMax**2, 2*xC**2*np.tan(phi)/mMax**2, 0], 
                            [2*xC**2*np.tan(phi)/mMax**2, 2*xC**2/mMax**2 + 2/fMax**2, -xC/np.cos(phi)**2]])
            gm.addConstr(pusherV + dVdx@xBar[:,i] + dVdu@uBar[:,i] >= pusherVLB)
            gm.addConstr(pusherV + dVdx@xBar[:,i] + dVdu@uBar[:,i] <= pusherVUB)
            # friction cone
            gm.addConstr(ft + uBar[1,i] <= mu*(fn + uBar[0,i]))
            gm.addConstr(ft + uBar[1,i] >= -mu*(fn + uBar[0,i]))
            # STC
            H1 = np.array([mu*pdot, -pdot, fn*mu - ft])
            H2 = np.array([mu*pdot, pdot, fn*mu + ft])
            if pdot != 0:
                gm.addConstr(np.min(-pdot, 0)*(ft-mu*fn) + H1@uBar[:,i] == 0)
                gm.addConstr(np.min(pdot, 0)*(ft+mu*fn) + H2@uBar[:,i] == 0)
            
        gm.optimize()
        nextControl = uBar.X[:,0] + controlNomi[:,0]
        J = np.array([[1,0,-yC], [0,1,xC]])
        B = J.transpose()
        Gc = np.zeros([2,3])
        Gc[:,0:2] = J@L@B
        Gc[:,2] = np.array([0, -xC/(np.cos(currentPhi)**2)])
        vPusher = Gc@nextControl
        pusherVel = Twist()
        pusherVel.linear.x = np.cos(sliderPose.theta) * vPusher[0] - np.sin(sliderPose.theta) * vPusher[1]
        pusherVel.linear.y = np.sin(sliderPose.theta) * vPusher[0] + np.cos(sliderPose.theta) * vPusher[1];
        self.pusherMovePub.publish(pusherVel)
        stepCounter += 1
        print('------')
        error = np.sqrt(xBarLB[0,0]**2 + xBarUB[1,0]**2)
        print(error)
        print(vPusher[0],', ',vPusher[1])


if __name__ == "__main__":
    rospy.init_node('control_gurobi')
    try:
        control = Control()
        rospy.spin()
    except rospy.ROSInterruptException:  pass
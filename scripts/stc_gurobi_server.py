#!/usr/bin/env python3
import gurobipy as gp 
from gurobipy import GRB 
import numpy as np

from pusher.srv import CalculateControl, CalculateControlResponse
import rospy
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


def handle_calculate_control(req):
    # Get nominal trajectory
    stateNomi = np.array(req.state_nominal)
    stateNomi = np.reshape(stateNomi, (4, nStep), 'F')
    controlNomi = np.array(req.control_nominal)
    controlNomi = np.reshape(controlNomi, (3, nStep-1), 'F')
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
    xBarLB[:, 0] = np.array(req.state_error)
    xBarUB[:, 0] = np.array(req.state_error)
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
    nextControl = tuple(uBar.X[:,0] + controlNomi[:,0])
    return CalculateControlResponse(nextControl)

def calculate_control_server():
    rospy.init_node('calculate_control_server')
    s = rospy.Service('/calculate_control', CalculateControl, handle_calculate_control)
    print("Service is ready!")
    rospy.spin()

if __name__ == "__main__":
    calculate_control_server()

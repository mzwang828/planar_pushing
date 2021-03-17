import gurobipy as gp 
from gurobipy import GRB 
import numpy as np

# parameters
mu = 0.60
muG = 0.35
m = 1.05
length = 0.09 # along x direction 
width = 0.09 # along y direction
mMax = 0.360 # numerical integrated
fMax = mu*m*9.81
xC = - length/2
# params for MPC
nStep = 20
tStep = 0.05
QFinal = 2000.0 * np.diag(np.array([3,3,0.1,0]))
Q = 10.0 * np.diag(np.array([3,3,0.1,0]))
R = 0.5 * np.diag(np.array([1,1,0.01]))

# define reference trajectory
targetLV = 0.05
stateNomi = np.zeros((4, nStep))
controlNomi = np.zeros((3, nStep-1))
for i in range(nStep):
    stateNomi[0, i] = targetLV*tStep*i
controlNomi[0, :] = 0.3


xBarLB = np.empty((4, nStep))
xBarLB[0, :] = -GRB.INFINITY
xBarLB[1, :] = -GRB.INFINITY
xBarLB[2, :] = -GRB.INFINITY
xBarLB[3, :] = -0.70 - stateNomi[3, :] # so the pusher not exceeds the slider's edge
xBarUB = np.empty((4, nStep))
xBarUB[0, :] = GRB.INFINITY
xBarUB[1, :] = GRB.INFINITY
xBarUB[2, :] = GRB.INFINITY
xBarUB[3, :] = 0.70 - stateNomi[3, :]
uBarLB = np.empty((3, nStep-1))
uBarLB[0, :] = -controlNomi[0, :]
uBarLB[1, :] = -GRB.INFINITY
uBarLB[2, :] = -GRB.INFINITY
uBarUB = np.empty((3, nStep-1))
uBarUB[0, :] = GRB.INFINITY
uBarUB[1, :] = GRB.INFINITY
uBarUB[2, :] = GRB.INFINITY

# set current error state x_bar_0
xBarLB[:, 0] = np.array([0,0.03,0,-0])
xBarUB[:, 0] = np.array([0,0.03,0,-0])

# pusher velocity bounds
pusherVLB = np.array([0, -0.3])
pusherVUB = np.array([0.3, 0.3])

gm = gp.Model("qp")
xBar = gm.addMVar((4, nStep), lb = xBarLB, ub = xBarUB, name="state")
uBar = gm.addMVar((3, nStep-1), lb = uBarLB, ub = uBarUB, name="control")
z = gm.addMVar((3, nStep-1), vtype = GRB.BINARY, name="mode")
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
    # mode selection
    M = 1.0
    # stick
    gm.addConstr(pdot+uBar[2,i] >= M*(z[0,i] - 1))
    gm.addConstr(pdot+uBar[2,i] <= M*(-z[0,i] + 1))
    # up
    gm.addConstr(pdot+uBar[2,i] >= (5*M*(z[1,i] - 1)))
    gm.addConstr(mu*(fn+uBar[0,i]) - (ft+uBar[1,i]) >= M*(z[1,i]-1))
    gm.addConstr(mu*(fn+uBar[0,i]) - (ft+uBar[1,i]) <= M*(-z[1,i]+1))
    # down
    gm.addConstr(pdot+uBar[2,i] <= (5*M*(-z[2,i] + 1)))
    gm.addConstr(mu*(fn+uBar[0,i]) + (ft+uBar[1,i]) >= M*(z[2,i]-1))
    gm.addConstr(mu*(fn+uBar[0,i]) + (ft+uBar[1,i]) <= M*(-z[2,i]+1))
    
    gm.addConstr(z[0,i] + z[1,i] + z[2,i] == 1)

gm.write("py.lp")
gm.optimize()
print('runtime:', gm.Runtime)

gm.printQuality()
print(gm.Status)
print('state:')
# print(xBar.X)
print(xBar.X + stateNomi)
print('control:')
# print(uBar.X)
print(uBar.X + controlNomi)
print('mode: ')
print(z.X)

fn = uBar.X[0,0] + controlNomi[0,0]
ft = uBar.X[1,0] + controlNomi[1,0]
pdot = uBar.X[2,0] + controlNomi[2,0]

pusherV = np.array([fn*(2*xC**2*np.tan(phi)**2/mMax**2 + 2/fMax**2) + 2*ft*xC**2*np.tan(phi)/mMax**2, 
                        2*fn*xC**2*np.tan(phi)/mMax**2 + ft*(2*xC**2/mMax**2 + 2/fMax**2) - pdot*xC/np.cos(phi)**2])

print(pusherV)

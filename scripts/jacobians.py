from sympy import sin, cos, tan, Matrix
from sympy import symbols
import numpy as np

# motion model
x, y, theta, phi, fn, ft, pdot, xC, fMax, mMax, mu = symbols('x y theta phi fn ft pdot xC fMax mMax mu')
state = Matrix([x, y, theta, phi])
control = Matrix([fn, ft, pdot])
L = Matrix(3,3,[2/(fMax*fMax), 0, 0, 0, 2/(fMax*fMax), 0, 0, 0, 2/(mMax*mMax)])
R = Matrix(3,3,[cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 1])
J = Matrix(2,3,[1, 0, xC * tan(phi), 0, 1, xC])
B = J.T
dynamics = Matrix.zeros(4,3)
dynamics[0:3, 0:2] = R*L*B
dynamics[3,2] = 1
dynamics = dynamics * control
# dfdx
# print(dynamics.jacobian(state))
# dfdu
# print(dynamics.jacobian(control))



# velocity mapping
Gc = Matrix.zeros(2,3)
Gc[:, 0:2] = J*L*B
Gc[:, 2] = Matrix([0, -xC/cos(phi)**2])
pusherV = Gc * control
# print(pusherV)
# dGcdx
print(pusherV.jacobian(state))
# dGcdu
# print(pusherV.jacobian(control))

# stc
stc1 = Matrix([-pdot*(ft - mu*fn)])
stc2 = Matrix([pdot*(ft + mu*fn)])
# print(stc1.jacobian(control))
# print(stc2.jacobian(control))
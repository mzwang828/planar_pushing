#!/usr/bin/env python
import yaml
import numpy as np
import rospy
import roslaunch
from gazebo_msgs.msg import ModelState
from gazebo_msgs.srv import SetModelState
from std_srvs.srv import Empty
from tf.transformations import quaternion_from_euler
from geometry_msgs.msg import Quaternion

with open ('/home/mzwang/qsp_ws/src/pusher/Config/initials.yaml') as pfile:
    params = yaml.safe_load(pfile)
initials = params.get('init')

set_state = rospy.ServiceProxy('/gazebo/set_model_state', SetModelState)
reset_world = rospy.ServiceProxy('/gazebo/reset_world', Empty)
count = 1

pkg = 'pusher'
executable = 'control_snopt'
node = roslaunch.core.Node(pkg,executable)
launch = roslaunch.scriptapi.ROSLaunch()
launch.start()

for init in initials:
    x = init[0]
    y = init[1]
    theta = init[2]
    phi = init[3]

    pusher_x = np.cos(theta) * (-0.0551) - np.sin(theta) * 0.0551 * np.tan(phi)
    pusher_y = np.sin(theta) * (-0.0551) + np.cos(theta) * 0.0551 * np.tan(phi)
    pusher_x += x
    pusher_y += y

    pusher_state = ModelState()
    pusher_state.model_name = 'pusher'
    pusher_state.pose.position.x = pusher_x
    pusher_state.pose.position.y = pusher_y
    pusher_state.pose.position.z = 0.05
    pusher_state.pose.orientation.x = 0
    pusher_state.pose.orientation.y = 0
    pusher_state.pose.orientation.z = 0
    pusher_state.pose.orientation.w = 1

    slider_state = ModelState()
    slider_state.model_name = 'slider'
    slider_state.pose.position.x = x
    slider_state.pose.position.y = y
    slider_state.pose.position.z = 0.05
    quat = quaternion_from_euler(0.0, 0.0, theta)
    slider_state.pose.orientation = Quaternion(*quat)

    try:
        print(init)
        reset_resp = reset_world()
        pusher_resp = set_state(pusher_state)
        slider_resp = set_state(slider_state)
    except rospy.ServiceException as e:
        print(count, " Set state failed: %s"%e)
    
    process = launch.launch(node)
    rospy.sleep(15)
    process.stop()
    



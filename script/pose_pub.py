#!/usr/bin/env python
# license removed for brevity
import rospy
import tf
from std_msgs.msg import String
from geometry_msgs.msg import Pose2D

def posePub():
    sliderPub = rospy.Publisher('/pose/slider', Pose2D, queue_size=10)
    pusherPub = rospy.Publisher('/pose/pusher', Pose2D, queue_size=10)
    rospy.init_node('odomPub', anonymous=True)
    rate = rospy.Rate(50) # 10hz
    sliderPose = Pose2D();
    pusherPose = Pose2D();
    listener = tf.TransformListener() 
    while not rospy.is_shutdown():
        try:
            (trans,rot) = listener.lookupTransform('/base_link', '/tag1', rospy.Time(0))
        except (tf.LookupException, tf.ConnectivityException, tf.ExtrapolationException):
            continue
        quaternion = (rot[0], rot[1], rot[2], rot[3])
        euler = tf.transformations.euler_from_quaternion(quaternion)
        sliderPose.x = trans[0]
        sliderPose.y = trans[1]
        sliderPose.theta = euler[2]
        
        try:
            (trans1,rot1) = listener.lookupTransform('/base_link', '/finger', rospy.Time(0))
        except (tf.LookupException, tf.ConnectivityException, tf.ExtrapolationException):
            continue
        pusherPose.x = trans1[0]
        pusherPose.y = trans1[1]
        pusherPose.theta = 0

        sliderPub.publish(sliderPose)
        pusherPub.publish(pusherPose)
        rate.sleep()

if __name__ == '__main__':
    try:
        posePub()
    except rospy.ROSInterruptException:
        pass
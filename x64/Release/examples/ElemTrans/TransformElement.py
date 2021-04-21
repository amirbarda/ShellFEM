import numpy as np
import random 
from math import pi, cos, sin, sqrt

N = 20 #how many copies

def get_x_roatation(theta):
	return np.matrix([[1, 0, 0],
                      [0, cos(theta), -sin(theta)],
                      [0, sin(theta), cos(theta)]])
    
def get_y_roatation(theta):
	return np.matrix([[cos(theta), 0, sin(theta)],
                      [0, 1, 0],
                      [-sin(theta), 0, cos(theta)]])

def get_z_roatation(theta):
	return np.matrix([[cos(theta), -sin(theta), 0],
                      [sin(theta), cos(theta), 0],
                      [0, 0, 1]])

def get_random_rotation():
	alpha = random.uniform(0, 2*pi)
	beta = random.uniform(0, 2*pi)
	gamma = random.uniform(0, 2*pi)

	mat_x = get_x_roatation(alpha)
	mat_y = get_y_roatation(beta)
	mat_z = get_z_roatation(gamma)

	rotation = mat_x.dot(mat_y.dot(mat_z))
	return rotation


v_lst = [np.array([7.878, 0.0, 8.09]), np.array([5.866, 0.369, 8.09]), np.array([5.832, 0.737, 8.09])]

for i in range(N):
	rotation = get_random_rotation()
	transformation = np.random.rand(3,1)
	print(transformation)
	for v in v_lst:
		tmp = rotation.dot(v.transpose())
		tmp += transformation.transpose()
		print(f"v {tmp.item(0):.9f} {tmp.item(1):.9f} {tmp.item(2):.9f}")

print()
for i in range(1, N + 1):
    print(f"f {i*3 -2} {i*3 -1} {i*3 -0}")

import numpy as np
import random 
from math import pi, cos, sin, sqrt

N = 20 #how many copies
A = 20 #magnitude for radial translation

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

f = open("mesh.obj", "r")
tmp = f.readlines()
v_lst = []
f_lst = []

for i in range(len(tmp)):
	if tmp[i].startswith("v "):
		tmp[i] = tmp[i][len("v "):]
		v_lst.append(np.fromstring(tmp[i], dtype=float, sep=' '))
	if tmp[i].startswith("f "):
		tmp[i] = tmp[i][len("f "):]
		f_lst.append(np.fromstring(tmp[i], dtype=int, sep=' '))

for i in range(N):
	rotation = get_random_rotation()
	transformation = np.random.rand(3,1)
	for v in v_lst:
		tmp = rotation.dot(v.transpose())
		tmp += A * transformation.transpose()
		print(f"v {tmp.item(0):.9f} {tmp.item(1):.9f} {tmp.item(2):.9f}")

print()

for i in range(N):
	k = i * len(v_lst)
	for f in f_lst:
		print(f"f {k + f.item(0)} {k + f.item(1)} {k + f.item(2)}")


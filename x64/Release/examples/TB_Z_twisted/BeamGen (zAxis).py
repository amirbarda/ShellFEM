import numpy as np

pi = np.pi

h = 1.1		#width [diameter]
r = h/2
L = 12		#length of beam [on z axis]
N = 12			#amount of plates for beam
Nr = 1
ang = pi/2

for m in range(N + 1):
	z = (m/N) * L
	print("v", 0, 0, f"{z/Nr:.9f}")

for m in range(N + 1):
	z = (m/N) * L
	phi = ang * (z / L)
	x = r * np.cos(phi)
	y = r * np.sin(phi)
	print("v", f"{x/Nr:.9f}", f"{y/Nr:.9f}", f"{z/Nr:.9f}")

for m in range(N + 1):
	z = (m/N) * L
	phi = pi + ang * (z / L)
	x = r * np.cos(phi)
	y = r * np.sin(phi)
	print("v", f"{x/Nr:.9f}", f"{y/Nr:.9f}", f"{z/Nr:.9f}")

print()
for m in range(1, N + 1):
	print("f", m, m + 1, m + N + 1)
	print("f", m + 1, m + N + 2, m + N + 1)

for m in range(1, N + 1):
	print("f", m + (2 * N) + 2, m + (2 * N) + 3, m)
	print("f", m + (2 * N) + 3, m + 1, m)

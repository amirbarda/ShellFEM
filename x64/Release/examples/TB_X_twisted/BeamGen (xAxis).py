import numpy as np

pi = np.pi

h = 1.1e-2		#width [diameter]
r = h/2
L = 12e-2		#length of beam [on x axis]
N = 12			#amount of plates for beam
ang = pi/2

for m in range(N + 1):
	x = (m/N) * L
	phi = pi + ang * (x / L)
	y = r + r * np.cos(phi)
	z = r * np.sin(phi)
	print("v", f"{x:.9f}", f"{y:.9f}", f"{z:.9f}")

for m in range(N + 1):
	x = (m/N) * L
	print("v", f"{x:.9f}", f"{r:.9f}", 0)

for m in range(N + 1):
	x = (m/N) * L
	phi = ang * (x / L)
	y = r + r * np.cos(phi)
	z = r * np.sin(phi)
	print("v", f"{x:.9f}", f"{y:.9f}", f"{z:.9f}")

print()
for m in range(1, N + 1):
	print("f", m, m + 1, m + N + 1)
	print("f", m + 1, m + N + 2, m + N + 1)

for m in range(1, N + 1):
	print("f", m + N + 1, m + N + 2, m + 2*N + 2)
	print("f", m + N + 2, m + 2*N + 3, m + 2*N + 2)

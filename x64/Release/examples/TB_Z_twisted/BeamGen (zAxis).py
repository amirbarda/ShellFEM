import numpy as np

h = 1.1				#width [diameter]
r = 0.55
L = 12			#length [on x axis]
pi = np.pi
ang = pi / 2
#ang = 0

for m in range(L + 1):
	print("v", 0, 0, m)

for m in range(L + 1):
	z = m
	phi = ang * (z / L)
	x = r * np.cos(phi)
	y = r * np.sin(phi)
	if (abs(x) < 10**(-7)): x = 0
	print("v", x, y, z)

for m in range(L + 1):
	z = m
	phi = pi + ang * (z / L)
	x = r * np.cos(phi)
	y = r * np.sin(phi)
	if (abs(y) < 10**(-7)): y = 0
	print("v", x, y, z)

print()
for m in range(1, L + 1):
	print("f", m, m + 1, m + L + 1)
	print("f", m + 1, m + L + 2, m + L + 1)

for m in range(1, L + 1):
	print("f", m + (2 * L) + 2, m + (2 * L) + 3, m)
	print("f", m + (2 * L) + 3, m + 1, m)

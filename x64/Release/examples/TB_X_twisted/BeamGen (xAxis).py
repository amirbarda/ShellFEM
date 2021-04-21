import numpy as np

h = 1.1				#width [diameter]
r = h/2
L = 12				#length [on x axis]
pi = np.pi
ang = pi / 2
Nr = 1000
#ang = 0

for m in range(L + 1):
	x = m
	phi = pi + ang * (x / L)
	y = r + r * np.cos(phi)
	z = r * np.sin(phi)
	if (abs(z) < 10**(-7)): z = 0
	print("v", f"{x/Nr:.9f}", f"{y/Nr:.9f}", f"{z/Nr:.9f}")

for m in range(L + 1):
	print("v", f"{m/Nr:.9f}", f"{r/Nr:.9f}", 0)

for m in range(L + 1):
	x = m
	phi = ang * (x / L)
	y = r + r * np.cos(phi)
	z = r * np.sin(phi)
	if (abs(y) < 10**(-7)): y = 0
	print("v", f"{x/Nr:.9f}", f"{y/Nr:.9f}", f"{z/Nr:.9f}")

print()
for m in range(1, L + 1):
	print("f", m, m + 1, m + L + 1)
	print("f", m + 1, m + L + 2, m + L + 1)

for m in range(1, L + 1):
	print("f", m + L + 1, m + 1 + L + 1, m + L + 1 + L + 1)
	print("f", m + 1 + L + 1, m + L + 2 + L + 1, m + L + 1 + L + 1)

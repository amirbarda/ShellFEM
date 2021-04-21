import numpy as np


M = 10  #latitude lines (horizontal).	[even though we compute a quarter sphere, we assign the amount of lines for a whole sphere] [M = 10 in article]
N = 16  #longitude lines (vertical).	[even though we compute a quarter sphere, we assign the amount of lines for a whole sphere] [N = 16 in article]
r = 10
pi = np.pi 

for m in range(M + 1):
	for n in range(N):
		x = r * np.sin(pi * m/M) * np.cos(2 * pi * n/N)
		y = r * np.sin(pi * m/M) * np.sin(2 * pi * n/N)
		z = r * np.cos(pi * m/M)
		if (0 <= z < r):		
			print("v", round(x, 3), round(y, 3), round(z ,3))

print()

for i in range(M//2 -1):
	k = i * N
	for j in range(1, N + 1):
		#print("\n(i,j) = (", i, j, ")")
		if (((j -1) % (N//4)) < (N//8)):
			print("f", j + k, j + k + N, j + k + 1)
			print("f", j + k + 1, j + k + N, j + k + N + 1)
		else:
		 	print("f", j + k, j + k + N, (j % N) + k + N + 1)
		 	print("f", j + k, (j % N) + k + N + 1, (j % N) + k + 1)

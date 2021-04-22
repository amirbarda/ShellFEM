import numpy as np
pi = np.pi 

# [r - radius] 
# [M - #latitude lines (horizontal)] 
# [N - #longitude lines (vertical)] 
# [MM - how many faces to draw by latitude lines]
def GenSphere(r, M, N, MM, top, bottom):
	LAST = 0	#count how many vertices where generated.
	for m in range(M + 1):
		for n in range(N):
			x = r * np.sin(pi * m/M) * np.cos(2 * pi * n/N)
			y = r * np.sin(pi * m/M) * np.sin(2 * pi * n/N)
			z = r * np.cos(pi * m/M)
			#if (0 <= z < r):		
			if (r * np.cos(pi * (MM +1)/M) < z < r):
				LAST += 1
				print("v", round(x, 3), round(y, 3), round(z ,3))

	if (top == True):
		print("v", round(0, 3), round(0, 3), round(r ,3))
	if (bottom == True):
		print("v", round(0, 3), round(0, 3), round(-r ,3))

	print()

	for i in range(MM -1):
		k = i * N
		for j in range(1, N + 1):
			#print("\n(i,j) = (", i, j, ")")
			if (((j -1) % (N//4)) < (N//8)):
				print("f", j + k, j + k + N, j + k + 1)
				print("f", j + k + 1, j + k + N, j + k + N + 1)
			else:
			 	print("f", j + k, j + k + N, (j % N) + k + N + 1)
			 	print("f", j + k, (j % N) + k + N + 1, (j % N) + k + 1)

	if (top == True):
		for i in range(1, N +1):
			print("f", LAST +1, i, (i % N) +1)
	if (bottom == True):
		for i in range(1, N +1):
			print("f", LAST +2, LAST +1 - i, LAST - (i % N))

#GenSphere(10, 10, 16, 10 -1, True, True) #sphere
#GenSphere(10, 10, 16, 10//2, True, True) #dreidel
#GenSphere(10, 10, 16, 10//2, True, False) #hemisphere
#GenSphere(10, 10, 16, 10//2, False, False) #punctured hemisphere
GenSphere(10, 100, 16, 78, False, False)

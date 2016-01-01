import numpy as np 
import matplotlib.pyplot as plt

f = lambda x: 10*(x-1)**2 -1
x = np.linspace(1,2,10)
x_0 = 1 + 1./3
x_1 = 1 + 2./3
x_2 = 1 + 3./4

psi = [lambda x: 1,lambda x: x, lambda x: x*x]
points = [x_0, x_1, x_2]

A = np.zeros((len(points),len(psi)))
b = np.zeros(len(points))

for i in range(len(points)):
	for j in range(len(psi)):
		A[i][j] = psi[j](points[i]) 
	b[i] = f(points[i])

c = np.linalg.solve(A,b)

plt.figure(1)
plt.plot(c[0]+c[1]*psi[1](x)+c[2]*psi[2](x))
plt.plot(f(x))
#plt.axis([0, 8, -1, 9])
plt.show()
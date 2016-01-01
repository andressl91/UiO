import numpy as np
import matplotlib.pyplot as plt 



N = 10
x = np.linspace(0,1,N+1)
dx = x[1]-x[0]
A = np.identity(N)
b = np.zeros(N)

def f(x):
	return x*(1-x)

def psi(x,i):
	return x**i

b[2] = f(2*dx); b[4] = f(4*dx)

for i in range(N):
	A[i][i] = psi(i*dx,1)

print A
#c = np.linalg.solve(A,b)

plt.figure(1)
#plt.plot(c)
plt.plot(f(x))
#plt.show()
import numpy as numpy

N = 10
x = np.zeros(N+1)
t = np.zeros(N+1)
y = np.zeros(N+1)

dt = 0.1
dx = 0.2
c = (dx/dt)**2
for i in range(N):
	q1 = 0.5*(q[i-1]+q[i])
	u[i+1] = -u[i-1]+2*u[i] + c*2*q1*(u[i-1]-u[i]) + c*f(i*dx)


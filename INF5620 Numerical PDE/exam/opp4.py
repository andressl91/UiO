import numpy as np
import matplotlib.pyplot as plt

N = 80
x = np.linspace(0,1,N+1)
dx = x[1]-x[0]
eps = 0.5
alpha = 1-2*eps/dx
beta = eps/dx - 1

#Centered
alpha1 = (1+2*eps/dx) ;beta1 = (2*eps/dx -1)

mat = np.zeros((N+1,N+1)) ;mat1 = np.zeros((N+1,N+1))
b = np.zeros((N+1,1)) ;b1 = np.zeros((N+1,1))
for i in range(1,N):
	#Upwind
	mat[i][i-1] = eps/dx
	mat[i][i] = alpha
	mat[i][i+1] = beta
	
	#centered
	mat1[i][i-1] = alpha1
	mat1[i][i] = -4*eps/dx
	mat1[i][i+1] = beta1
	

mat[0][0] = 1;mat[-1][-1] = 1;

mat1[0][0] = 1 ;mat1[-1][-1] = 1;

b[-1] = 1 ;b1[-1] = 1

u = np.linalg.solve(mat, b)
u1 = np.linalg.solve(mat1, b)


exact = lambda x: (1-np.exp(x/eps))/(1-np.exp(1./eps))
error = np.max(u[:][0]-exact(x))
error1 = np.max(u1[:][0]-exact(x))
print error, error1



plt.figure(1)
plt.plot(exact(x),'g')
plt.plot(u)
plt.plot(u1, 'r')
plt.show()
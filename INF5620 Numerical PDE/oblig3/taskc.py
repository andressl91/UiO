import numpy as np
import matplotlib.pyplot as plt

def exact(x):
	return x*(0.5*x-1)

def approx(x,i):
	return -16./(np.pi*(2*i+1))**3 * np.sin((2*i+1)*np.pi*x/2.)

N = [0,1,20]
Nx = 50
x = np.linspace(0,1,Nx+1)
u = np.zeros(Nx+1)


plt.figure(1)
for i in N:
	for j in range(i+1):
		u[:] += approx(x[:],j)
	plt.plot(x,u,label='$N={i}$'.format(i=i))
	
	u[:]=0
plt.plot(x,exact(x),'black',label='$Exact solution$')
plt.legend(loc='best')
plt.title("Solve $u''=1$ with basis $sin(2i+1)x\pi/2$ ")
plt.xlabel("$x$")
plt.ylabel("$u(x)$")
plt.savefig("Taskc.png")
plt.show()



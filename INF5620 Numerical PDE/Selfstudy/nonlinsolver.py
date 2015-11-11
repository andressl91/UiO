import numpy as np

def F(u,a,b,c):
	return a*u**2 + b*u + c

T = 5 #Time
Nt = 500 #Number of time steps
k = 1 #Iteration counter
t = np.linspace(0, T, Nt+1)
dt = t[1]-t[0]
u = np.zeros(Nt+1)
u_exact = np.zeros(Nt+1)
#Initiate start value
u0 = 5
u[0] = u_exact
u0[0] = u0

#second order coefficients

for i in range(1,Nt+1):
	a = dt; b = (1.0-dt); c = -u[i-1] #c is the last calculated term
	u_ = u[i-1] #This is the guess for the iteration, good guess last iteration result
	while F(u_,a,b,c) > 1E-10 and k < 500:
		print "here"
		u_ = -c/(a*u_ + b) 
		k += 1

	u[i] = u_
	u_exact[i] = -b + np.sqrt(b**2 - 4*a*c)/(2.0*a)
import matplotlib.pyplot as plt 	

plt.figure(1)
plt.plot(t, u_exact, 'g')
plt.plot(t, u, 'r')
plt.show()
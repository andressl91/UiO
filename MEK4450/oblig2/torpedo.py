rho = 1000 #Density
S = 1.5 #Reference area 
Cd = 0.2 #Drag Coefficient
m = 50*10**3 #Mass of Torpedo anchors
ma = 5*10**3 #Added mass
g = 9.81 #Gravitation acceleration

import numpy as np
N = 100
T = 20
t = np.linspace(0,T,N+1)
dt = t[1]-t[0]
v = np.zeros(N+1)
x = np.zeros(N+1)
v[0] = 0
x[0] = 0
crit = np.zeros(N+1)
crittime = np.zeros(N+1)

for i in range(N):
	v[i+1] = v[i] - dt*(0.5*rho*Cd*(v[i])**2/(m+ma) - m*g/(m+ma)) 
	x[i+1] = x[i] + dt*v[i+1]


crit = np.abs(v[:]-50).argmin()

import matplotlib.pyplot as plt

fig = plt.figure(1)
fig.suptitle('Torpedo-anchor velocity')
ax = fig.add_subplot(111)
ax.set_title("Critical velocity reached at %.2f meters after %.1f seconds" % (x[crit],t[crit]))
ax.set_xlabel('Time seconds')
ax.set_ylabel('Veloity m/s')
plt.axhline(y=v[crit],color="red")
plt.plot(t,v)
plt.savefig("velocity.jpg")
#plt.show()


import numpy as np



def f(L,m):
	return 1/(2*np.pi)*np.sqrt(200000*10**3/L/m)

L = 30
m = 50000

x = np.linspace(0,200,201)
print x
import matplotlib.pyplot as plt 
plt.figure()
plt.plot(x, 1./f(x,m))
plt.show()
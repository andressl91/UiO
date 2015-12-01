import numpy as np

b = 2.05
f = 0.3
ma = 2.65

rt = 80 #True formation resistivity
rw = 0.2 #Formation water resistivity
		

def por(b,j,ma):
	return (b-ma)/(f-ma)

def sat(po):
	return 1.0/po*np.sqrt(rw/rt)

p = por(b,f,ma)
s = sat(p)

print p, s
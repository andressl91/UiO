import numpy as np

b = 1.95
f = 0.8
ma = 2.7

rt = 2*10**1 #Formation water resistivity
rw = 190*10**-3 #True formation resistivity


def por(b,j,ma):
	return (b-ma)/(f-ma)

def sat(po):
	return 1.0/po*np.sqrt(rw/rt)

p = por(b,f,ma)
s = sat(p)

print p, s
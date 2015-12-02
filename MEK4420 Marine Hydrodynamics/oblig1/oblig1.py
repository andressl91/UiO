import numpy as np
import matplotlib.pyplot as plt

def points(N, a, b):	
	if a == b:
		o = 2.0*np.pi*a
	else:
		o = 2*np.pi
	
	dx = o/(N)
	x = np.zeros(N)
	y = np.zeros(N)
	for i in range(N):
		x[i] = a*np.cos(dx*i)
		y[i] = b*np.sin(dx*i)
 	return x, y
"""
x, y = points(100,2,2)
plt.figure()
plt.axis([-3,3,-3,3])
plt.plot(x,y)
plt.show()
"""
def theta(N, a, b, j):
	x, y = points(N, a, b)
	theta = np.zeros(N) 
	initx = (x[j+1] + x[j])*0.5
	inity = (y[j+1] + y[j])*0.5
	theta = np.zeros(N)

	theta[:-1] = np.arccos(((x[:-1]-initx)*(x[1:]-initx) + (y[:-1]-inity)*(y[1:]-inity))\
	/(np.sqrt((x[:-1]-initx)**2+(y[:-1]-inity)**2)\
	*np.sqrt((x[1:]-initx)**2+(y[1:]-inity)**2)))

	xa = x[len(theta)-1] - initx; xb = x[0] - initx
	ya = y[len(theta)-1] - inity; yb = y[0] - inity
	theta[len(theta)-1] = np.arccos((xa*xb + ya*yb)/(np.sqrt(xa**2 + ya**2) * np.sqrt(xb**2 + yb**2)))
	return theta

def integrate(N, eq, direction,a,b):
	area = 0
	x, y = eq
	integral = np.zeros(N)
	if direction == 11 or 66:
		num = x
		n1_a = -((x[:-2]+x[1:-1])/(2*a**2))/np.sqrt(((x[1:-1]+x[:-2])/(2*a**2))**2 + ((y[1:-1]+y[:-2])/(2*b**2))**2)
		n2_a = -((y[:-2]+y[1:-1])/(2*b**2))/np.sqrt(((x[1:-1]+x[:-2])/(2*a**2))**2 + ((y[1:-1]+y[:-2])/(2*b**2))**2)
		rx_a = 0.5*(x[:-2]+x[1:-1]); ry_a = 0.5*(y[:-2]+y[1:-1]); 

		n1_b = -((x[1:-1]+x[2:])/(2*a**2))/np.sqrt(((x[2:]+x[1:-1])/(2*a**2))**2 + ((y[2:]+y[1:-1])/(2*b**2))**2)
		n2_b = -((y[1:-1]+y[2:])/(2*b**2))/np.sqrt(((x[2:]+x[1:-1])/(2*a**2))**2 + ((y[2:]+y[1:-1])/(2*b**2))**2)
		rx_b = 0.5*(x[1:-1]+x[2:]); ry_b = 0.5*(y[1:-1]+y[2:]); 
	
	if direction == 22:
		num = y

	x_0 = np.zeros(N-1); y_0 = np.zeros(N-1)
	x_0[:] = 0.5*(x[:-1]+x[1:])
	y_0[:] = 0.5*(y[:-1]+y[1:])

	ds = np.sqrt((x[1:-1]-x[:-2])**2 + (y[1:-1]-y[:-2])**2)
	na = -(num[:-2]+num[1:-1])/(2*a**2)/np.sqrt(((x[:-2]+x[1:-1])/(2*a**2))**2 + ((y[:-2]+y[1:-1])/(2*b**2))**2)
	nb = -(num[1:-1]+num[2:])/(2*a**2)/np.sqrt(((x[1:-1]+x[2:])/(2*a**2))**2 + ((y[1:-1]+y[2:])/(2*b**2))**2)


	for i in range(N-1): #N-1
		radius_a = np.sqrt((x_0[i] -x[:-2])**2 + (y_0[i] - y[:-2])**2)
		radius_b = np.sqrt((x_0[i] -x[1:-1])**2 + (y_0[i] - y[1:-1])**2)
		
		if direction == 66:
			
			n_a = np.cross([n1_a[i],n2_a[i]],[rx_a[i],ry_a[i]])
			n_b = np.cross([n1_b[i],n2_b[i]],[rx_b[i],ry_b[i]])

			f_a = np.log(radius_a[:]) * na
			f_b = np.log(radius_b[:]) * nb
			integral[i] = np.sum(0.5*(f_a[:] + f_b[:])*ds[:])

	
		else:
			f_a = np.log(radius_a[:]) * na[:]
			f_b = np.log(radius_b[:]) * nb[:]
			integral[i] = np.sum(0.5*(f_a[:] + f_b[:])*ds[:])
	return integral

a = 2.; b=2.; N = 6; direction = 22
#print points(N,a,b)
print integrate(N, points(N,a,b),direction,a,b)

def matrix(N, a, b):
	mat = np.zeros((N,N))
	for i in range(N):
		if i == N-1:
			mat[i][:] = -theta(N, a, b, -1)
		else:
			mat[i][:] = -theta(N, a, b, i)
	return mat

def added(sol, N, a,b, direction):
	x, y = points(N,a,b)
	area = 0
	if direction == 11:
		num = x*a**-2
		for i in range(N-2):
			na = -0.5*(num[i]+num[i+1])/np.sqrt((((x[i+1]+x[i])/(2*a**2))**2+((y[i+1]+y[i])/(2*b**2))**2))
			nb = -0.5*(num[i+1]+num[i+2])/np.sqrt((((x[i+2]+x[i+1])/(2*a**2))**2+((y[i+2]+y[i+1])/(2*b**2))**2))
			ds = np.sqrt((x[i+1]-x[i])**2 + (y[i+1]-y[i])**2)
			area = area + 0.5*(sol[i]*na+sol[i+1]*nb)*ds
		return area
	if direction == 22:
		numerator = y*b**-2
		for i in range(N-1):
			n1 = -numerator[i]/np.sqrt((x[i+1]**2+y[i+1]**2))
			ds = np.sqrt((x[i+1]-x[i])**2 + (y[i+1]-y[i])**2)
			area = area + 0.5*(sol[i]*n1+sol[i+1]*n1)*ds
		return area
	
	if direction == 66:
		for i in range(N-2):
			n1a = -((x[i]+x[i+1])/(2*a**2))/np.sqrt(((x[i+1]+x[i])/(2*a**2))**2 + ((y[i+1]+y[i])/(2*b**2))**2)
			n2a = -((y[i]+y[i+1])/(2*b**2))/np.sqrt(((x[i+1]+x[i])/(2*a**2))**2 + ((y[i+1]+y[i])/(2*b**2))**2)
			rxa = 0.5*(x[i]+x[i+1]); rya = 0.5*(y[i]+y[i+1]); 
			na = np.cross([n1a,n2a],[rxa,rya])

			n1b = -((x[i+1]+x[i+2])/(2*a**2))/np.sqrt(((x[i+2]+x[i+1])/(2*a**2))**2 + ((y[i+2]+y[i+1])/(2*b**2))**2)
			n2b = -((y[i+1]+y[i+2])/(2*b**2))/np.sqrt(((x[i+2]+x[i+1])/(2*a**2))**2 + ((y[i+2]+y[i+1])/(2*b**2))**2)
			rxb = 0.5*(x[i+1]+x[i+2]); ryb = 0.5*(y[i+1]+y[i+2]); 
			nb = np.cross([n1b,n2b],[rxb,ryb])


			ds1 = np.sqrt((x[i+1]-x[i])**2 + (y[i+1]-y[i])**2)
			ds2 = np.sqrt((x[i+2]-x[i+1])**2 + (y[i+2]-y[i+1])**2)
			area = area + 0.5*(sol[i]*na*ds1 + sol[i+1]*nb*ds2)
		return area


def result(N, a, b, direction):
	 return np.linalg.solve(matrix(N,a,b), integrate(N, points(N,a,b),direction,a,b))

def main():
	N = 500
	a = 2.0
	b = 2.0
	angle = np.zeros(N)
	direction = 11
	
	solve = result(N, a, b, direction)
	mass =  added(solve, N,a,b, direction)
	#exact = 4.753 square
	
	if a == b:
		if direction == 66:
			exact = 0
			print "Exact value: %.5f, Calculated mass: %.5f" % (exact,mass)
		else:
			for i in range(N):
				angle[i] = i*2*np.pi/N 
			error = max(-np.cos(angle)-solve)/max(-np.cos(angle))*100
			masse = float((a**2*np.pi-mass)/np.pi * 100)
			plt.plot(solve, 'r')
			plt.plot(-np.cos(angle))
			plt.title("Numerical solution agains analytical circle solution, error %.3f %% \n Added mass Calculated %.2f  error = %1.3f %% " %(error, mass, masse)) 
			plt.show()
	
	else:
		if direction == 66:
			exact = (np.pi/8.0)*(a**2-b**2)**2
			masse = float(exact - mass)/(exact) * 100.0
			print "Exact value: %.5f, Calculated mass: %.5f, error: %.5f %%" % (exact,mass, masse)
		else :
			exact = np.pi*b**2
			masse = float((np.pi*b**2 - mass)/(np.pi*b**2) * 100.0)
			print "Exact value: %.5f, Calculated mass: %.5f, error: %.5f %%" % (exact,mass, masse)

if __name__ == "__main__":
	main()
	print 

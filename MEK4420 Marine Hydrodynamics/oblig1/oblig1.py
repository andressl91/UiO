import numpy as np
import matplotlib.pyplot as plt

class Figure:
	def __init__(self, r_a, r_b, N):
		self.r_a = r_a ;self.r_b = r_b
		self.N = N; 
		self.x = None; self.y = None
		self.mat = None
	def figure_points(self):
		N = self.N; r_a = self.r_a ; r_b = self.r_b 
		dx = 2*np.pi/(N)
		x = np.zeros(N+1)
		y = np.zeros(N+1)
		for i in range(N+1):
			x[i] = r_a*np.cos(dx*i)
			y[i] = r_a*np.sin(dx*i)
			self.x = x; self.y = y
		return x,y

	def angles(self):
		N = self.N; 
		x, y = self.figure_points()
		mat = np.zeros((N,N))
		theta = np.zeros(N)
		for i in range(N):
			initx = (x[i+1] + x[i])*0.5
			inity = (y[i+1] + y[i])*0.5
			test = ((x[:-1]-initx)*(x[1:]-initx) + (y[:-1]-inity)*(y[1:]-inity))\
			/(np.sqrt((x[:-1]-initx)**2+(y[:-1]-inity)**2)\
			*np.sqrt((x[1:]-initx)**2+(y[1:]-inity)**2)) 
			#print i
			theta = np.arccos(test[:])
	
			mat[i][:i+1] = theta[:(i+1)]
			mat[i][i] = -np.pi
			mat[i][i+1:] = theta[i+1:]
		self.mat = mat
		return mat
	def integrate(self, direction):
		area = 0
		N = self.N
		x = self.x; y = self.y
		a = self.r_a; b = self.r_b
		integral = np.zeros(N)
		x_0 = 0.5*(x[:-1]+x[1:])
		y_0 = 0.5*(y[:-1]+y[1:])
		ds = np.sqrt((x[1:]-x[:-1])**2 + (y[1:]-y[:-1])**2)
		if direction == 11:
			num = x #PASS PAA RETNING
			n = -(num[:-1]+num[1:])/(2*a**2)/np.sqrt(((x[:-1]+x[1:])/(2*a**2))**2 + ((y[:-1]+y[1:])/(2*b**2))**2)
		if direction == 22:
			num = y
			n = -(num[:-1]+num[1:])/(2*a**2)/np.sqrt(((x[:-1]+x[1:])/(2*a**2))**2 + ((y[:-1]+y[1:])/(2*b**2))**2)
		if direction == 66:
			nx = (x_0/(a**2))/np.sqrt((x_0**2/(a**4)) + ((y_0**2)/(b**4)) )
			ny = (y_0/(a**2))/np.sqrt((x_0**2/(a**4)) + ((y_0**2)/(b**4)) )
			rx = x_0
			ry = y_0
			n = x_0*ny - y_0*nx
			for i in range(N):
				radius_a = np.sqrt((x_0[i] -x[:-1])**2 + (y_0[i] - y[:-1])**2)
				radius_b = np.sqrt((x_0[i] -x[1:])**2 + (y_0[i] - y[1:])**2)

				f_a = np.log(radius_a[:]) * nx#na
				f_b = np.log(radius_b[:]) * nx#nb
				integral[i] = np.sum(0.5*(f_a[:] + f_b[:])*ds[:])
			return integral

		for i in range(N): #N-1
			radius_a = np.sqrt((x_0[i] -x[:-1])**2 + (y_0[i] - y[:-1])**2)
			radius_b = np.sqrt((x_0[i] -x[1:])**2 + (y_0[i] - y[1:])**2)
			f_a = np.log(radius_a[:]) * n[:]
			f_b = np.log(radius_b[:]) * n[:]
			integral[i] = np.sum(0.5*(f_a[:] + f_b[:])*ds[:])
		return integral

	def added(self, sol, direction):
		x = self.x; y = self.y
		a = self.r_a; b = self.r_b
		area = 0
		ds = np.sqrt((x[1:]-x[:-1])**2 + (y[1:]-y[:-1])**2)
		if direction == 11:
			num = x*a**-2
			n = -0.5*(num[:-1]+num[1:])/np.sqrt((((x[:-1]+x[1:])/(2*a**2))**2+((y[:-1]+y[1:])/(2*b**2))**2))
			
		if direction == 22:
			num = y*a**-2
			n = -0.5*(num[:-1]+num[1:])/np.sqrt((((x[:-1]+x[1:])/(2*a**2))**2+((y[:-1]+y[1:])/(2*b**2))**2))

		if direction == 66:
			x_0 = 0.5*(x[:-1]+x[1:])
			y_0 = 0.5*(y[:-1]+y[1:])
			nx = -(x_0/(b**2)) / np.sqrt((x_0**2/(a**4)) + ((y_0**2)/(b**4)) )
			ny = -(y_0/(a**2)) / np.sqrt((x_0**2/(a**4)) + ((y_0**2)/(b**4)))
			rx = x_0
			ry = y_0
			n = ny*rx - nx*ry

		area = np.sum((sol)*ds*n)	
		return area

fig = Figure(3,2,1000)
mat = fig.angles()
direction = 11
integ = fig.integrate(direction)
test = np.linalg.solve(mat, integ)
print fig.added(test,direction)
print np.pi*2**2


def main():
	N = 500
	a = 2.0
	b = 2.0
	angle = np.zeros(N)
	direction = 22
	
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
			plt.title("Numerical solution agains analytical circle solution, error %.3f %% \n Added mass Calculated %.5f  error = %1.3f %% " %(error, mass, masse)) 
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
	#main()
	print 
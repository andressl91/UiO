import numpy as np
import matplotlib.pyplot as plt

def simulate(beta, Theta, epsilon, num_periods, time_steps_per_period, plot):
	t = num_periods * 2*np.pi
	dt = 2*np.pi/time_steps_per_period
	x = np.zeros(num_periods * time_steps_per_period+1)
	y = np.zeros(num_periods * time_steps_per_period+1)
	time = np.zeros(num_periods *  time_steps_per_period+1)
	angle = np.zeros(num_periods * time_steps_per_period +1)


	x[0] = (1 + epsilon)* np.sin(Theta)
	y[0] = 1 - (1 - epsilon)*np.cos(Theta)
	dimL0 = np.sqrt(x[0]**2 + (y[0]-1)**2) 
	angle[0] = np.arctan(x[0]/(1 - y[0]))
	x[1] = 0.5*(x[0]*(-dt**2*beta/(1 - beta)*(1 - beta/dimL0) + 2))
	y[1] = 0.5*(-dt**2*beta/(1 - beta)*(1 - beta/dimL0)*(y[0]-1) - dt**2*beta + 2*y[0])

	for i in range(1, num_periods * time_steps_per_period):
		angle[i+1]  = np.arctan(x[i]/(1 - y[i]))
		dimL = np.sqrt(x[i]**2 + (y[i]-1)**2)
		x[i+1] = x[i]*(-dt**2*beta/(1 - beta)*(1 - beta/dimL) + 2) - x[i-1]
		y[i+1] = -dt**2*beta/(1 - beta)*(1 - beta/dimL)*(y[i]-1) - dt**2*beta + 2*y[i] - y[i-1]
		time[i+1] = i*dt

	if plot == True:
		plt.figure(1)
		plt.title("Pendulum with  beta = %1.0f and epsilon %1.0f" % (beta, epsilon))
		plt.gca().set_aspect('equal')
		plt.plot(x, y )
		plt.show()

	return x, y, time

def vertical(plot, beta, epsilon):
	#Vertical checks if the simulate function is accurate for 
	#only vertical movement, against analytical solution for 
	#vibration equation. 
	x_0, y_0, time = simulate(beta, Theta = 0, epsilon = 0.1, 
									num_periods = 2, time_steps_per_period = 400, 
									plot = False)
	ex = epsilon * np.cos(np.sqrt(beta/(1-beta))*time) 
	if plot == True:
		plt.figure(2)
		plt.xlabel("Time")
		plt.plot(time, ex)
		plt.plot(time, y_0, '*')
		plt.show()

def demo():
	x_0, y_0, time = simulate(beta = 0.9, Theta = 30, epsilon =0, 
									num_periods = 3, time_steps_per_period = 600, 
									plot = False)
	
	plt.figure(1)
	plt.title("Demo of pendulum, beta = 0.9, Theta = 30")
	plt.gca().set_aspect('equal')
	plt.plot(x_0, y_0)
	plt.show()


def validate():
	#Validate checks if simulate returns a still pendulum,
	#for Theta, epsilon = 0, aka initial angle and stretch of
	#pendulum
	x_0, y_0, time = simulate(beta = 0.99, Theta = 0, epsilon = 0, 
									num_periods = 4, time_steps_per_period = 60, 
									plot = False)

	if int(np.sum(x_0)) == 0 and int(np.sum(y_0)) == 0:
		print "Pendulum is not moving"
	else:
		print "ERROR, Pendulum should not be moving"


def main():
	#Beta value close to 1, gives ordinary pendulum without elastisy.. 
	
	simulate(beta = 0.99, Theta = 30, epsilon = 0, 
			num_periods = 6, time_steps_per_period = 60, 
			plot = True)
	validate()

	#vertical(plot = True, beta = 0.9, epsilon = 0.1)
	demo()

if __name__ == '__main__':
	main()

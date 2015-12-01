import numpy as np

def density(pressure):
	return 0.73*pressure + 0.123

def diff_press(p,del_h,f,L,D,u):
	hydro = -density(p)*9.81*del_h 
	fric = f*L/D*0.5*density(p)*u*u
	return hydro + fric
def Re(u,l,my):
	return u*l/my

def drop(p_init):
	#pressure drop from Reservoir to well
	del_h = 4000. #Elevation change of control volume
	L = 4000. #Length of control volume
	D = 5*0.0254 #Diameter of controlvolume from inch to m
	e = 5E-5 #wall roughtness
	u = 9. #velocity
	p = p_init#Pressure Bar #Pascal ? 10**6
	my = 2E-5 #Dynamic viscosity
	f = (-1.8*np.log10((e/D*3.7)**1.11 + 6.9/Re(u,D,my)))**-2 #Friction
	p_well = p + diff_press(p,del_h,2*f,L,D,u)/10**6
	#print p_well

	#pressure drop from well to topside
	del_h = 1000. #Elevation change of control volume
	L = 50.*10**3 #Length of control volume
	D = 12*0.0254 #Diameter of controlvolume from inch to m
	e = 5E-5 #wall roughtness
	#Conservation of mass calc
	# = density(p_init)/density(p_well)\
	#	* 5*0.0254/(12*0.0254)*u #2. #velocity
	u = 3.
	p = p_well	#Pressure Bar
	my = 2E-5 #Dynamic viscosity
	f = (-1.8*np.log10((e/D*3.7)**1.11 + 6.9/Re(u,D,my)))**-2 #Friction
	p_top = p + diff_press(p,del_h,f,L,D,u)/10**6
	return p_top

T = 23
#T = 27
pressure = [360-15*i for i in range(T)]
time = [i for i in range(T)]
pressure = np.asarray(pressure)
time = np.asarray(time)
print drop(pressure)

##################################################
"""
Q = np.asarray([30, 25, 20, 17, 14])
Density = np.zeros(len(Q))
u_vel = np.zeros(len(Q))
D = 5*0.0254 #Diameter of controlvolume from inch to m
well, top = drop(pressure)
Density[:] = density(pressure[:-4:4])
u_vel[:] = Q[:]/Density[:]*(4./(np.pi*D**2))
print u_vel
#print pressure
#print pressure[:-4:4]

well, top = drop(pressure)
#Q = np.asarray([60, 50, 40, 34, 28])
Q = np.asarray([60, 60, 60, 60, 60])
D = 12*0.0254 #Diameter of controlvolume from inch to m
Density[:] = density(well[:-4:4])
u_vel[:] = Q[:]/Density[:]*(4./(np.pi*D**2))
print u_vel
"""
##################################################
import matplotlib.pyplot as plt 
fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.set_title("Plotting pressure at reciever")
ax.set_xlabel('Time   Years')
ax.set_ylabel('Pressure  Bar')
#plt.plot(time, pressure)
plt.plot(time, drop(pressure))
plt.axhline(y=50,color="red")
#plt.axis([1, 22, 0, 360])
plt.savefig("Pressure.png")
plt.show()
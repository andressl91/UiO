"""
Exercise 13

This code is based on wave1D_dn_vc.py by HP Langtangen. The following is
his documentation for wave1D_dn_vc.py:
"""
"""
1D wave equation with Dirichlet or Neumann conditions
and variable wave velocity::
 u, x, t, cpu = solver(I, V, f, c, U_0, U_L, L, dt, C, T,
                       user_action=None, version='scalar',
                       stability_safety_factor=1.0)
Solve the wave equation u_tt = (c**2*u_x)_x + f(x,t) on (0,L) with
u=U_0 or du/dn=0 on x=0, and u=u_L or du/dn=0
on x = L. If U_0 or U_L equals None, the du/dn=0 condition
is used, otherwise U_0(t) and/or U_L(t) are used for Dirichlet cond.
Initial conditions: u=I(x), u_t=V(x).
T is the stop time for the simulation.
dt is the desired time step.
C is the Courant number (=max(c)*dt/dx).
stability_safety_factor enters the stability criterion:
C <= stability_safety_factor (<=1).
I, f, U_0, U_L, and c are functions: I(x), f(x,t), U_0(t),
U_L(t), c(x).
U_0 and U_L can also be 0, or None, where None implies
du/dn=0 boundary condition. f and V can also be 0 or None
(equivalent to 0). c can be a number or a function c(x).
user_action is a function of (u, x, t, n) where the calling code
can add visualization, error computations, data analysis,
store solutions, etc.
"""
import numpy as np
import matplotlib.pyplot as plt

def solver(I, V, f, c, U_0, U_L, L, dt, C, T, u_exact, solve,
           user_action=None, version='scalar',
           stability_safety_factor=1.0):
    """Solve u_tt=(c^2*u_x)_x + f on (0,L)x(0,T]."""
    Nt = int(round(T/dt))
    t = np.linspace(0, Nt*dt, Nt+1)      # Mesh points in time

    # Find max(c) using a fake mesh and adapt dx to C and dt
    if isinstance(c, (float,int)):
        c_max = c
    elif callable(c):
        c_max = max([c(x_) for x_ in np.linspace(0, L, 101)])
    dx = dt*c_max/(stability_safety_factor*C)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)          # Mesh points in space

    # Treat c(x) as array
    if isinstance(c, (float,int)):
        c = np.zeros(x.shape) + c
    elif callable(c):
        # Call c(x) and fill array c
        c_ = np.zeros(x.shape)
        for i in range(Nx+1):
            c_[i] = c(x[i])
        c = c_

    q = c**2
    C2 = (dt/dx)**2; dt2 = dt*dt    # Help variables in the scheme

    # Wrap user-given f, I, V, U_0, U_L if None or 0
    if f is None or f == 0:
        f = (lambda x, t: 0) if version == 'scalar' else \
            lambda x, t: np.zeros(x.shape)
    if I is None or I == 0:
        I = (lambda x: 0) if version == 'scalar' else \
            lambda x: np.zeros(x.shape)
    if V is None or V == 0:
        V = (lambda x: 0) if version == 'scalar' else \
            lambda x: np.zeros(x.shape)
    if U_0 is not None:
        if isinstance(U_0, (float,int)) and U_0 == 0:
            U_0 = lambda t: 0
    if U_L is not None:
        if isinstance(U_L, (float,int)) and U_L == 0:
            U_L = lambda t: 0


    u   = np.zeros(Nx+1)   # Solution array at new time level
    u_1 = np.zeros(Nx+1)   # Solution at 1 time level back
    u_2 = np.zeros(Nx+1)   # Solution at 2 time levels back

    Ix = range(0, Nx+1)
    It = range(0, Nt+1)

    # Load initial condition into u_1
    for i in range(0,Nx+1):
        u_1[i] = I(x[i])

    """
    For the first time step we use special formulas
    """
    for i in Ix[1:-1]:
        u[i] = u_1[i] + dt*V(x[i]) + \
        0.5*C2*(0.5*(q[i] + q[i+1])*(u_1[i+1] - u_1[i]) - \
                0.5*(q[i] + q[i-1])*(u_1[i] - u_1[i-1])) + \
        0.5*dt2*f(x[i], t[0])
        
    """
    Setting the boundary values for the first time step. If solve = s1, 
    the scheme from a) and b) is used. If solve = s2, the scheme from c) and
    d) is used.
    """

    if solve=='s1':
        i = Ix[0]
        if U_0 is None:
            # Set boundary values (x=0: i-1 -> i+1 since u[i-1]=u[i+1]
            # when du/dn = 0, on x=L: i+1 -> i-1 since u[i+1]=u[i-1])
            ip1 = i+1
            im1 = ip1  # i-1 -> i+1
            u[i] = u_1[i] + dt*V(x[i]) + \
                   0.5*C2*(0.5*(q[i] + q[ip1])*(u_1[ip1] - u_1[i])  - \
                           0.5*(q[i] + q[im1])*(u_1[i] - u_1[im1])) + \
            0.5*dt2*f(x[i], t[0])
        else:
            u[i] = U_0(dt)
    
        i = Ix[-1]
        if U_L is None:
            im1 = i-1
            ip1 = im1  # i+1 -> i-1
            u[i] = u_1[i] + dt*V(x[i]) + \
                   0.5*C2*(0.5*(q[i] + q[ip1])*(u_1[ip1] - u_1[i])  - \
                           0.5*(q[i] + q[im1])*(u_1[i] - u_1[im1])) + \
            0.5*dt2*f(x[i], t[0])
        else:
            u[i] = U_L(dt)
            
    if solve=='s2':
        #Boundary condition u[i+1] - u[i] = 0 for i = 0
        i = Ix[0]
        if U_0 is None:      
            u[i] = u[i+1]
        else:
            u[i] = U_0(dt)
    
        #Boundary condition u[i] - u[i-1] = 0 for i = Nx
        i = Ix[-1]
        if U_L is None:
            u[i] = u[i-1]
        else:
            u[i] = U_L(dt)

    # Computing the first error, later used to find the convergence rate
    error = abs(u_exact(x,0)-u)
    e2_sum = sum(error**2) 

    # Update data structures for next step
    u_2, u_1, u = u_1, u, u_2

    for n in It[1:-1]:
        # Update all inner points
        if version == 'scalar':
            for i in Ix[1:-1]:
                u[i] = - u_2[i] + 2*u_1[i] + \
                    C2*(0.5*(q[i] + q[i+1])*(u_1[i+1] - u_1[i])  - \
                        0.5*(q[i] + q[i-1])*(u_1[i] - u_1[i-1])) + \
                dt2*f(x[i], t[n])
        else:
            raise ValueError('version=%s' % version)

        # Insert boundary conditions
        """
        Setting the boundary values for the first time step. If solve = s1, 
        the scheme from a) and b) is used. If solve = s2, the scheme from c) 
        and d) is used.
        """
        if solve=='s1':
            i = Ix[0]
            if U_0 is None:
                # Set boundary values
                # x=0: i-1 -> i+1 since u[i-1]=u[i+1] when du/dn=0
                # x=L: i+1 -> i-1 since u[i+1]=u[i-1] when du/dn=0
                ip1 = i+1
                im1 = ip1
                u[i] = - u_2[i] + 2*u_1[i] + \
                       C2*(0.5*(q[i] + q[ip1])*(u_1[ip1] - u_1[i])  - \
                           0.5*(q[i] + q[im1])*(u_1[i] - u_1[im1])) + \
                dt2*f(x[i], t[n])
            else:
                u[i] = U_0(t[n+1])
    
            i = Ix[-1]
            if U_L is None:
                im1 = i-1
                ip1 = im1
                u[i] = - u_2[i] + 2*u_1[i] + \
                       C2*(0.5*(q[i] + q[ip1])*(u_1[ip1] - u_1[i])  - \
                           0.5*(q[i] + q[im1])*(u_1[i] - u_1[im1])) + \
                dt2*f(x[i], t[n])
            else:
                u[i] = U_L(t[n+1])
                
        if solve=='s2':
            # Insert boundary conditions
            i = Ix[0]
            if U_0 is None:
                u[i] = u[i+1]
            else:
                u[i] = U_0(t[n+1])
    
            i = Ix[-1]
            if U_L is None:
                u[i] = u[i-1]
            else:
                u[i] = U_L(t[n+1])
                
        # Computing the error, later used to find the convergence rate
        error = abs(u_exact(x,t[n])-u)
        e2_sum += sum(error**2) 

        # Update data structures for next step

        u_2, u_1, u = u_1, u, u_2

    # Important to correct the mathematically wrong u=u_2 above
    # before returning u
    u = u_1
    
    # Total error:
    E = np.sqrt(dt*dx*e2_sum)
    
    return  u, x, E



def test_q(q, solve):
    L = 14.  # Length
    Nx = 100 # Number of grid points 
    
    # The exact u we have fitted f to:
    u_exact = lambda x, t: np.cos(np.pi*x/L)*np.cos(t)
    
    # Initial conditions:
    I = lambda x: u_exact(x, 0)
    """
    du(x, t)/dt = -sin(t)*cos(pi*x/L)
    V = du(x, 0)/dt = -sin(0)*cos(pi*x/L) = 0
    """
    V = lambda x: 0
    
    """
    Different q values for a) and b) give different f and c.
    F is found using sympy in a separate program.
    """
    if q == 'q1': # q1 = 1 + (x- L/2)**4
        f = lambda x, t: (-16*L**2*np.cos(np.pi*x/L) - \
                        8*np.pi*L*(L - 2*x)**3*np.sin(np.pi*x/L) + \
                        np.pi**2*((L - 2*x)**4 + \
                        16)*np.cos(np.pi*x/L))*np.cos(t)/(16*L**2)
        c = lambda x: np.sqrt(1 + (x- L/2)**4)
        
    if q == 'q2': # q2 = 1 + cos(pi*x/L)
        f = lambda x, t: (-L**2*np.cos(np.pi*x/L) - \
                      2*np.pi**2*np.sin(np.pi*x/L)**2 + \
                      np.pi**2*np.cos(np.pi*x/L) + np.pi**2)*np.cos(t)/L**2
        c = lambda x: np.sqrt(1 + np.cos(np.pi*x/L))
                
    """
    Boundary condtions. Use 'None' for Neumann conditions, give a value for
    Dirichlet condtions. In this exercise we will only need 'None'.
    """            
                
    U_0 = None #lambda t: u_exact(0, t)
    U_L = None #lambda t: u_exact(L, t)
    
    
    x_e = np.linspace(0, L, Nx) #x needed to find an estimated c_max
    c_max = max(c(x_e))         #c_max needed to make dt fit
    C = 0.75                    #Courant number, C <= 1
    dt = C*((L)/Nx)/c_max       #Increase Nx (above) to make dt smaller
    T =  10                     #Number of periods for simulation
   

    u, x, E = solver(I, V, f, c, U_0, U_L, L, dt, C, T, u_exact, 
                     solve, version='scalar', stability_safety_factor=1)

    """
    Plotting the exact solution and the numerical solution to compare, using 
    the values for T and Nx given above.
    """
    plt.plot(x, u)
    plt.hold('on')
    plt.plot(x_e, u_exact(x_e, T))
    plt.title('Nx = %.0f, T = %.1f, Solver = %s, q = %s' % (Nx, T, solve, q))
    plt.legend(['Numerical solution', 'Exact solution'], loc=0)
    plt.xlabel('x')
    plt.ylabel('u')
    plt.show()
    
    
    # Some empty lists to store values in
    dt_list = []
    E_list = []
    """
    List of Nx values too be looped through to see what effect increasing Nx 
    (= decreasing dt) has on the error.
    """
    Nx_list = [4, 8, 16, 32, 64, 128, 256, 512]
    
    for Nx in Nx_list:
        #Calculating everything that depends on Nx again, the rest is as before
        x_e = np.linspace(0, L, Nx)
        c_max = max(c(x_e))                  
        dt = C*((L)/Nx)/c_max
        u, x, E = solver(I, V, f, c, U_0, U_L, L, dt, C, T, u_exact, 
                          solve, version='scalar', stability_safety_factor=1)
        dt_list.append(dt)
        E_list.append(E)
        
    #Convergence rate. Not sure if this is correct, please let me know!    
    r = np.zeros(len(Nx_list)-1)        
    for i in range(len(E_list)-1):
        r[i] = (np.log(E_list[i+1]/E_list[i]))/(np.log(dt_list[i+1]/dt_list[i]))
    
    """
    Plotting the error as a function of Nx and the convergence rate (r).
    """    
    plt.plot(Nx_list, E_list)
    plt.title('Solver = %s, q = %s' % (solve, q))
    plt.xlabel('Nx')
    plt.ylabel('Error')
    plt.show()
    plt.plot(Nx_list[:-1], r)
    plt.title('Solver = %s, q = %s' % (solve, q))
    plt.xlabel('Nx')
    plt.ylabel('r')
    plt.show()
    
    

if __name__ == '__main__':
    #Try one (or more) of the following:
    
    test_q('q1', 's1')      #Exercise 13 a)
    test_q('q2', 's1')     #Exercise 13 b)
    test_q('q1', 's2')     #Exercise 13 c) and d)
    test_q('q2', 's2')     #Exercise 13 c) and d)
import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
V, t, I, w, dt = sym.symbols('V t I w dt')  # global symbols
b = sym.symbols('b') #global symbols for quadtratic
c = sym.symbols('c') #global symbols for poly
f = None  # global variable for the source term in the ODE

def ode_source_term(u):
    """Return the terms in the ODE that the source term
    must balance, here u'' + w**2*u.
    u is symbolic Python function of t."""
    return sym.diff(u(t), t, t) + w**2*u(t) 

def residual_discrete_eq(u):
    """Return the residual of the discrete eq. with u inserted."""
 
    diff = DtDt(u, dt) + w**2*u(t)
    R = ode_source_term(u) - diff
    return sym.simplify(R)

def residual_discrete_eq_step1(u):
    """Return the residual of the discrete eq. at the first
    step with u inserted."""
    f_0 = w**2*I + DtDt(u, dt).subs(t,0)
    first = 0.5*(dt**2*f_0 - dt**2*w**2*I + 2*I + 2*dt*V)
    R = (ode_source_term(u).subs(t, dt) - sym.diff(u(t), t, t).subs(t, dt))*w**-2 - first
    return sym.simplify(R)

def DtDt(u, dt):
    """Return 2nd-order finite difference for u_tt.
    u is a symbolic Python function of t.
    """
    return (u(t+dt)-2*u(t) + u(t-dt) )/dt**2 

def solver(I, w, dt, T, V, eq):
    N = int(round(T/dt)) +1 

    u = np.zeros(N+1)
    time = np.zeros(N+1)
    #f_0 = 0
    f_0 = w**2*eq(0.0) + DtDt(eq, dt).subs(t, 0)
    u[0] = I
    u[1] = 0.5*(dt**2*f_0 - dt**2*w**2*I + 2*I + 2*dt*V)

    for i in range(1, len(u)-1):
        #f_last = 0 
        f_last = w**2*u[i] + DtDt(eq, dt).subs(t, i*dt)
        u[i+1] = dt**2*f_last + u[i]*(2 - dt**2*w**2) - u[i-1]
        time[i+1] = i*dt
 
    plt.figure(1)
    plt.title("Numerical against Exact solution")
    plt.plot(time, u, 'r*')
    #Testing for correct numerical scheme for f(t) = 0
    #plt.plot(time, I*np.cos(w*time) + V/w*np.sin(w*time))
    plt.plot(time, eq(time))
    plt.show()

    return time, u

def exact(eq, t):
    return eq(t)


def typ(dim, I, V):
    if dim == "linear":
        return  lambda t: V*t + I
    if dim == "quadratic":
        return  lambda t: t**2 + V*t + I
    if dim ==  "poly":
        return lambda t: t**3 + t**2 + V*t + I

def convergence_rates(eq, I, V, w, m, solver, num_periods=8):
    """
    Return m-1 empirical estimates of the convergence rate
    based on m simulations, where the time step is halved
    for each simulation.
    """
    dt = 2*np.pi/w/30 # 30 time step per period 2*pi/w
    T = 2*np.pi/w*num_periods
    dt_values = []
    E_values = []
    for i in range(m):
        time, u = solver(I, w, dt, T, V, eq)
        #u_e = exact(I, w, t)
        u_e = eq(time)
        E = np.sqrt(dt*sum((u_e-u)**2))
        dt_values.append(dt)
        E_values.append(E)
        dt = dt/2

    r = [np.log(E_values[i-1]/E_values[i])/
        np.log(dt_values[i-1]/dt_values[i])
        for i in range(1, m, 1)]
    return r
    

def main(u, ty):
    """
    Given some chosen solution u (as a function of t, implemented
    as a Python function), use the method of manufactured solutions
    to compute the source term f, and check if u also solves
    the discrete equations.
    """
    print '=== Testing exact solution: %s ===' % u
    print "Initial conditions u(0)=%s, u'(0)=%s:" % \
          (u(t).subs(t, 0), sym.diff(u(t), t).subs(t, 0))

    # Method of manufactured solution requires fitting f
    global f  # source term in the ODE
    f = sym.simplify(ode_source_term(u))
    # Residual in discrete equations (should be 0)
    print 'residual step1:', residual_discrete_eq_step1(u)
    print 'residual:', residual_discrete_eq(u)

    I = 1.0; T = 5.0; dt = 0.01; w = 1.0; V = 1.0
    eq = typ(ty, I, V)
    solver(I, w, dt, T, V, eq)

def linear():
    ty = "linear"
    print ty
    main(lambda t: V*t + I, ty)
   
def quadratic():
    ty = "quadratic"
    print ty
    main(lambda t: b*t**2 + V*t + I, ty)

def poly():
    ty = "poly"
    print ty
    main(lambda t: c*t**3 + b*t**2 + V*t + I, ty)

def calc():
    return null

#This ensures that the program will NOT exectute when if it is 
#imported to another program.
#Say if this was executed with linear() only, it would execute every
#Time I imported the file
if __name__ == '__main__':
    linear()
    quadratic()
    poly()
    #r = convergence_rates(5, solver)


    #a u[1] = dt**2*f(t[0]) - dt**2*u[0] + 2*u[0]- u[-1]
    # (u[1]-u[-1])*dt**-2 = V; u[-1] = u[1] - v*dt**-2
    #u[1] = 0.5*(dt**2*f(t[0]) - dt**2*u[0] + 2*u[0] + v*dt**-2)

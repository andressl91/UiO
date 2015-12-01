from dolfin import *
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

set_log_active(False) # Removes callinf FFC...
class solver(object):

    def __init__(self, **kwargs):
    	for key in kwargs:
    		setattr(self, key, kwargs[key])
        self.u0 = None
        self.u = None
        self.f = None
        self.alpha = None
        self.V = None #Last CHANGE
    def define_mesh(self, mesh):
        self.mesh = mesh

    def define_init(self, exp):
        self.u0 = exp #expression

    def define_dt(self, dt):
        self.dt = dt

    def define_f(self, f):
        self.f = f

    def define_alpha(self, alpha):
        self.alpha = alpha

    def alpha_func(self, u=0):
        #return 1+ u*u
        return self.alpha(u)

    def f_func(self,t=0):
        return self.f

    def get_solution(self):
        return self.u

    def Piccard(self, iterations = 1, getstep = None, P_element = 1):
    	#Define Functionspace
        mesh = self.mesh
    	V = FunctionSpace(mesh, 'CG', P_element); self.V = V
    	u = TrialFunction(V)
    	v = TestFunction(V)

		#Initial Condition
        self.u0.t = 0 #Fuckup Fenics ?
        self.f.t = 0 #Fuckup Fenics ?
        u0 = self.u0  	
    	u_1 = interpolate(u0, V)
    	#Define variatonal problem
        dt = self.dt; q = self.q; T = self.T
    	a = u*v*dx + dt/q*self.alpha_func(u_1)*inner(nabla_grad(u), nabla_grad(v))*dx
    	L = (u_1 + dt/q*self.f_func() )*v*dx

    	#u = Function(V)
        Nt = int(round(T/dt)+1)
        E_list = np.zeros(Nt+1)
    	error = 1
        u_ = Function(V)

        for i in range(1, Nt+1):
            k = 0
            self.f.t = (i-1)*dt
            A = assemble(a)

            while error > 1E-5 or k < iterations:   
                b = assemble(L)
                #solve(A, u.vector(), b)
                solve(A, u_.vector(), b)
                k += 1
                u_1.assign(u_) 
                e = u_1.vector().array() - u_.vector().array()
                error = np.sqrt(np.sum(e*e)/u_.vector().array().size)

            

            if getstep == i*dt:
                u0.t = i*dt
                ue = interpolate(u0, V)
                e = ue.vector().array() - \
                    u_1.vector().array()
                E = np.sqrt(np.sum(e*e)/u_1.vector().array().size)
                self.u = u_
                return E

        self.u = u_1 #This is only the last
    	return self.u
		
class Verify(object):

        def __init__(self, solverclass):
            self.sol = solverclass

        def constant(self):
            dim = 1
            Dim1 = IntervalMesh(10,0,6)
            Dim2 = UnitSquareMesh(6, 4)
            Dim3 = UnitCubeMesh(6, 10, 5)
            Dim = [Dim1, Dim2, Dim3]
            #sol.define_alpha(Constant(1))
            sol.define_alpha(lambda u: Constant(1))
            sol.define_f(Constant(0))
            sol.define_init(Expression("1", t=0))
            print "Exercise D"
            print "Checking solver for 1-3 dimentions for P1 and P2 elements, be patient"
            for i in range(1,3): #Testing P1 and P2 elements
                while dim <= 3: #Increasing dimention of solutionspace
                    sol.define_mesh(Dim[dim-1])
                    u = self.sol.Piccard(P_element = i)
                    u_vector = u.vector().array()
                    test = all(x >= u_vector[0]-1E-5 for x in u_vector)
                    assert(test == True) #Test if all values are equal
                    dim += 1
                dim = 1
            print "All values of the vector are constant"
            print

        #Checking convergencerate for 
        #u_e = exp(-pi**2*t)cos(pi*x)
        def analytic(self):
            ue = Expression("exp(-pi*pi*t)*cos(pi*x[0])", t=0)
            Dim2 = UnitSquareMesh(10, 10)
            sol.define_mesh(Dim2)
            sol.define_init(ue)
            sol.define_alpha(lambda u: Constant(1))
            sol.define_f(Constant(0))
 
            Nx = int(1./np.sqrt(sol.dt))
            getstep = sol.dt*10
            print "Exercise E"
            print "Getting solution at time %.3f seconds" % getstep
            for i in range(5): 
                sol.define_mesh(UnitSquareMesh(Nx, Nx))
                error = sol.Piccard(getstep = getstep)
                print "E/h = %.4f for t=%.5f seconds" % (error/sol.dt, sol.dt)
                sol.define_dt(sol.dt/2.)
                Nx = int(1./np.sqrt(sol.dt))
            print "E/h is approximately constant"
            print 

        def manufactured(self, plot):
            Dim1 = IntervalMesh(1000,0,1)
            sol.define_mesh(Dim1)
            ue = Expression("t*x[0]*x[0]*(0.5-x[0]/3.0)", t=0)
            sol.define_init(ue)
            sol.define_f(Expression("-q*x[0]*x[0]*x[0]/3\
                            + q*x[0]*x[0]/2 + \
            8*t*t*t *x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]/9 -\
            28*t*t*t *x[0]*x[0]*x[0]*x[0]*x[0]*x[0]/9 + \
            7*t*t*t *x[0]*x[0]*x[0]*x[0]*x[0]/2 - \
            5*t*t*t *x[0]*x[0]*x[0]*x[0]/4 + \
            2*t*x[0]-t", t=0,q=sol.q))
            sol.define_alpha(lambda u: 1+u*u)
            getstep = [i*sol.dt for i in [1,10,20,30]]
            count = 0 #for filemaking
            print "Exercise F"
            print "Check solution plots in the provided PDF"
            print
            for i in getstep:
                sol.Piccard(getstep = i)
                ue.t = i
                u_e = interpolate(ue, sol.V)
                if plot == True:
                    fig = plt.figure(1)
                    ax = fig.add_subplot(111)
                    ax.set_title("Time %.4f seconds" % i)
                    ax.set_xlim([0, 1000]); ax.set_ylim([0, .07])
                    ax.plot(sol.u.vector().array(), 'r')
                    ax.plot(u_e.vector().array())
                    plt.savefig('frame%.1d.png' % count)
                    count += 1
                    plt.close(fig)
                
        def convergence(self):
            ue = Expression("t*x[0]*x[0]*(0.5-x[0]/3.0)", t=0)
            sol.define_init(ue)
            sol.define_f(Expression("q*pow(x[0],2)*(-2*x[0]+3)/6 -\
                ( -12*t*x[0] + 3*t*(-2*x[0]+3) ) *\
                ( pow(x[0],4)*pow(-dt+t,2)*pow(-2*x[0] + 3,2) + 36 )/324 -\
                ( -6*t*pow(x[0],2) + 6*t*x[0]*(-2*x[0]+3) ) *\
                ( 36*pow(x[0],4)*(-dt+t)*(-dt+t)*(2*x[0]-3) +\
                36*x[0]*x[0]*x[0]*(-dt+t)*(-dt+t) *\
                (-2*x[0]+3)*(-2*x[0]+3) )/5832", t=0.0,dt=sol.dt,q=sol.q))
            sol.define_alpha(lambda u: 1+u*u)
            cut = 6
            E = np.zeros(cut)
            h = np.zeros(cut)
            Nx = int(1./np.sqrt(sol.dt))
            get = sol.dt*10
            for i in range(cut):
                sol.u0.dt = sol.dt 
                sol.define_mesh(IntervalMesh(Nx, 0, 1))
                h[i] = sol.dt
                E[i] = sol.Piccard(getstep=get)
                sol.define_dt(sol.dt/2.)
                Nx = int(1./np.sqrt(sol.dt))
            print "Exercise H"
            r = np.zeros(len(E)-1)
            r[:] = np.log(E[:-1]/E[1:])/np.log(h[:-1]/h[1:])
            print r
            print "r stays around 1"
            print

    
        def gaussian(self, sigma = 0.09, beta = 100.0, Plot = False):
            mesh = UnitSquareMesh(20, 20)
            sol.define_mesh(mesh)
            sol.define_init(Expression("10*exp(-1/(2*sigma*sigma)*\
                        (x[0]*x[0]+x[1]*x[1]))", sigma = 0.09))
            sol.alpha = lambda u: 1+beta*u*u
            sol.f = Expression("t",t=0)
            Nt = int(round(sol.T/sol.dt)+1)
            t = [i*sol.dt for i in range(50)]
            import matplotlib.pyplot as plt
            u_file = File("velocity/Gaussian.pvd")
            print "Exercise I"
            print "Solution presented in the GIF file provided"
            print

            if Plot == True:
                count = 0
                for i in t:
                    sol.Piccard(getstep = i)
                    #u_file << sol.u
                    wiz = plot(sol.get_solution(),
                    wireframe = True,           # use wireframe rendering
                    interactive = False,           # do not hold plot on screen
                    scalarbar = True,           # hide the color mapping bar
                    hardcopy_prefix = "myplot",  # default plotfile name
                    scale = 2.0    ,               # scale the warping/glyphs
                    title = "Diffusion plot ",
                    range_max = 0.4)
                    wiz.write_png("result/myplot%04d" % count)
                    count += 1
                
if __name__ == "__main__":
    sol = solver(dt = 0.01, T = 1.0, q = 1.0)
 
    v = Verify(sol)
    v.constant() #Exercise D
    v.analytic() #Exercise E
    v.manufactured(plot = False) #Exercise F
    v.convergence() #Exercise H
    v.gaussian(Plot = False) #Exercise I

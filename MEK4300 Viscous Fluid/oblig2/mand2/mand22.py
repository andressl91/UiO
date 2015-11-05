from dolfin import *

mesh = Mesh("tube2.xml")
mesh = refine(mesh) #To get a smoother solution, original was too rough 

V = VectorFunctionSpace(mesh, 'CG', 2)
Q = FunctionSpace(mesh, 'CG', 1)
VQ = V*Q
u, p= TrialFunctions(VQ)
v, q = TestFunctions(VQ)
#These are used for the initial guess in Falkner-Skan
ut, pt= TrialFunctions(VQ)
vt, qt = TestFunctions(VQ)


#Condition
rho = 1.0
nu = 10**-3
D =0.1
R = D/2
Re = 20
Um = 0.3 
H = 0.41

#Inflow
u_in = Expression(("4*Um*x[1]*(H-x[1])/(H*H)", "0.0"), H=H, Um = Um)

def top(x, on_boundary): return x[1] > (H-DOLFIN_EPS)
def bottom(x, on_boundary): return x[1] < DOLFIN_EPS
def circle(x, on_boundary): return sqrt((x[0] - 0.20)**2 + (x[1] - 0.20)**2) < (R+DOLFIN_EPS)

bc0 = DirichletBC(VQ.sub(0), Constant((0, 0)), top)
bc1 = DirichletBC(VQ.sub(0), Constant((0, 0)), bottom)
bc2 = DirichletBC(VQ.sub(0), Constant((0, 0)), circle)
inflow = DirichletBC(VQ.sub(0), u_in, "x[0] < DOLFIN_EPS")
bcs = [bc0, bc1, bc2, inflow]


F1 = nu*inner(grad(ut), grad(vt))*dx - q*div(ut)*dx - div(vt)*pt*dx + Constant(0)*qt*dx 
uguess = Function(VQ)

solve(lhs(F1) == rhs(F1), uguess, bcs=bcs)

k = 0
error = 1
uh_1 = Function(VQ) #The new unknown function

while k < 10 and error > 1e-12:
	F = nu*inner(grad(u), grad(v))*dx + inner(grad(u)*uguess.sub(0), v)*dx - div(v)*p*dx-q*div(u)*dx+Constant(0)*q*dx
	solve(lhs(F) == rhs(F), uh_1 , bcs=bcs) 
	error = errornorm(uguess, uh_1)
	uguess.assign(uh_1) #updates the unknown for next it.
	k += 1

us, ps = uguess.split()

#Making parameters to calculate forces and coefficients
n = -FacetNormal(mesh) #Fetching normal 
x1 = as_vector((1,0)) #unit vectors
x2 = as_vector((0,1))
nx = dot(n, x1)#Fetching the x and y components of the normal
ny = dot(n, x2)
nn = as_vector((ny, -nx))

cir = AutoSubDomain(circle) #Making a domain from our circle function
ff = FacetFunction("size_t", mesh) #Function with a value for each facet
ff.set_all(0)
cir.mark(ff, 1)
ds = ds[ff]
#ds is in Fenics, representing in use for integration over facet at boundary
un = dot(nn, us)

#Now we calculate the unknowns
rud = (2*Um)/3.0

#Drag
Fd = assemble((rho*nu*dot(grad(un), n)*ny - ps*nx)*ds(1), exterior_facet_domains=ff)
#Lift
Fl = assemble(-(rho*nu*dot(grad(un), n)*nx + ps*ny)*ds(1), exterior_facet_domains=ff)
#Drag coefficient
Cd = (2*Fd)/(rho*rud*rud*D)
#Lift coefficient
Cl = (2*Fl)/(rho*rud*rud*D)
#Pressure difference
pressure = ps.compute_vertex_values()
pdiff = pressure[5] - pressure[7] #
mc = mesh.coordinates()
print mc[8][0]

print('Calulated values: Cd = %f, Cl = %f pdiff = %f' % ( Cd, Cl, pdiff))




import numpy as np 
from math import sin
import sympy as sym
from sympy.utilities.lambdify import lambdify, implemented_function

x,y,t,A,B,w,kx,ky,q,b= sym.symbols('x,y,t,A,B,w,kx,ky,q,b')  # global symbols
#q = 1+x
#w = sym.sqrt(q*(kx**2+ky**2))
#ue = (A*sym.cos(w*t)+B*sym.sin(w*t))*sym.exp(-sym.sqrt(q)*t)*sym.cos(kx*x)*sym.cos(ky*y)

q = sym.Function('q')(x,y)
q_y = sym.diff(q,y)
print q_y
q = q.subs(q,x+y)
q_y = sym.diff(q,y)
print q_y

x,y,t,A,B,w,kx,ky,q,b= sym.symbols('x,y,t,A,B,w,kx,ky,q,b')  # global symbols
q = sym.Function('q')(x,y)
w = sym.sqrt(q*(kx**2+ky**2))
ue = (A*sym.cos(w*t)+B*sym.sin(w*t))*sym.exp(-sym.sqrt(q)*t)*sym.cos(kx*x)*sym.cos(ky*y)
u = ue
u_x = sym.diff(u,x)
u_y = sym.diff(u,y)
q_x = sym.diff(q,x)
q_y = sym.diff(q,y)
rhs = q_x*u_x + q*sym.diff(u,x,x) + q_y*u_y + q*sym.diff(u,y,y)
lhs = sym.diff(u,t,t) + b*sym.diff(u,t)
q = q.subs(q,y)
q_x = sym.diff(q,x)
q_y = sym.diff(q,y)
print q_x

rhs = q_x*u_x + q*sym.diff(u,x,x) + q_y*u_y + q*sym.diff(u,y,y)
print rhs
#f = np.vectorize(sym.lambdify((x,y,t), f))



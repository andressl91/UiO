import numpy as np

N = 100
x = np.linspace(0, 1, N)
y = x
dx = x[1]- x[0]
def u(x):
    return np.sin(x)

def v(x):
    return np.cos(x)

grad_uv = [(u(x[i+1])*v(x[i+1]) - u(x[i])*v(x[i]))/dx for i in range(20, 80)]
gradu_v = [(u(x[i+1]) - u(x[i]))/dx * v(x[i+1]) for i in range(20, 80)]
u_gradv = [u(x[i+1]) * (v(x[i+1]) - v(x[i]))/dx for i in range(20, 80)]

diff = [grad_uv[i] - gradu_v[i] - u_gradv[i] for i in range(len(grad_uv))]
L2_diff = np.linalg.norm(diff)
print(L2_diff)
#print(grad_uv)
#print(gradu_v + u_gradv)

#u_vec = [u(x), 0 for i in range(len(x))]
#print(u_vec)

import numpy as np
import matplotlib.pyplot as plt

def func(x):
    if (1./3. <= x <= 2./3.):
        return 1.0
    else:
        return 0.0

def integrate(N, k):
    #x = np.linspace(0,1,N+1)
    #h = x[1]-x[0]
    #area = h/2.0* (f(x[:-1])*np.sin(k*np.pi*x[:-1]) + f(x[1:])*np.sin(k*np.pi*x[1:]))
    #return 2*np.sum(area) #2* from inner product
    area = 1./(np.pi*k)*(- np.cos(k*np.pi*2./3) + np.cos(k*np.pi*1./3))
    return area*2


def sinefunc(N):
    k = [i for i in range(1,90)]
    x = np.linspace(0, 1, N+1)

    inte = np.vectorize(integrate)
    c_k = inte(N, k[:])
    y = 0

    for i in range(len(c_k)):
        y = y + c_k[i]*np.sin(k[i]*np.pi*x)

    return y

def plotsolution(N):
    x = np.linspace(0, 1, N+1)
    f = np.vectorize(func)
    y = sinefunc(N)

    plt.figure(1)
    plt.axis([0,1,-1,1.5])
    plt.plot(x, f(x))
    plt.plot(x, y)
    plt.show()

plotsolution(100)

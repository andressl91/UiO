import numpy as np
import matplotlib.pyplot as plt

def func(x):
    if (1./3. <= x <= 2./3.):
        return 1.0
    else:
        return 0.0

def integrate(N, k):
    area = 1./(np.pi*k)*(- np.cos(k*np.pi*2./3) + np.cos(k*np.pi*1./3))
    return area*2


def fourier(N, K):
    k = [i for i in range(1,K+1)]
    x = np.linspace(0, 1, N+1)

    inte = np.vectorize(integrate)
    c_k = inte(N, k[:])
    y = 0

    for i in range(len(c_k)):
        y = y + c_k[i]*np.sin(k[i]*np.pi*x)

    return y

def norms(N, K):
    x = np.linspace(0, 1, N+1)
    h = x[1]-x[0]
    k = [i for i in range(1,K+1)]
    u = fourier(N, K)

    f = np.vectorize(func)
    e = f(x)-u

    print "For K = %d, N = %d \n" % (K, N)

    #Make L2 norm
    L2 = np.sqrt( np.sum(e*e) / len(e) ) #better answer outside

    #Make max_norm
    arg = np.argmax(e)
    L_max = np.amax(e)

    #Make H1 norm
    H1 = np.sqrt ( np.sum(e[:-1]*e[:-1]) + np.sum( ( (e[1:]-e[:-1])/h )**2 ) ) #/ len(e)

    print "L2 = %.5f  L_max = %.3f " % (L2, L_max)
    print "H1 = %.5f\n" % H1

def loop():
    for k in [10,50,100, 200]:
        print "Number of fourier sine functions are %d" % k
        for j in [100]:
            norms(j, k)

loop()

def plotsolution(N, K):
    x = np.linspace(0, 1, N+1)
    f = np.vectorize(func)
    y = fourier(N, K)

    plt.figure(1)
    plt.axis([0,1,-1,1.5])
    plt.plot(x, f(x))
    plt.plot(x, y)
    plt.show()

#norms(100, 10)
#plotsolution(100, 100)

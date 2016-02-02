import numpy as np
import matplotlib.pyplot as plt


#func = lambda x: climb*x if 0.2 < x < 0.3 else 0
def func(x):
    if 0.2 < x <= 0.3:
        return (1-0)/(0.3-0.2)*(x-0.2)

    if 0.4 > x >= 0.3:
        return (0-1.)/(0.4-0.3)*(x-0.3) + 1

    else:
        return 0*x

def hat(x, h):
    if (0.5 - h) <= x <= 0.5:
        return 1./h * x  - 1./h*(0.5 - h)
    if  0.5 <= x < (0.5 + h):
        return 1./h * x + 1./h*(0.5 - h)
    else:
        return x*0

def norms(x):
    N = len(x)
    f = np.vectorize(func)

    L2 = np.sqrt( 1./N*np.sum(f(x)*f(x)) )

    der = f(x[:-1])-f(x[1:])/(x[:-1]-x[1:])
    H1 = np.sqrt( 1./(N-1) * ( np.sum(f(x)*f(x)) + np.sum(der*der) ) )
    print; print '-----------------------------------'
    print 'For %d points ' % N
    print 'L2 = %f, H1 =%f ' % (L2, H1)
    print '-----------------------------------'; print


#for i in [50, 100, 200, 400, 1000, 5000, 10000]:
for i in [10]:
    x = np.linspace(0, 1, i)
    h = x[2]-x[1]
    #norms(x)

    y = np.vectorize(hat)
    f = np.vectorize(func)
    plt.figure(1)
    plt.plot(x, y(x, h),label = 'N = %d' % i )
    plt.legend()
plt.show()

import numpy as np
import matplotlib.pyplot as plt

p = 300  # size of the image for sampling, m=p*p
t = np.linspace(0, 1, p)
[V, U] = np.meshgrid(t, t)
Y = np.concatenate((U.flatten()[None, :], V.flatten()[None, :]))

n = 30
X = .5+.5j + np.exp(1j*np.pi/4) * 1 * \
    (.1*(np.random.rand(1, n)-.5)+1j*(np.random.rand(1, n)-.5))
X = np.concatenate((np.real(X), np.imag(X)))
a = np.ones(n)/n

def Gauss(mx, my, s):
    return np.exp((-(U-mx)**2-(V-my)**2)/(2*s**2))


Mx = [.6, .4]  # means
My = [.9, .1]
S = [.07, .09]  # variance
W = [.5, .5]  # weights
b = W[0]*Gauss(Mx[0], My[0], S[0]) + W[1]*Gauss(Mx[1], My[1], S[1])
b = b/np.sum(b.flatten())

Col = np.random.rand(n, 3)
plt.imshow(-b[::-1, :], extent=[0, 1, 0, 1], cmap='gray')
plt.scatter(X[1, :], X[0, :], s=30, c=.8*Col)

f = np.zeros(n)
def distmat(x, y): return np.sum(
    x**2, 0)[:, None] + np.sum(y**2, 0)[None, :] - 2*x.transpose().dot(y)


D = distmat(Y, X) - f[:].transpose()
fC = np.min(D, axis=1)
I = np.reshape(np.argmin(D, axis=1), [p, p])

OT = np.sum(f*a) + np.sum(fC*b.flatten())
print(OT)

plt.imshow(I[::-1, :], extent=[0, 1, 0, 1])
plt.scatter(X[1, :], X[0, :], s=20, c='k')
plt.contour(t, t, I, np.linspace(-.5, n-.5, n), colors='k')
plt.axis('off')

plt.show()

def exercise1():
    tau = .02  # step size
    niter = 200  # iteration for the descent
    q = 6  # number of displays
    ndisp = np.unique(np.round(1 + (niter/4-1)*np.linspace(0, 1, q)**2))
    kdisp = 0
    f = np.zeros(n)
    E = np.zeros(niter)
    for it in range(niter):
        # compute Laguerre cells and c-transform
        D = distmat(Y, X) - f[:].transpose()
        fC = np.min(D, axis=1)
        I = np.reshape(np.argmin(D, axis=1), [p, p])
        E[it] = np.sum(f*a) + np.sum(fC*b.flatten())
        # display
        if (kdisp < len(ndisp)) and (ndisp[kdisp] == it):
            plt.subplot(2, 3, kdisp+1)
            plt.imshow(I[::-1, :], extent=[0, 1, 0, 1])
            plt.scatter(X[1, :], X[0, :], s=20, c='k')
            plt.contour(t, t, I, np.linspace(-.5, n-.5, n), colors='k')
            plt.axis('off')
            kdisp = kdisp+1
        # gradient
        R = (I[:, :, None] == np.arange(0, n)[None, None, :]) * b[:, :, None]
        nablaE = a-np.sum(R, axis=(0, 1)).flatten()
        f = f+tau*nablaE

        print(E[it])

    plt.show();

def exercise2():
    f = np.zeros(n)
    niter = 300
    q = 6
    ndisp = np.unique(np.round(1 + (niter/2-1)*np.linspace(0, 1, q)**2))
    kdisp = 0
    E = np.zeros(niter)
    for it in range(niter):
        # sample
        k = np.int(np.random.rand(1) < W[1])  # select one of the two Gaussian
        y = np.array((S[k] * np.random.randn(1) + Mx[k],
                      S[k] * np.random.randn(1) + My[k]))
        # detect Laguerre cell where y is
        R = np.sum(y**2) + np.sum(X**2, axis=0) - 2*y.transpose().dot(X) - f[:]
        i = np.argmin(R)
        # gradient
        nablaEy = a.copy()
        nablaEy[i] = nablaEy[i] - 1
        # gradient ascent
        l0 = 10  # warmup phase.
        tau = .1/(1 + it/l0)
        f = f + tau*nablaEy
        # compute Laguerre cells and c-transform
        D = distmat(Y, X) - f[:].transpose()
        fC = np.min(D, axis=1)
        I = np.reshape(np.argmin(D, axis=1), [p, p])
        E[it] = np.sum(f*a) + np.sum(fC*b.flatten())
        # display
        if (kdisp < len(ndisp)) and (ndisp[kdisp] == it):
            plt.subplot(2, 3, kdisp+1)
            plt.imshow(I[::-1, :], extent=[0, 1, 0, 1])
            plt.scatter(X[1, :], X[0, :], s=20, c='k')
            plt.contour(t, t, I, np.linspace(-.5, n-.5, n), colors='k')
            plt.axis('off')
            kdisp = kdisp+1
        print(E[it])
    plt.show()

def exercise3():
    niter = 60
    q = 6
    ndisp = np.unique(np.round(1 + (niter/4-1)*np.linspace(0, 1, q)**2))
    kdisp = 0
    E = np.zeros(niter)
    X1 = X.copy()
    for it in range(niter):
        # compute Voronoi cells
        D = D = distmat(Y, X1)
        fC = np.min(D, axis=1)
        I = np.reshape(np.argmin(D, axis=1), [p, p])
        E[it] = np.sum(fC*b.flatten())
        # display
        if (kdisp < len(ndisp)) and (ndisp[kdisp] == it):
            plt.subplot(2, 3, kdisp+1)
            plt.imshow(I[::-1, :], extent=[0, 1, 0, 1])
            plt.scatter(X[1, :], X[0, :], s=20, c='k')
            plt.contour(t, t, I, np.linspace(-.5, n-.5, n), colors='k')
            plt.axis('off')
            kdisp = kdisp+1
        # update barycenter
        A = (I[:, :, None] == np.arange(0, n)[None, None, :]) * b[:, :, None]
        B = (I[:, :, None] == np.arange(0, n)[None, None, :]) * \
            b[:, :, None] * (U[:, :, None] + 1j*V[:, :, None])
        X1 = np.sum(B, axis=(0, 1)) / np.sum(A, axis=(0, 1))
        X1 = np.concatenate((np.real(X1)[None, :], np.imag(X1)[None, :]))
        print(E[it])

    plt.show()


exercise2();

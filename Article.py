import math
from matplotlib.colors import BoundaryNorm 
import numpy as np
from numpy import pi
import matplotlib.pyplot as mlt

h_bar,m,w    = 1, 1, 1

def R(N):
    BoundaryR, N = 10, N
    R = np.linspace(BoundaryR,-BoundaryR,N)
    return R

def V(R):
    return np.diag(0.5 * m * w**2 * R**2 )

def T(R):
    n = len(R)
    h2m = h_bar**2/2*m
    k = np.pi/(R[1]-R[0])
    T = np.zeros((n,n))
    for i in range(n):
        T[i,i] = h2m * ((k**2)/3)*(1+(2/(n**2)))
        for j in range(i+1,n):
            T[i,j] = h2m * (2*k**2*(-1)**(j-i))/(n**2 * (np.sin(np.pi*(j-i)/n))**2)
            T[j,i] = T[i,j]
    return T   


def H(N):
    H = T(R(N)) + V(R(N))
    E, psi = np.linalg.eigh(H)
    return E, psi
    

def plot(n,i):
    psi = H(n)[1]
    for j in range(i+1):
        mlt.plot(R(n),psi[:,j])

plot(1000,5)


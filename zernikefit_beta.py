import numpy as np
import pandas as pd
from math import factorial,cos,sin
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def radial(r,n,m):
    R = 0
    for k in range(int((n-m)/2)):
        R += ((-1)**k)*factorial(n-k)*(r**(n - (2*k)))
    
    return R

def zernike(r,theta,n,m):
    return radial(r,n,m)*cos(m*theta)

def zernike_poly(r,theta,n,m,params):
    z = 0
    for i in range(len(n)):
        z += params[i]*zernike(r,theta,n[i],m[i])
    
    return z

def zernike_vec(r,theta,n,m,params):
    z = np.zeros((len(r),))
    for i in range(len(r)):
        z[i] = zernike_poly(r[i],theta[i],n,m,params)
        
    return z

def coeff_matrix(r,theta,n,m):
    if len(r) != len(theta):
        raise ValueError('Length of r and theta vectors does not match')
    if len(n) != len(m):
        raise ValueError('Length of zernike indices are not the same')
        
    coeff = np.zeros((len(r),len(n)))
    for i in range(len(r)):
        for j in range(len(n)):
            coeff[i,j] = zernike(r[i],theta[i],n[j],m[j])
    
    return coeff


data = pd.read_csv('ANT.txt', sep = '\t')
data = data.values
r = data[:,0]/6
theta = data[:,1]*np.pi/180
z = data[:,2]

n_poly = 8
n = np.zeros((n_poly))
m = n.copy()
k = 0
for i in range(1,n_poly+1):
    n[k] = i
    m[k] = -i
    k += 1
    
coeff = coeff_matrix(r,theta,n,m)
np.linalg.lstsq(coeff,z)
"""
disc = 100
r = np.linspace(0,1,disc)
theta = np.linspace(0,2*np.pi,disc)
r,theta = np.meshgrid(r,theta)
r = np.reshape(r, (1,disc**2))[0]
theta = np.reshape(theta, (1,disc**2))[0]
z = zernike_vec(r,theta,[1],[1],[1])

x = r*np.cos(theta)
y = r*np.sin(theta)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x,y,z,s=0.8)
plt.show()
"""

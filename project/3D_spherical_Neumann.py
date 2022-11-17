import numpy as np
import scipy as sp
from scipy import integrate
import matplotlib.pyplot as plt



parameters = [] #list with the data

doc = open("data.txt", "r") #opens a file and reads the data
for line in doc:
        point = line.split()[-1]
        element = float(point)
        parameters = np.append(parameters,element) 
lambda_t, lambda_f, v, tau, nu  = parameters
doc.close()

eta=v*(nu-1)/lambda_f
mu=lambda_t*v/3

#in order to compute the critical radius we need to find the 
#zeros of this function
def equation(R):
    eq_r=-1+(R*np.sqrt(eta/mu)*(np.tan(R*np.sqrt(eta/mu)))**(-1))\
        +((3/2)*(R/lambda_t))   
    return eq_r

print("the critical radius is ", sp.optimize.fsolve(equation,0.1) , "m")
r1=(sp.optimize.fsolve(equation,0.1)+0.002) 
#radius bigger than the critical one

def kappa(k): #we now substitute r1 in the previous equation and solve for k
    eq_k=-1+k*r1*(1/np.tan(k*r1))+(3/2)*(r1/lambda_t)
    return eq_k
k=sp.optimize.fsolve(kappa,10)

alpha=mu*k**2-eta 
A=r1/np.sin(k*r1)


def n(r,t): #define the neutron density
    n = A*np.exp(-alpha*t)*np.sin(k*r)/r
    return n

#plot
R_arr=np.linspace(-r1, r1, 100)#range on the r axis
T_arr=np.linspace(0, 3E-6, 100)#range on the t axis

R, T = np.meshgrid(R_arr, T_arr)
fig = plt.figure()

from mpl_toolkits.mplot3d import Axes3D
ax = Axes3D(fig)

ax.plot_surface(R, T, n(R, T), cmap ='viridis',\
                edgecolor ='black', color = "white")

ax.set_xlabel('r')
ax.set_ylabel('t')
ax.set_zlabel('n(r,t)')

plt.show()

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

print("the critical radius is ", np.pi*np.sqrt(mu/eta) , "m")
r1=np.pi*np.sqrt(mu/eta)+0.002
#radius bigger than the critical one




integrand = lambda R: (2/r1)*R*(1-((R/r1)**2))*np.sin(p*np.pi*R/r1)
#defines the function that has to be integrated

a_p_value=[] #vector with only the values of the coefficients
a_p_error=[] #vector with only the errors of the coefficients

threshold = 1E-6 
#under this value we neglect the coefficients and we stop computing them

p=1 #labels the coefficients a_p. Inizialized to 1

coefficients = open("coefficients_spherical_D", "w+") 
#file where all the coefficients will be written

valid = True
while valid:
    value, error = integrate.quad(integrand, 0, r1) #compute the integral
    if np.abs(value) < threshold: #verify if it is bigger than the treshold
        valid = False
    else:
        a_p_value.append(value)#add the value of the integral to the list
        a_p_error.append(error)#add the error of the integral to the list
        print("a_{}".format(p), "=", value,"+-", error)
        coefficients.write("a_{} = {} +- {}  \n".format(p, value, error))
    p = p+1  #when p is even the integral is always zero

coefficients.close()
    
def n(r,t): #define the neutron density
    n = 0
    for p in range(len(a_p_value)):
        n += a_p_value[p]/r * np.exp(eta*t-mu*(p**2*np.pi**2*t/r1**2))\
            *np.sin(p*np.pi*r/r1)
    return abs(n)    

#plot
R_arr=np.linspace(-r1, r1, 100)#range on the r axis
T_arr=np.linspace(0, 2E-6, 100)#range on the t axis

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






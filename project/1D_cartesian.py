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

A = 1        #parameter setted for the gaussian initial condition
Lambda = 100 #parameter setted for the gaussian initial condition

print("the critical lenght is ", np.pi*np.sqrt(mu/eta), "m" )
L = np.pi*np.sqrt(mu/eta)+0.002    #lenght bigger than the critical one

def f(X):  # gaussian initial conditions
    function_IC = A * np.exp((-4 * Lambda * (X - 0.5 * L) ** 2 / (L ** 2)))
    return function_IC

integrand = lambda X: (2 / L) * f(X) * np.sin(p * np.pi * X / L) 
#defines the function that has to be integrated

a_p_value = [] #vector with only the values of the coefficients
a_p_error = [] #vector with only the errors of the coefficients 

threshold = 1E-6 
#under this value we neglect the coefficients and we stop computing them

p = 1 #labels the coefficients a_p. Inizialized to 1

coefficients = open("coefficients_1D", "w+") 
#file where all the coefficients will be written

valid = True
while valid:
    value, error = integrate.quad(integrand, 0, L) #compute the integral
    if np.abs(value) < threshold: #verify if it is bigger than the threshold
        valid = False
    else:
        a_p_value.append(value) #add the value of the integral to the list
        a_p_error.append(error) #add the error of the integral to the list 
        print("a_{}".format(p), "=", value,"+-", error)
        coefficients.write("a_{} = {} +- {}  \n".format(p, value, error)) 
        #writes on the file the coefficients
    p = p+2  #when p is even the integral is always zero

coefficients.close()

def n(x, t): #define the neutron density
    n = 0
    for i in range(len(a_p_value)):
        P = 2 * i + 1
        n += a_p_value[i] * np.exp(eta*t - mu*(P**2)*(np.pi**2)*t/(L**2))\
            * np.sin(P * np.pi * x/L)
    return n

#plot
X_arr=np.linspace(0, L, 100) #range on the x axis 
T_arr=np.linspace(0, 0.00002, 100) #range on the t axis

X, T = np.meshgrid(X_arr, T_arr)

fig = plt.figure()

from mpl_toolkits.mplot3d import Axes3D
ax = Axes3D(fig)

ax.plot_surface(X, T, n(X, T), cmap ='viridis',\
                edgecolor ='black', color = "white")

ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('n(x,t)')

plt.show()












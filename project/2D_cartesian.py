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

L = np.pi*np.sqrt(2*mu/eta)+0.002 #lenght bigger than the critical one
print("the critical lenght is ", np.pi*np.sqrt(2*mu/eta), "m" )

def f_2D(X,Y):  #initial conditions
    function_IC = (16*X*Y/(L**2))*(1-(X/L))*(1-(Y/L))
    return function_IC

integrand = lambda X,Y:\
            (4/(L**2))*f_2D(X,Y)*np.sin(p*np.pi*X/L)*np.sin(q*np.pi*Y/L) 
#defines the function that has to be integrated

a_pq_value=[] #vector with only the values of the coefficients
a_pq_error=[] #vector with only the errors of the coefficients

a_q_value=[] 
#vector with only the values of the coefficients with a certain fixed p
a_q_error=[] 
#vector with only the errors of the coefficients with a certain fixed p

threshold = 1E-6
#under this value we neglect the coefficients and we stop computing them

p=1
q=1
#labels of the coefficients a_pq. Inizialized to 1

coefficients = open("coefficients_2D", "w+")
#file where all the coefficients will be written

valid=True
while valid:
    while valid:
        value, error = integrate.dblquad(integrand, 0, L, 0, L) 
        #compute the integral
        if np.abs(value) < threshold: #verify if it is bigger than the threshold
            valid = False
        else:
            a_q_value.append(value)#add the value of the integral to the list
            a_q_error.append(error)#add the error of the integral to the list 
            #these are coefficients with a fixed p
            print("a_{},{}".format(p,q), "=", value,"+-", error)
            coefficients.write("a_{},{} = {} +- {} \
                               \n".format(p, q, value, error)) 
            #writes on the file the coefficients
            q += 2  #when q is even the integral is always zero

    a_pq_value.append(a_q_value) 
    a_pq_error.append(a_q_error) 
    #add the coefficients with a fixed p to this list and their relative errors

    q = 1
    p += 2
    #reset the q lable and increase the p lable. 
    #When p is even the integral is always zero
    
    if len(a_q_value) == 0: break
    else: 
        a_q_value = [] 
        a_q_error = []
        #reset of the two lists
        valid=True

coefficients.close()

a_pq_value = a_pq_value[:-1] 
a_pq_error = a_pq_error[:-1]
#the last element of these two lists is an empty list, so we delete it

t=1E-7 #we fix a particular time 

def n(x, y): #define the neutron density
    n = 0
    for i in range(len(a_pq_value)):
        P = 2 * i + 1
        for j in range(len(a_pq_value[i])):
            Q = 2 * j + 1
            n += a_pq_value[i][j] * np.exp(eta*t -mu*np.pi**2*((P**2)/(L**2)\
                + (Q**2)/(L**2))*t) * np.sin(P*np.pi*x/L)* np.sin(Q*np.pi*y/L)
    return n


#plot
X_arr=np.linspace(0, L, 1000)#range on the x axis
Y_arr=np.linspace(0, L, 1000)#range on the y axis

X, Y = np.meshgrid(X_arr, Y_arr)

fig = plt.figure()

from mpl_toolkits.mplot3d import Axes3D
ax = Axes3D(fig)

ax.plot_surface(X, Y, n(X, Y), cmap ='viridis',\
                edgecolor ='black', color = "white")

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('n(x,y)')

plt.show()


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

L = np.pi*np.sqrt(3*mu/eta)+0.002    #lenght bigger than the critical one
print("the critical lenght is ",np.pi*np.sqrt(3*mu/eta) , "m" )

def f_3D(X,Y,Z):  # initial conditions
    function_IC = (8*X*Y*Z/(L**3))*(1-(X/L))*(1-(Y/L))*(1-(Z/L))
    return function_IC

integrand = lambda X,Y,Z: (8/(L**3))*f_3D(X,Y,Z)*np.sin(p*np.pi*X/L)\
                          *np.sin(q*np.pi*Y/L)*np.sin(r*np.pi*Z/L)
#defines the function that has to be integrated

a_pqr_value=[] #vector with only the values of the coefficients
a_pqr_error=[] #vector with only the errors of the coefficients

a_qr_value=[] 
#vector with only the values of the coefficients with a certain fixed p
a_qr_error=[]
#vector with only the errors of the coefficients with a certain fixed p

a_r_value=[] 
#vector with only the values of the coefficients with p and q fixed
a_r_error=[] 
#vector with only the errors of the coefficients with p and q fixed

threshold = 1E-6 
#under this value we neglect the coefficients and we stop computing them

p=1
q=1
r=1
#labels the coefficients a_pqr. Inizialized to 1

coefficients = open("coefficients_3D", "w+")
#file where all the coefficients will be written

valid=True
while valid:
    while valid:
        while valid:
            value, error = integrate.tplquad(integrand, 0, L, 0, L, 0, L) 
            #compute the integral
            if np.abs(value) < threshold: 
                #verify if it is bigger than the threshold
                valid = False
            else:
                a_r_value.append(value)
                #add the value of the integral to the list
                a_r_error.append(error)
                #add the error of the integral to the list 
                print("a_{},{},{}".format(p,q,r), "=", value,"+-", error)
                coefficients.write("a_{},{},{} = {} +- {} \
                               \n".format(p, q, r, value, error))
                r += 2  #when r is even the integral is always zero
        
        r=1 
        q += 2
        #reset the r lable and increase the q lable. 
        #When q is even the integral is always zero
        if len(a_r_value) == 0: break 
#at a certain point, all the coefficients will be smaller than the threshold, 
#so the last element of the list will be an empty list
#this would be a break condition because it implies that 
#all the next lists would be empty lists
        else:
            a_qr_value.append(a_r_value)
            a_qr_error.append(a_r_error)
#add the coefficients with fixed p and q to this list and their relative errors
            a_r_value = []
            a_r_error = []
            #reset of the two lists
            valid=True
        
    q=1
    p += 2
    #reset the q lable and increase the p lable. 
    #When p is even the integral is always zero

    if len(a_qr_value) == 0: break
#at a certain point, all the coefficients will be smaller than the threshold, 
#so the last element of the list will be an empty list
#this would be a break condition because it implies that 
#all the next lists would be empty lists
    else:
        a_pqr_value.append(a_qr_value)
        a_pqr_error.append(a_qr_error)
#add the coefficients with fixed p and q to this list and their relative errors
        a_qr_value = []
        a_qr_error = []
        #reset of the two lists
        valid=True
        
coefficients.close()

t=2E-7 #we fix a particular time 

def n(x, y, z): #define the neutron density
    n = 0
    for i in range(len(a_pqr_value)):
        P = 2 * i + 1
        for j in range(len(a_pqr_value[i])):
            Q = 2 * j + 1
            for k in range(len(a_pqr_value[i][j])):
                R = 2 * k +1            
                n += a_pqr_value[i][j][k] *\
                    np.exp(eta*t -mu*np.pi**2*((P**2)/(L**2) + \
                    (Q**2)/(L**2) + (R**2)/(L**2))*t) * \
                    np.sin(P*np.pi*x/L)* np.sin(Q*np.pi*y/L)* \
                    np.sin(R*np.pi*z/L)
    return n

Z=L/2 #we fix a particular value of Z

#plot
X_arr=np.linspace(0, L, 1000)#range on the x axis
Y_arr=np.linspace(0, L, 1000)#range on the y axis

X, Y = np.meshgrid(X_arr, Y_arr)

fig = plt.figure()

from mpl_toolkits.mplot3d import Axes3D
ax = Axes3D(fig)

ax.plot_surface(X, Y, n(X, Y, Z), cmap ='viridis',\
                edgecolor ='black', color = "white")

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('n(x,y)')

plt.show()


        















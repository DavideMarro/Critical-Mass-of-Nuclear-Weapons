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

L = np.pi*np.sqrt(3*mu/eta)+0.002 #lenght bigger than the critical one
r1=np.sqrt((eta*L**2-np.pi**2*mu)*mu)*L*sp.special.jn_zeros(0, 1)/\
           (eta*L**2 - np.pi**2*mu)+0.002 #radius bigger than the critical one
print("the critical length is ", np.pi*np.sqrt(3*mu/eta), "m")
print("the critical radius is ", np.sqrt((eta*L**2-np.pi**2*mu)*mu)*L\
           *sp.special.jn_zeros(0, 1)/\
           (eta*L**2 - np.pi**2*mu) ,"m")


alpha=sp.special.jn_zeros(0, 100) 
#list with the zeros of the bessel function of the first tipe

integrand = lambda Z,R: (4/(L*(r1**2)*(sp.special.jv(1, alpha[q-1])**2)))\
                        *sp.special.jv(0, alpha[q-1]*R/r1)*R*\
                        (1-((R**2)/(r1**2)))*(np.sin(np.pi*Z/L))**2     
#defines the function that has to be integrated                           

a_1q_value=[] #vector with only the values of the coefficients
a_1q_error=[] #vector with only the errors of the coefficients

threshold = 0.0001 
#under this value we neglect the coefficients and we stop computing them

q=1
#labels the coefficients a_1q. Inizialized to 1

coefficients = open("coefficients_cilindrical", "w+") 
#file where all the coefficients will be written

valid = True
while valid:
    value, error = integrate.dblquad(integrand, 0, r1, 0, L) 
    #compute the integral
    if np.abs(value) < threshold: #verify if it is bigger than the threshold
        valid = False
    else:
        a_1q_value.append(value) #add the value of the integral to the list
        a_1q_error.append(error) #add the error of the integral to the list
        print("a_{},{}".format(1,q), "=", value,"+-", error)
        coefficients.write("a_{},{} = {} +- {}  \n".format(1, q, value, error)) 
    q = q+1  #when q is even the integral is always zero
    
coefficients.close()
    
t=1E-5 #we fix a particular time 
def n(r, z): #define the neutron density
    n = 0
    for q in range(len(a_1q_value)):
        n += a_1q_value[q] * sp.special.jv(0, alpha[q-1]*r/r1)*\
        np.sin(np.pi*z/L)*np.exp(((eta*(r1**2)*(L**2)-mu*((alpha[q-1]**2)*\
        (L**2)+(np.pi**2)*(r1**2))) / ((r1**2)*(L**2)))*t )
    return abs(n)
    
#plot
R_arr=np.linspace(-r1, r1, 100) #range on the r axis
Z_arr=np.linspace(0, L, 100)    #range on the z axis

R, Z = np.meshgrid(R_arr, Z_arr)
fig = plt.figure()

from mpl_toolkits.mplot3d import Axes3D
ax = Axes3D(fig)

ax.plot_surface(R, Z, n(R, Z), cmap ='viridis',     
                edgecolor ='black', color = "white")

ax.set_xlabel('r')
ax.set_ylabel('z')
ax.set_zlabel('n(r,z)')













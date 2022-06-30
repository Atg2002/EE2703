# -*- coding: utf-8 -*-
"""
Created on Thu May 12 13:00:04 2022

@author: Krithikavvgreat
"""
#Importing Libraries
import numpy as np
import math 
from pylab import*

# Defining Matrix_M of Hphi = M*J
def Matrix_M(N):
    return (1/(2*pi*a))*np.identity(2*N-2)

#Dataset used in the assignment 
pi = np.pi                                #Defining pi for ease of use
Im = 1                                    #Defining Max amplitude I              
l = 0.5                                   #Half length of antenna
c=2.9979e8                                #speed of light
mu0=4e-7*pi                               #permeability of free space
a = 0.01                                  #Radius of wire used in transmission
lamda = l*4.0                             #Wavelength 
f=c/lamda                                 #frequency 
k=2*pi/lamda                              #wavenumber
N =100                                    #Number of sections in each half section of the antenna
dz = l/N                                  #spacing of current sample

#Insert function defined to add the known values of I's @indexes 0,N,2N 
def Insert( J,val2, n1 =  N if N%2 ==0 else N+1 , val1 = 0,val3= 0):
    lst = J.tolist()
    lst.insert(0, val1)
    lst.insert(n1, val2)
    lst.append(val3)
    return np.array(lst)

z = np.linspace(-N,N,2*N+1)*dz            #Define points along the antenna as an array
I = Im*np.sin(k*(l - abs(z)))             #Defining current in the wire use modulus of z
u = np.delete(z,[0,N if N%2 ==0 else N+1,2*N])#removing the end(-l & l) and the middle(0)
J = Im*np.sin( k*(l - abs(u)) )           #Finding the value of current at locations in u

M = Matrix_M(N)                           #Finding M matrix corresponding to given N
c1, c2 =  np.meshgrid(z,z)                #Using Meshgrid to map all the point in z
Rz = np.sqrt((c1 - c2)**2 + a**2)         #Finding all the possible distances and defining Rij to be used for Pij
c1, c2 =  np.meshgrid(u,u)                #Using Meshgrid to map all the point in u
Ru = np.sqrt((c1 - c2)**2 + a**2)         #Finding all the possible distances and defining Rij to be used for Pij
RiN = np.sqrt(a**2 + u**2)               #Definging RiN the distance form center(z=0) to all the current sources

#Finding Pb anf Pij as defined by equations in the question
Pb = (mu0/(4*pi)) * (np.exp(-1j*RiN*k)) *dz/ RiN      
P = (mu0/(4*pi)) * ((np.exp(-1j*Ru*k))/ Ru) *dz

#Finding Vector potential ont he given values of z using P and Pb
A = np.dot(P,J) + Pb*Im

#Finding Qbi anf Q as defined by equations in the question
Qbi = Pb*a/mu0*(1j*k/RiN + 1/RiN**2)
Q = P*a/mu0*(1j*k/Ru +1/Ru**2)

#Finding Hphi using Q and Qbi
Hphi = np.dot(Q,J)+ Qbi*Im

#Finding Unknown current using (M-Q)*J_new = Qbi*Im; Other methods of sloving commnted out
'''
J_new = np.linalg.solve(M-Q, Qbi*Im)
J_new = np.dot(inverse(M-Q), Qbi*Im)'''
J_new = np.linalg.lstsq(M-Q, Qbi*Im)[0]
J_new = Insert(J_new,val2 = Im)                         #inserting the known values of 0(@ -l and +l) and 1 (@0)
plot(z,I, label = "True Current Vector")                #Plotting true Current vector obtained in 1. using the sin equation
plot(z,real(J_new), label = "Calculated Current vector")#Plotting J_found via acalculation of vector potential and Magnetic field
grid(True)
legend()
title("I vs z for N ={}".format(N))
xlabel(r"$z(meter)\rightarrow$")
ylabel(r"I(A)")
savefig("fig9-N-{}.png".format(N))
show()

#Printing all the asked matrices:
#Please uncomment to print
'''
print("z:\n",z,"\n")
print("u:\n",u,"\n")
I.shape = (-1,1)
print("\n I:\n", I.round(2))
J.shape = (-1,1)
print("\n J:\n", J.round(2))
J_new.shape = (-1,1)
print("\n J_calculated*1e5:\n", (J_new*1e5).round(2))    
print("\n Rz:\n",Rz.round(2))
print("\n Ru: \n",Ru.round(2))
print("\n P*1e8:\n", (P*1e8).round(2))
Pb = ((Pb*1e8).round(2))
Pb.shape = (-1,1)
print("\n Pb*1e8:\n", Pb )
A = ((A* 1e8).round(2))
A.shape = (-1,1)
print("\n A*1e8 :\n", A )
print("\n Q :\n",(Q.round(2) ))
Qbi.shape = (-1,1)
print("\n Qb:\n", Qbi.round(2) )
Hphi = (Hphi.round(2))
Hphi.shape =  (-1,1)
print("\n Hphi:\n" , Hphi )'''
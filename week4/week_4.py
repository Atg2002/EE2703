import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import scipy
from scipy.integrate import quad
from math import *

PI =np.pi

def exp(x):                                                   # function to find exponention of an array
    return np.exp(x)

def coscos(x):                                                # Function to find cos(cos()) of an array
    return np.cos(np.cos(x))

x = np.arange(-2*PI, 4*PI, 0.01)                              # range of x -2pi to 4pi
x_2 = np.arange(0, 2*PI, 0.01)                                # Restricting the limit to 0-2pi to avoid periodic repetitions

fcosine_1 = lambda x : exp(x) * cos(i*x)                      # getting cos values via a short function for fourier coefficients of exp 
fsine_1 = lambda x :   exp(x) * sin(i*x)                      # getting sin values via a short function for fourier coefficients of exp

fcosine_2 = lambda x : coscos(x) * cos(i*x)                   # getting cos values via a short function for fourier coefficients of exp 
fsine_2 = lambda x :   coscos(x) * sin(i*x)                   # getting sin values via a short function for fourier coefficients of exp

An_exp = []                                                   # Storing an coefficients of exp                                                   
Bn_exp = []                                                   # Storing bn coefficients of exp
exp_FT =0                                                     #Making the approximate function via an and bn upto 26 each(b0 =0)

An_coscos = []                                                # Storing an coefficients of exp                                                   
Bn_coscos = []                                                # Storing bn coefficients of exp
coscos_FT =0                                                  #Making the approximate function via an and bn upto 26 each(b0 =0)

n = 26                                                        #no. of coeffiecients is 26(total will be 52 with 26 an and 26 bn and b0 =0)

for i in range(n):                                            # finding an's and bn's of coscos and exp via quad 
    an = (quad(fcosine_1 , 0.0, 2*PI)[0]) / PI
    bn = (quad(fsine_1 , 0.0, 2*PI)[0]) / PI
    an_1 = (quad(fcosine_2 , 0.0, 2*PI)[0]) / PI
    bn_1 = (quad(fsine_2 , 0.0, 2*PI)[0]) / PI

    An_exp.append(an)
    Bn_exp.append(bn)
    An_coscos.append(an_1)
    Bn_coscos.append(bn_1)

for i in range(n):
    if i ==0 :
        exp_FT +=An_exp[i]/2
        coscos_FT += An_coscos[i]/2
        An_exp[i] =  An_exp[i]/2
        An_coscos[i] = An_coscos[i]/2
    else :
       exp_FT+= ( An_exp[i]*np.cos(i*x_2) + Bn_exp[i]*np.sin(i*x_2) )
       coscos_FT+= ( An_coscos[i]*np.cos(i*x_2) + Bn_coscos[i]*np.sin(i*x_2) )

#plotting expx true and fourier approx
plt.figure(1)
plt.title(r'$e^x$ on a semilogy plot')
plt.xlabel("x")
plt.ylabel("exp(x)")
plt.grid(True)
plt.semilogy(x, exp(x), label ="True")
plt.semilogy(x_2,exp_FT, "og" ,label = "Predicted" )
plt.legend()
plt.savefig('fig1.png')
plt.show()

#plotting coscos true and fourier approx 
plt.figure(2)
plt.grid(True)
plt.title(r'Plot of $cos(cos(x))$')
plt.xlabel("x")
plt.ylabel("cos(cos(x))")
plt.plot(x, coscos(x), label ="True")
plt.plot(x_2,coscos_FT, "og" ,label = "Predicted" )
plt.legend()
plt.savefig('fig2.png')
plt.show()

n_graph = np.arange(0, n, 1)                                        #forming array from 0-n-1

#plot of coeffiecicients of exponential x semilogy
#Note- Rather than plotting coeffiecients as one array I have plotted an's and bn's seperately thats why my n varies to 26 rather than till 52
plt.figure(3)
plt.xlabel("n")
plt.ylabel("fourier coefficients")
plt.title(r"Coefficients of fourier series of $e^x$ on a semilogy scale")
plt.grid(True)
plt.semilogy(n_graph, An_exp, "ro")
plt.semilogy(n_graph, np.abs(Bn_exp), "ro")
plt.savefig('fig3.png')
plt.show()

#plot of coeffiecicients of exponential x loglog
#Note- Rather than plotting coeffiecients as one array I have plotted an's and bn's seperately thats why my n varies to 26 rather than till 52
plt.figure(4)
plt.xlabel("n")
plt.ylabel("fourier coefficients")
plt.title(r"Coefficients of fourier series of $e^x$ on a loglog scale")
plt.grid(True)
plt.loglog(n_graph, An_exp, "ro")
plt.loglog(n_graph, np.abs(Bn_exp), "ro")
plt.savefig('fig4.png')
plt.show()

#plot of coeffiecicients of coscos x
#Note- Rather than plotting coeffiecients as one array I have plotted an's and bn's seperately thats why my n varies to 26 rather than till 52
plt.figure(5)
plt.grid(True)
plt.xlabel("n")
plt.ylabel("fourier coefficients")
plt.title(r"Coefficients of fourier series of $cos(cos(x))$ on a semilogy scale")
plt.semilogy(n_graph, np.abs(An_coscos), "ro")
plt.semilogy(n_graph, np.abs(Bn_coscos), "ro")
plt.savefig('fig5.png')
plt.show()

#plot of coeffiecicients of coscos x
#Note- Rather than plotting coeffiecients as one array I have plotted an's and bn's seperately thats why my n varies to 26 rather than till 52
plt.figure(6)
plt.xlabel("n")
plt.ylabel("fourier coefficients")
plt.title(r"Coefficients of fourier series of $cos(cos(x))$ on a loglog scale")
plt.grid(True)
plt.loglog(n_graph, np.abs(An_coscos), "ro")
plt.loglog(n_graph, np.abs(Bn_coscos), "ro")
plt.savefig('fig6.png')
plt.show()

#fourier coefficients via matrix solving and lstsq 
x_3=np.linspace(0,2*pi,401)
x_3=x_3[:-1]                                    # drop last term to have a proper periodic integral
b=exp(x_3)                                      # f has been written to take a vector
A=np.zeros((400,51))                            # allocate space for A
A[:,0]=1                                        # col 1 is all ones
for k in range(1,26):
    A[:,2*k-1]=np.cos(k*x_3)                    # cos(kx) column
    A[:,2*k]=np.sin(k*x_3)                      # sin(kx) column
#endfor

c1=scipy.linalg.lstsq(A,b)[0]                   # the ’[0]’ is to pull out the
                                                # best fit vector. lstsq returns a list.
A_temp =[]                                      #for storing an's
B_temp = []                                     #for sroring bn's

for i in range(len(c1) ):
    if i ==0 :
         A_temp.append(c1[i])
         B_temp.append(0.0)
    else:
        if i%2 ==0 : B_temp.append(c1[i])
        else : A_temp.append(c1[i])

#plotting coefficients obtained via above method and earlier ones on semilogy
plt.figure(7)
plt.grid(True)
plt.xlabel("n")
plt.ylabel("fourier coefficients")
plt.title(r"Coefficients of fourier series of $e^x$ on a semilogy scale")
plt.semilogy(n_graph, An_exp, "ro",label ="an")
plt.semilogy(n_graph, np.abs(Bn_exp), "ro", label = "bn")
plt.semilogy(n_graph, np.abs(A_temp), "go", label = "an_lstsq")
plt.semilogy(n_graph, np.abs(B_temp), "go", label = "bn_lstsq")
plt.legend()
plt.savefig('fig7.png')
plt.show()

#plotting coefficients obtained via above method and earlier ones on loglog
plt.figure(8)
plt.grid(True)
plt.xlabel("n")
plt.ylabel("fourier coefficients")
plt.title(r"Coefficients of fourier series of $e^x$ on a loglog scale")
plt.loglog(n_graph, An_exp, "ro", label = "an")
plt.loglog(n_graph, np.abs(Bn_exp), "ro", label ="bn")
plt.loglog(n_graph, np.abs(A_temp), "go", label="an_lstsq")
plt.loglog(n_graph, np.abs(B_temp), "go", label ="bn_lstsq")
plt.legend()
plt.savefig('fig8.png')
plt.show()

A_diff = np.abs(A_temp - np.array(An_exp))                           #finding the difference in an's 
print("The max deviation is found for a", np.where( A_diff == np.amax(A_diff))[0] , " = ", np.amax(A_diff)) #printing maximum an

B_diff = np.abs(B_temp - np.array(Bn_exp))                           #finding the difference in bn's
print("The max deviation is found for b", np.where( B_diff == np.amax(B_diff))[0] , " = ", np.amax(B_diff))#printing maximum bn

b_new_exp = np.dot(A,c1)

b=coscos(x_3) # f has been written to take a vector

A[:,0]=1 # col 1 is all ones
for k in range(1,26):
    A[:,2*k-1]=np.cos(k*x_3) # cos(kx) column
    A[:,2*k]=np.sin(k*x_3) # sin(kx) column
#endfor

c1=scipy.linalg.lstsq(A,b)[0] # the ’[0]’ is to pull out the
# best fit vector. lstsq returns a list.
#repeating same for coscos
A_temp =[]
B_temp = []

for i in range(len(c1) ):
    if i ==0 :
         A_temp.append(c1[i])
         B_temp.append(0.0)
    else:
        if i%2 ==0 : B_temp.append(c1[i])
        else : A_temp.append(c1[i])

#plotting coefficients obtained via above method and earlier ones on semilogy
plt.figure(9)
plt.grid(True)
plt.xlabel("n")
plt.ylabel("fourier coefficients")
plt.title(r"Coefficients of fourier series of $cos(cos(x))$ on a semilogy scale")
plt.semilogy(n_graph, np.abs(An_coscos), "ro", label= "an")
plt.semilogy(n_graph, np.abs(Bn_coscos), "ro", label ="bn")
plt.semilogy(n_graph, np.abs(A_temp), "go", label ="an_lstsq")
plt.semilogy(n_graph, np.abs(B_temp), "go", label = "bn_lstsq")
plt.legend()
plt.savefig('fig9.png')
plt.show()

#plotting coefficients obtained via above method and earlier ones on loglog
plt.figure(10)
plt.grid(True)
plt.xlabel("n")
plt.ylabel("fourier coefficients")
plt.title(r"Coefficients of fourier series of $cos(cos(x))$ on a loglog scale")
plt.loglog(n_graph, np.abs(An_coscos), "ro", label= "an")
plt.loglog(n_graph, np.abs(Bn_coscos), "ro", label = "bn")
plt.loglog(n_graph, np.abs(A_temp), "go", label ="an_lstsq")
plt.loglog(n_graph, np.abs(B_temp), "go", label = "bn_lstsq")
plt.legend()
plt.savefig('fig10.png')
plt.show()

A_diff = np.abs(A_temp - np.array(An_coscos))                    # Max difference an's
print("\nThe max deviation is found for a", np.where( A_diff == np.amax(A_diff))[0] , " = ", np.amax(A_diff))#max difference in bn's

B_diff = np.abs(B_temp - np.array(Bn_coscos))                    #Max difference bn's
print("The max deviation is found for b", np.where( B_diff == np.amax(B_diff))[0] , " = ", np.amax(B_diff))#max difference in bn's

b_new_coscos = np.dot(A,c1)

#plotting graph of exp x :true and via coefficients derived through lstsq
plt.figure(11)
plt.grid(True)
plt.title(r'$e^x$ on a semilogy plot')
plt.xlabel("x")
plt.ylabel("exp(x)")
plt.semilogy(x, exp(x), label ="True")
plt.semilogy(x_3,b_new_exp, "go" ,label = "new " )
plt.legend()
plt.savefig('fig11.png')
plt.show()

#pltting graph of coscos x : true and via coefficients derivied through lstsq
plt.figure(12)
plt.title(r'plot of $cos(cos(x))$ ')
plt.xlabel("x")
plt.ylabel("exp(x)")
plt.grid(True)
plt.plot(x, coscos(x), label ="True")
plt.plot(x_3,b_new_coscos, "go" ,label = "new" )
plt.legend()
plt.savefig('fig12.png')
plt.show()
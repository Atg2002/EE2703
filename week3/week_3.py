#importing module : - numpy, scipy, matplotlib
import numpy as np      
import matplotlib.pyplot as plt
import scipy
import scipy.special as sp
#defining A and B
A = 1.05
B = -0.105
sigma = np.logspace(-1, -3, 9)

#defining bassel fn for ease of use
def jn(num, var):
    return sp.jn(num, var)

#defining modulus fn for ease of use
def mod(val):
    return np.abs(val)

#Qn 2 :reading data from fitting.dat
def read_data():
    data = np.loadtxt("fitting.dat", dtype= float,delimiter=" ")
    time = data[0:,0]
    yy = data[0:, 1:]
    return time, yy

#Qn 3: Generating g(t; A,B)
def g(t, A_, B_) :
    return A_*jn(2,t) + B_*t

#Qn 2 : Plotting various noise filled data
def plotFig0_ErrorLabels(t, yy,i):
    plt.figure(0)
    plt.plot(t,yy)
    plt.plot(t, g(t, A, B), color = "black", linewidth =3)

    plt.xlabel("$t\\rightarrow$", size = 13)
    plt.ylabel("$f(t) + noise\\rightarrow$", size = 13)
    label = []
    for i in range(len(sigma)):
        label.append( "sigma" + "$_" + str(i + 1) + "$" + "=" + str(round(sigma[i], 3)) )
    label.append("True Value")
    plt.legend(label)
    plt.title("Qn4 : Plot of f(t) vs t for various noises")
    plt.grid(True)
    plt.savefig('fig{}.png'.format(i))
    plt.show()

#Qn5 : Plotting Errorbars
def plotErrorbar(t, data,i):
    stdv = np.std(data[:,0] - g(t, A, B))
    
    plt.plot(t, g(t, A, B), color = "black", label = "f(t)" )
    plt.errorbar(t[::5],data[::5,0],stdv,fmt='ro', label = "Errorbar")
    plt.legend()
    plt.title("Q5: Data points for σ = 0.10 along with exact function ")
    plt.grid(True)
    plt.savefig('fig{}.png'.format(i))
    plt.show()

#definign c_ fn for combining two couloumn vectors  
def c_(x,y):
    M = np.zeros((len(x),2))

    for i in range(len(y)):
        M[i,0] = x[i]
        M[i,1] = y[i]
    
    return M

# Qn 6 : getting M mtrix
def get_Matrix(t):
    x = [jn(2, t[i]) for i in range(len(t))]
    y= t
    M = c_(x,y)
    return M

#Qn 7 : Forming the error matrix epsilon
def error_matrix(t,data,a,b):
    error = np.zeros((len(a), len(b)))
    for i in range(len(a)):
        for j in range(len(b)): 
            error[i][j] = ( (data[:,0] - g(t, a[i],b[j]))**2 ).mean()

    return error

#Qn 9 : Estimating A and B
def estimate_A_B(M, g):
    return scipy.linalg.lstsq(M,g)

#Qn 10 : Plotting error in A and B vs noise
def plot_error_noise(A_error, B_error,i):
    plt.plot(sigma ,A_error, linestyle = ":", marker ="o", label = "A_error", color = "red"  )
    plt.plot(sigma ,B_error, linestyle = ":", marker ="o", label = "B_error", color = "green")
    plt.legend()
    plt.xlabel("$Noise standard deviation\\rightarrow$", size =13)
    plt.ylabel("$MS error\\rightarrow$", size =13)
    plt.title("Variation of error with noise")
    plt.grid(True)
    plt.savefig('fig{}.png'.format(i))
    plt.show()

#Qn 11 half part 
def plot_error_noise_log(A_error,B_error,i):
    plt.loglog(sigma ,A_error, linestyle = ":", marker ="o", label = "A_error", color = "r" )
    plt.loglog(sigma ,B_error, linestyle = ":", marker ="o", label = "B_error", color = "b")
    plt.xlabel("$Noise standard deviation\\rightarrow$", size =13)
    plt.ylabel("$MS error\\rightarrow$", size =13)
    plt.grid(True)
    plt.legend()
    plt.title("Variation of error with noise in logarithmic scale")
    plt.savefig('fig{}.png'.format(i))
    plt.show()

# Qn 11: Using stem to get proper graph in logarithim
def plotStemPlot(A_error,B_error,i):
    plt.loglog(sigma ,A_error, linestyle = "", marker ="o", label = "A_error", color = "r" )
    plt.loglog(sigma ,B_error, linestyle = "", marker ="o", label = "B_error", color = "b")
    plt.stem(sigma ,A_error, "ro", label="Aerr")
    plt.stem(sigma ,B_error, "b",  label= "Berr")
    plt.xlabel("$Noise standard deviation\\rightarrow$", size =13)
    plt.ylabel("$MS error\\rightarrow$", size =13)
    plt.title("Qn11: Variation of error with noise(in stem)")
    plt.legend
    plt.savefig('fig{}.png'.format(i))
    plt.show()

#storing time and data in t and yy
t, yy = read_data()
plotFig0_ErrorLabels(t,yy,1)
plotErrorbar(t, yy,2)
M = get_Matrix(t)

p = np.array((A,B))

result = np.allclose( g(t, A, B),np.dot(M,p), 1e-05, 1e-08 )
print("The result of g(t, A, B) == np.dot(M,p) is : ", result)

A_new = np.arange(0,2.1,0.1)
B_new = np.linspace(-0.2,0,21)

mse = error_matrix(t,yy,A_new,B_new)

#Qn 8 : Plotting Contours
X,Y =np.meshgrid(A_new,B_new)
contour_val = np.linspace(0.025,0.5, 20)
temp = plt.contour(X,Y,mse, contour_val)
plt.clabel(temp, temp.levels[:4], inline=1, fontsize=8)
plt.plot(1.05, -0.105, 'ro', color ="red", linewidth = 3)
plt.annotate("Exact location", (A,B), xytext =(0.9,-0.1),weight = "bold", color = "black", size =10)
plt.xlabel("$A\\rightarrow$", size =13)
plt.ylabel("$B\\rightarrow$", size =13)
plt.title("Qn8: contour plot of ε $_i$$_j$")
plt.savefig('fig{}.png'.format(3))
plt.show()

AB_array = estimate_A_B(M, yy[:,0])[0]
print("estimated value of A = ", AB_array[0],"\nestimated value of B = ", AB_array[1])

A_error = [mod(AB_array[0] - A)]
B_error = [mod(AB_array[1] - B)]

for i in range(1,len(sigma)):
    AB_array = estimate_A_B(M, yy[:,i])[0]

    A_error.append(mod(AB_array[0]-A))
    B_error.append(mod(AB_array[1]-B))

plot_error_noise(A_error,B_error,4)
plot_error_noise_log(A_error,B_error,5)
plotStemPlot(A_error,B_error,6)
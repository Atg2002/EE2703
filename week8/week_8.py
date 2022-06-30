import numpy as np
from pylab import*
import pandas as pd
pi = np.pi                                             #Defining pi for ease of use

#Defined to obtain true fourier transfrom of gaussian fn
def true_gauss_fft(w):
    return np.exp(-w**2/2)*sqrt(2*pi)

# graph spectrum fn to plot spectrum for a given function
def graph_specturm(x,y,x_range,title, magnitude = -1, return_error = False ):
    w = linspace(-32,32, len(x)+1)[:-1]               #Defining omega
    Y = fftshift(fft(y))/len(x)                       #Obtaing fft(y)
    if return_error == True:                          #returning error if asked
        Y = fftshift(abs(fft(y)))/len(x)         
        Y = Y*((2*pi)**0.5)/(Y.max())
        error = abs(true_gauss_fft(w) -Y)           
        return error.max()                            
    fig, (ax1,ax2) = subplots(2,1)                    #Plotting magnitude and phase subplots
    ax1.set_title(title)
    ax1.plot(w, abs(Y))
    ax1.set_ylabel(r"$|Y|$")
    ax1.grid(True)
    ax1.set_xlim(x_range)
    ax2.set_xlabel(r"$k$")
    ax2.plot(w, angle(Y),"ro")
    ax2.set_ylabel(r"Phase of $Y$")
    ax2.grid(True)
    if magnitude>=0:
        ii=where(abs(Y)>magnitude)
        ax2.plot(w[ii],angle(Y[ii]),'go')
    ax2.set_xlim(x_range)
    graph_specturm.count+=1
    savefig("fig8-{}.png".format(graph_specturm.count))
    show()
graph_specturm.count = 0                              #Defined to save the graphs

#Various sampling types for better visualisation
x = np.linspace(0, 2*pi, 129)
x = x[:-1]

x1 = linspace(-4*pi,4*pi, 513)
x1 = x1[:-1]

#Qn1: a. fft of sin(5t)
graph_specturm(x = x, y = sin(5*x),x_range = [-10,10],title =  r"Spectrum of $\sin(5t)$", magnitude=1e-3)
#Qn1: b. fft of (1+0.1cos(t)cos10t
graph_specturm(x = x,y = cos(10*x)*(1+0.1*cos(x)),x_range = [-15,15],title = r"Spectrum of $(1+0.1\cos(t)*\cos10t$")
#Qn1 : c. fft of (1+cos(t)cos10t
graph_specturm(x = x1,y = cos(10*x1)*(1+0.1*cos(x1)),x_range = [-15,15],title = r"Spectrum of $(1+0.1\cos(t)*\cos10t$")
#Qn2 : a. fft of sin^3(t)
graph_specturm(x = x,y = (sin(x))**3,x_range = [-5,5],title = r"Spectrum of $(\sin(t))^{3}$", magnitude=1e-2)
#Qn2 : b. fft of cos^3(t)
graph_specturm(x = x,y = (cos(x))**3,x_range = [-5,5],title = r"Spectrum of $(\cos(t))^{3}$", magnitude=1e-2)
#Qn3 : fft of cos(20t +5cos(t))
graph_specturm(x = x1,y = cos(20*x1 +5*cos(x1))  ,x_range = [-20,20],title = r"Spectrum of $\cos(20t +5\cos(t))$", magnitude=1e-3)
#Qn4 : fft of e^(-t^2 / 2)

Data = pd.DataFrame(columns = ["Interval","Samples Taken","Frequency", "Max Error"])
#Iterations to store each interval and sampling valuse corresponding error in a dataframe
for i in range(0,16):
    if i<4:
        x = np.linspace(-2*pi,2*pi,2**(i+7)+1)
        Data.loc[i] = [ [x[0],x[-1]], 2**(i+7),len(x)/( 2*pi ), graph_specturm(x = x,y = exp(-(x**2)/2),x_range = [-25,25],title = r"Spectrum of $e^{-t^{2}/2}$", return_error=True) ] 
    elif i<8:
        x = np.linspace(-4*pi,4*pi,2**(i+3)+1)
        Data.loc[i] = [ [x[0],x[-1]], 2**(i+3),len(x)/( 4*pi ), graph_specturm(x = x,y = exp(-(x**2)/2),x_range = [-25,25],title = r"Spectrum of $e^{-t^{2}/2}$", return_error=True) ] 
    elif i<12:
        x = np.linspace(-8*pi,8*pi,2**(i-1)+1)[:-1]
        Data.loc[i] = [ [x[0],x[-1]], 2**(i-1),len(x)/( 8*pi ), graph_specturm(x = x,y = exp(-(x**2)/2),x_range = [-25,25],title = r"Spectrum of $e^{-t^{2}/2}$", return_error=True) ] 
    elif i<16:
         x = np.linspace(-12*pi,12*pi,2**(i-5)+1)
         Data.loc[i] = [ [x[0],x[-1]],2**(i-5) ,len(x)/( 8*pi ), graph_specturm(x = x,y = exp(-(x**2)/2),x_range = [-25,25],title = r"Spectrum of $e^{-t^{2}/2}$", return_error=True) ] 
    i+=1
Data.style
#Plotting the graph which has least deviation from true fft
error = Data["Max Error"].min()
row = Data.loc[Data["Max Error"] == error]
interval = row.iloc[0]['Interval']
samp = row.iloc[0]['Samples Taken']
print("The minimum error from all the above iterations occurs for {} interval with {} samples taken with max error being = {}\n The corresponding Graph is:-\n".format(interval, samp, error) )   
x = linspace(interval[0], interval[1], samp+1)[:-1]
graph_specturm(x = x,y = exp(-(x**2)/2),x_range = [-25,25],title = r"Spectrum of $e^{-t^{2}/2}$", magnitude=1e-3)
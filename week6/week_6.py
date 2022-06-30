import scipy.signal as sp
import numpy as np
import matplotlib.pyplot as plt

#Qn 1 : Time response by calculating lapalce tranform and using impulse for decay =0.5
X = sp.lti([1, 0.5], np.polymul([1,0,2.25], [1,1,2.5]))       #X(s) calculated
t = np.linspace(0, 50, 501)                                   #Defining time interval

t, x = sp.impulse(X, None,t )                                 #Finding x(t) with impulse and plotting
plt.figure(1)
plt.plot(t,x)
plt.xlabel(r"Time(s)$\rightarrow$")
plt.ylabel(r"x(t)$\rightarrow$")
plt.title("Forced oscillator with decay = 0.5")
plt.savefig('fig1.png')
plt.show()

#Qn2 :Time response by calculating lapalce tranform and using impulse for decay =0.05 (smaller)
X = sp.lti([1, 0.05], np.polymul([1,0,2.25], [1,0.1,2.2525])) #X(s) calculated

t, x = sp.impulse(X, None,t )                                 #x(t) with impulse fn and plotting
plt.figure(2)
plt.xlabel(r"Time(s)$\rightarrow$")
plt.ylabel(r"x(t)$\rightarrow$")
plt.title("Forced oscillator with decay = 0.05")
plt.plot(t,x)
plt.savefig('fig2.png')
plt.show()

#Qn 3:  Solving for decay of 0.05, and for frequencies 1.4,1.45,1.5,1.55,1.6
H = sp.lti([1], [1,0,2.25])                                  #Impulse response of system
frequencies = np.linspace(1.4,1.6,5)                         #Frequency range to see response of system at each
t = np.linspace(0, 150, 1000)                                #time interval
i=1
for freq in frequencies:                                     #Plotting seperately
    f = (np.cos(freq*t)) * (np.exp(-0.05*t))
    t, x, svec = sp.lsim(H, f,t )
    plt.xlabel(r"Time(s)$\rightarrow$")
    plt.ylabel(r"x(t)$\rightarrow$")
    string = "Forced oscillator with frquency = {}".format(freq)
    plt.title(string)
    plt.plot(t,x)
    plt.savefig('fig3.{}.png'.format(i))
    i+=1
    plt.show()
i=1
for freq in frequencies:                                     #Plotting the responses toghether to visualize the difference
    f = (np.cos(freq*t)) * (np.exp(-0.05*t))
    t, x, svec = sp.lsim(H, f,t )

    string = "frquency = {}".format(freq)
    plt.plot(t,x, label = string )
    plt.savefig("fig4.{}.png".format(i))
    i+=1
plt.title("Combined Visualization for different frequencies")
plt.xlabel(r"Time(s)$\rightarrow$")
plt.ylabel(r"x(t)$\rightarrow$")
plt.legend()
plt.show()


#Qn4: Coupled oscillations of x and y
X = sp.lti ( [1,0,2], [1,0,3,0] )                           #X(s)
Y = sp.lti ( [2], [1,0,3,0] )                               #Y(s)
t = np.linspace(0,20, 501)

tx,x = sp.impulse(X, None, t)                               #Plotting x(t) and y(t)
ty, y = sp.impulse(Y, None, t)
plt.plot(tx,x, label ="x(t)")
plt.plot(ty, y, label = "y(t)")
plt.xlabel(r"Time(s)$\rightarrow$")
plt.ylabel(r"x(t) and y(t)$\rightarrow$")
plt.title("Coupled oscillations x and y")
plt.legend()
plt.savefig('fig5.png')
plt.show()

#Qn 5 : Bode plot of Given transfer function
H = sp.lti([1],[10**-12, 10**-4, 1])

w,S,phi=H.bode()
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.semilogx(w,S)
ax2 = fig.add_subplot(212)
ax2.semilogx(w,phi)
ax1.title.set_text('Magnitude response')
ax2.title.set_text('Phase response')
ax1.set_xlabel(r'$\omega$')
ax1.set_ylabel('|H(jw)|l')
ax2.set_ylabel("/_H(jw)")
ax2.set_xlabel(r'$\omega$')
plt.savefig('fig6.png')
plt.show()


#Qn6 : Calculation of vout for given vin for t in micro and milli seconds
#For microseconf timescale
t= np.linspace(0, 30e-6,10000)
vi = np.cos( (1e3) *t) -np.cos( (1e6) *t )
t, x, svec = sp.lsim(H, vi,t )
plt.xlabel(r"Time(s)$\rightarrow$")
plt.ylabel(r"$v_{o}(t)\rightarrow$")
plt.title("RLC response for time in microseconds")
plt.plot(t,x)
plt.savefig('fig7.png')
plt.show()

#for milli second timescale
t= np.linspace(0, 30e-3,10000)
vi = np.cos( (1e3) *t) -np.cos( (1e6) *t )
t, x, svec = sp.lsim(H, vi,t )
plt.xlabel(r"Time(s)$\rightarrow$")
plt.ylabel(r"$v_{o}(t)\rightarrow$")
plt.title("RLC response for time in milliseconds")
plt.plot(t,x)
plt.savefig('fig8.png')
plt.show()
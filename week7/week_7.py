from sympy import*
import scipy.signal as sp
import numpy as np
import matplotlib.pyplot as plt

PI = np.pi                                                                     #Defining PI constant
s = symbols('s')                                                               #Defining s as a symbol
#High pass filter circuit
def highpass(R1,R2,C1,C2,G,Vi):
    A = Matrix([[0,0,1,-1/G],[-(s*C2*R2)/(1+s*R2*C2),1,0,0],[0,-G,G,1],[(s*C1 + s*C2+ 1/R1),-s*C2,0,-1/R1]])
    b = Matrix([0,0,0,Vi*s*C1])
    V = A.inv() *b
    return V[3]

#Low pass filter circuit
def Lowpass(R1,R2,C1,C2,G,Vi):
    A = Matrix([[0,0,1,-1/G],[-1/(1+s*R2*C2),1,0,0],[0,-G,G,1],[-1/R1-1/R2-s*C1,1/R2,0,s*C1]])
    b=Matrix([0,0,0,-Vi/R1])
    V = A.inv() *b
    return V[3]

def get_rational_coeffs(expr):
    num, denom = expr.as_numer_denom()
    num = [float(i) for i in Poly(num, s).all_coeffs()]
    denom = [float(i) for i in Poly(denom, s).all_coeffs()]
    return num,denom

def vi(t):
    return np.heaviside(t, 1)*(np.sin(2000*PI*t) + np.cos(2e+6 *PI *t))

def vi_decay(t, freq, decay):
    return np.heaviside(t, 1)*np.cos(freq*t)*np.exp(-decay*t) 

#Various timescales used for plotting; this is done to better visualise the graphs
t0 = np.linspace(0,1,1000) 
t1 = np.linspace(0,0.1,100000)
t2 = np.linspace(0,0.001,100000)
t3 = np.linspace(0, 0.00001, 1000000)
t4 = np.linspace(0,0.01, 100000)
#Finding H for low pass filter by passing Vi = 1(s domanin)
Vo=Lowpass(10000,10000,1e-9,1e-9,1.586,1)

#Defining H for low pass filter
num, denom = get_rational_coeffs(Vo)
H = sp.lti(num, denom) 
w,S,phi=H.bode()
#Plotting magnitude response of low pass filter
plt.title("Magnitude plot of LowPass Filter")
plt.xlabel(r"$\omega \rightarrow$")
plt.ylabel(r"$20log(|H(j\omega)|)$")
plt.semilogx(w,S)
plt.savefig('fig1.png')
plt.show()
#Plotting pahse response of low pass filter 
plt.title("Phase plot of LowPass Filter")
plt.xlabel(r"$\omega \rightarrow$")
plt.ylabel(r"$\angle H(j\omega)$")
plt.semilogx(w,phi)
plt.savefig('fig2.png')
plt.show()


#defining u(t)
u = [1 for i in t2 if i>=0 ]   

# Qn1: Plotting step response for low pass filter 
t2,y,svec = sp.lsim(H,u,t2)
plt.title("Step response of given Low pass filter")   
plt.xlabel(r"$t \rightarrow$")    
plt.ylabel(r"$V_{out}$") 
plt.plot(t2,y)
plt.savefig('fig3.png')
plt.show()

#t = symbols('t', real = True)
#Qn2: Output response for given sinosuidal input Low pass filter
plt.title("$V_{in}$ = [$sin(2000\pi t)$ + $cos(2x10^6\pi t)$]u(t)(Low pass filter)")   
plt.xlabel(r"$t \rightarrow$")    
plt.ylabel(r"$V_{out}$") 
plt.plot(t1,vi(t1))
plt.savefig('fig4.png')
plt.show()
fig, (ax1, ax2) = plt.subplots(2,1) 
t1, vo, svec = sp.lsim(H, vi(t1), t1)
ax1.set_title("Vout for Low pass filter$")   
ax1.set_xlabel(r"$t \rightarrow$")    
ax1.set_ylabel(r"$V_{out}$") 
ax1.plot(t1,vo)
#zoomed in version ie t value in 1e-5
t2, vo, svec = sp.lsim(H, vi(t2), t2)  
ax2.set_xlabel(r"$t \rightarrow$ (zoomed in)")    
ax2.set_ylabel(r"$V_{out}$") 
ax2.plot(t2,vo)
plt.savefig('fig6.png')
plt.show()

#Qn 3: Output response for given sinosuidal input High pass filter
Vo=highpass(10000,10000,1e-9,1e-9,1.586,1)                                     #Vo for Vi = 1(in s domain)
num, denom = get_rational_coeffs(Vo)
H1 = sp.lti(num, denom)                                                        #H1 is impulse reponse for high pass filter in s domain

w,S,phi=H1.bode()
#Plotting magnitude response of low pass filter
plt.title("Magnitude plot of HighPass Filter")
plt.xlabel(r"$\omega \rightarrow$")
plt.ylabel(r"$20log(|H(j\omega)|)$")
plt.semilogx(w,S)
plt.savefig('fig7.png')
plt.show()
#Plotting pahse response of low pass filter 
plt.title("Phase plot of HighPass Filter")
plt.xlabel(r"$\omega \rightarrow$")
plt.ylabel(r"$\angle H(j\omega)$")
plt.semilogx(w,phi)
plt.savefig('fig8.png')
plt.show()

t3, vo, svec = sp.lsim(H1, vi(t3), t3)                                         #vo for given sinouidal input
plt.plot(t3,vo) 
plt.title("Vout for High pass filter$")   
plt.xlabel(r"$t \rightarrow$")    
plt.ylabel(r"$V_{out}$")    
plt.savefig('fig9.png')                                                
plt.show()

#Qn 4:
#Vout for high frequency input:-
#1.low pass filter
plt.title(r" Given high frequency $V_{in}$ = $cos(10^5t)e^{-5t}$")   
plt.xlabel(r"$t \rightarrow$")    
plt.ylabel(r"$V_{in}$") 
plt.plot(t0,vi_decay(t0, 1e5,5))
plt.savefig('fig10.png')
plt.show()
t0, vo, svec = sp.lsim(H, vi_decay(t0, 1e5,5), t0)
plt.title("Vout for low pass filter for given decaying sinosouidal input(high frequency)$")   
plt.xlabel(r"$t \rightarrow$")    
plt.ylabel(r"$V_{out}$") 
plt.plot(t0,vo)
plt.savefig('fig11.png')
plt.show()
#2.High pass filter

fig, (ax1, ax2) = plt.subplots(2,1) 
t4, vo, svec = sp.lsim(H1, vi_decay(t4, 1e5,5), t4)
ax1.set_title("Vout for High pass filter for an decaying sinosouidal input(high frequency)")   
ax1.set_xlabel(r"$t \rightarrow$")    
ax1.set_ylabel(r"$V_{out}$") 
ax1.plot(t4,vo)
#zoomed in version ie t value in 1e-5
t2, vo, svec = sp.lsim(H1, vi_decay(t2, 1e5,5), t2)
ax2.set_xlabel(r"$t \rightarrow$ zoom in")    
ax2.set_ylabel(r"$V_{out}$") 
ax2.plot(t2,vo)
plt.savefig('fig12.png')
plt.show()


#Vout for low frquency input
#1.Low pass filter
#Plotting Vi
t1 = np.linspace(0,10, 10000)
plt.title("given low frequency $V_{in}$ = $cos(6\pi t)e^{-0.5t}$")   
plt.xlabel(r"$t \rightarrow$")    
plt.ylabel(r"$V_{in}$") 
plt.plot(t1,vi_decay(t1, 6*PI,0.5))
plt.savefig('fig13.png')
plt.show()

#Plotting Vout for LPF
t1, vo, svec = sp.lsim(H, vi_decay(t1, 6*PI,0.5), t1)
plt.title("$V_{out}$ for given low frequency input(low pass filter)")   
plt.xlabel(r"$t \rightarrow$")    
plt.ylabel(r"$V_{out}$") 
plt.plot(t1,vo)
plt.savefig('fig14.png')
plt.show()

#Plotting Vout for HPF
fig, (ax1,ax2) = plt.subplots(2,1)
t2, vo, svec = sp.lsim(H1, vi_decay(t2, 6*PI,0.5), t2)
ax1.set_title("$V_{out}$ for given low frequency input(High pass filter)")   
ax1.set_xlabel(r"$t \rightarrow$")    
ax1.set_ylabel(r"$V_{out}$") 
ax1.plot(t2,vo)
#zoomed in version
t3, vo, svec = sp.lsim(H1,  vi_decay(t3, 6*PI,0.5), t3)
ax2.set_xlabel(r"$t \rightarrow$ zoomed in")    
ax2.set_ylabel(r"$V_{out}$") 
ax2.plot(t3,vo)
plt.savefig('fig15.png')
plt.show()

#Qn 5:
t1 = np.linspace(0, 0.0001,1000)
Vo = highpass(10000,10000,1e-9,1e-9,1.586,1/s)

num , denom = get_rational_coeffs(Vo)
Vo = sp.lti(num, denom)
t1,vo = sp.impulse(Vo, None, t1)
plt.title("Step response of the given High Pass filter")   
plt.xlabel(r"$t \rightarrow$")    
plt.ylabel(r"$V_{out}$") 
plt.plot(t1, vo)
plt.savefig('fig16.png')
plt.show()
from pylab import*
import numpy as np
#Qn 1:
# Running the example codes in the file
#Spectrum of sin(root2 t)
t = linspace(-pi,pi,65)[:-1]
dt=t[1]-t[0];fmax=1/dt
y=sin(sqrt(2)*t)
y[0]=0 # the sample corresponding to -tmax should be set zeroo
y=fftshift(y) # make y start with y(t=0)
Y=fftshift(fft(y))/64.0
w=linspace(-pi*fmax,pi*fmax,65);w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-10,10])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\sin\left(\sqrt{2}t\right)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-10,10])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
savefig("fig0.png")
show()
#Plotting y for -3pi to 3pi
t1=linspace(-pi,pi,65);t1=t1[:-1]
t2=linspace(-3*pi,-pi,65);t2=t2[:-1]
t3=linspace(pi,3*pi,65);t3=t3[:-1]
figure(2)
plot(t1,sin(sqrt(2)*t1),'b',lw=2)
plot(t2,sin(sqrt(2)*t2),'r',lw=2)
plot(t3,sin(sqrt(2)*t3),'r',lw=2)
ylabel(r"$y$",size=16)
xlabel(r"$t$",size=16)
title(r"$\sin\left(\sqrt{2}t\right)$")
grid(True)
savefig("fig9-2.png")
show()
#Plotting the periodic extension of sin(root2 t) obtained in the interval (-pi,pi)
y=sin(sqrt(2)*t1)
figure(3)
plot(t1,y,'bo',lw=2)
plot(t2,y,'ro',lw=2)
plot(t3,y,'ro',lw=2)
ylabel(r"$y$",size=16)
xlabel(r"$t$",size=16)
title(r"$\sin\left(\sqrt{2}t\right)$ with $t$ wrapping every $2\pi$ ")
grid(True)
savefig("fig9-3.png")
show()
#Plotting the magnitude response of ramp function and seeing that it decays as 1/w
y=t
y[0]=0 # the sample corresponding to -tmax should be set zeroo
y=fftshift(y) # make y start with y(t=0)
Y=fftshift(fft(y))/64.0
w=linspace(-pi*fmax,pi*fmax,65);w=w[:-1]
figure()
semilogx(abs(w),20*log10(abs(Y)),lw=2)
xlim([1,10])
ylim([-20,0])
xticks([1,2,5,10],["1","2","5","10"],size=16)
ylabel(r"$|Y|$ (dB)",size=16)
title(r"Spectrum of a digital ramp")
xlabel(r"$\omega$",size=16)
grid(True)
savefig("fig9-4.png")
show()
#Windowing and plotting time response
n=arange(64)
wnd=fftshift(0.54+0.46*cos(2*pi*n/63))
y=sin(sqrt(2)*t1)*wnd
figure(3)
plot(t1,y,'bo',lw=2)
plot(t2,y,'ro',lw=2)
plot(t3,y,'ro',lw=2)
ylabel(r"$y$",size=16)
xlabel(r"$t$",size=16)
title(r"$\sin\left(\sqrt{2}t\right)\times w(t)$ with $t$ wrapping every $2\pi$ ")
grid(True)
savefig("fig9-5.png")
show()
#Windowed functions DFT
y=sin(sqrt(2)*t)*wnd
y[0]=0 # the sample corresponding to -tmax should be set zeroo
y=fftshift(y) # make y start with y(t=0)
Y=fftshift(fft(y))/64.0
w=linspace(-pi*fmax,pi*fmax,65);w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-8,8])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\sin\left(\sqrt{2}t\right)\times w(t)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-8,8])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
savefig("fig9-6.png")
show()
#Using more samples to get better accuracy in the above DFT
t4=linspace(-4*pi,4*pi,257);t4=t4[:-1]
dt4=t4[1]-t4[0];fmax1=1/dt4
n1=arange(256)
wnd1=fftshift(0.54+0.46*cos(2*pi*n1/256))
y=sin(sqrt(2)*t4)
y=y*wnd1
y[0]=0 # the sample corresponding to -tmax should be set zeroo
y=fftshift(y) # make y start with y(t=0)
Y=fftshift(fft(y))/256.0
w1=linspace(-pi*fmax1,pi*fmax1,257);w1=w1[:-1]
figure()
subplot(2,1,1)
plot(w1,abs(Y),'b',w1,abs(Y),'bo',lw=2)
xlim([-4,4])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\sin\left(\sqrt{2}t\right)\times w(t)$")
grid(True)
subplot(2,1,2)
plot(w1,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
savefig("fig9-7.png")
show()
#Qn 2 :
#Spectrum of cos^3(w0t) without and with a Hamming window
t = linspace(-4*pi,4*pi,257); t = t[:-1]                #defining t
w0 = 0.86                                               #Given w0
dt = t[1]-t[0]; fmax = 1/dt                             #getting fmax
n = arange(256)                                         #defining n for wnd
wnd = fftshift(0.54+0.46*cos(2*pi*n/256))               #wnd the hamming function
y = (cos(w0*t))**3                                      #deyfining 
y[0] = 0
y = fftshift(y)
Y=fftshift(fft(y))/256
w = linspace(-pi*fmax,pi*fmax,257); w = w[:-1]
figure(2)
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-4,4])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\cos^{3}\left(\omega_{0}t\right)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
savefig("fig9-8.png")
show()
#Spectrum of cos^3(w0t)*wnd i.e. with windowing
y = cos(w0*t)**3
y = y*wnd
y[0] = 0
y = fftshift(y)
Y = fftshift(fft(y))/256.0
figure(3)
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-4,4])
ylabel(r"$|Y|\rightarrow$",size=16)
title(r"Spectrum of $\cos^{3}(\omega_{0}t)$ with Hamming window")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"Phase of $Y\rightarrow$",size=16)
xlabel(r"$\omega\rightarrow$",size=16)
grid(True)
savefig("fig9-9.png")
show()
#Qn3
#Estimating w0 and delta using weighted average agianst |Y(w)|
w0 =1.5
d =0.5
t = linspace(-pi,pi,129)[:-1]
dt = t[1]-t[0]; fmax = 1/dt
n = arange(128)
wnd = fftshift(0.54+0.46*cos(2*pi*n/128))
y = cos(w0*t +d)*wnd
y[0]=0
y = fftshift(y)
Y = fftshift(fft(y))/128.0
w = linspace(-pi*fmax,pi*fmax,129)[:-1]
figure(4)
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-4,4])
ylabel(r"$|Y|\rightarrow$",size=16)
title(r"Spectrum of $\cos(w_0t+\delta)$ with Hamming window")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"Phase of $Y\rightarrow$",size=16)
xlabel(r"$\omega\rightarrow$",size=16)
grid(True)
savefig("fig9-10.png")
show()
#Calculating w0 with weighted averages
indexes = where(w>=0)
w_cal = sum(abs(Y[indexes])**2*w[indexes])/sum(abs(Y[indexes])**2)
#We will calculate delta by using phase. We will find this using the w closest to w0
i = abs(w-w_cal).argmin()
delta = angle(Y[i])
print("Calculated value of w0 without noise: ",w_cal)
print("Calculated value of delta without noise: ",delta)

#Qn4:
#Adding "white gaussian noise"
wgn = 0.1*randn(128)             #Generating the "white gaussian noise"
y = cos(w0*t +d)*wnd + wgn*wnd   #Calculating y by adding to original y and then multiplying with Hamming fn
y[0]=0                           #Setting y(-pi) = 0
y = fftshift(y)                  #Calculating freq response of this y
Y = fftshift(fft(y))/128.0   
figure(5)                        #Plotting
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-4,4])
ylabel(r"$|Y|\rightarrow$",size=16)
title(r"Spectrum of $\cos(w_0t+\delta)$ with Hamming window and white gaussian noise")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"Phase of $Y\rightarrow$",size=16)
xlabel(r"$\omega\rightarrow$",size=16)
grid(True)
savefig("fig9-11.png")
show()
#Finding w0 and delta by the method mentioned in qn3 
indexes = where(w>=0)
w_cal = sum(abs(Y[indexes])**2*w[indexes])/sum(abs(Y[indexes])**2)
i = abs(w-w_cal).argmin()
delta = angle(Y[i])
print("Calculated value of w0 without noise: ",w_cal)
print("Calculated value of delta without noise: ",delta)

#Qn 5
#plotting DFT of chirped signal
t = linspace(-pi,pi,1025)[:-1]
dt = t[1] - t[0]; fmax = 1/dt
n = arange(1024)
wnd = fftshift(0.54 + 0.46*cos(2*pi*n/1024))
y = cos(16*(1.5+ t/(2*pi))*t)*wnd
y[0] =0
y = fftshift(y)
Y = fftshift(fft(y)) /1024
w = linspace(-pi*fmax,pi*fmax,1025)[:-1]
figure(6)
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-100,100])
ylabel(r"$|Y|\rightarrow$",size=16)
title(r"Spectrum of $\cos(16t(1.5+\frac{t}{2\pi}))$ with Hamming window")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-100,100])
ylabel(r"Phase of $Y\rightarrow$",size=16)
xlabel(r"$\omega\rightarrow$",size=16)
grid(True)
savefig("fig9-12.png")
show()

#Qn6:
#Finding time-freq variation of DFT using surface plot 
n = arange(64)
wnd = fftshift(0.54 + 0.46*cos(2*pi*n/64))
Wnd = meshgrid(wnd, wnd)[0][:16]                    #Forming Wnd from the meshgrid of wnd ; to be used to find
t_arrays = np.array(split(t,16))                    # Splitting t into 64 rows
y_arrays = cos(16*(1.5+ t_arrays/(2*pi))*t_arrays)*Wnd#Find y_arrays 
y_arrays[:0] = 0                                    #putting y(-pi) =0 in all rows
y_arrays = fftshift(y_arrays)                       #computing DFT
Y_arrays = fftshift(fft(y_arrays))

t = t[::64]	
w = linspace(-fmax*pi,fmax*pi,64+1); w = w[:-1]
t,w = meshgrid(t,w)
Y_mag = abs(Y_arrays)                               #Finding Y_magnitude
fig1 = figure(7)
ax = fig1.add_subplot(111, projection='3d')
surf=ax.plot_surface(w,t,Y_mag.T,cmap='viridis',linewidth=0, antialiased=False)
fig1.colorbar(surf, shrink=0.5, aspect=5)
ax.set_title('surface plot');
ylabel(r"$\omega\rightarrow$")
xlabel(r"$t\rightarrow$")
savefig("fig9-13.png")
show()

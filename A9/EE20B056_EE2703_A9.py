from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3

def dftplot(x,w,xl,tit):
    X=fft(x)
    X=fftshift(X)/len(X)

    figure()

    subplot(2,1,1)
    title(tit)
    ylabel(r'$|X|$')
    xlim(-xl,xl)
    plot(w,abs(X),'g')

    subplot(2,1,2)
    xlabel(r'$\omega$')
    ylabel(r'Phase of $X$')
    xlim(-xl,xl)
    plot(w,unwrap(angle(X)),'b')

    show()

#1 DFT of sin(sqrt(2)t) without windowing

t=linspace(-pi,pi,65)
t=t[:-1]
fm1=1/(t[1]-t[0])
y=sin(2**0.5*t)
y[0]=0
w=pi*linspace(-fm1,fm1,65)[:-1]
dftplot(y,w,10,r'sin($\sqrt{2}t)$')

#1 DFT of sin(sqrt(2)t) with windowing

t=linspace(-4*pi,4*pi,257)
t=t[:-1]
fm=1/(t[1]-t[0])
n = arange(256)
wind = fftshift(0.54+0.46*cos(2*pi*n/256))
y=sin(2**0.5*t)*wind
y[0]=0
w=pi*linspace(-fm,fm,257)[:-1]
dftplot(y,w,8,r'sin($\sqrt{2}t)$ with windowing')

#2 DFT of cos^3(w0t) with and without windowing

y = cos(0.86*t)**3
y1 = y*wind
y[0]=0
y1[0]=0
Y = fftshift(y)
Y1 = fftshift(y1)

dftplot(Y,w,4,r"Spectrum of $\cos^{3}(0.86t)$ without Hamming window")
dftplot(Y1,w,4,r"Spectrum of $\cos^{3}(0.86t)$ with Hamming window")

#3 w0 and d determination from a given cos(w0t+d) sample vector
# w0=1.234 d=0.456

w0 = 1.456
d = 0.456

t = linspace(-pi,pi,129)[:-1]
dt = t[1]-t[0]; fmax = 1/dt
n = arange(128)
wnd = fftshift(0.54+0.46*cos(2*pi*n/128))
y = cos(w0*t + d)*wnd
y[0]=0
y = fftshift(y)
Y = fftshift(fft(y))/128.0
w = linspace(-pi*fmax,pi*fmax,129); w = w[:-1]
dftplot(y,w,4,r"Spectrum of $\cos(w_0t+\delta)$ with Hamming window")

# w0 is calculated by finding the weighted average of all w>0. Delta is found by calculating the phase at w closest to w0.
ii = where(w>=0)
w_cal = sum(abs(Y[ii])**2*w[ii])/sum(abs(Y[ii])**2)
i = abs(w-w_cal).argmin()
delta = angle(Y[i])
print("Calculated value of w0 without noise: ",w_cal)
print("Calculated value of delta without noise: ",delta)

#4 w0 and d determination from a given noised cos(w0t+d) sample vector
y = (cos(w0*t + d) + 0.1*randn(128))*wnd
y[0]=0
y = fftshift(y)
Y = fftshift(fft(y))/128.0
dftplot(y,w,4,r"Spectrum of a noisy $\cos(w_0t+\delta)$ with Hamming window")

# w0 is calculated by finding the weighted average of all w>0. Delta is found by calculating the phase at w closest to w0.
ii = where(w>=0)
w_cal = sum(abs(Y[ii])**2*w[ii])/sum(abs(Y[ii])**2)
i = abs(w-w_cal).argmin()
delta = angle(Y[i])
print("\nCalculated value of w0 with noise: ",w_cal)
print("Calculated value of delta with noise: ",delta)

#5 DFT of cos(16(1.5+t/2pi)t)

t = linspace(-pi,pi,1025); t = t[:-1]
dt = t[1]-t[0]; fmax = 1/dt
n = arange(1024)
wnd = fftshift(0.54+0.46*cos(2*pi*n/1024))
y = cos(16*t*(1.5 + t/(2*pi)))*wnd
y[0]=0
y = fftshift(y)
Y = fftshift(fft(y))/1024.0
w = linspace(-pi*fmax,pi*fmax,1025); w = w[:-1]
dftplot(y,w,100,r"Spectrum of chirped function")

#6 surface plot of 64 samples wide DFTs vs t
t_array = split(t,16)
Y_mag = zeros((16,64))
Y_phase = zeros((16,64))

for i in range(len(t_array)):
	n = arange(64)
	wnd = fftshift(0.54+0.46*cos(2*pi*n/64))
	y = cos(16*t_array[i]*(1.5 + t_array[i]/(2*pi)))*wnd
	y[0]=0
	y = fftshift(y)
	Y = fftshift(fft(y))/64.0
	Y_mag[i] = abs(Y)
	Y_phase[i] = angle(Y)

t = t[::64]	
w = linspace(-fmax*pi,fmax*pi,64+1); w = w[:-1]
t,w = meshgrid(t,w)

fig1 = figure(7)
ax = fig1.add_subplot(111, projection='3d')
surf=ax.plot_surface(w,t,Y_mag.T,cmap=cm.jet,linewidth=0, antialiased=False)
fig1.colorbar(surf, shrink=0.5, aspect=5)
ax.set_title('surface plot');
ylabel(r"$\omega\rightarrow$")
xlabel(r"$t\rightarrow$")

show()


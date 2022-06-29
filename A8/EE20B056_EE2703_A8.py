from pylab import *

#1

#DFT of sin(5t)
x=linspace(0,2*pi,129)
x=x[:-1]
y=sin(5*x)
Y=fftshift(fft(y))/128.0
w=linspace(-64,63,128)
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-10,10])
ylabel(r'$|Y|$',size=16)
title(r'Spectrum of $\sin(5t)$')
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-10,10])
ylabel(r"Phase of $Y$",size=16)
xlabel(r'$k$',size=16)
grid(True)
show()

#DFT of (1+0.1 cos(t)) cos(10t)
x=linspace(-4*pi,4*pi,513)
x=x[:-1]
y=(1+0.1*cos(x))*cos(10*x)
Y=fftshift(fft(y))/512.0
w=linspace(-64,64,513)
w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-15,15])
ylabel(r'$|Y|$',size=16)
title(r'Spectrum of $\left(1+0.1\cos\left(t\right)\right)\cos\left(10t\right)$')
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-15,15])
ylabel(r"Phase of $Y$",size=16)
xlabel(r'$\omega$',size=16)
grid(True)
show()

#DFT of sin^3t and cos^3t
x=linspace(-4*pi,4*pi,513)
x=x[:-1]
y=cos(x)**3
Y=fftshift(fft(y))/512.0
w=linspace(-64,64,513)
w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-15,15])
ylabel(r'$|Y|$',size=16)
title(r'Spectrum of $cos^{3}(t)$')
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-15,15])
ylabel(r"Phase of $Y$",size=16)
xlabel(r'$\omega$',size=16)
grid(True)
show()

x=linspace(-4*pi,4*pi,513)
x=x[:-1]
y=sin(x)**3
Y=fftshift(fft(y))/512.0
w=linspace(-64,64,513)
w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-15,15])
ylabel(r'$|Y|$',size=16)
title(r'Spectrum of $sin^{3}(t)$')
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-15,15])
ylabel(r"Phase of $Y$",size=16)
xlabel(r'$\omega$',size=16)
grid(True)
show()

#DFT of cos(20t +5 cos(t))
x=linspace(-4*pi,4*pi,513)
x=x[:-1]
y=cos(20*x+5*cos(x))
Y=fftshift(fft(y))/512.0
w=linspace(-64,64,513)
w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-40,40])
ylabel(r'$|Y|$',size=16)
title(r'Spectrum of $cos(20t+5cos(t))$')
grid(True)
subplot(2,1,2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-40,40])
ylabel(r"Phase of $Y$",size=16)
xlabel(r'$\omega$',size=16)
grid(True)
show()

#DFT of a gaussian exp(t**2/2)
#FT{exp(t^2/2)}=(1/sqrt(2*pi))*exp(-0.5*w**2)

T = 2*pi
N = 128
iter = 0
tolerance = 1e-6

while True:

	t = linspace(-T/2,T/2,N+1)[:-1]
	w = N/T * linspace(-pi,pi,N+1)[:-1] 
	y = exp(-0.5*t**2)
	iter = iter + 1

	Y = fftshift(fft(y))*T/(2*pi*N)
	Y_actual = (1/sqrt(2*pi))*exp(-0.5*w**2)
	error = mean(abs(abs(Y)-Y_actual))

	if error < tolerance:
		break
	
	T = T*2
	N = N*2

print("Best values of N and T for which the DFT with actual FT of gaussian:")
print("N:"+str(N))
print("T:"+str(T))
print("error:"+str(error))

plot(w,abs(Y),label='DFT')
plot(w,abs(Y_actual),label='FT',alpha=0.5)
title("DFT and FT of gaussian")
legend()
grid(True)
xlim([-10,10])
show()
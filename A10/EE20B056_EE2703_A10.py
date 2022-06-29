import csv
from pylab import *
from scipy import signal

def readsequence(name):
    sequence=[]
    with open(name,'r') as file:
        f=csv.reader(file)
        for row in f:
            sequence.append(float(row[0]))
    return(sequence)

def responseplot(w,h):
    fig, ax1 = plt.subplots()
    ax1.set_title('Digital filter frequency response')
    ax1.plot(w, 20 * np.log10(abs(h)), 'b')
    ax1.set_ylabel('Amplitude [dB]', color='b')
    ax1.set_xlabel('Frequency [rad/sample]')
    ax2 = ax1.twinx()
    angles = np.unwrap(np.angle(h))
    ax2.plot(w, angles, 'g')
    ax2.set_ylabel('Angle (radians)', color='g')
    ax2.grid()
    ax2.axis('tight')
    plt.show()

#1
h=readsequence("h.csv")

#2
wh, magh = signal.freqz(h)
responseplot(wh,magh)

#3
n = linspace(1,2**10,2**10)
x = cos(0.2*pi*n) + cos(0.85*pi*n)

plot(n,x)
xlabel(r'$n\rightarrow$')
ylabel=(r'$x\rightarrow$')
xlim([1,100])
title(r'Plot of $\cos\left(0.2\pi t\right)+\cos\left(0.85\pi t\right)$')
show()

#4 linear convolution of x and h
linconv=zeros(len(h)+len(x)-1)
for N in range(len(h)+len(x)-1):
    for k in range(len(x)):
        if N>k and k>0 and k<len(x) and N-k<len(h):
            linconv[N]+=x[k]*h[N-k]
plot(linconv)
title("Linear Convolution")
show()

#5 circular convolution of x and h
h_=concatenate((h,zeros(len(x)-len(h))))
y=ifft(fft(x)*fft(h_))

plot(y)
xlabel(r'$n\rightarrow$')
ylabel=(r'$y\rightarrow$')
title("Circular convolution")
show()

#6 linear convolution using circular convolution


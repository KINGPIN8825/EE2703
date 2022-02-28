import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy

ll=-2*np.pi
ul=4*np.pi
n=10000
nfourier=51

def f1(x):
    return(np.exp(x))

def f2(x):
    return(np.cos(np.cos(x)))

def plotf1():
    x = np.linspace(ll,ul,n)
    plt.title(r'$e^x$ on a semilogy plot')
    plt.xlabel('x')
    plt.ylabel(r'$log(e^x)$')
    plt.grid(True)
    plt.semilogy(x,f1(x),color='r',label='actual')
    plt.semilogy(x,f1(np.remainder(x,2*np.pi)),color='g',label='fourier')
    plt.legend()
    plt.show()

def plotf2():
    x = np.linspace(ll,ul,n)
    plt.title(r'Plot of $cos(cos(x))$')
    plt.xlabel('x')
    plt.ylabel(r'$cos(cos(x))$')
    plt.grid(True)
    plt.plot(x,f2(x),color='r',label='actual')
    plt.plot(x,f2(np.remainder(x,2*np.pi)),color='g',label='fourier')
    plt.legend()
    plt.show()

def FT(fn,function):
    a = np.zeros(fn)
    def fcos(x,k,f):
        return f(x)*np.cos(k*x)/np.pi
    def fsin(x,k,f):
        return f(x)*np.sin(k*x)/np.pi

    a[0] = integrate.quad(function,0,2*np.pi)[0]/(2*np.pi)
    for i in range(1,fn):
        if(i%2==1):
            a[i] = integrate.quad(fcos,0,2*np.pi,args=(int(i/2)+1,function))[0]
        else:
            a[i] = integrate.quad(fsin,0,2*np.pi,args=(int(i/2),function))[0]
    return a

def fplot(f1FT,f2FT,color = 'ro'):
    f1FT = np.abs(f1FT)
    f2FT = np.abs(f2FT)

    plt.title(r"Coefficients of fourier series of $e^x$ on a semilogy scale")
    plt.xlabel(r'$n$')
    plt.ylabel(r'$log(coeff)$')
    plt.semilogy(f1FT,color)
    plt.grid(True)
    plt.show()

    plt.title(r"Coefficients of fourier series of $e^x$ on a loglog scale")
    plt.xlabel(r'$log(n)$')
    plt.ylabel(r'$log(coeff)$')
    plt.loglog(f1FT,color)
    plt.grid(True)
    plt.show()


    plt.title(r"Coefficients of fourier series of $cos(cos(x))$ on a semilogy scale")
    plt.xlabel(r'$n$')
    plt.ylabel(r'$log(coeff)$')
    plt.semilogy(f2FT,color)
    plt.show()

    plt.title(r"Coefficients of fourier series of $cos(cos(x))$ on a loglog scale")
    plt.xlabel(r'$log(n)$')
    plt.ylabel(r'$log(coeff)$')
    plt.loglog(f2FT,color)
    plt.grid(True)
    plt.show()

plotf1()#1a
plotf2()#1b

f1FT = FT(nfourier,f1)
f2FT = FT(nfourier,f2)
fplot(f1FT,f2FT)#3

#--------------------------------------------------------------------------------------------------

def generateAb(x,f):
    A = np.zeros((x.shape[0],nfourier))
    A[:,0] = 1
    for i in range(1,int((nfourier+1)/2)):
        A[:,2*i-1]=np.cos(i*x)
        A[:,2*i]=np.sin(i*x)
    return A,f(x)

x=np.linspace(0,2*np.pi,401)
x=x[:-1]

Af1,bf1=generateAb(x,f1);cf1=scipy.linalg.lstsq(Af1,bf1)[0]
Af2,bf2=generateAb(x,f2);cf2=scipy.linalg.lstsq(Af2,bf2)[0]

def plotlstsq(f1FT,f2FT,cf1,cf2,color = 'go'):
    f1FT = np.abs(f1FT)
    f2FT = np.abs(f2FT)
    cf1 = np.abs(cf1)
    cf2 = np.abs(cf2)
    plt.title(r"Coefficients of fourier series of $e^x$ on a semilogy scale")
    plt.xlabel(r'$n$')
    plt.ylabel(r'$log(coeff)$')
    plt.semilogy(f1FT,'ro')
    plt.semilogy(cf1,color)
    plt.legend(["true","pred"])
    plt.grid(True)
    plt.show()
    plt.title(r"Coefficients of fourier series of $e^x$ on a loglog scale")
    plt.xlabel(r'$log(n)$')
    plt.ylabel(r'$log(coeff)$')
    plt.loglog(f1FT,'ro')
    plt.semilogy(cf1,color)
    plt.legend(["true","pred"])
    plt.grid(True)
    plt.show()

    plt.title(r"Coefficients of fourier series of $cos(cos(x))$ on a semilogy scale")
    plt.xlabel(r'$n$')
    plt.ylabel(r'$log(coeff)$')
    plt.semilogy(f2FT,'ro')
    plt.semilogy(cf2,color)
    plt.legend(["true","pred"])
    plt.grid(True)
    plt.show()
    plt.title(r"Coefficients of fourier series of $cos(cos(x))$ on a loglog scale")
    plt.xlabel(r'$log(n)$')
    plt.ylabel(r'$log(coeff)$')
    plt.loglog(f2FT,'ro')
    plt.semilogy(cf2,color)
    plt.legend(["true","pred"])
    plt.grid(True)
    plt.show()

plotlstsq(f1FT,f2FT,cf1,cf2,'go')

print("The error in Coefficients of e^x =",np.amax(np.abs(f1FT -cf1)))
print("The error in Coefficients of cos(cos(x)) =",np.amax(np.abs(f2FT -cf2)))

cf1 = np.reshape(cf1,(nfourier,1))
#Finding values of the function from the Coefficients obtained using lstsq
TTT = np.matmul(Af1,cf1)
#plotting results
x = np.linspace(0,2*np.pi,400,endpoint=True)
plt.title(r"Plot of $e^x$")
t = np.linspace(-2*np.pi,4*np.pi,n,endpoint=True)
plt.semilogy(t,f1(t))
plt.semilogy(t,f1(np.remainder(t,2*np.pi)))
plt.semilogy(x,TTT,'go')
plt.xlabel('x')
plt.ylabel(r'$e^x$')
plt.legend(["true","expected","pred"])
plt.grid(True)
plt.show()

cf2 = np.reshape(cf2,(nfourier,1))
#Finding values of the function from the Coefficients obtained using lstsq
TTT = np.matmul(Af2,cf2)
#plotting results
x = np.linspace(0,2*np.pi,400,endpoint=True)
plt.title(r"Plot of $cos(cos(x))$")
t = np.linspace(-2*np.pi,4*np.pi,n,endpoint=True)
plt.plot(x,TTT,'go')
plt.plot(t,f2(t))
plt.xlabel('x')
plt.ylabel(r'$cos(cos(x))$')
plt.legend(["prediction","true"])
plt.grid(True)
plt.show()
from pylab import *

l=0.5           #quarter wavelength
c=2.9979e8      #speed of light
mu0=4e-7*pi     #permeability of free space
Im=1            #current injected into the antenna
a=0.01          #radius of wire
lamda=l*4       #wavelength
f=c/lamda       #frequency
k=2*pi/lamda    #wavenumber


def M(N):   #function to create M matrix  
    return identity(2*N-2)/(2*pi*a)

def R(N):   #function returning distance between points
    t=linspace(0,2*N,2*N+1)
    x,y=meshgrid(t,t)
    Rz=sqrt(((l/N)*(x-y))**2+a**2)    #distance between all sets of points
    Rn=concatenate( ( (Rz[:,N])[1:N], (Rz[:,N])[N+1:2*N] ) )    #distance between Im current element to all points
    Ru=delete(Rz,[0,N,-1],0);Ru=delete(Ru,[0,N,-1],1)   #distance between all sets of points except Im element and end elements
    return(Rz,Ru,Rn)

def consmat(N): #function to CONStruct the MATrices Q and Qb

    Rz,Ru,Rn=R(N)
    dz=l/N

    Pb=mu0*dz/(4*pi)*exp(-1j*k*Rn)/Rn
    Qb=a*(1j*k*divide(Pb,Rn)+divide(Pb,square(Rn)))/mu0

    P=mu0*dz*divide(exp(complex(0,-k)*Ru),Ru)/(4*pi)
    Q=a*(complex(0,k)*divide(P,Ru)+divide(P,square(Ru)))/mu0
    return(Q,Qb)


N=100
Q,Qb=consmat(N)
t=linspace(0,2*N,2*N+1)-(N)*ones(2*N+1)
Iexp=Im*sin(k*(l-(l/N)*abs(t)))
J=(matmul(inv(M(N)-Q),Qb)*Im)   #J vector representing unknown currents
I=zeros(2*N+1)  #'I' vector representing all currents
I[0]=I[2*N]=0;I[N]=Im   #inserting boundary values in 'I' vector
I[1:N],I[N+1:2*N]=real(J[0:N-1]),real(J[N-1:])  #inserting the unknown currents in the 'I' vector

plot(t,I,'r',label='Simulated current')
plot(t,Iexp,'g',label='Expected current')
grid()
title('Simulated current in half-wave dipole antenna for N='+str(N),size=12)
xlabel(r'$Position(x)\rightarrow$',size=10)
ylabel(r'$Current\rightarrow$',size=10)
legend()
show()




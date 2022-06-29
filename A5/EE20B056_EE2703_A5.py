from pylab import *
import sys
import mpl_toolkits.mplot3d.axes3d as p3

#size of plate is 1*1 fixed
Nx=25; # size along x
Ny=25; # size along y
radius=0.35;# radius of central lead
Niter=1500; # number of iterations to perform

if(len(sys.argv)==3):
	Nx = int(sys.argv[1])
	Ny = int(sys.argv[2])
	Niter = int(sys.argv[3])
elif(len(sys.argv)==1):
    print()
else:
    print("Invalid argument input!")
    exit()

x = linspace(-0.5,0.5,Nx)    
y = linspace(0.5,-0.5,Ny)    
X,Y = meshgrid(x,y)            
phi = zeros((Nx,Ny))     
ii = where(X*X + Y*Y <= radius*radius)  
phi[ii] = 1.0                  

errors=zeros(Niter)
for k in range(Niter):
    oldphi=phi.copy()
    phi[1:-1,1:-1] = 0.25*(phi[1:-1,0:-2] + phi[1:-1,2:] + phi[0:-2,1:-1] + phi[2:,1:-1])

    phi[1:-1,0] = phi[1:-1,1]
    phi[1:-1,-1] = phi[1:-1,-2]
    phi[0,1:-1] = phi[1,1:-1]
    
    phi[ii] = 1.0

    errors[k]=(abs(phi-oldphi)).max()

c_approx_500 = lstsq(c_[ones(Niter-500),arange(500,Niter)],log(errors[500:]),rcond=None)
a_500,b_500 = exp(c_approx_500[0][0]),c_approx_500[0][1]

c_approx = lstsq(c_[ones(Niter),arange(Niter)],log(errors),rcond=None)
a, b = exp(c_approx[0][0]), c_approx[0][1]

title('The semilog plots of error function')
semilogy(arange(500,Niter),a_500*exp(b_500*arange(500,Niter)),color='r')
semilogy(arange(0,Niter),a*exp(b*arange(0,Niter)),color='b')
plot(arange(0,Niter,50),errors[::50],'go')
legend(('fit1','fit2','errors'))
show()

fig1 = figure(4)
ax = p3.Axes3D(fig1)
title('The 3-D surface plot of the potential')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('Potential $(\phi)$')
surf = ax.plot_surface(X, Y, phi, rstride=1, cstride=1, cmap=cm.jet,linewidth=0, antialiased=True)
show()

contour(x,y,phi)
plot(x[ii[0]],y[ii[1]],'ro')
xlabel('x')
ylabel('y')
title('The Contour plot of final potential')
grid()
show()

Jx = np.zeros((Nx,Ny))
Jy = np.zeros((Nx,Ny))
Jy[1:-1,1:-1] = 0.5*(phi[1:-1,2:] - phi[1:-1,0:-2])
Jx[1:-1,1:-1] = 0.5*(phi[2:,1:-1] - phi[0:-2,1:-1])
quiver(y,x,Jy[::-1,:],Jx[::-1,:])
contour(x,y,phi)
plot(x[ii[0]],y[ii[1]],'ro')
show()

T = np.zeros((Nx,Ny))
T[:,:] = 300
sigma = 6*(10**7)
kappa = 386
#in the following calculation J is actually dV and hence it is divided with the distance between the grid points and multiplied with sigma
for i in range(Niter):
    T[1:-1,1:-1] = 0.25*(T[1:-1,0:-2] + T[1:-1,2:] + T[0:-2,1:-1] + T[2:,1:-1] + ((((Jx*Nx)**2)[1:-1,1:-1] + ((Jy*Ny)**2)[1:-1,1:-1])*sigma*(10**-4))/kappa)
    T[1:-1,0]=T[1:-1,1]
    T[1:-1,Nx-1]=T[1:-1,Nx-2]
    T[0,1:-1]=T[1,1:-1]
    T[ii] = 300.0
fig1=figure(4)
ax=p3.Axes3D(fig1)
title('The 3-D surface plot of the temperature')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('Temperature')
ax.plot_surface(X, Y, T, rstride=1, cstride=1, cmap=cm.jet,linewidth=0, antialiased=True)
show()
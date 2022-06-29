import scipy.signal as sp
import numpy as np
import matplotlib.pyplot as plt
#1
F1=sp.lti([1,0.5],np.polymul([1,1,2.5],[1,0,2.25]))
t1,x1=sp.impulse(F1,None,np.linspace(0,100,1001))
plt.plot(t1,x1)
plt.xlabel(r'$t\rightarrow$')
plt.ylabel(r'$x(t)\rightarrow$')
plt.show()

#2
F2=sp.lti([1,0.05],np.polymul([1,0.1,2.2525],[1,0,2.25]))
t2,x2=sp.impulse(F2,None,np.linspace(0,100,1001))
plt.plot(t2,x2)
plt.xlabel(r'$t\rightarrow$')
plt.ylabel(r'$x(t)\rightarrow$')
plt.show()

#3
H3=sp.lti([1],[1,0,2.25])
for w in np.linspace(1.4,1.6,5):
    t3=np.linspace(0,50,1001)
    f3=np.cos(w*t3)*np.exp(-0.05*t3)
    t3,x3,svec=sp.lsim(H3,f3,t3)
    plt.plot(t3,x3,label='w='+str(w))
plt.legend()
plt.xlabel(r'$t\rightarrow$')
plt.ylabel(r'$x(t)\rightarrow$')
plt.show()

#4
t4=np.linspace(0,20,201)
x4=sp.lti([1,0,2],[1,0,3,0])
y4=sp.lti([2],[1,0,3,0])
t4,x4 = sp.impulse(x4,None,t4)
t4,y4 = sp.impulse(y4,None,t4)
plt.plot(t4,x4,label='x(t)')
plt.plot(t4,y4,label='y(t)')
plt.xlabel(r'$t\rightarrow$')
plt.ylabel(r'$f(t)\rightarrow$')
plt.legend()
plt.show()

#5
H5=sp.lti([1],[1e-12,1e-4,1])
w,S,phi=H5.bode()
plt.semilogx(w,S,label=r'$20\log|H(j\omega)|$')
plt.semilogx(w,phi,label=r'$\angle H(j\omega)$')
plt.xlabel(r'$\omega(log)\rightarrow$')
plt.legend()
plt.show()

#6
t6=np.linspace(0,10e-3,int(1e5))
v6=np.cos(1e3*t6) - np.cos(1e6*t6)
t6,y6,svec6=sp.lsim(H5,v6,t6)
plt.plot(t6,y6)
plt.xlabel(r'$t\rightarrow$')
plt.ylabel(r'$V_o(t)\rightarrow$')
plt.show()
plt.plot(t6[0:300],y6[0:300])
plt.xlabel(r'$t\rightarrow$')
plt.ylabel(r'$V_o(t)\rightarrow$')
plt.show()
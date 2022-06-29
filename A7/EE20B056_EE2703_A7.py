from sympy import *
import pylab as p
import scipy.signal as sp
def lowpass(R1,R2,C1,C2,G,Vi):
    s=symbols('s')
    A=Matrix([[0,0,1,-1/G],[-1/(1+s*R2*C2),1,0,0],[0,-G,G,1],[-1/R1-1/R2-s*C1,1/R2,0,s*C2]])
    b=Matrix([0,0,0,-Vi/R1])
    V=A.inv()*b
    return(A,b,V)

def highpass(R1,R3,C1,C2,G,Vi):
    s = symbols('s')
    A = Matrix([[0,-1,0,1/G],[s*C2*R3/(s*C2*R3+1),0,-1,0],[0,G,-G,1],[-s*C2-1/R1-s*C1,0,s*C2,1/R1]])
    b = Matrix([0,0,0,-Vi*s*C1])
    V = A.inv()*b
    return A,b,V

def coeffextract(Vo):
    s=symbols('s')
    top, bot = [[float(i) for i in Poly(j, s).all_coeffs()] for j in Vo.as_numer_denom()]
    return(sp.lti(top,bot))

s=symbols('s')
ww=p.logspace(0,8,801)
ss=1j*ww


A,b,V=lowpass(10000,10000,1e-9,1e-9,1.586,1)
Vo=V[3]
H1=coeffextract(Vo)
h1func=lambdify(s,Vo,'numpy')
p.loglog(ww,abs(h1func(ss)),lw=2)
p.xlabel(r'$\omega(log)\rightarrow$')
p.title('Plot of magnitude response of low-pass filter')
p.grid(True)
p.show()

A,b,V=highpass(10000,10000,1e-9,1e-9,1.586,1)
Vo=V[3]
H2=coeffextract(Vo)
h2func=lambdify(s,Vo,'numpy')
p.loglog(ww,abs(h2func(ss)),lw=2)
p.xlabel(r'$\omega(log)\rightarrow$')
p.title('Plot of magnitude response of high-pass filter')
p.grid(True)
p.show()


#1
t=p.linspace(0,1e-3,1001)
t,y1=sp.step(H1,None,t)
p.title('Plot of unit step response of low-pass filter')
p.plot(t,y1,label='Vo')
p.plot(t,p.ones(len(t)),label='Vi')
p.xlabel(r'$t\rightarrow$')
p.legend()
p.show()

#2
t=p.linspace(0,1e-2,100001)
Vi=p.sin(2000*p.pi*t) + p.cos(2e6*p.pi*t)
t,y2,svec = sp.lsim(H1,Vi,t)
p.title('Plot of low-pass filter response to multiple sinusoid')
p.plot(t,y2,label='Vo',color='r')
p.plot(t,Vi,label='Vi',color='b',alpha=0.5)
p.xlabel(r'$t\rightarrow$')
p.legend()
p.show()

damp=1000
#4
t=p.linspace(0,1e-2,1001)
Vi=p.sin(1e3*p.pi*t)*p.exp(-damp*t)
t,y3,svec=sp.lsim(H2,Vi,t)
p.title('Plot of high-pass filter response to low frequency damped sinusoid')
p.plot(t,y3,label='Vo')
p.plot(t,Vi,label='Vi')
p.xlabel(r'$t\rightarrow$')
p.legend()
p.show()


t=p.linspace(0,1e-5,1001)
Vi=p.sin(1e6*p.pi*t)*p.exp(-damp*t)
t,y4,svec=sp.lsim(H2,Vi,t)
p.title('Plot of high-pass filter response to high frequency damped sinusoid')
p.plot(t,y4,label='Vo')
p.plot(t,Vi,label='Vi')
p.xlabel(r'$t\rightarrow$')
p.legend()
p.show()

#5
t=p.linspace(0,1e-3,1001)
t,y1=sp.step(H2,None,t)
p.title('Plot of unit step response of high-pass filter')
p.plot(t,y1,label='Vo')
p.plot(t,p.ones(len(t)),label='Vi')
p.xlabel(r'$t\rightarrow$')
p.legend()
p.show()
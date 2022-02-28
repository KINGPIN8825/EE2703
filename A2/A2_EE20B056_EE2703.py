import numpy as np
from sys import argv
import math as m

def f(self, freq):
    return float(self.tokens[4]) if freq==0 else complex(float(self.tokens[4])*m.cos(float(self.tokens[5])),float(self.tokens[4])*m.sin(float(self.tokens[5])))
class component():
    def __init__(self,line) -> None:
        self.line=line
        self.tokens=line.split()
        self.from_node = self.tokens[1]
        self.to_node = self.tokens[2]
    def value(self,f=0) -> complex:
        if self.tokens[0][0]=='R':
            return float(self.tokens[3])
        elif self.tokens[0][0]=='L':
            return complex(0,2*np.pi*f*float(self.tokens[3]))
        elif self.tokens[0][0]=='C':
            return complex(0,(-1)/(2*np.pi*f*float(self.tokens[3])))
        elif self.tokens[0][0]=='V' or self.tokens[0][0]=='I':
            if f==0:
                return float(self.tokens[4])
            else:
                return complex(float(self.tokens[4])*m.cos(float(self.tokens[5])),float(self.tokens[4])*m.sin(float(self.tokens[5])))
        elif self.tokens[0][0]=='E' or self.tokens[0][0]=='G' or self.tokens[0][0]=='H' or self.tokens[0][0]=='F':
            return float(self.tokens[5])
    def type(self)-> str:
        if self.tokens[0][0]=='V' or self.tokens[0][0]=='I':
            return 'a'
        else:
            return 'p'

if len(argv)!=2:	#CHECKING IF THE NUMBER OF ARGUMENTS SUPPLIED IS CORRECT
	print('\nUsage :python3 %s <input-file>'% argv[0])
	exit()

CIRCUIT='.circuit'
END = '.end'
AC='.ac'


#------------------------------------------------------------------------------------

try:
    with open(argv[1]) as f:
        lines=f.readlines()
        netlist=[]
        start=-1;end=0;srcdict=dict()#dictionary when put key as the source gives its frequency as output
        for line in lines:
            line=line.split('#')[0]
            if CIRCUIT==line[:len(CIRCUIT)]:
                start=lines.index(line)
            if END==line[:len(END)]:
                end=lines.index(line)
            if AC==line[:len(AC)]:
                srcdict[line.split()[1]]=line.split()[2]
            netlist.append(line)
        lines=netlist
        netlist=[]
        for line in lines:#changing 'dc' in dc sources with 0 and 'ac' in ac sources with its frequency
            temp=line.split()
            if len(temp)>=4:
                if temp[3]=='dc':
                    temp[3]='0'
                    srcdict[temp[0]]='0'
                elif temp[3]=='ac':
                    temp[3]=srcdict[temp[0]]
            line=' '.join(temp)
            netlist.append(line)
        if end<start:
            print("Invalid circuit definition")
            exit()
        else:
            complist=[]
            for i in range(start+1,end):
                complist.append(component(netlist[i]))#list of class instances holding the definition of each component(i.e., each line)

            for i in range(len(complist)):
                if len(complist[i].tokens[-1].split('e'))>1:# converting string aeb to the float a*10^b
                    complist[i].tokens[-1]=str(float(complist[i].tokens[-1].split('e')[0])*10**float(complist[i].tokens[-1].split('e')[1]))
                print(complist[i].tokens)
#-------------------------------------------------------------------------------------------------------------------
#the part between comments works perfectly fine, now only work on 'complist' list
 
#MX=b ;  X construction
        Xdict=dict()
        k=0 
        for comp in complist:
            for i in range(2):
                if comp.tokens[i+1] not in Xdict:
                    Xdict[comp.tokens[i+1]]=k
                    k=k+1
        for comp in complist:
            if comp.tokens[0][0]=='V' or comp.tokens[0][0]=='E' or comp.tokens[0][0]=='H':
                if 'I'+comp.tokens[0] not in Xdict:
                    Xdict['I'+comp.tokens[0]]=k
                    k=k+1
        print(Xdict)
#MX=b ; M and b construction
        Xf=dict()

        for src in complist:
            if src.type()=='a':
                #print(src.tokens)
                M=np.zeros((len(Xdict),len(Xdict)),dtype='complex_')
                b=np.zeros(len(Xdict),dtype='complex_')
                if src.tokens[0][0]=='I' :
                    b[Xdict[src.tokens[1]]]-=src.value(float(src.tokens[3]))
                    b[Xdict[src.tokens[2]]]+=src.value(float(src.tokens[3]))
                elif src.tokens[0][0]=='V':
                    M[Xdict[src.tokens[1]],Xdict['I'+src.tokens[0]]]+=1
                    M[Xdict[src.tokens[2]],Xdict['I'+src.tokens[0]]]-=1
                    M[Xdict['I'+src.tokens[0]],Xdict[src.tokens[1]]]+=1
                    M[Xdict['I'+src.tokens[0]],Xdict[src.tokens[2]]]-=1
                    #print(src.tokens[3])
                    b[Xdict['I'+src.tokens[0]]]+=src.value(float(src.tokens[3]))
                    #print(src.value(float(src.tokens[3])))
                    #print(b)
                for comp in complist:
                    if comp.type()=='p':
                        if comp.tokens[0][0]=='R' or comp.tokens[0][0]=='L' or comp.tokens[0][0]=='C':
                            M[Xdict[comp.tokens[1]],Xdict[comp.tokens[1]]]+=1/comp.value()
                            M[Xdict[comp.tokens[1]],Xdict[comp.tokens[2]]]+=(-1)/comp.value()
                            M[Xdict[comp.tokens[2]],Xdict[comp.tokens[1]]]+=(-1)/comp.value()
                            M[Xdict[comp.tokens[2]],Xdict[comp.tokens[2]]]+=1/comp.value()
                            print(comp.tokens[0])
                            print(M)
                        elif comp.tokens[0][0]=='G':
                            M[Xdict[comp.tokens[1]],Xdict[comp.tokens[3]]]+=float(comp.tokens[5])
                            M[Xdict[comp.tokens[1]],Xdict[comp.tokens[4]]]-=float(comp.tokens[5])
                            M[Xdict[comp.tokens[2]],Xdict[comp.tokens[3]]]-=float(comp.tokens[5])
                            M[Xdict[comp.tokens[2]],Xdict[comp.tokens[4]]]+=float(comp.tokens[5])
                        elif comp.tokens[0][0]=='E':
                            M[Xdict[comp.tokens[1]],Xdict['I'+comp.tokens[0]]]+=1
                            M[Xdict[comp.tokens[2]],Xdict['I'+comp.tokens[0]]]-=1
                            M[Xdict['I'+comp.tokens[0]],Xdict[comp.tokens[1]]]+=1
                            M[Xdict['I'+comp.tokens[0]],Xdict[comp.tokens[2]]]-=1
                            M[Xdict['I'+comp.tokens[0]],Xdict[comp.tokens[3]]]-=float(comp.tokens[5])
                            M[Xdict['I'+comp.tokens[0]],Xdict[comp.tokens[4]]]+=float(comp.tokens[5])
                        elif comp.tokens[0][0]=='F':
                            M[Xdict[comp.tokens[1]],Xdict['I'+comp.tokens[3]]]+=float(comp.tokens[3])
                            M[Xdict[comp.tokens[2]],Xdict['I'+comp.tokens[3]]]-=float(comp.tokens[3])
                            for comp1 in complist:
                                if comp1.tokens[0]==comp.tokens[3]:
                                    M[Xdict[comp1.tokens[1]],Xdict['I'+comp.tokens[3]]]+=1
                                    M[Xdict[comp1.tokens[2]],Xdict['I'+comp.tokens[3]]]-=1
                                    M[Xdict['I'+comp.tokens[3]],Xdict[comp1.tokens[1]]]+=1
                                    M[Xdict['I'+comp.tokens[3]],Xdict[comp1.tokens[2]]]-=1
                        elif comp.tokens[0][0]=='H':
                            for comp1 in complist:
                                if comp1.tokens[0]==comp.tokens[3]:
                                    M[Xdict[comp.tokens[1]],Xdict['I'+comp.tokens[3]]]+=1
                                    M[Xdict[comp.tokens[2]],Xdict['I'+comp.tokens[3]]]-=1
                                    M[Xdict[comp1.tokens[1]],Xdict['I'+comp1.tokens[3]]]+=1
                                    M[Xdict[comp1.tokens[2]],Xdict['I'+comp1.tokens[3]]]-=1
                                    M[Xdict['I'+comp1.tokens[3]],Xdict[comp1.tokens[1]]]+=1
                                    M[Xdict['I'+comp1.tokens[3]],Xdict[comp1.tokens[2]]]-=1
                                    M[Xdict['I'+comp.tokens[3]],Xdict[comp.tokens[1]]]+=1
                                    M[Xdict['I'+comp.tokens[3]],Xdict[comp.tokens[2]]]-=1
                                    M[Xdict['I'+comp.tokens[3]],Xdict[comp1.tokens[3]]]-=float(comp.tokens[3])
                M=np.delete(M,Xdict['GND'],0)
                M=np.delete(M,Xdict['GND'],1)
                b=np.delete(b,Xdict['GND'],0)
                print(M)
                print(b)
                Xf[src.tokens[0]]=np.linalg.solve(M,b)
        print(Xf)
except FileNotFoundError:
    print("Invalid File.")
    exit()
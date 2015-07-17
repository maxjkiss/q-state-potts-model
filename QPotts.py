from numpy import *
import numpy.random
from matplotlib import pyplot as plt

q = 4

def initializeSpin(L):
    C = random.randint(q,size=(L,L))
    return C

def kronecker(C1, C2):
    if(C1 == C2):
        return 1
    else:
        return 0

def initializeEnergy(C):
    E,L = 0,len(C)
    for i in range(L):
        for j in range(L):
            Eright  = kronecker(C[i,j], C[i,(j+1)%L])
            Ebottom = kronecker(C[i,j], C[(i+1)%L,j])
            E -= Eright+Ebottom
    return E

def calcM(C):
    L,m = len(C),0
    for i in range(L):
        for j in range(L):
            if(C[1,1] == C[i,j]):
                m += 1
    m /= (L**2)
    return abs(m)

def calcEDiff(C,(i,j)):
    L = len(C)

    Eright1, Eleft1 = kronecker(C[i,j], C[i, (j+1)%L]), kronecker(C[i,j], C[i, (j-1)%L])
    Etop1, Ebot1 = kronecker(C[i,j], C[(i-1)%L, j]), kronecker(C[i,j], C[(i+1)%L, j])
    Sum1 = -(Eright1+Eleft1+Etop1+Ebot1)

    Crand = random.randint(q)

    Eright2, Eleft2 = kronecker(Crand, C[i, (j+1)%L]), kronecker(Crand, C[i, (j-1)%L])
    Etop2, Ebot2 = kronecker(Crand, C[(i-1)%L, j]), kronecker(Crand, C[(i+1)%L, j])
    Sum2 = -(Eright2+Eleft2+Etop2+Ebot2)

    return Sum2-Sum1,Crand

def update(C,beta):
    L = len(C)
    result = 0
    for k in range(L*L):
        i,j = random.randint(L),random.randint(L)
        ediff,Crand = calcEDiff(C,(i,j))
        if(random.random() < exp(-beta*ediff)):
            C[i,j] = Crand
            result += ediff
    return result

def main(L):
    L, N, Ntherm,index = 8, 10000, 1000,0
    Tc = 2./(log(1+sqrt(2)))
    Marr, Tarr, Earr = [],[],[]
        
    for T in linspace(0.1,2,num = 50):
        beta = 1./T
        
        C = initializeSpin(L)
        E = initializeEnergy(C)
        
        ESum = 0
        magDensity = 0
            
        for i in range(Ntherm):
            E += update(C,beta)
                
        for i in range(N):
            E += update(C,beta)
            ESum += E/(L*L)
            m = calcM(C)
            magDensity += m
                
        Marr.append(1.*magDensity/N), Tarr.append(1.*T)
        index += 1
    return Marr, Tarr

M4, T4 = main(4)
M6, T6 = main(6)
M8, T8 = main(8)

plt.legend(['L=4','L=6','L=8'])
plt.plot(T4,M4)
plt.plot(T6,M6)
plt.plot(T8,M8)
plt.xlabel('T')
plt.ylabel('M')
plt.title('Magnetization for Various Lattices for q = 4')
plt.show()
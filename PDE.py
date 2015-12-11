import numpy as np
from scipy.linalg import solve
import matplotlib.pyplot as plt

def f(x):
    y=2*np.exp(x) + x*np.exp(x)
    return y
    
for n in range(1,6):
    h=2**-n
    x = np.arange(0,1+h,h)
    N = len(x)
    al=0
    be=np.exp(1)

    P=np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            if i==j:
                P[i,j]=-2
            elif i+1==j:
                P[i,j]=1
            elif i==j+1:
                P[i,j]=1

    A=P/(h**2)
    b=np.zeros((N,1))
    b[0]=f(x[0])-(al/(h**2))
    b[N-1]=f(x[N-1])-(be/(h**2))
    for i in range(1,N-1):
        b[i]=f(x[i])
    #X=np.zeros((5,N))
    X =solve(A,b)
    y=np.zeros(N)
    for i in range (N):
        y[i]=x[i]*np.exp(x[i])

    plt.plot(x,X,label='n='+str(n))
    plt.plot(x,y,label='real',linestyle='--')
    
plt.legend()
plt.show()

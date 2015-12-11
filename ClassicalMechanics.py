import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.figure

#Exercise 1
h=10**-1
t = np.arange(0, 10+h,h)
N = len(t)

k1 = 0.5
k2 = 3
m1 = 2
m2 = 3
a = k2/float(m1)
b = (k1+k2)/float(m1)

def f(s,t):
    dsdt = np.zeros(len(s))
    dsdt[0] = (a * s[3]) - (b * s[1])
    dsdt[1] = s[0]
    dsdt[2] = s[1] - s[3]
    dsdt[3] = s[2]
    return dsdt
    

X = np.zeros((N,4))
X[0,:] = (0,0,0,3)
Xode = odeint(f, X[0,:], t)

E = np.zeros(N)
for i in range(N):
    E[i] = 0.5 * ((m1*(Xode[i,0])**2) + (m2*(Xode[i,2])**2) + (k1*(Xode[i,1])**2) + (k2*(Xode[i,3]-Xode[i,1])**2))

fig = plt.figure()
fig.suptitle("Ex. 1: Motion for a system of two springs")
plt.plot(t, Xode[:,1], label="$x_1(t)$")
plt.plot(t, Xode[:,3], label="$x_2(t)$")
plt.plot(t, E, label="$E(t)$")
plt.legend()
#We can see that energy is conserved is E(t) remains constant every time.

#Exercise 2

g = 9.8

def ftab(s,t):
    dsdt = np.zeros(len(s))
    dsdt[0] = s[2]
    dsdt[1] = s[3]
    dsdt[2] = -2 *(s[2] * s[3])/ s[1]
    dsdt[3] = (0.75* s[1] * (s[2]**2)) - (0.25*g) 
    return dsdt
    

Y = np.zeros((N,4))
Y[0,:] = (0,4,1,0)
Yode = odeint(ftab, Y[0,:], t)

Z = np.zeros((N,4))
Z[0,:] = (0,4,np.sqrt(5/6.),0)
Zode = odeint(ftab, Z[0,:], t)


fig = plt.figure()
fig.suptitle("Ex. 3f: Motion for $\dot{\phi}_0=1$")
plt.plot(t, Yode[:,1], label="$r(t)$ for $\dot{\phi}_0=1$")
plt.plot(t, Yode[:,0], label="$\phi(t)$ for $\dot{\phi}_0=1$")
plt.legend()

fig = plt.figure()
fig.suptitle("Ex. 3g: Motion for $\dot{\phi}_0=\sqrt{5/6}$")
plt.plot(t, Zode[:,1], label="$r(t)$ for $\dot{\phi}_0=\sqrt{5/6}$")
plt.plot(t, Zode[:,0], label="$\phi(t)$ for $\dot{\phi}_0=\sqrt{5/6}$")
plt.legend()
plt.show()

#We notice that there is minimal displacement per time when initial velocity in the phi direction is sqrt(5/6) than when it is 1. 
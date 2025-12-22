import matplotlib.pyplot as plt
import numpy as np 
import matplotlib.animation as animation
import IPython

ur=2.;ul=-2.

U0  = lambda x: ul*(x<0.5) + ur*(x>=0.5)
#U0 = lambda x: 1*(x<0.6)*(x>0.4)
#U0 = lambda x: np.exp(-(x-0.5)**2 /(2*0.01)) 
f  = lambda u: u**3
f_p= lambda u: 3*u**2

def g(u,v):
    if(u<v): return min(f(u),f(v))
    else :   return max(f(u),f(v))

N=1000; T=1; CFL=1.9


dx = 1./N
X= np.linspace(0,1,N)
FG = np.empty_like(X) 
FD = np.empty_like(X)
U  = np.empty_like(X)

Film=[]; F=[[],[]]

U = U0(X)
t=0.; n =0
dt=0.

while t<T and n<10000:
    
    c=0
    
    #calcul de dt
    vitesse = max(abs(f_p(U)))
    if(vitesse !=0) : dt = (min(CFL * dx /(2* vitesse), T-t))
    else :            dt = CFL*dx  
    
    if(dt<0): dt = 10e-5

    ## Calcul du flux
    for j in range(len(X)):
        
        if(j==0):   FG[j]= g(U[j],U[j+1])
        else :      FG[j]= g(U[j-1],U[j])
        
        if(j==N-1): FD[j]= g(U[j-1],U[j])
        else :      FD[j]= g(U[j],U[j+1]) 
        
        # if(j==0):     
        #     FD[j] = 0.5*(f(U[j] + f(U[j+1]))  +0.5*c*(U[j] - U[j+1]))
        #     FG[j] = 0.5*(f(U[-1] + f(U[j]))  +0.5*c*(U[-1] - U[j])) 
        # elif(j==N-1):     
        #     FD[j] = 0.5*(f(U[j] + f(U[0]))  +0.5*c*(U[j] - U[0]))
        #     FG[j] = 0.5*(f(U[j-1] + f(U[j]))  +0.5*c*(U[j-1] - U[j])) 
        # else : 
        #     FD[j] = 0.5*(f(U[j] + f(U[j+1]))  +0.5*c*(U[j] - U[j+1]))
        #     FG[j] = 0.5*(f(U[j-1] + f(U[j]))  +0.5*c*(U[j-1] - U[j])) 
                    
    ## Calcul de la solution
    #U = U - dt/dx*(FD-FG)
    
    for j in range(len(X)):
        if((j==0 or j==1) or ((j==N-2 or j==N-1))):
            print(j,U[j],-(dt/dx)* (FD[j]-FG[j]))
        U[j] = U[j] -(dt/dx)* (FD[j]-FG[j])
        
    print(dt)
    print('----------',t,'------------')
    
    Film.append(U.copy()); F[0].append(FD.copy()); F[1].append(FG.copy())
    n+=1; t+= dt

print(t,'bloop')

U_xt = np.zeros((len(Film),len(U)))
for i in range(len(Film)) : U_xt[i,:] = Film[-i]

fig, ax = plt.subplots()
ax.imshow(U_xt,extent=[0,1,0,1])
plt.show()
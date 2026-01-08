import matplotlib.pyplot as plt
import numpy as np 

ul=-1.;ur=1.
a=0.2;b=1

#u= lambda x,t: (a*x+b)/(a*t+1)
#U0 = lambda x: a*x+b

U0  = lambda x: ul*(x<0.5) + ur*(x>=0.5)
#U0 = lambda x: 1*(x<0.6)*(x>0.4)
#U0 = lambda x: np.exp(-(x-0.5)**2 /(2*0.01)) 

#U0 = lambda x: np.cos(2*np.pi*x)

f  = lambda u: 0.5* u**2
f_p= lambda u: u

N=1000; T=1; CFL=0.9


dx = 1./N
X= np.linspace(0,1,N)
FG = np.empty_like(X) 
FD = np.empty_like(X)
U  = np.empty_like(X)

Film=[]; F=[[],[]]

U = U0(X)

def g(u,v):
    if(u<v):    return min(f(u),f(v))
    else :      return max(f(u),f(v))

    

# g = lambda u,v: 0*u + 0*v   

flux_methode="Godunov"; condition_limite="Neumann"; ghost_nodes=[0,0]; film_bool=True

# if(flux_methode == "Godunov"): 
#     g= lambda u,v : min(f(u),f(v))*(u<v) + max(f(u),f(v))*(u>=v)
# if(flux_methode == "Rusanov"):
#     g= lambda u,v : 0.5*(f(u)+f(v))-0.5* max(np.abs(f_p(u)), np.abs(f_p(v))) *(v-u)
# if(flux_methode == "Roe"):
#     g= lambda u,v : 0.5*(f(u)+f(v))-0.5*f_p(u)*(u-v)*(u!=v) + f(u)*(u==v)

FG = np.empty_like(X) 
FD = np.empty_like(X)

Film=[]; F=[[],[]]  

t=0.; n =0; dt=0.
while t<T and n<10000:
    #calcul de dt
    vitesse = max(abs(f_p(U)))
    if(vitesse !=0) : dt = (min(CFL * dx /(2* vitesse), T-t))
    else :            dt = CFL*dx  

    if(dt<0): dt = 10e-12; print("\n !!!! dt<0 !!! \n")

    ## Calcul du flux
    for j in range(len(X)):
        
        if(j==0):   
            if(condition_limite =="Dirichlet"):    FG[j] = g(ghost_nodes[0],       U[j])
            elif(condition_limite =="Neumann"):    FG[j] = g(U[j]+ghost_nodes[0],U[j])
            elif(condition_limite =="Periodique"): FG[j] = g(U[-1], U[j])

        else :      FG[j]= g(U[j-1],U[j])
        
        if(j==N-1):   
            if(condition_limite =="Dirichlet"):    FD[j] = g(U[j], ghost_nodes[1])
            elif(condition_limite =="Neumann"):    FD[j] = g(U[j], U[j]+ghost_nodes[1])
            elif(condition_limite =="Periodique"): FD[j] = g(U[j], U[0])

        else :      FD[j]= g(U[j],U[j+1]) 

        
    ## eval
    if(n%50==0):
        plt.plot(X,U,'b')
        #plt.plot(X,u(X,t),'r')
        #plt.plot(X,FD,'r')
        #plt.plot(X,FG,'g')
        plt.title('temps = '+str(t))
        plt.show()
        print(sum(U*U/2))
        print(dt)
        
    ## Calcul de la solution
    for j in range(len(X)):
        U[j] = U[j] -(dt/dx)* (FD[j]-FG[j])

    if(film_bool): Film.append(U.copy()); F[0].append(FD.copy()); F[1].append(FG.copy())
    n+=1; t+= dt
    

if(film_bool): 
    U_xt = np.zeros((len(Film),len(U)))
    for i in range(len(Film)) : U_xt[i,:] = Film[-i]

print(t,'bloop')

U_xt = np.zeros((len(Film),len(U)))
for i in range(len(Film)) : U_xt[i,:] = Film[-i]

fig, ax = plt.subplots()
im = ax.imshow(U_xt,extent=((0,1,0,1)))

X = np.linspace(0,1,N); I = np.linspace(0, 1, 50)

for i in range(50): plt.plot(X*f_p(U0(I[i]))+I[i],X,'--r')

plt.xlim(0,1)

fig.colorbar(im, ax=ax, label='Interactive colorbar')
plt.show()
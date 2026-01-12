import matplotlib.pyplot as plt
import numpy as np

function = lambda u: 0.5*u**2
derivee_de_la_function = lambda u: u

ul=-1.;ur=1.
u0  = lambda x: ul*(x<0.5) + ur*(x>=0.5)

J=500; T=1; CFL = 0.9;
condition_limite = "Neumann"; ghost_nodes =[0,0]; film_bool= True


dx = 1./J
X= np.linspace(0,1,J)
FG_u = np.empty_like(X);        FG_v = np.empty_like(X)      
FD_u = np.empty_like(X);        FD_v = np.empty_like(X)

Film_u=[]; Film_v = []; F=[[],[]]

U = u0(X); V = u0(X)


def g_godunov(u,v):
    U = np.linspace(u, v, J)
    if(u<v):    return np.min(function(U))
    else :      return np.max(function(U))

g_rusanov   = lambda u,v : 0.5*(function(u)+function(v))-0.5*max(np.abs(derivee_de_la_function(u)), np.abs(derivee_de_la_function(v)))*(v-u)

t=0.; n =0; dt=0.
while t<T and n<10000:
    #calcul de dt
    vitesse_u = max(abs(derivee_de_la_function(U)));            vitesse_v = max(abs(derivee_de_la_function(V)))
    vitesse = max(vitesse_u,vitesse_v)
    if(vitesse !=0) : dt = (min(CFL * dx /(2* vitesse), T-t))
    else :            dt = CFL*dx  

    if(dt<0): dt = 10e-12; print("\n !!!! dt<0 !!! \n")

    ## Calcul du flux
    for j in range(len(X)):
        
        if(j==0):   
            if(condition_limite =="Dirichlet"):    FG_u[j] = g_godunov(ghost_nodes[0],       U[j])
            elif(condition_limite =="Neumann"):    FG_u[j] = g_godunov(U[j]+ghost_nodes[0],U[j])
            elif(condition_limite =="Periodique"): FG_u[j] = g_godunov(U[-1], U[j])

            if(condition_limite =="Dirichlet"):    FG_v[j] = g_rusanov(ghost_nodes[0],       V[j])
            elif(condition_limite =="Neumann"):    FG_v[j] = g_rusanov(V[j]+ghost_nodes[0],V[j])
            elif(condition_limite =="Periodique"): FG_v[j] = g_rusanov(V[-1], V[j])

        else :      FG_u[j]= g_godunov(U[j-1],U[j]);    FG_v[j]= g_rusanov(V[j-1],V[j])
        
        if(j==J-1):   
            if(condition_limite =="Dirichlet"):    FD_u[j] = g_godunov(U[j], ghost_nodes[1])
            elif(condition_limite =="Neumann"):    FD_u[j] = g_godunov(U[j], U[j]+ghost_nodes[1])
            elif(condition_limite =="Periodique"): FD_u[j] = g_godunov(U[j], U[0])

            if(condition_limite =="Dirichlet"):    FD_v[j] = g_rusanov(V[j], ghost_nodes[1])
            elif(condition_limite =="Neumann"):    FD_v[j] = g_rusanov(V[j], V[j]+ghost_nodes[1])
            elif(condition_limite =="Periodique"): FD_v[j] = g_rusanov(V[j], V[0])

        else :      FD_u[j]= g_godunov(U[j],U[j+1]);    FD_v[j] =g_rusanov(V[j],V[j+1])

    ## eval
    if(n%1==0):
        plt.plot(X,U-V,'r')
        plt.plot(X,U,'b')
        plt.plot(X,V,'g')
        #plt.plot(X,u(X,t),'r')
        #plt.plot(X,FD,'r')
        #plt.plot(X,FG,'g')
        plt.title('temps = '+str(t))
        plt.show()
        print(dt)
        
    ## Calcul de la solution
    for j in range(len(X)):
        U[j] = U[j] -(dt/dx)* (FD_u[j]-FG_u[j])
        V[j] = V[j] -(dt/dx)* (FD_v[j]-FG_v[j])

    if(film_bool): Film_u.append(U.copy());Film_v.append(V.copy())
    n+=1; t+= dt
    

print(t,'bloop')    

if(film_bool): 
    U_xt = np.zeros((len(Film_u),len(U)))
    V_xt = np.zeros((len(Film_v),len(V)))
    for i in range(len(Film_u)) : U_xt[i,:] = Film_u[-i]
    for i in range(len(Film_v)) : V_xt[i,:] = Film_v[-i]
    
fig, ax = plt.subplots()
im = ax.imshow(V_xt-U_xt,extent=((0,1,0,1)))
fig.colorbar(im, ax=ax, label='Interactive colorbar')
plt.show()
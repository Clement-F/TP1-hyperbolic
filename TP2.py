import numpy as np
from scipy import linalg


gamma = 1.4

prim_to_cons = lambda rho,u,p:      (rho, rho*u,rho*(0.5*u**2 +p/((gamma-1)*rho))       )
cons_to_prim = lambda rho,rhou,rhoE:(rho, rhou/rho, (gamma-1)*(rhoE-0.5*(rhou**2)/rho)) 

f = lambda rho,u,p : 0*rho+ 0*u + 0*p
u = lambda rho,u,p : u
c = lambda rho,u,p : np.sqrt(gamma * p/rho)

l1 = lambda rho,u,p: u(rho,u,p)+c(rho,u,p)
l2 = lambda rho,u,p: u(rho,u,p)
l3 = lambda rho,u,p: u(rho,u,p)-c(rho,u,p)

u0 = lambda x: (rho1,u1,p1) if x<0.5 

J = 100; dx= 1/J
X= np.linspace(dx,1,J)-dx/2

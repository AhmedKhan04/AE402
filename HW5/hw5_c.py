import sympy as sp 
import matplotlib.pyplot as plt 
# Dimensional values, 

M_sun = 1.989e30  
M_earth = 5.972e24
M_moon = 7.347e22
M_jupiter = 1.898e27
M_europa = 4.799e22
M_apophis = 2.699e10

d_sun_apophis = 149e6 
d_sun_jupiter = 750e6 
d_earth_moon = 384e3 
d_jupiter_europa = 670e3 

# CD3RP Set UP 

mu_1 = M_apophis/(M_apophis + M_sun)

m1 = 1 - mu_1
m2 = mu_1 
R = 1
G = 1 

# Solve for the sun-apophis system L points 

x = sp.symbols('x')
L1_EQ = sp.Eq(x - (1- mu_1)/(mu_1 + x)**2 + (mu_1)/(x - 1 + mu_1)**2, 0) 
L1_x = sp.nsolve(L1_EQ, x, 0.9, prec=50)
L1_y = 0 

Sun_x = -mu_1
Sun_y = 0 

Apophis_x = 1-mu_1
Apophis_y = 0


DU = d_sun_apophis
L1 = sp.Matrix([sp.N(L1_x), sp.N(L1_y)]) * DU 
Apophis = sp.Matrix([sp.N(Apophis_x), sp.N(Apophis_y)]) * DU
print(f"L1 : {(L1-Apophis).norm()} km")


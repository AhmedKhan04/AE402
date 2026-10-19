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

mu_1 = M_europa/(M_europa + M_jupiter)

m1 = 1 - mu_1
m2 = mu_1 
R = 1
G = 1 

# Solve for the Jupiter-Europa system L points 

x = sp.symbols('x')
L1_EQ = sp.Eq(x - (1- mu_1)/(mu_1 + x)**2 + (mu_1)/(x - 1 + mu_1)**2, 0) 
L1_x = list(sp.solveset(L1_EQ, x, domain=sp.S.Reals))[0]
L1_y = 0 



L2_EQ = sp.Eq(x - (1-mu_1)/(mu_1 + x)**2 - (mu_1)/(x - 1 + mu_1)**2, 0) 
L2_x = list(sp.solveset(L2_EQ, x, domain=sp.S.Reals))[0]
L2_y = 0

L3_EQ = sp.Eq(x + (1-mu_1)/(mu_1 + x)**2 + (mu_1)/(x - 1 + mu_1)**2, 0) 
L3_x = list(sp.solveset(L3_EQ, x, domain=sp.S.Reals))[0]
L3_y = 0 

L4_x = 0.5 - mu_1
L4_y = sp.sqrt(3)/2

L5_x = 0.5 - mu_1
L5_y = -sp.sqrt(3)/2

Jupiter_x = -mu_1
Jupiter_y = 0

Europa_x = 1-mu_1
Europa_y = 0 


DU = d_jupiter_europa
L1 = sp.Matrix([sp.N(L1_x), sp.N(L1_y)]) * DU 
Europa = sp.Matrix([sp.N(Europa_x), sp.N(Europa_y)]) * DU
print(f"L1 : {(L1-Europa).norm()} km")


L2 = sp.Matrix([sp.N(L2_x), sp.N(L2_y)]) * DU 
print(f"L2 : {(L2-Europa).norm()} km")

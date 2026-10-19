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

mu_1 = M_jupiter/(M_jupiter + M_sun)
m1 = 1 - mu_1
m2 = mu_1 
R = 1
G = 1 

# Solve for the sun-jupiter system L points 

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

Sun_x = -mu_1
Sun_y = 0 

Jupiter_x = 1-mu_1
Jupiter_y = 0

plt.figure()

plt.scatter(L1_x, L1_y, label=f"L1: x = {float(L1_x):.4f}, y = {float(L1_y):.4f}")
plt.scatter(L2_x, L2_y, label=f"L2: x = {float(L2_x):.4f}, y = {float(L2_y):.4f}")
plt.scatter(L3_x, L3_y, label=f"L3: x = {float(L3_x):.4f}, y = {float(L3_y):.4f}")
plt.scatter(L4_x, L4_y, label=f"L4: x = {float(L4_x):.4f}, y = {float(L4_y):.4f}")
plt.scatter(L5_x, L5_y, label=f"L5: x = {float(L5_x):.4f}, y = {float(L5_y):.4f}")
plt.scatter(Sun_x, Sun_y, s=200, color='orange', label=f"Sun: x = {float(Sun_x):.4f}, y = {float(Sun_y):.4f}")
plt.scatter(Jupiter_x, Jupiter_y, s=100, color='brown', label=f"Jupiter: x = {float(Jupiter_x):.4f}, y = {float(Jupiter_y):.4f}")
plt.xlabel('x (DU)')
plt.ylabel('y (DU)')
plt.title("Sun Jupiter System 5 Lagrange Points")
plt.xlim(-1.5,1.5)
plt.ylim(-1.5,1.5)
plt.legend()
plt.tight_layout()
plt.show()
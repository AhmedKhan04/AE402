import numpy as np 

p_closest = 132252 
R_jupiter = 71492 
mu_jupiter = 1.2669e8
r_circ_jupiter = 5.2 
r_circ_earth = 1
mu_sun = 1.327e11


def get_a(r_a, r_p):
    return (r_a  + r_p)/2 

def get_period(a, mu): 
    return 2*np.pi * np.sqrt(a**3/mu)

a_hohmann = get_a(r_circ_jupiter, r_circ_earth)
time = get_period(a_hohmann, 4*np.pi**2)/2 # since in au the mu is 4 * pi**2 
print(f"ideal hohmann transfer time: {time}") # years

# From horizons: 
v_jupiter = np.array([1.300549460195539e00, -9.883165034312524e00, 4.813490050639437e-01])

r_jupiter = np.array([-1.925734664929881e06, 9.137231909342932e06, -6.294401986104711e05])

v_inf = np.sqrt( np.linalg.norm(v_jupiter)**2 -  2 * mu_jupiter / (np.linalg.norm(r_jupiter)))
print(f"v_inf: km/s: {v_inf}")

radius_pariapsis_jupiter = p_closest + R_jupiter
print(f"radius of pariapsis: {radius_pariapsis_jupiter}")
def get_v_p(v_inf, mu, r_p): 
    l = v_inf**2 + (2*mu)/r_p
    return np.sqrt(l)

v_p = get_v_p(v_inf, mu_jupiter, radius_pariapsis_jupiter)
print(f"v_p (km/s) {v_p}")
def get_turn_angle(r_p, v_inf, mu): 
    e = 1 + (r_p * v_inf**2)/mu
    return np.degrees(2*np.arcsin(e**-1))

print(f"turn angle: degrees: {get_turn_angle(radius_pariapsis_jupiter, v_inf, mu_jupiter)}")


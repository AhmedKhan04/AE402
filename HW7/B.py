import numpy as np 

v_initial = 2.5146
r_initial = 160.57
mu = 4902.8 
R_moon = 1737.4
Isp = 314.5


def get_v_inf(v_current, mu, r_current): 
    return np.sqrt(v_current**2  - 2*mu/r_current)
v_inf = get_v_inf(v_initial, mu, r_initial + R_moon)
print(f"v_inf {v_inf}")

R_p = R_moon + 111.3 
def get_turn_angle(v_inf, R_p, mu):
    e = 1 + (R_p*v_inf**2)/mu 

    return 2*np.arcsin(e**-1)
turn_angle = get_turn_angle(v_inf, R_p, mu)
print(f"turn angle {(turn_angle)}")
def get_delta_v(Isp, mf, m0):
    g0 = 9.81
    return Isp * g0 * np.log(m0/mf)

delta_v_rocket = get_delta_v(Isp, 32675.7, 43572.8)
print(f"delta v rocket {delta_v_rocket/1000}")

def get_v_circular(mu, r_mag):
    return np.sqrt(mu/r_mag)

r_circ = R_moon + 111.31 

v_circ = get_v_circular(mu, r_circ)
print(f"v_circ {v_circ}")
# convervation of energy 
# 0.5 * m * v**2 -  GMm/r = C 
# 0.5 * v**2  - mu/r = C 

def get_free_fall_v (v_1, r_1, r_2, mu): 
    C = 0.5 * v_1**2  - mu/r_1
    v_2 = np.sqrt((C + mu/r_2)/0.5)
    return v_2 

v_2 = get_free_fall_v(0, r_circ, R_moon, mu)

print(f"v free fall {v_2}")

r_p = R_moon 
r_a = R_moon + 111.3 
a = (r_p + r_a)/2

def get_v_from_vis_viva(mu, a, r_mag): 
    return np.sqrt(mu * (2/r_mag - 1/a))

v_initial = v_circ
v_final = get_v_from_vis_viva(mu, a, r_a)
delta_v = v_final - v_initial
print(f"delta v hohmann {delta_v * 1000}")

v_at_landing = get_v_from_vis_viva(mu, a, r_p)
print(f"delta v at landing {v_at_landing}")

def get_period(mu,a):
    return 2 * np.pi * np.sqrt(a**3/mu)

time_to_landing = get_period(mu,a)/2 
print(f"time to landing {time_to_landing}")

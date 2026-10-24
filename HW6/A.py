import numpy as np 


AU = 149597871 

r_earth  = 1  * AU # AU -> km
r_venus = 0.723  * AU # AU -> km 
r_mercury = 0.387 * AU # AU -> km 



mu_sun = 1.327e11
mu_earth = 398600
mu_venus = 324859 

p_R_earth  = 6378 
p_R_venus = 6052

p_closest = 5768 

# First lets do a hohmann transfer 

def get_a(r_a, r_p):
    return (r_a  + r_p)/2 

def get_velocity(a, r, mu):
    v_s = mu * (2/r - 1/a)
    return np.sqrt(v_s)

a_transfer_e_v = get_a(r_earth, r_venus)
v_departure = get_velocity(a_transfer_e_v, r_earth, mu_sun)
v_earth = get_velocity(r_earth, r_earth, mu_sun)

v_inf = v_earth - v_departure
# print(v_departure)
# print(v_earth)
print(f"v_inf for Earth to Venus: {v_inf}") # km/s

v_venus = get_velocity(r_venus, r_venus, mu_sun)

print(f"v_venus in circular {v_venus}")


a_transfer_v_merc = get_a(r_venus, r_mercury)
v_departure_2 = get_velocity(a_transfer_v_merc, r_venus, mu_sun)

v_inf_2 = v_venus - v_departure_2
# print(v_departure)
# print(v_earth)
print(f"v_inf for Venus to Mercury: {v_inf_2}") # km/s

v_arrival = get_velocity(a_transfer_e_v, r_venus, mu_sun)
v_inf_3 = v_arrival - v_venus
print(f"v_inf for arrival Venus: {v_inf_3}") # km/s

radius_pariapsis_venus = p_closest + p_R_venus

print(f"r_p for arrival Venus: {radius_pariapsis_venus}") # km

def get_turn_angle(r_p , v_inf, mu):
    sin_delt = (1 + (r_p * v_inf**2)/mu)**-1
    sin_delt = np.arcsin(sin_delt)
    sin_delt *= 2 
    return np.degrees(sin_delt) 

v_inf_real = 7.9466 
angle_turn = get_turn_angle(radius_pariapsis_venus, v_inf_real, mu_venus)
print(f"turn_angle: {angle_turn} degrees")

def get_e(angle):
    return (np.sin(np.radians(angle)/2))**-1

e_hyperbolic = get_e(angle_turn)
print(f"eccentricity: {e_hyperbolic}")
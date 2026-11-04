# delta_v apollo
import numpy as np  

mu = 398600 
r_earth = 6378 

orbit_height = 187.376 ## apogee
v_apogee = 6.91883

r_orbit = r_earth + orbit_height 

def get_a_from_vis_viva(v_mag,r_mag): 
    v=  v_mag 
    r = r_mag 
    return ((-(v**2)/mu) + (2/r))**-1

a = get_a_from_vis_viva(v_apogee,r_orbit)

def get_eccentricity_from_ra(a, r_mag): 
    return  r_mag/a - 1 

e = get_eccentricity_from_ra(a,r_orbit)

def get_f_at_r(a, e, r_mag, other_solution=False): 
    if other_solution: 
        return 2*np.pi -  np.arccos((a * (1-e**2) - r_mag)/(e*r_mag))
    return np.arccos((a * (1-e**2) - r_mag)/(e*r_mag))

f_at_apogee = np.pi #get_f_at_r(a,e,r_orbit)
f_at_ground = get_f_at_r(a,e,r_earth, other_solution=True)

def f_to_E(f,e):
    E = 2 * np.arctan2(np.sqrt(1-e) * np.sin(f/2), np.sqrt(1+e) * np.cos(f/2))
    return E 

E_at_Apogee = f_to_E(f_at_apogee,e)
E_at_Ground = f_to_E(f_at_ground,e)
# print(e)
# print(np.degrees(f_at_ground))
# print(np.degrees(E_at_Apogee))
# print(np.degrees(E_at_Ground))

def get_n (mu, a):
    return np.sqrt(mu/a**3)

n = get_n(mu,a)
delta_time_from_perigee_to_apogee = (E_at_Apogee - e*np.sin(E_at_Apogee))/n 
delta_time_from_perigee_to_ground = (E_at_Ground - e*np.sin(E_at_Ground))/n
delta_time_from_apogee_to_ground = delta_time_from_perigee_to_ground - delta_time_from_perigee_to_apogee
#print(delta_time_from_apogee_to_ground)

print(f"semi_major axis {a}")
print(f"eccentricity {e}")
print(f"time: {delta_time_from_apogee_to_ground}")

a_circle = r_earth + orbit_height 

def get_v_circular(mu, r_mag):
    return np.sqrt(mu/r_mag)

def get_v_elliptical(mu, a, r_mag):
    print ("a", a)
    print("r", r_mag)
    print("mu", mu)
    return np.sqrt(mu * (2/r_mag - 1/a))

v_circ = get_v_circular(mu, a_circle)
delta_v = v_circ - v_apogee 
print(f"delta v {delta_v}")

Isp = 421 

def get_mass_ratio(delta_v, Isp): 
    g = 9.81 
    return np.exp(delta_v/(Isp*g))

mass_ratio = get_mass_ratio(delta_v * 1000, Isp)

print(f"mass ratio {mass_ratio}")

# mass ratio = m0/mf 
# mass ratio ** -1 = mf/m0 
m0 = 166571
mf = mass_ratio**-1 * m0

print(f"m_final {mf}")
m0_r = 134046
mf_r = 63116 
mass_real = m0_r/mf_r

def get_delta_v(ratio, Isp): 
    g = 9.81 
    return Isp * g * np.log(ratio)

delta_v_real = get_delta_v(mass_real, Isp)
print(f"delta v real {delta_v_real/1000}")
r_earth_moon = np.array([-3.923817963190890E+05, 4.118232536241712E+04, -6.892403934792474E+02])
r_earth_moon_mag = np.linalg.norm(r_earth_moon)
print(f"distance to moon {r_earth_moon_mag}")
def get_r_SOI(mass_ratio, a):
    return a * (mass_ratio/3)**(1/2)

EarthMoonratio = 81.3005690769
moonEarthratio = EarthMoonratio**-1 
r_SOI_moon = get_r_SOI(moonEarthratio, r_earth_moon_mag)
print(f"sphere of influence {r_SOI_moon}")

rp =  r_earth + orbit_height 
ra = r_earth_moon_mag 
a = (rp + ra) / 2

v2 = get_v_elliptical(mu, a, rp)
v1 = get_v_circular(mu, rp)
dv1 = v2 - v1

v2 = get_v_elliptical(mu, a, ra)
v1 = get_v_circular(mu, ra)
dv2 = v1 - v2

print(f"v1 {v1}")
print(f"v2 {v2}")
print(f"theoretical delta v {dv1}")
print(f"difference {delta_v_real/1000 - np.abs( dv1)}")
# 2440422.208333333 = A.D. 1969-Jul-19 17:00:00.0000 TDB 
#  X =-3.923817963190890E+05 Y = 4.118232536241712E+04 Z =-6.892403934792474E+02
#  VX=-5.615334915052516E-02 VY=-9.891802933004656E-01 VZ=-8.852424885555377E-02
#  LT= 1.316035840094653E+00 RG= 3.945376193180709E+05 RR=-4.725069871887726E-02
#  81.3005690769   mass ratio earth ot moon from horizons
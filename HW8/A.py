# Earth 
# 2444907.500000000 = A.D. 1981-Oct-30 00:00:00.0000 TDB [del_T=     52.182506 s]
#  X = 1.190386338685324E+08 Y = 8.886577991548955E+07 Z = 3.845452046155930E+03
#  VX=-1.831566510772723E+01 VY= 2.376758260846223E+01 VZ= 1.497952092988086E-03
#  LT= 4.955119413738054E+02 RG= 1.485507428728050E+08 RR=-4.587454800452966E-01

# Venus
# 2445029.500000000 = A.D. 1982-Mar-01 00:00:00.0000 TDB [del_T=     52.185403 s]
#  X =-1.073335973692154E+08 Y =-7.393965232027010E+06 Z = 6.096223037918333E+06
#  VX= 2.193935361550072E+00 VY=-3.509250894549594E+01 VZ=-6.045300912325420E-01
#  LT= 3.594504980951503E+02 RG= 1.077605483532694E+08 RR= 1.884220741067862E-01

# Output units    : KM-S

import numpy as np 
from scipy.optimize import fsolve
np.set_printoptions(precision=3, suppress=True)

AU = 1.496e+8 # km 
TU = 365/(2*np.pi)

pos_earth = np.array([1.190386338685324e8,8.886577991548955e7,3.845452046155930e3])/AU 
pos_venus = np.array([-1.073335973692154e8,-7.393965232027010e6,6.096223037918333e6])/AU
v_earth = np.array([-1.831566510772723e1,  2.376758260846223e1,  1.497952092988086e-3]) * (TU*86400/AU)
v_venus = np.array([2.193935361550072e0,-3.509250894549594e1,-6.045300912325420e-1]) * (TU*86400/AU)

mag_earth_pos = np.linalg.norm(pos_earth)
mag_venus_pos = np.linalg.norm(pos_venus)

print(f"earth mag: {mag_earth_pos:.3f}")
print(f"venus mag: {mag_venus_pos:.3f}")



def get_theta(r1, r2, v1, v2):
    r1_mag = np.linalg.norm(r1)
    r2_mag = np.linalg.norm(r2)
    theta = np.arccos(np.dot(r1, r2)/(r1_mag * r2_mag))
    if ((np.cross(r1, v1)[2]< 0 and np.cross(r2, v2)[2]< 0 )  ):
        print("retrograde")
        theta = np.pi - theta
    elif((np.cross(r1, v1)[2] > 0 and np.cross(r2, v2)[2] > 0)):
        print("prograde") 
    else: 
        print("mismatch!!")

    return r1_mag, r2_mag,  theta
    

def get_chord(r1,r2, v1, v2):
    r1_mag, r2_mag, theta = get_theta(r1,r2, v1, v2)
    return np.sqrt(r1_mag**2 + r2_mag**2 - 2*r1_mag*r2_mag*np.cos(theta))

chord = get_chord(pos_earth, pos_venus, v_earth, v_venus)
print(f"chord: {chord:.3f}")

def get_semi_perimeter(r1_m, r2_m, c): 
    return (r1_m + r2_m + c)/2

s = get_semi_perimeter(mag_earth_pos, mag_venus_pos, chord)

print(f"semi_perimeter: {s:.3f}")

def get_am(semi_perimeter): 
    return semi_perimeter/2 

am = get_am(s)
print(f"minimum semi major {am:.3f}")


_, _, theta = get_theta(pos_earth, pos_venus, v_earth, v_venus)
print(f"angle between vectors: {theta:.3f} radians")

def get_beta(s, c, a): 
    rhs = np.sqrt((s-c)/(2*a))
    return 2*np.arcsin(rhs)

def get_alpha(s, a): 
    rhs = np.sqrt((s)/(2*a))
    return 2*np.arcsin(rhs)

def get_tm(a, mu, alpha, beta):
    x = np.sqrt(a**3/mu)
    x *= (alpha - beta - (np.sin(alpha) - np.sin(beta)))
    return x 

mu = 1
tm = get_tm(am, 1, get_alpha(s, am), get_beta(s, chord, am)) # TU 

print(f"tm: {tm:.3f}")
target_time = 122 / TU


def cost_lambert(a_guess):
    alpha = get_alpha(s, a_guess)
    beta = get_beta(s, chord, a_guess)
    TOF = 2*a_guess**(3/2)*((1/2)*(alpha - beta) - (1/2)*(np.sin(alpha) - np.sin(beta))) # sqrt mu assumed = 1 
    cost = 0.5 * (TOF - target_time)**2
    # print(f"cost!! {cost}")
    return cost 

a_final = fsolve(cost_lambert, [am])
print(f"a-final solved for: {a_final}")

u_r1 = pos_earth/mag_earth_pos
u_r2 = pos_venus/mag_venus_pos

print(f"u_r1: {u_r1}")
print(f"u_r2: {u_r2}")

uc = (pos_venus-pos_earth)/chord
print(f"unit_vec_chord: {uc}")

alpha_final = get_alpha(s, a_final)
beta_final = get_beta(s, chord, a_final)

def get_A(a, alpha):
    return np.sqrt(1/(4*a))  * (np.tan(alpha/2))**-1


def get_B(a, beta):
    return np.sqrt(1/(4*a))  * (np.tan(beta/2))**-1
    
A = get_A(a_final, alpha_final)
B = get_B(a_final, beta_final)

print(f"alpha_final: {alpha_final}")
print(f"beta_final: {beta_final}")
print(f"A_final: {A}")
print(f"B_final: {B}")

v1 = (B+A)*uc + (B-A)*u_r1
print(f"v_1: {v1}")
print(f"v_earth: {v_earth}")

delta_v = np.linalg.norm(v1 - v_earth)
delta_v *= (TU*86400/AU)**-1

print(f"delta_v {delta_v:.3f} km/s")

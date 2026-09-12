import numpy as np 

# Question 1: 

# To obtain the semi major axis of TelStar-1, we apply the vis viva equation 

# Givens: 
Earth_Mu = 398600 # km^3/s^2
R_earth = 6378 # km 
r_p = 952 # km 
v_p = 8.256 #km/s
# v = sqrt(Earth_Mu * (2/r - 1/a))

def a_from_vis_vis(v, r):
    a = (-(v**2)/Earth_Mu + 2/r) ** -1 
    return a 

# to calculate eccentricity we Earth_Must use the definition of perigee radius 

def e_from_rp(rp, a):
    e = 1 - rp/a
    return e 

# we now use angular momentum definition 

def get_angular_momentum(a, e): 
    h = np.sqrt(Earth_Mu * a * (1-e**2))
    return h 

# we now use radius definition 

def altitude_at_f(a,e,f): 
    r = (a*(1-e**2))/(1+e*np.cos(np.radians(f))) # to center
    r -= R_earth # altitude 
    return r

# Flight path angle 

def get_flight_path_angle(h, e, f): 
    vr = (Earth_Mu/h) * e*np.sin(np.radians(f))
    v_perpindicular = (Earth_Mu/h) * e*np.sin(np.radians(f))
    flight_path_angle = np.atan2(vr, v_perpindicular)
    return flight_path_angle

Telstar_a = a_from_vis_vis(v_p, r_p)
Telstar_e = e_from_rp(r_p, Telstar_a)
Telstar_h = get_angular_momentum(Telstar_a, Telstar_e)
Telstar_alt = altitude_at_f(Telstar_a,Telstar_e, f = 120)
Telstar_flight_angle = get_flight_path_angle(Telstar_h, Telstar_e, f=120)

# Question 2 

molniya_a = 26500
molniya_e = 0.73


# f = 90 degrees by definition of semi latus rectum 

molniya_f = 90 

def get_ecentric_anomaly(e, f):
    tan_E_2 = np.sqrt(((1-e)/(1+e)) * np.tan(np.radians(f)/2))
    E = 2 * np.arctan2(tan_E_2, 1)
    return E 


def get_Mean_anomaly(E, e): 
    M = E - e*np.sin(E)
    return M

def get_mean_angular_motion(a): 
    n = np.sqrt(Earth_Mu/a**3)
    return n 

def get_time_change(M, n): 
    return  M/n

molniya_E = get_ecentric_anomaly(molniya_e, molniya_f)
molniya_M = get_Mean_anomaly(molniya_E, molniya_e)
molniya_n = get_mean_angular_motion(molniya_a)
molniya_time = get_time_change(molniya_M, molniya_n)


# Question 3

Mars_Mu = 42828
R_mars = 3389 
r_p_2 = 1650  + R_mars
r_a_2 = 16680  + R_mars
target_alt = r_a_2 


def get_eccentricity_mars(r_a, r_p):
    e = r_a - r_p
    e /= r_a + r_p 
    return e 

def get_semi_major_axis_mars(r_a, e): 
    a = r_a/(1+e)
    return a 


def get_angular_momentum_mars(a, e):
    h = np.sqrt(Mars_Mu * a * (1-e**2))
    return h 


def get_period_mars(a): 
    P = 2*np.pi * np.sqrt(a**3/Mars_Mu)
    return P 

def get_mean_angular_motion_mars(a): 
    n = np.sqrt(Mars_Mu/a**3)
    return n  

def get_anomoly_from_alt_mars(altitude, e, a): 
    r = altitude + R_mars
    r = r**-1 
    r *= a*(1-e**2)
    r -= 1 
    r /= e 
    return np.acos(r)

def get_ecentric_anomoly_mars(e, true_anomoly): 
    ec = np.sqrt(((1-e)/(1+e)) * np.tan(true_anomoly/2))
    ec = np.atan2(ec,1) * 2 
    return ec 

def get_mean_anomoly_mars(e, E):
    M = E - e*np.sin(E)
    return M 

def get_time_on_true_anomoly_mars(n, M):
    dt = M/n
    return dt 

def get_total_above_alt_time_mars(P, time_to_get_at_f): 
    # we consider the orbit symettircal 
    # thus the time to get to a certain altittude is t, and the total time spent above it is P - 2t, where 2 accounts for the second half 

    return P - 2*time_to_get_at_f


e_mars = get_eccentricity_mars(r_a=r_a_2, r_p=r_p_2)
a_mars = get_semi_major_axis_mars(r_a_2, e_mars)
h_mars = get_angular_momentum_mars(a_mars, e_mars)
P_mars = get_period_mars(a_mars)
n_mars = get_mean_angular_motion_mars(a_mars)
f_mars_4000 = get_anomoly_from_alt_mars(4000, e_mars, a_mars)
E_mars_4000 = get_ecentric_anomoly_mars(e_mars, f_mars_4000)
M_mars_4000 = get_mean_anomoly_mars(e_mars, E_mars_4000)
dt_mars_4000 = get_time_on_true_anomoly_mars(n_mars, M_mars_4000)
total_time_spent = get_total_above_alt_time_mars(P_mars, dt_mars_4000)

Q1 = [Telstar_a, Telstar_e, Telstar_h, Telstar_alt, Telstar_flight_angle]
Q2 = [molniya_E , molniya_M ,molniya_n , molniya_time]
Q3 = [e_mars, a_mars, h_mars , P_mars , n_mars , f_mars_4000 , E_mars_4000, M_mars_4000, dt_mars_4000, total_time_spent]

print("%" + "*"*50 + "%")
print('Question 1 :')
for val in Q1: 
    print(f"{val}\n")

print("%" + "*"*50 + "%")
print('Question 2 :')
for val in Q2: 
    print(f"{val}\n")

print("%" + "*"*50 + "%")
print('Question 3 :')
for val in Q3: 
    print(f"{val}\n")
print("%" + "*"*50 + "%")


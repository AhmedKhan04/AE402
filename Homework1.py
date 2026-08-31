import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sp 

R = 6378 # km
mu = 398600 # km^3/s^2
time = 3 * 60 * 60 # 3 hours -> seconds

def main():
    r_0 = np.array([-307.8439, -4517.4816, -5747.6257]) # km
    v_0 = np.array([5.1593, 3.795, -3.2591]) #km
    state_0 = np.array([*v_0, *r_0])
    
    sol  = sp.solve_ivp(ODE, (0, time), state_0, rtol = 1e-13)
    time_pts = sol['t']
    state_vects = sol['y']

    velocity_vecs = state_vects[:3].T
    position_vecs = state_vects[3:].T


    
    x, y, z = position_vecs.T[0], position_vecs.T[1], position_vecs.T[2]
    fig = plt.figure() 
    ax = fig.add_subplot(111, projection ="3d")
    ax.plot(x, y, z)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    plt.show()

    

    



def ODE(time, state):
    v = state[:3]
    r = state[3:]
    r_dd = -mu*r/(R**3)
    r_d = v 
    return np.array([*r_dd, *r_d]) 



main() 
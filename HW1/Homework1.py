import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sp 

R = 6378 # km
mu = 398600 # km^3/s^2
time = 3 * 60 * 60 # 3 hours -> seconds

def main():
    # initial conditions 
    r_0 = np.array([-307.8439, -4517.4816, -5747.6257]) # km
    v_0 = np.array([5.1593, 3.795, -3.2591]) #km
    state_0 = np.array([*v_0, *r_0]) # initial state

    # numerical integration 
    sol  = sp.solve_ivp(ODE, (0, time), state_0, method='RK45', rtol = 1e-13, atol = 1e-13) # solving the ODE 
    time_pts = sol['t'] # time points 
    state_vects = sol['y'] # state vectors 

    velocity_vecs = state_vects[:3].T # velocity vecs
    position_vecs = state_vects[3:].T # position vecs 

    angular_vel_vecs = [] 
    h_t0 = np.linalg.cross(position_vecs[0], velocity_vecs[0]) # initial specific angular momentum 
    for vel_vec, pos_vec in zip(velocity_vecs, position_vecs): 
        h_t = np.linalg.cross(pos_vec,vel_vec)   # angular momentum 
        angular_vel_vecs.append(np.linalg.norm(h_t - h_t0)/np.linalg.norm(h_t0)) # conservation error


    sol2  = sp.solve_ivp(ODE, (0, time), state_0,method='RK45', rtol = 1e-7, atol = 1e-7) # solving the ODE 
    time_pts2 = sol2['t'] # time points 
    state_vects2 = sol2['y'] # state vectors 

    velocity_vecs2 = state_vects2[:3].T # velocity vecs
    position_vecs2 = state_vects2[3:].T # position vecs 

    angular_vel_vecs2 = [] 
    h_t02 = np.linalg.cross(position_vecs2[0], velocity_vecs2[0]) # initial specific angular momentum 
    for vel_vec, pos_vec in zip(velocity_vecs2, position_vecs2): 
        h_t2 = np.linalg.cross(pos_vec,vel_vec)  # angular momentum 
        angular_vel_vecs2.append(np.linalg.norm(h_t2 - h_t02)/np.linalg.norm(h_t02)) # conservation error

    # Plotting the results
    x, y, z = position_vecs.T[0], position_vecs.T[1], position_vecs.T[2]
    fig = plt.figure(dpi=300) 
    ax = fig.add_subplot(111, projection ="3d")
    ax.plot(x, y, z)
    ax.set_title("Sputnik Orbit (tol = 1e-13)")
    ax.set_xlabel("X (km)")
    ax.set_ylabel("Y (km)")
    ax.set_zlabel("Z (km)")
    
    fig2 = plt.figure(dpi=300) 
    ax2 = fig2.add_subplot(111)
    ax2.plot(time_pts, angular_vel_vecs)
    ax2.set_title("Specific Angular Momentum Conservation Error v Time (tol = 1e-13)")
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("Specific Angular Momentum Conservation Error")
    ax2.set_yscale("log")

    x2, y2, z2 = position_vecs2.T[0], position_vecs2.T[1], position_vecs2.T[2]
    fig3 = plt.figure(dpi=300) 
    ax3 = fig3.add_subplot(111, projection ="3d")
    ax3.plot(x2, y2, z2)
    ax3.set_title("Sputnik Orbit (tol = 1e-7)")
    ax3.set_xlabel("X (km)")
    ax3.set_ylabel("Y (km)")
    ax3.set_zlabel("Z (km)")
    
    fig4 = plt.figure(dpi=300) 
    ax4 = fig4.add_subplot(111)
    ax4.plot(time_pts2, angular_vel_vecs2)
    ax4.set_title("Specific Angular Momentum Conservation Error v Time (tol = 1e-7)")
    ax4.set_xlabel("Time (s)")
    ax4.set_ylabel("Specific Angular Momentum Conservation Error")
    ax4.set_yscale("log")

    fig.savefig("./HW1/orbit_13")
    fig2.savefig("./HW1/conv_13")
    fig3.savefig("./HW1/orbit_7")
    fig4.savefig("./HW1/conv_7")
    plt.show()

    

    



def ODE(time, state):
    v = state[:3]
    r = state[3:]
    r_dd = -mu*r/(np.linalg.norm(r)**3)
    r_d = v 
    return np.array([*r_dd, *r_d]) 


# https://chatgpt.com/share/6a9607e9-fbe0-83ea-9a53-09eb333a7289
main() 
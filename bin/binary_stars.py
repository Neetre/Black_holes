'''
Simulation of a system of binary stars.

Neetre 2024
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Const
G = 6.674e-11
M_sun = 1.989e30

# Starting values
m1 = M_sun * 2.5
m2 = M_sun * 2

r1 = 1.2e9
r2 = 0.8e9

x1, y1 = -1e11, 0
x2, y2 = 1e11, 0

vx1, vy1 = 0, -2e4
vx2, vy2 = 0, 2e4


def derivates(state, t):
    '''
    This function calculates the derivatives of the state vector.
    '''
    x1, y1, x2, y2, vx1, vy1, vx2, vy2 = state
    
    # Distance between the stars
    dx, dy = x2 - x1, y2 - y1
    r = np.sqrt(dx**2 + dy**2)
    fx, fy = G * m1 * m2 * dx / r**3, G * m1 * m2 * dy / r**3
    
    # Gravitational acceleration
    ax1, ay1 = fx / m1, fy / m1
    ax2, ay2 = -fx / m2, -fy / m1
    
    # Tider forces
    k = 1e-7
    fx_tide = k * (m1 * r1**5 / r**6 + m2 * r2**5 / r**6) * dx
    fy_tide = k * (m1 * r1**5 / r**6 + m2 * r2**5 / r**6) * dy
    
    # Tider acceleration
    ax1_tide, ay1_tide = -fx_tide / m1, -fy_tide / m1
    ax2_tide, ay2_tide = fx_tide / m2, fy_tide / m2
    
    # Total acceleration
    ax1 += ax1_tide
    ay1 += ay1_tide
    ax2 += ax2_tide
    ay2 += ay2_tide
    
    return [vx1, vy1, vx2, vy2, ax1, ay1, ax2, ay2]


def plot_binary_stars(x1, y1, x2, y2):
    plt.figure(figsize = (10, 6))
    plt.plot(x1, y1, label="Stella 1")
    plt.plot(x2, y2, label="Stella 2")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.title("Evoluzione sistema di Stelle Binarie")
    plt.grid()
    plt.legend()
    plt.axis("equal")
    plt.show()


state0 = [x1, y1, x2, y2, vx1, vy1, vx2, vy2] # Initial state

def binary_stars():
    global state0
    
    t_span = np.linspace(0, 5 *365 * 24 * 3600, 1000) # 5 years in seconds
    solution = odeint(derivates, state0, t_span)  # Solution of the ODE

    # Results extraction
    x1, y1 = solution[:, 0], solution[:, 1]
    x2 , y2 = solution[:, 2], solution[:, 3]
    
    plot_binary_stars(x1, y1, x2, y2)


if __name__ == "__main__":
    binary_stars()
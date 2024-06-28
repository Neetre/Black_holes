'''
This script is used to calculate the BH mass and the event horizon radius of a black hole.

Neetre 2024
'''

import numpy as np
import matplotlib.pyplot as plt


# Const
G = 6.67430e-11
c = 299792458
M_s = 1.98847e30


def eh_rotating_bh(M, a, theta, plotting=False):
    '''
    Calculate the event horizon of a black hole given its mass M, its spin a
    and the angle theta.
    '''
    M = M * M_s
    a = a * G * M / (c**2)   # https://en.wikipedia.org/wiki/Schwarzschild_radius
    r = G * M / c**2

    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    
    # rho = np.sqrt(r**2 + a **2)
    r_plus = r + np.sqrt(r**2 - a**2)
    print(f"Radius of the outer event horizon: {r_plus} m")
    r_minus = r - np.sqrt(r**2 - a**2)
    print(f"Radius of the inner event horizon: {r_minus} m")
    
    if plotting:
        x_plus = r_plus * sin_theta
        y_plus = r_plus * cos_theta
        x_minus = r_minus * sin_theta
        y_minus = r_minus * cos_theta
        
        coord = x_plus, y_plus, x_minus, y_minus
        
        return coord
    
    
def plot_ev_rotating_bh(coord: tuple):
    x_plus, y_plus, x_minus, y_minus = coord
    plt.plot(x_plus, y_plus, 'r', label='Outer EV')
    plt.plot(x_minus, y_minus, 'r', label='Inner EV')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Event Horizon of a Rotating Black Hole')
    plt.legend()
    plt.show()


def ev_non_rotating_bh(M, a):
    return


def plot_ev_bh(coord: tuple):
    return


def ev():
    M = 4.297e6  # This is the 'M' of Sagittarius A* in the Milky Way
    a = 0.99616  # This is the 'a' of Sagittarius A* in the Milky Way
    theta = np.linspace(0, 2 * np.pi, 1000)
    coord = eh_rotating_bh(M, a, theta, True)
    plot_ev_rotating_bh(coord)
    

if __name__ == "__main__":
    ev()

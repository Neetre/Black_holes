'''
This program simulates the effect of a gravitational lensing by a black hole.
The black hole is assumed to be rotating.
The event horizon of the black hole is calculated and plotted.
The event horizon of a non-rotating black hole is also calculated and plotted.
The effect of gravitational lensing is then simulated by plotting the path of a light ray that passes close to the black hole.
The light ray is bent by the gravitational field of the black hole, resulting in a distortion of the image of a background object.

Neetre 2024
'''

import numpy as np
import matplotlib.pyplot as plt
from ev_bh import eh_rotating_bh, plot_ev_rotating_bh, ev_non_rotating_bh

M = 10e14
DL = 1 # Gpc
DS = 2 # Gpc

G = 6.67430e-11
c = 299792458
M_s = 1.98847e30


def einstein_radius(M, DL, DS):
    '''
    Calculate the Einstein radius of a black hole given its mass M, the distance
    from the observer to the lens DL and the distance from the observer to the source DS.
    '''
    
    M = M * M_s
    DL = DL * 3.086*10**16 # m
    DS = DS * 3.086*10**16 # m
    DLS = DS - DL
    theta_E = np.sqrt(4 * G * M / (c**2) * DLS / (DL * DS))
    print(f"Einstein radius: {theta_E} rad")
    return theta_E


def position_image(theta_E, theta, a):
    '''
    Calculate the position of the image of a background object given the Einstein radius theta_E,
    the angle theta and the spin of the black hole a.
    '''
    r_plus, r_minus = eh_rotating_bh(M, a, theta, rad=True)
    x_plus = r_plus * np.sin(theta)
    y_plus = r_plus * np.cos(theta)
    
    x_minus = r_minus * np.sin(theta)
    y_minus = r_minus * np.cos(theta)
    
    x_image = x_plus + theta_E
    y_image = y_plus + theta_E
    return x_image, y_image, x_plus, y_plus, x_minus, y_minus


def plot_image(coord):
    '''
    Plot the image of a background object that has been gravitationally lensed
    by a rotating black hole
    '''
    
    x_image, y_image, r_plus, y_plus, r_minus, y_minus = coord
    
    # Event Horizon plots
    x_plus = r_plus * np.sin(np.linspace(0, 2 * np.pi, 1000))
    y_plus = r_plus * np.cos(np.linspace(0, 2 * np.pi, 1000))
    x_minus = r_minus * np.sin(np.linspace(0, 2 * np.pi, 1000))
    y_minus = r_minus * np.cos(np.linspace(0, 2 * np.pi, 1000))
    
    plt.plot(x_plus, y_plus, 'r', label='Outer EV')
    plt.plot(x_minus, y_minus, 'k', label='Inner EV')
    plt.plot(x_image, y_image, 'bo', label='Image')
    plt.plot(0, 0, 'yo', label='Black Hole')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Gravitational Lensing by a Rotating Black Hole')
    plt.legend()
    plt.show()



def lensing_bh():
    theta_E = einstein_radius(M, DL, DS)
    theta = np.pi / 4
    a = 0.5
    coord = position_image(theta_E, theta, a)
    plot_image(coord)
    
    return


if __name__ == "__main__":
    lensing_bh()

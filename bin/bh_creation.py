'''
This script simulate the creation of a black hole from a star.

Neetre 2024
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import odeint



# Const
G = 6.67430e-11
c = 299792458
M_s = 1.989e30
R_s = 6.957e8
sigma_SB = 5.670e-8  # Stefan-Boltzmann constant
a_rad = 4 * sigma_SB / c
G_F = 1.436e-62  # Fermi constant

initial_mass = 40 * M_s  # kg
initial_rad = 20 * R_s  # m
initial_temp = 1e7  # K
omega = 1e-4
B_0 = 1e8


# Simultaion param
t_max = 1000
dt = 0.1
N = 200


t = np.arange(0, t_max, dt)
mass = np.linspace(0, initial_mass, N)
radius = np.zeros((len(t), N))
temp = np.zeros((len(t), N))
rho = np.zeros((len(t), N))
P = np.zeros((len(t), N))
Y_e = np.ones((len(t), N)) * 0.5
X = np.ones((len(t), N)) * 0.7


radius[0] = np.linspace(1e-10, initial_rad, N)
rho[0] = initial_mass / (4/3 * np.pi * initial_rad**3) * (1 - (radius[0]/initial_rad)**2)
temp[0] = initial_temp * (1 - 0.9 * (radius[0]/initial_rad))
P[0] = G * initial_mass**2 / (8 * np.pi * initial_rad**4) * (1 - (radius[0]/initial_rad)**2)


def TOV_equations(y, m_r, M, R):
    r, P = y
    rho = P / (G * M**2 / (8 * np.pi * R**4))  # Simplified EOS
    drdm = 1 / (4 * np.pi * rho * r**2)
    dPdm = -(G * m_r * rho / r**2) * (1 + P / (rho * c**2)) * (1 + 4 * np.pi * r**3 * P / (m_r * c**2)) / (1 - 2 * G * m_r / (r * c**2))
    return [drdm, dPdm]


def update_rad(rho, T, r, dm):
    kappa = 0.2 * (1 + X[0]) + 0.1  # Thomson + heavy elements
    P_rad = a_rad * T**4 / 3
    L_rad = -64 * np.pi * sigma_SB * r**4 * T**3 * np.gradient(T, dm) / (3 * kappa * rho)
    return P_rad, L_rad


def update_neu(rho, T, Y_e, r):
    # Stephen Hawking my father
    Gamma_ec = 1e-10 * rho * Y_e * T**5
    Gamma_pc = 1e-10 * rho * (1 - Y_e) * T**5
    Gamma_pp = 1e-10 * rho * T**4
    
    dY_edt = (Gamma_pc - Gamma_ec) / rho
    Y_e += dY_edt * dt
    
    L_nu = 4 * np.pi * r**2 * (Gamma_ec + Gamma_pc + 2 * Gamma_pp) * (4 * G_F * 1.67e-27 * T**2)

    return Y_e, L_nu


def update_rotation_magnetic(r, rho, omega, R):
    omega_r = omega * R**2 / r**2
    P_rot = 1/3 * rho * (omega_r * r)**2
    
    B = B_0 * (R / r)**2
    P_mag = B**2 / (8 * np.pi)
    
    return P_rot, P_mag


bh_formed = False
bh_time = 0



for i in range(1, len(t)):
    M = np.sum(mass * (radius[i-1, 1] - radius[i-1, 0]))
    R = radius[i-1, 1]
    
    for j in range(1, N):
        y_0 = [radius[i-1, j], P[i-1, j]]
        sol = odeint(TOV_equations, y_0, [mass[j-1], mass[j]], args=(M, R))
        radius[i, j], P[i, j] = sol[-1]
        
    rho[i] = P[i] / (G * M**2 / (8 * np.pi * R**4))
    temp[i] = (3 * P[i] / a_rad) ** (1/4)


if not bh_formed:
    print("Black hole did not form within the simulation time.")
else:
    print(f"Final black hole mass: {mass[-1]/M_s:.2f} solar masses")
    print(f"Final black hole radius: {radius[-1]:.2e} meters")

# Print final state regardless of black hole formation
print(f"Final mass: {mass[-1]/M_s:.2f} solar masses")
print(f"Final radius: {radius[-1]/R_s:.2f} solar radii")
print(f"Final temperature: {temp[-1]:.2f} K")


# plotting time
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 15))

line1, = ax1.plot([], [], 'r-', label='Mass')
ax1.set_xlim(0, t_max)
ax1.set_ylim(0, 1.1 * initial_mass)
ax1.set_ylabel('Mass (Kg)')
ax1.legend()


line2, = ax2.plot([], [], 'b-', label='Radius')
ax2.set_xlim(0, t_max)
ax2.set_ylim(0, 1.1 * initial_rad)
ax2.set_ylabel('Radius (m)')
ax2.legend()


line3, = ax3.plot([], [], 'g-', label='Temperature')
ax3.set_xlim(0, t_max)
ax3.set_ylim(0, 1.1 * initial_temp)
ax3.set_xlabel("Time")
ax3.set_ylabel("Temperature (k)")
ax3.legend()


def animate(i):
    line1.set_data(t[:i], mass[:i])
    line2.set_data(t[:i], radius[:i])
    line3.set_data(t[:i], temp[:i])
    return line1, line2, line3


anim = FuncAnimation(fig, animate, frames=len(t), interval=50, blit=True)
plt.tight_layout()
plt.show()

anim.save('../data/black_hole_formation.mp4', writer='ffmpeg', fps=30)

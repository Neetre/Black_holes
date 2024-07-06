'''
This script simulate the creation of a black hole from a star.

Neetre 2024
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import solve_ivp
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)



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

count = 0
def TOV_equations(m_r, y, M, R):
    global count
    r, P = y
    if r <= 0 or P <= 0:
        return [0, 0]
    rho = max(P / (G * M**2 / (8 * np.pi * R**4)), 1e-10)
    drdm = 1 / (4 * np.pi * rho * r**2)
    dPdm = -(G * m_r * rho / r**2) * (1 + P / (rho * c**2)) * (1 + 4 * np.pi * r**3 * P / (m_r * c**2)) / (1 - 2 * G * m_r / (r * c**2))
    print(count)
    count += 1
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

count = 0
for i in range(1, len(t)):
    M = np.sum(mass * (radius[i-1, 1] - radius[i-1, 0]))
    R = radius[i-1, -1]
    
    for j in range(1, N):
        sol = solve_ivp(TOV_equations, [mass[j-1], mass[j]], [radius[i-1, j], P[i-1, j]], args=(M, R))
        radius[i, j], P[i, j] = sol.y[:, -1]
    
    rho[i] = np.maximum(P[i] / (G * M**2 / (8 * np.pi * R**4)), 1e-10)
    temp[i] = np.maximum((3 * P[i] / a_rad) ** (1/4), 1)
    
    P_rad, L_rad = update_rad(rho[i], temp[i], radius[i], mass[1]-mass[0])
    Y_e[i], L_nu = update_neu(rho[i], temp[i], Y_e[i-1], radius[i])
    P_rot, P_mag = update_rotation_magnetic(radius[i], rho[i], omega, R)
    
    P[i] += P_rad + P_rot + P_mag
    
    if np.min(radius[i]) <= 2 * G * M / c**2:
        print(f"Black Hole formed at t= {t[i]:.2f}")
        break
    print(count)
    count += 1


# plotting time
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))


def animate(frame):
    ax1.clear()
    ax2.clear()
    
    r = radius[frame]
    theta = np.linspace(0, 2*np.pi, 100)
    R, Theta = np.meshgrid(r, theta)
    X, Y = R * np.cos(Theta), R * np.sin(Theta)
    
    im = ax1.pcolormesh(X, Y, np.log10(rho[frame][np.newaxis, :]), cmap='viridis', shading='auto')
    ax1.set_xlim(-R_s, R_s)
    ax1.set_ylim(-R_s, R_s)
    ax1.set_aspect('equal')
    ax1.set_title(f't = {t[frame]:.2f} s')
    plt.colorbar(im, ax=ax1, label='log(ρ) [kg/m³]')
    
    ax2.plot(r / R_s, P[frame], label='Total')
    ax2.plot(r / R_s, P_rad, label='Radiation')
    ax2.plot(r / R_s, P_rot, label='Rotation')
    ax2.plot(r / R_s, P_mag, label='Magnetic')
    ax2.set_xlabel('r/R_sun')
    ax2.set_ylabel('Pressure [Pa]')
    ax2.set_yscale('log')
    ax2.set_title('Pressure Components')
    ax2.legend()
    
    return ax1, ax2


def create_bh():
    anim = FuncAnimation(fig, animate, frames=len(t), interval=50, blit=True)
    plt.tight_layout()
    plt.show()

    anim.save('../data/black_hole_formation.mp4', writer='ffmpeg', fps=30)
    

if __name__ == "__main__":
    create_bh()

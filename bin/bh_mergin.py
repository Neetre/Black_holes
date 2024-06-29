'''
Merge of two Black Holes and
the resulting gravitational waves

Neetre 2024
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# Const
G = 6.67430e-11
c = 299792458
M_s = 1.989e30

# BH param
M1 = 30 * M_s
M2 = 20 * M_s
initial_separation = 1e9 # m

# Tide param
t_max = 100
dt = 0.1

M_total = M1 + M2
Mu = (M1 * M2) / M_total


def orbital_freq(r):
    return np.sqrt(G * M_total / (r**3))


def gw_strain(t, r):
    omega = orbital_freq(r)
    h_plus = (4 * G**2 * Mu * M_total) / (c**4 * r) * (2 * np.cos(2 * omega * t))
    h_cross = (4 * G**2 * Mu * M_total) / (c**4 * r) * (2 * np.sin(2 * omega * t))
    return h_plus, h_cross


t = np.arange(0, t_max, dt)
r = np.zeros_like(t)
h_plus = np.zeros_like(t)
h_cross = np.zeros_like(t)

for i in range(len(t)):
    r[i] = max(initial_separation * (1 - t[i]/t_max)**0.25, 2*G*M_total/c**2)
    h_plus[i], h_cross[i] = gw_strain(t[i], r[i])
    

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

# Merge
scatter1, = ax1.plot([], [], 'bo', markersize=10, label='Black Hole 1')
scatter2, = ax1.plot([], [], 'ro', markersize=8, label='Black Hole 2')
ax1.set_xlim(-1.2*initial_separation, 1.2*initial_separation)
ax1.set_ylim(-1.2*initial_separation, 1.2*initial_separation)
ax1.set_aspect('equal')
ax1.set_title('Black Hole Merger')
ax1.legend()


# Waves
line_plus, = ax2.plot([], [], 'b-', label='h_plus')
line_cross, = ax2.plot([], [], 'r-', label='h_cross')
ax2.set_xlim(0, t_max)
ax2.set_ylim(-1.1*max(abs(h_plus.max()), abs(h_cross.max())), 1.1*max(abs(h_plus.max()), abs(h_cross.max())))
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Strain')
ax2.set_title('Gravitational Waves')
ax2.legend()


def animate(i):
    theta = orbital_freq(r[i]) * t[i]
    x1 = r[i] * np.cos(theta) * M2/ M_total
    y1 = r[i] * np.sin(theta) * M2/ M_total
    x2 = -r[i] * np.cos(theta) * M1/ M_total
    y2 = -r[i] * np.sin(theta) * M1/ M_total
    
    scatter1.set_data([x1], [y1])
    scatter2.set_data([x2], [y2])
    
    line_plus.set_data(t[:i], h_plus[:i])
    line_cross.set_data(t[:i], h_cross[:i])
    
    return scatter1, scatter2, line_plus, line_cross


anim = FuncAnimation(fig, animate, frames=len(t), interval=50, blit=True)
plt.tight_layout()
plt.show()

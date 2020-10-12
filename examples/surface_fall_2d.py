import numpy as np
from tqdm import tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt
from gransim import GranSim, animate_matplotlib

plt.style.use('dark_background')
Nparticles = 1000
radii = np.array([.01]*Nparticles)
mass = np.array([.01]*Nparticles)

initial = np.empty([Nparticles,2])
for i in range(Nparticles):
    x = .03*(i%10) + np.random.normal(0, 1e-3)
    y = .1 + .03*(i//10)
    initial[i] = (x,y)

Nsteps = 15000
sim = GranSim(position=initial,
              radii=radii,
              mass=mass,
              young_mod=1.0e6,
              damp_normal=0.01,
              damp_tangent=1.0,
              friction=0.5,
              dt=1e-4)

traj = np.empty([Nsteps,Nparticles,2], dtype=float)

for i in tqdm(range(Nsteps)):
    sim.step()
    traj[i] = sim.position

fig, ax = plt.subplots()
colors = mpl.colors.TABLEAU_COLORS
anim = animate_matplotlib(traj[::100], radii, colors)

ax.axhspan(-.5, 0.8, color='black', hatch='//', zorder=0)
ax.axhline(0, color='gray')
ax.set_ylim([-.5,.8])
ax.axis('off')

plt.show()

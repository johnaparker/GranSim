import numpy as np
from tqdm import tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt
from gransim import GranSim, animate_matplotlib
import vdynamics

plt.style.use('dark_background')
Nparticles = 1000
radii = np.array([.01]*Nparticles)
mass = np.array([.01]*Nparticles)

initial = np.empty([Nparticles,2])
for i in range(Nparticles):
    x = .03*(i%40) + np.random.normal(0, 1e-3)
    y = .1 + .03*(i//40)
    initial[i] = (x,y)

Nsteps = 10000
sim = GranSim(position=initial,
              radii=radii,
              mass=mass,
              young_mod=1.0e6,
              damp_normal=0.01,
              damp_tangent=3.0,
              friction=0.5,
              dt=1e-4)

traj = np.empty([Nsteps,Nparticles,2], dtype=float)

for i in tqdm(range(Nsteps)):
    sim.step()
    traj[i] = sim.position

colors = mpl.colors.TABLEAU_COLORS
vdynamics.animate_2d(traj[::50], radii, colors)

plt.show()

import numpy as np
from tqdm import tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt
from gransim import granular_media_2d, animate_matplotlib
import vdynamics

Nparticles = 3000
radii = np.array([.01]*Nparticles)
mass = np.array([.01]*Nparticles)

initial = np.empty([Nparticles,2])
for i in range(Nparticles):
    x = .03*(i%40) + np.random.normal(0, 1e-3)
    y = .1 + .03*(i//40)
    initial[i] = (x,y)

Nsteps = 10000

sim = granular_media_2d(dt=1e-4)
sim.add_grains(position=initial,
              radii=radii,
              mass=mass,
              young_mod=1.0e6*np.ones(Nparticles),
              damp_normal=0.01*np.ones(Nparticles),
              damp_tangent=3.0*np.ones(Nparticles),
              friction=0.5*np.ones(Nparticles))

sim.add_wall(point=[.6,-.5],
             normal=[1,1])
sim.add_wall(point=[.6,-.5],
             normal=[-1,1])
traj = np.empty([Nsteps,Nparticles,2], dtype=float)

for i in tqdm(range(Nsteps)):
    sim.step()
    traj[i] = sim.position

colors = mpl.colors.TABLEAU_COLORS
vdynamics.animate_2d(traj[::50], radii, colors='gray', background_color='white', edge_color='k', linewidth=.1)

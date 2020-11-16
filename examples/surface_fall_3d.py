import numpy as np
from tqdm import tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt
from gransim import granular_media, animate_matplotlib
import vdynamics
from itertools import cycle

Nparticles = 10000
radii = np.array([.01]*Nparticles)
mass = np.array([.01]*Nparticles)

initial = np.empty([Nparticles,3])
for zp in range(100):
    z = zp*.03
    for i in range(100):
        x = .03*(i%10) + np.random.normal(0, 1e-3)
        y = .3 + .03*(i//10) + np.random.normal(0, 1e-3)

        initial[zp*100 + i] = (x,y,z)

Nsteps = 10000

sim = granular_media(dt=1e-4)
sim.add_grains(position=initial,
              radii=radii,
              mass=mass,
              young_mod=1.0e6*np.ones(Nparticles),
              damp_normal=0.01*np.ones(Nparticles),
              damp_tangent=3*np.ones(Nparticles),
              friction=0.5*np.ones(Nparticles))

sim.add_wall(point=[0,0,0],
             normal=[0,0,1])
traj = np.empty([Nsteps,Nparticles,3], dtype=float)

for i in tqdm(range(Nsteps)):
    sim.step()
    traj[i] = sim.position

colors = mpl.colors.TABLEAU_COLORS
vdynamics.animate_3d(traj[::50], .01*np.ones(traj.shape[1]), colors=colors, background_color='white')

import numpy as np
from tqdm import tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt
from gransim import granular_media_2d, animate_matplotlib
import vdynamics
from itertools import cycle

Nparticles = 1000
radii = np.array([.01]*Nparticles)
mass = np.array([.01]*Nparticles)

initial = np.empty([Nparticles,2])
for i in range(Nparticles):
    x = .03*(i%40) + np.random.normal(0, 1e-3)
    y = .3 + .03*(i//40)
    initial[i] = (x,y)

Nsteps = 20000

sim = granular_media_2d(dt=1e-4)
sim.add_grains(position=initial,
              radii=radii,
              mass=mass,
              young_mod=1.0e6*np.ones(Nparticles),
              damp_normal=0.01*np.ones(Nparticles),
              damp_tangent=3.0*np.ones(Nparticles),
              friction=0.5*np.ones(Nparticles))

pos = [[x,.2] for x in np.linspace(0,1,60)]
Vparticles = len(pos)
sim.add_static_grains(position=np.array(pos),
              radii=np.array([.01]*Vparticles),
              mass=np.array([.01]*Vparticles),
              young_mod=1.0e6*np.ones(Vparticles),
              damp_normal=0.01*np.ones(Vparticles),
              damp_tangent=3.0*np.ones(Vparticles),
              friction=0.5*np.ones(Vparticles))

sim.add_wall(point=[0,0],
             normal=[0,1])
traj = np.empty([Nsteps,Nparticles+Vparticles,2], dtype=float)

for i in tqdm(range(Nsteps)):
    sim.step()
    traj[i] = sim.position

wheel = cycle(mpl.colors.TABLEAU_COLORS)
colors = [next(wheel) for i in range(Nparticles)]
colors.extend(['k']*Vparticles)
vdynamics.animate_2d(traj[::50], .01*np.ones(traj.shape[1]), colors=colors, background_color='white', edge_color='k', linewidth=.1)

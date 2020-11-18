import numpy as np
from tqdm import tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt
from gransim import granular_media_2d, hopper_2d
import vdynamics
from itertools import cycle

Nparticles = 3000
radii = np.array([.01]*Nparticles)
mass = np.array([.01]*Nparticles)

initial = np.empty([Nparticles,2])
for i in range(Nparticles):
    x = -.3 + .03*(i%20) + np.random.normal(0, 1e-3)
    y = 1 + .03*(i//20)
    initial[i] = (x,y)

Nsteps = 20000

sim = granular_media_2d(dt=1e-4)
sim.add_grains(position=initial,
              radii=radii,
              mass=mass,
              young_mod=1.0e6*np.ones(Nparticles),
              damp_normal=0.01*np.ones(Nparticles),
              damp_tangent=10.0*np.ones(Nparticles),
              friction=0.5*np.ones(Nparticles))

pos, _ = hopper_2d([0,.5], .22, 1.7, np.pi/7, .01, closed=True)
Vparticles = len(pos)
sim.add_static_grains(position=np.array(pos),
              radii=np.array([.01]*Vparticles),
              mass=np.array([.01]*Vparticles),
              young_mod=1.0e6*np.ones(Vparticles),
              damp_normal=0.01*np.ones(Vparticles),
              damp_tangent=10.0*np.ones(Vparticles),
              friction=0.5*np.ones(Vparticles))

sim.add_wall(point=[0,0],
             normal=[0,1])
traj = np.empty([Nsteps,Nparticles+Vparticles,2], dtype=float)

for i in tqdm(range(Nsteps)):
    sim.step()
    traj[i] = sim.position

sim2 = granular_media_2d(dt=1e-4)
sim2.add_grains(position=np.copy(sim.position[:Nparticles]),
              radii=radii,
              mass=mass,
              young_mod=1.0e6*np.ones(Nparticles),
              damp_normal=0.01*np.ones(Nparticles),
              damp_tangent=10.0*np.ones(Nparticles),
              friction=0.5*np.ones(Nparticles))

pos, _ = hopper_2d([0,.5], .22, 1.7, np.pi/7, .01, closed=False)
Vparticles = len(pos)
sim2.add_static_grains(position=np.array(pos),
              radii=np.array([.01]*Vparticles),
              mass=np.array([.01]*Vparticles),
              young_mod=1.0e6*np.ones(Vparticles),
              damp_normal=0.01*np.ones(Vparticles),
              damp_tangent=10.0*np.ones(Vparticles),
              friction=0.5*np.ones(Vparticles))

sim2.add_wall(point=[0,-1],
             normal=[0,1])
Nsteps = 40000
traj = np.empty([Nsteps,Nparticles+Vparticles,2], dtype=float)

for i in tqdm(range(Nsteps)):
    sim2.step()
    traj[i] = sim2.position

yi = traj[0,:Nparticles,1]

colors = ['C0' if yi[i] % .24 < .12 else 'C1' for i in range(Nparticles)]
colors.extend(['k']*Vparticles)
vdynamics.animate_2d(traj[::50], .01*np.ones(traj.shape[1]), colors=colors, background_color='white', edge_color='k', linewidth=.1)

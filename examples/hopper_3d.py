import numpy as np
from tqdm import tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt
from gransim import granular_media, hopper
import vdynamics
from itertools import cycle

Nparticles = 10000
radii = np.array([.01]*Nparticles)
mass = np.array([.01]*Nparticles)

initial = np.empty([Nparticles,3])
for zp in range(100):
    z = zp*.03 + 1
    for i in range(100):
        x = -.15 + .03*(i%10) + np.random.normal(0, 1e-3)
        y = -.15 + .03*(i//10) + np.random.normal(0, 1e-3)

        initial[zp*100 + i] = (x,y,z)

Nsteps = 10000

sim = granular_media(dt=1e-4)
sim.add_grains(position=initial,
              radii=radii,
              mass=mass,
              young_mod=1.0e6*np.ones(Nparticles),
              damp_normal=0.001*np.ones(Nparticles),
              damp_tangent=1.0*np.ones(Nparticles),
              friction=0.5*np.ones(Nparticles))

pos, _ = hopper([0,0,.5], .2, 0.6, np.pi/7, .01, closed=True)
Vparticles = len(pos)
sim.add_static_grains(position=np.array(pos),
              radii=np.array([.01]*Vparticles),
              mass=np.array([.01]*Vparticles),
              young_mod=1.0e6*np.ones(Vparticles),
              damp_normal=0.001*np.ones(Vparticles),
              damp_tangent=1.0*np.ones(Vparticles),
              friction=0.5*np.ones(Vparticles))

sim.add_wall(point=[0,0,0],
             normal=[0,0,1])
traj = np.empty([Nsteps,Nparticles+Vparticles,3], dtype=float)

for i in tqdm(range(Nsteps)):
    sim.step()
    traj[i] = sim.position

sim2 = granular_media(dt=1e-4)
sim2.add_grains(position=np.copy(traj[-1,:Nparticles]),
              radii=np.array([.01]*Nparticles),
              mass=np.array([.01]*Nparticles),
              young_mod=1.0e6*np.ones(Nparticles),
              damp_normal=0.001*np.ones(Nparticles),
              damp_tangent=1.0*np.ones(Nparticles),
              friction=0.5*np.ones(Nparticles))

pos, _ = hopper([0,0,.5], .2, 0.6, np.pi/7, .01, closed=False)
Vparticles = len(pos)
sim2.add_static_grains(position=pos,
              radii=np.array([.01]*Vparticles),
              mass=np.array([.01]*Vparticles),
              young_mod=1.0e6*np.ones(Vparticles),
              damp_normal=0.001*np.ones(Vparticles),
              damp_tangent=1.0*np.ones(Vparticles),
              friction=0.5*np.ones(Vparticles))

Nsteps = 40000
sim2.add_wall(point=[0,0,0],
             normal=[0,0,1])
traj = np.empty([Nsteps,Nparticles+Vparticles,3], dtype=float)

for i in tqdm(range(Nsteps)):
    traj[i] = sim2.position
    sim2.step()

zi = traj[0,:Nparticles,2]

# color_wheel = cycle(mpl.colors.TABLEAU_COLORS)
# colors = [next(color_wheel) for i in range(Nparticles)]
colors = ['C0' if zi[i] % .18 < .09 else 'C3' for i in range(Nparticles)]
colors.extend([(0,0,0,.2)]*Vparticles)

vdynamics.animate_3d(traj[::50], .01*np.ones(traj.shape[1]), colors=colors, background_color='white')

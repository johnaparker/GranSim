from cpp import GranSim
import numpy as np
from tqdm import tqdm
from stoked import trajectory, trajectory_animation, circle_patches
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.collections import PatchCollection, CircleCollection
import matplotlib as mpl
from itertools import cycle
from my_pytools.my_matplotlib.animation import save_animation

class UpdatablePatchCollection(PatchCollection):
    def __init__(self, patches, *args, **kwargs):
        self.patches = patches
        PatchCollection.__init__(self, patches, *args, **kwargs)

    def get_paths(self):
        self.set_paths(self.patches)
        return self._paths

Nparticles = 250
radii = np.array([.01]*Nparticles)
mass = np.array([.01]*Nparticles)

initial = np.empty([Nparticles,2])
for i in range(Nparticles):
    x = .03*(i%10) + np.random.normal(0, 1e-3)
    y = .1 + .03*(i//10)
    initial[i] = (x,y)

Nsteps = 10000
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

# fig, ax = plt.subplots()
# anim = trajectory_animation(trajectory(traj)[::10], patches=circle_patches(radii), trail=0)

def update(frame):
    collection.set_offsets(traj[frame])
    for i, circle in enumerate(circles):
        circle.center = traj[frame,i]

fig, ax = plt.subplots()
circles = []
for i,pos in enumerate(traj[0]):
    circles.append(plt.Circle(pos, radii[i]))

collection = UpdatablePatchCollection(circles)
ax.add_collection(collection)
ax.autoscale()
ax.set_aspect('equal')
anim = FuncAnimation(fig, update, range(0,len(traj),100), interval=30)

color_wheel = cycle([f'C{i}' for i in range(10)])
colors = [next(color_wheel) for i in range(Nparticles)]
collection.set_facecolor(colors)

xmax = np.max(traj[...,0]) + radii[0]
xmin = np.min(traj[...,0]) - radii[0]
ymax = np.max(traj[...,1]) + radii[0]
ymin = np.min(traj[...,1]) - radii[0]
ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])

plt.show()

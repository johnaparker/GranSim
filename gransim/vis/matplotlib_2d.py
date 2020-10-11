from itertools import cycle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.collections import PatchCollection, CircleCollection

class UpdatablePatchCollection(PatchCollection):
    def __init__(self, patches, *args, **kwargs):
        self.patches = patches
        PatchCollection.__init__(self, patches, *args, **kwargs)

    def get_paths(self):
        self.set_paths(self.patches)
        return self._paths

def animate_matplotlib(traj, radii, colors=None, ax=None):
    if ax is None:
        ax = plt.gca()
    if colors is None:
        colors = ['C0']

    fig = ax.figure
    Nparticles = traj.shape[1]

    def update(frame):
        collection.set_offsets(traj[frame])
        for i, circle in enumerate(circles):
            circle.center = traj[frame,i]

        return collection,

    circles = []
    for i,pos in enumerate(traj[0]):
        circles.append(plt.Circle(pos, radii[i]))

    collection = UpdatablePatchCollection(circles)
    ax.add_collection(collection)

    ax.set_aspect('equal')
    anim = FuncAnimation(fig, update, len(traj), interval=30)

    color_wheel = cycle(colors)
    colors = [next(color_wheel) for i in range(Nparticles)]
    collection.set_facecolor(colors)

    xmax = np.max(traj[...,0]) + radii[0]
    xmin = np.min(traj[...,0]) - radii[0]
    ymax = np.max(traj[...,1]) + radii[0]
    ymin = np.min(traj[...,1]) - radii[0]
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

    return anim

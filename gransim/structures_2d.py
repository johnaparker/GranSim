import numpy as np

def hopper_2d(position, width, height, angle, grain_radius, closed=False):
    """A hopper container \_/

    Arguments:
        position       (x,y) position of the bottom (center) of the hopper
        width          size of the bottom
        height         height of the hopper
        angle          angle of the side walls relative to vertical
        grain_radius   size of grain positions making up hopper
        closed         if not closed, remove the particles along the bottom edge (default: False)"""
    D = 2*grain_radius

    Nx = int(np.ceil(width/D))
    x1 = np.linspace(-width/2, width/2, Nx) + position[0]
    y1 = np.zeros_like(x1) + position[1]

    L = height/np.cos(angle)
    Ny = int(np.ceil(L/D))
    x2 = np.linspace(0, L, Ny)*np.cos(np.pi/2-angle) + position[0] + width/2
    y2 = np.linspace(0, L, Ny)*np.sin(np.pi/2-angle) + position[1]
    x3 = -np.linspace(0, L, Ny)*np.cos(np.pi/2-angle) + position[0] - width/2
    y3 = np.linspace(0, L, Ny)*np.sin(np.pi/2-angle) + position[1]

    if closed:
        x = np.concatenate((x1, x2, x3))
        y = np.concatenate((y1, y2, y3))
    else:
        x = np.concatenate((x2, x3))
        y = np.concatenate((y2, y3))

    radius = grain_radius*np.ones(len(x))
    return np.c_[x,y], radius

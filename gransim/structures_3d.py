import numpy as np

def hopper(position, width, height, angle, grain_radius, closed=False):
    """A hopper container \_/

    Arguments:
        position       (x,y,z) position of the bottom (center) of the hopper
        width          size of the bottom
        height         height of the hopper
        angle          angle of the side walls relative to vertical
        grain_radius   size of grain positions making up hopper
        closed         if not closed, remove the particles along the bottom edge (default: False)"""
    D = 2*grain_radius

    Nx = int(np.ceil(width/D/2)) + 1
    p1 = [[0,0,0]] 

    for i in range(1,Nx):
        radius = (width - D)/2*i/(Nx-1)
        Np = int(np.ceil(2*np.pi*radius/D))
        for j in range(Np):
            theta = j/Np*2*np.pi
            xval = radius*np.cos(theta)
            yval = radius*np.sin(theta)

            p1.append([xval,yval,0])

    p1 = np.array(p1) + position

    L = height/np.cos(angle)
    Ny = int(np.ceil(L/D))
    p2 = [] 

    for i in range(Ny):
        Li = (L - D/2)*i/(Ny-1)
        z = np.cos(angle)*Li
        radius = np.sin(angle)*Li + width/2
        Np = int(np.ceil(2*np.pi*radius/D))
        for j in range(Np):
            theta = j/Np*2*np.pi
            xval = radius*np.cos(theta)
            yval = radius*np.sin(theta)

            p2.append([xval,yval,z])

    p2 = np.array(p2) + position

    if closed:
        particles = np.concatenate((p1, p2), axis=0)
    else:
        particles = p2

    radii = grain_radius*np.ones(len(position))
    return particles, radii

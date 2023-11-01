import numpy as np
import matplotlib.pyplot as plt


path = "Write the path here"

# plate size, mm
w = h = 85.
# intervals in x-, y- directions, mm
dx = dy = 0.1
# Diffusion constant, mm2.s-1
D = 1.1772e-5

minall, c0 = 0, 5

nx, ny = int(w/dx), int(h/dy)

dx2, dy2 = dx*dx, dy*dy
dt = dx2 * dy2 / (2 * D * (dx2 + dy2))

u0 = minall * np.ones((nx, ny))
u = u0.copy()

# Initial conditions - circle of radius r centred at (cx,cy) (mm)
r, cx, cy = 3.5, 42.5, 42.5
r2 = r**2
for i in range(nx):
    for j in range(ny):
        p2 = (i*dx-cx)**2 + (j*dy-cy)**2
        if p2 < r2:
            u0[i,j] = c0

def do_timestep(u0, u):
    # Propagate with forward-difference in time, central-difference in space
    u[1:-1, 1:-1] = u0[1:-1, 1:-1] + D * dt * (
          (u0[2:, 1:-1] - 2*u0[1:-1, 1:-1] + u0[:-2, 1:-1])/dx2
          + (u0[1:-1, 2:] - 2*u0[1:-1, 1:-1] + u0[1:-1, :-2])/dy2 )

    u0 = u.copy()
    return u0, u

# Number of timesteps
nsteps = int((48*3600)/dt)+2

# Output
mfig = np.linspace(0, nsteps, 49)
for i in range(len(mfig)):
    mfig[i]=int(mfig[i])
print(mfig)

fignum = 0

for m in range(nsteps):

    u0, u = do_timestep(u0, u)
    if int(m*dt/3600)==35 :
        print("here")
        flag=0
        rp = 10
        r2p = rp**2
        for i in range(nx):
            for j in range(ny):
                p2 = (i*dx-cx)**2 + (j*dy-cy)**2
                if int(p2) == r2p:
                    if u0[i,j] < 2 and u0[i,j]!=0:
                        flag=1


        fig, ax = plt.subplots()
        im = ax.imshow(u.copy(), cmap="jet", vmin=minall,vmax=c0)
        ax.set_xlabel("x(mm)")
        ax.set_ylabel("y(mm)")
        ax.set_title('{} hour'.format(int(m*dt/3600)+1))
        cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.set_label(r"Concentration($\frac{\mu g}{ml})$")
        
        plt.savefig(path + "{}.png".format(int(m*dt/3600)), dpi = 400)

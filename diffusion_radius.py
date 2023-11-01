import numpy as np
import matplotlib.pyplot as plt


path = "Write the path here"
# plate size, mm
w = h = 85.
# intervals in x-, y- directions, mm
dx = dy = 0.1
# Diffusion constant, mm2.s-1
D = 1.1772e-5

radius = np.linspace(0,14,140)

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
nsteps = int((48*3600)/dt)
print(nsteps)
# Output
mfig_list = np.linspace(0, nsteps, 48)

mfig = [int(mfig_list[0]), int(mfig_list[5])+1,  int(mfig_list[11])+1,  int(mfig_list[23])+1,  int(mfig_list[35])+1, int(mfig_list[46])+1]
fignum = 0
for m in range(nsteps):
    u0, u = do_timestep(u0, u)
    if m in mfig:
        fignum += 1
        print(m, fignum)
        fig, ax = plt.subplots()
        ax.plot(radius, u[425, 425:565])
        ax.set_xlabel("Radius(mm)", fontsize = 16)
        ax.set_ylabel(r"Concentration($\frac{\mu g}{ml})$" , fontsize = 16)
        ax.set_ylim([0, c0+0.05*c0])
        if int(m*dt/3600)+1==47:
            ax.set_title('{} hour'.format(48) , fontsize = 18)
        else:
            ax.set_title('{} hour'.format(int(m*dt/3600)+1) , fontsize = 18 )

        plt.savefig(path + "{}_his.png".format(int(m*dt/3600)), dpi = 400)

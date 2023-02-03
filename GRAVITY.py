import numpy as np
from numpy import array as arr
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.constants import G
from numpy.random import normal as norm
from time import perf_counter as timer

# Shorthand functions to access x & y coordinates of positions
def x(pos):
    return pos[:,0]

def y(pos):
    return pos[:,1]

# Initial parameters & general settings
N = 300
softening = 30000
velsigma = 1000
possigma = 1e5
umass = 1e22
randmass = False  # makes masses of objects vary according to a normal distribution of σ = umass, μ = 0
realtime = True  # set to False when calculation speed is slower than animation speed, to avoid lag

mode = 's'

logpos, logcom, logt = [], [], []

# Time Settings
time = 0
start = 0
stop = 100
dt = 0.075

# Velocity Modifier Parameters
ivm = 10e-1
svm = 1e-1
ovm = 0.25

# Calculation Time Estimation
baseline = 1.65394  # average time it takes in seconds to render 100 frames of 1000 objects
timeestimate = baseline * (N/1000)**2 * ((stop-start)/dt)/100

# Initialisation of mass & acceleration
posx, posy = norm(0, possigma, N), norm(0, possigma, N)
pos = np.hstack(([[i] for i in posx], [[i] for i in posy]))
acc = arr([arr([0, 0]) for i in range(N)])

if randmass:
    mass = abs(norm(0, umass, N))
else:
    mass = arr([umass for i in range(N)])
    
totmass = sum(mass)

# Data Plot Initialisation
fig, ax = plt.subplots()
ax.grid(True)
ax.set_aspect('equal')
ax.set_xlabel('x-pos (m)')
ax.set_ylabel('y-pos (m)')
ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
plot, = ax.plot(x(pos), y(pos) , 'o', markersize=2)
comassx, comassy = sum(x(pos)*mass)/totmass, sum(y(pos)*mass)/totmass
comassv = arr([comassx, comassy])
comassplot, = ax.plot(comassx, comassy, 'o', markersize=5)

sq1, sq2 = max(abs(x(pos))), max(abs(y(pos)))

if N <= 50:
    if sq1 > sq2:
        sq2 = sq1
    elif sq1 < sq2:
        sq1 = sq2
if N >= 750:
    sq1 *= 2
    sq2 *= 2
    
ax.set_xlim(-3*sq1, 3*sq1)
ax.set_ylim(-3*sq2, 3*sq2)

# Velocity Modes
if mode == 'r':
    vel = np.hstack(([[i] for i in norm(0, velsigma, N)], [[i] for i in norm(0, velsigma, N)]))      # Random 
elif mode == 'o':
    vel = np.hstack(([[-ovm*y(pos)[i]] for i in range(N)], [[ovm*x(pos)[i]] for i in range(N)]))     # Orbit
elif mode == 's':
    vel = np.hstack(([[svm*x(pos)[i]] for i in range(N)], [[svm*y(pos)[i]] for i in range(N)]))      # Supernova
elif mode == 'n':
    vel = np.hstack(([[0] for i in range(N)], [[0] for i in range(N)]))                              # Null
elif mode == 'i':
    vel = np.hstack(([[ivm*-x(pos)[i]] for i in range(N)], [[ivm*-y(pos)[i]] for i in range(N)]))    # Implosion

# Returns array of accelerations for each object, given positions, masses, and softening parameter (prevents division by 0)
def getAcc(pos, mass, softening):
    x = pos[:,0:1]
    y = pos[:,1:2]
    
    dx = x.T - x
    dy = y.T - y
    
    inv_r3 = (dx**2 + dy**2 + softening**2)**-1.5
    
    ax = G * (dx * inv_r3) @ mass
    ay = G * (dy * inv_r3) @ mass
    
    a = np.hstack(([[ax[i]] for i in range(N)], [[ay[i]] for i in range(N)]))
    return a

# Update animation function for real time simulation (fastest simulation time, lowest framerate)
# Recommended for lower values of N & (stop-start)/dt (number of frames)
def rtupdate(t):
    global pos, vel, acc, sq1, sq2, comassx, comassy
    vel += acc * dt/2.0
    pos += vel * dt
    acc = getAcc(pos, mass, softening)
    vel += acc * dt/2.0
    
    comassx, comassy = sum(x(pos)*mass)/totmass, sum(y(pos)*mass)/totmass
    comassplot.set_data(comassx, comassy)
    comassplot,
    
    plot.set_data(x(pos), y(pos))
    plot,

# Update animation function for beforehand calculation simulation (slowest simulation time, highest framerate)
# Recommended for higher values of N & (stop-start)/dt
def bcupdate(t):
    i = int(t)
    plot.set_data(logpos[i][:,0], logpos[i][:,1])
    plot,
    
    comassplot.set_data(logcom[i][0], logcom[i][1])
    comassplot,
    
    simtime.set_text('t = {}s'.format(round(logt[i], 2)))


if realtime:
    # Initialises real time animation
    ani = FuncAnimation(fig, rtupdate, frames=np.arange(0, (stop-start)/dt, 1), interval = 1)
else:
    # Calculates positions and velocities using main loop, then initialises animation.
    print("Time Estimate: {}s".format(timeestimate))
    t0 = timer()
    simtime = ax.text(0.85, 1.05, 't = 0.00s', transform=ax.transAxes)
    while time < stop:
        logpos.append(pos)
        logcom.append(comassv)
        logt.append(time)
        
        vel = vel + acc * dt/2.0
        pos = pos + vel * dt
        acc = getAcc(pos, mass, softening)
        vel = vel + acc * dt/2.0
        
        comassv = arr([sum(x(pos)*mass)/totmass, sum(y(pos)*mass)/totmass])
        
        time += dt
    
    calctime = ax.text(0, 1.1, "Calculation Time: {}s".format(round(timer()-t0, 4)), transform=ax.transAxes)
    Ntext = ax.text(0, 1.05, "N = {}".format(N), transform=ax.transAxes)
    ani = FuncAnimation(fig, bcupdate, frames=np.arange(0, (stop-start)/dt, 1), interval = 15)
         
plt.show()
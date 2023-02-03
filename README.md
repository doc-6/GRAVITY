# GRAVITY
Simulation of gravitationally bound objects (n-body problem)

This Python project can simulate thousands of objects based on the provided initial conditions, such as number of objects, masses, positions, and velocities.
All positions are simulated using a normal distribution function.
Velocities are simulated according to different 'modes'. There are 5 modes currently:
'r' - Random - Objects are given a random velocity each, according to a normal distribution function.
'o' - Orbit - Objects are given a velocity perpendicular to the line connecting them and the origin. The magnitude is proportional to the distance from the origin. Orbital motion is counter-clockwise, as a vector field V = [-y, x] is used for the velocities.
's' - Supernova - Objects are given a velocity pointing outwards from the origin, proportional to the distance from the origin. A vector field V = [x, y] is used to simulate these velocities.
'n' - Null - Objects are given zero initial velocity.
'i' - Implosion - Objects are given a velocity pointing inwards from the origin, proportional to the distance from the origin. A vector field V = [-x, -y] is used to simulate these velocities.

Modes 'o', 's', and 'i' are each given a Velocity Modifier Parameter, which is used to tune the magnitude of the velocities so that they're not too low or high. This also allows for experimentation.

The user can also choose to either run the simulation in real time, or calculate the positions preemptively, and then animate. Real time simulations of over 1000 objects will result in a severe drop in framerate. Hence, by calculating the positions before animating, the animation itself becomes much smoother and more pleasant to view, however do note that preemptive calculations may take a while, especially for larger numbers of objects & frames.

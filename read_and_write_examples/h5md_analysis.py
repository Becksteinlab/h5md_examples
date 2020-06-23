import numpy as np
import pyh5md
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import TRR

# initialize stock TRR example:
u = mda.Universe(TRR)

# open h5 file:
f = pyh5md.File('TRR3.h5', 'r')
# extract particles group
particles = f.particles_group('particles')

# extract particle postitions and velocities with pyh5md.element function 
part_pos = pyh5md.element(particles, 'position').value[:]
part_vel = pyh5md.element(particles, 'velocity').value[:]

# extract observables data from f group 
com = pyh5md.element(f, 'observables/center_of_mass').value[:]
mean_vel = pyh5md.element(f, 'observables/mean_velocity').value[:]
f.close()

#TEST POSITION
print('positions tests:')
for ts in range(len(u.trajectory)):
    diff = u.trajectory[ts].positions - part_pos[ts,:,:]
    if np.all(diff) == np.all(np.zeros(u.atoms.positions.shape)):
        print ('Timestep:',ts, 'Pass')
    else:
        print ('Fail')
        
#TEST CENTER OF MASS
print('center_of_mass tests:')
for ts in range(len(u.trajectory)):
    u.trajectory[ts]
    diff = u.atoms.center_of_geometry() - com[ts,:]
    if np.all(diff) == np.all(np.zeros(u.atoms.center_of_geometry().shape)):
        print ('Timestep:',ts, 'Pass')
    else:
        print ('Fail')
        
#TEST MEAN VELOCITY
print('mean velocity tests:')
for ts in range(len(u.trajectory)):
    u.trajectory[ts]
    diff = u.atoms.velocities.mean() - mean_vel[ts]
    
    if np.all(diff) == np.all(np.zeros(u.atoms.velocities.mean().shape)):
        print ('Timestep:',ts, 'Pass')
    else:
        print ('Fail')
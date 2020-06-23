import pyh5md
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import TRR

# Open an H5MD file to read
with pyh5md.File('TRR3.h5', 'r') as f:
    
    # extract particles group
    particles = f.particles_group('particles')

    # extract particle postitions and velocities with pyh5md.element function 
    part_pos = pyh5md.element(particles, 'position')
    part_vel = pyh5md.element(particles, 'velocity')

    # extract observables data from f group 
    com = pyh5md.element(f, 'observables/center_of_mass')
    mean_vel = pyh5md.element(f, 'observables/mean_velocity')
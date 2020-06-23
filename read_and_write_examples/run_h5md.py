import numpy as np
import pyh5md
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import TRR

# Create the trajectory data
u = mda.Universe(TRR)
positions = u.atoms.positions

# Open a H5MD file to write
with pyh5md.File('TRR3.h5', 'w', author='Edis') as f:
   
    ts = u.trajectory.ts
    step_from = mda.lib.util.Namespace(step=ts.frame, time=ts.time)

    # Add a trajectory group
    particles = f.particles_group('particles')

    # Add the position data element in the trajectory group
    part_pos = pyh5md.element(particles, 'position', step_from=step_from, store='time', 
                              shape=positions.shape, dtype=positions.dtype, time=True)
    # Add the velocity data element into the trajectory group
    part_vel = pyh5md.element(particles, 'velocity', step_from=step_from, store='time', 
                              shape=positions.shape, dtype=positions.dtype, time=True)
    
    # Add observable groups
    # Center of mass:
    com = u.atoms.center_of_geometry()
    obs_com = pyh5md.element(f, 'observables/center_of_mass', step_from=step_from, store='time',
                      shape=com.shape, dtype=com.dtype, time=True)
    # Mean velocity:
    mean_vel = u.atoms.velocities.mean()
    obs_mean_vel = pyh5md.element(f, 'observables/mean_velocity', step_from=step_from, store='time', 
                           shape=mean_vel.shape, dtype=mean_vel.dtype, time=True)
    
    # Define edges to be 3x3 matrix which means the unit cell is triclinic
    edges = u.trajectory.ts.triclinic_dimensions
    # Create the box
    particles.create_box(dimension=3, boundary=['none', 'none', 'none'], store='time',
                    data=edges, step_from=step_from)
    
    
    for ts in u.trajectory: 

        # Append the data to the H5MD file.
        part_pos.append(u.atoms.positions, ts.frame, time=ts.time)
        part_vel.append(u.atoms.velocities, ts.frame, time=ts.time)
        obs_com.append(u.atoms.center_of_geometry(), ts.frame, time=ts.time)
        obs_mean_vel.append(u.atoms.velocities.mean(), ts.frame, time=ts.time)
        particles.box.edges.append(edges, ts.frame, time=ts.time)
        
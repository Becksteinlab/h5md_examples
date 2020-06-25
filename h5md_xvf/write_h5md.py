import numpy as np
import pyh5md
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import TPR_xvf, TRR_xvf

u = mda.Universe(TPR_xvf, TRR_xvf)

# Open a H5MD file to write
with pyh5md.File('H5MD.h5', 'w', author='Edis') as f:
   
    # Add a trajectory group
    atoms = f.particles_group('atoms')

    ts = u.trajectory.ts
    
    # Add the position data element in the trajectory group
    atoms_positions = pyh5md.element(atoms,'positions', store='time', data=u.atoms.positions, time=True)
    atoms_velocities = pyh5md.element(atoms, 'velocities', data=u.atoms.velocities, step_from=atoms_positions, 
                                      store='time', time=True)
    atoms_forces = pyh5md.element(atoms, 'forces', data=u.atoms.forces, step_from=atoms_positions, 
                                  store='time', time=True)
    atoms_masses = pyh5md.element(atoms, 'masses', store='fixed', data=u.atoms.masses)
    atoms_n_atoms = pyh5md.element(atoms, 'n_atoms', store='fixed', data=u.atoms.n_atoms)
    
    # Define edges to be 3x3 matrix which means the unit cell is triclinic
    # Create the box
    atoms.create_box(dimension=3, boundary=['periodic', 'periodic', 'periodic'], store='time',
                    data=u.trajectory.ts.triclinic_dimensions, step_from=atoms_positions)
    
    # Append the data to the H5MD file
    for ts in u.trajectory: 
        atoms.box.edges.append(u.trajectory.ts.triclinic_dimensions, ts.frame, time=ts.time)
        atoms_positions.append(u.atoms.positions, ts.frame, time=ts.time)
        atoms_velocities.append(u.atoms.velocities, ts.frame, time=ts.time)
        atoms_forces.append(u.atoms.forces, ts.frame, time=ts.time)
        atoms_masses.append(u.atoms.masses, ts.frame, time=ts.time)
        atoms_n_atoms.append(u.atoms.n_atoms, ts.frame, time=ts.time)
   
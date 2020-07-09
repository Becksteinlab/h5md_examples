import numpy as np
import h5py
import pyh5md
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import TPR_xvf, TRR_xvf

u = mda.Universe(TPR_xvf, TRR_xvf)

# open h5 file:
f = pyh5md.File('cobrotoxin.h5md', 'r')
# extract particles group
atoms = f.particles_group('trajectory')

# extract datasets from .h5 file with pyh5md.element function 
box = pyh5md.element(atoms['box'], 'edges').value[:]
positions = pyh5md.element(atoms, 'positions').value[:]
velocities = pyh5md.element(atoms, 'velocities').value[:]
forces = pyh5md.element(atoms, 'forces').value[:]
n_atoms = pyh5md.element(atoms, 'n_atoms').value[()]



#TESTS
print('Box Dimensions test (True=Pass False=Fail):')
for ts, x in zip(u.trajectory, box):
    print('Timestep', u.trajectory.ts.frame, np.allclose(u.trajectory.ts.triclinic_dimensions, x))
    
print('Positions test (True=Pass False=Fail):')
for ts, x in zip(u.trajectory, positions):
    print('Timestep', u.trajectory.ts.frame, np.allclose(u.atoms.positions, x))
    
print('Velocities test (True=Pass False=Fail):')
for ts, x in zip(u.trajectory, velocities):
    print('Timestep', u.trajectory.ts.frame, np.allclose(u.atoms.velocities, x))
    
print('Forces test (True=Pass False=Fail):')
for ts, x in zip(u.trajectory, forces):
    print('Timestep', u.trajectory.ts.frame, np.allclose(u.atoms.forces, x))
    

print('N_atoms test:')
if u.atoms.n_atoms == n_atoms:
    print('Pass')
else:
    print('Fail')
    
f.close()
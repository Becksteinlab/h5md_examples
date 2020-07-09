import pyh5md
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import TPR_xvf, TRR_xvf

u = mda.Universe(TPR_xvf, TRR_xvf)

# Open a H5MD file to write
with pyh5md.File('cobrotoxin.h5md', 'w', creator='write_h5md.py') as f:

    # Add a trajectory group into particles group
    trajectory = f.particles_group('trajectory')

    # Add the positions, velocities, forces, masses, n_atoms groups into the trajectory group
    trajectory_positions = pyh5md.element(trajectory,'positions', store='time', data=u.trajectory.ts.positions, time=True)
    trajectory_velocities = pyh5md.element(trajectory, 'velocities', data=u.trajectory.ts.velocities, step_from=trajectory_positions,
                                      store='time', time=True)
    trajectory_forces = pyh5md.element(trajectory, 'forces', data=u.trajectory.ts.forces, step_from=trajectory_positions,
                                  store='time', time=True)
    trajectory_n_atoms = pyh5md.element(trajectory, 'n_atoms', store='fixed', data=u.atoms.n_atoms)
    data_step = pyh5md.element(trajectory, 'data/step', store='time', data=u.trajectory.ts.data['step'])
    data_lambda = pyh5md.element(trajectory, 'data/lambda', store='time', data=u.trajectory.ts.data['lambda'])
    data_dt = pyh5md.element(trajectory, 'data/dt', store='time', data=u.trajectory.ts.data['dt'])
    
    # Data entry is 3x3 matrix which means unitcell is triclinic
    trajectory.create_box(dimension=3, boundary=['periodic', 'periodic', 'periodic'], store='time',
                    data=u.trajectory.ts.triclinic_dimensions, step_from=trajectory_positions)

    # Append the value, step, and time datasets
    for ts in u.trajectory:
        trajectory.box.edges.append(u.trajectory.ts.triclinic_dimensions, ts.frame, time=ts.time)
        trajectory_positions.append(u.trajectory.ts.positions, ts.frame, time=ts.time)
        trajectory_velocities.append(u.trajectory.ts.velocities, ts.frame, time=ts.time)
        trajectory_forces.append(u.trajectory.ts.forces, ts.frame, time=ts.time)
        data_step.append(u.trajectory.ts.data['step'], ts.frame, time=ts.time)
        data_lambda.append(u.trajectory.ts.data['lambda'], ts.frame, time=ts.time)
        data_dt.append(u.trajectory.ts.data['dt'], ts.frame, time=ts.time)
        

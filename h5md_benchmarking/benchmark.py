import MDAnalysis as mda
from MDAnalysis.analysis.rms import rmsd
import time



def benchmark_rmsd(u):
    """Benchmarks rmsd calculation for a given universe"""

    CA = u.select_atoms("protein and name CA")
    x_ref = CA.positions.copy()

    total_io = 0
    total_rmsd = 0
    total_loop = -time.time()
    for i in range(len(u.trajectory)):

        start_io = time.time()
        ts = u.trajectory[i]
        total_io += time.time() - start_io

        start_rmsd = time.time()
        result = rmsd(CA.positions, x_ref, superposition=True)
        total_rmsd += time.time() - start_rmsd

    total_loop += time.time()

    return {"Loop": total_loop,
            "Loop_per_frame": (total_loop/u.trajectory.n_frames),
            "I/O": total_io,
            "I/O_per_frame": (total_io/u.trajectory.n_frames),
            "RMSD": total_rmsd,
            "RMSD_per_frame": (total_rmsd/u.trajectory.n_frames),
            "Overhead": (total_loop - (total_io + total_rmsd)),
            "Overhead_per_frame": ((total_loop - (total_io + total_rmsd))/u.trajectory.n_frames)}


def benchmark_rmsd_n_times(u, n):
    """benchmark a universe 'u' n times"""

    loop_time = []
    loop_time_per_frame = []
    io_time = []
    io_time_per_frame = []
    rmsd_time = []
    rmsd_time_per_frame = []
    overhead_time = []
    overhead_time_per_frame = []

    for i in range(n):
        total_times = benchmark_rmsd(u)
        loop_time.append(total_times["Loop"])
        loop_time_per_frame.append(total_times["Loop_per_frame"])
        io_time.append(total_times["I/O"])
        io_time_per_frame.append(total_times["I/O_per_frame"])
        rmsd_time.append(total_times["RMSD"])
        rmsd_time_per_frame.append(total_times["RMSD_per_frame"])
        overhead_time.append(total_times["Overhead"])
        overhead_time_per_frame.append(total_times["Overhead_per_frame"])

    return {"Loop": loop_time,
            "Loop/frame": loop_time_per_frame,
            "IO": io_time,
            "IO/frame": io_time_per_frame,
            "RMSD": rmsd_time,
            "RMSD/frame": rmsd_time_per_frame,
            "Overhead": overhead_time,
            "Overhead/frame": overhead_time_per_frame}


# H5MD vs pyh5md
import pyh5md
def benchmark_pyh5md_rmsd(filename):
    with pyh5md.File(filename, 'r') as f:

        x_ref = f['particles/trajectory/position/value'][0].copy()
        n_frames = f['particles/trajectory/position/value'].shape[0]

        total_io = 0
        total_rmsd = 0
        total_loop = -time.time()
        for frame in range(n_frames):
            x_ref = f['particles/trajectory/position/value'][0].copy()

            start_io = time.time()
            positions = f['particles/trajectory/position/value'][frame]
            total_io += time.time() - start_io

            start_rmsd = time.time()
            result = rmsd(positions, x_ref, superposition=True)
            total_rmsd += time.time() - start_rmsd

        total_loop += time.time()

    return {"Loop": total_loop,
            "Loop_per_frame": (total_loop/n_frames),
            "I/O": total_io,
            "I/O_per_frame": (total_io/n_frames),
            "RMSD": total_rmsd,
            "RMSD_per_frame": (total_rmsd/n_frames),
            "Overhead": (total_loop - (total_io + total_rmsd)),
            "Overhead_per_frame": ((total_loop - (total_io + total_rmsd))/n_frames)}

def benchmark_pyh5md_rmsd_n_times(filename, n):
    """benchmark an H5MD file n times"""

    loop_time = []
    loop_time_per_frame = []
    io_time = []
    io_time_per_frame = []
    rmsd_time = []
    rmsd_time_per_frame = []
    overhead_time = []
    overhead_time_per_frame = []

    for i in range(n):
        total_times = benchmark_pyh5md_rmsd(filename)
        loop_time.append(total_times["Loop"])
        loop_time_per_frame.append(total_times["Loop_per_frame"])
        io_time.append(total_times["I/O"])
        io_time_per_frame.append(total_times["I/O_per_frame"])
        rmsd_time.append(total_times["RMSD"])
        rmsd_time_per_frame.append(total_times["RMSD_per_frame"])
        overhead_time.append(total_times["Overhead"])
        overhead_time_per_frame.append(total_times["Overhead_per_frame"])

    return {"Loop": loop_time,
            "Loop/frame": loop_time_per_frame,
            "IO": io_time,
            "IO/frame": io_time_per_frame,
            "RMSD": rmsd_time,
            "RMSD/frame": rmsd_time_per_frame,
            "Overhead": overhead_time,
            "Overhead/frame": overhead_time_per_frame}

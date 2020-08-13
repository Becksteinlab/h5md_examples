import MDAnalysis as mda
from benchmark import benchmark_rmsd_n_times, benchmark_pyh5md_rmsd_n_times
import numpy as np
import matplotlib.pyplot as plt


# benchmark different file formats
u_DCD = mda.Universe("cobrotoxin.tpr", "cobrotoxin1000x.dcd")
DCDdict = benchmark_rmsd_n_times(u_DCD, 5)

u_TRR = mda.Universe("cobrotoxin.tpr", "cobrotoxin1000x.trr")
TRRdict = benchmark_rmsd_n_times(u_TRR, 5)

u_XTC = mda.Universe("cobrotoxin.tpr", "cobrotoxin1000x.xtc")
XTCdict = benchmark_rmsd_n_times(u_XTC, 5)

u_H5MD = mda.Universe("cobrotoxin.tpr", "cobrotoxin1000x.h5md")
H5MDdict = benchmark_rmsd_n_times(u_H5MD, 5)


Loop_means = [np.mean(DCDdict["Loop"]),
              np.mean(TRRdict["Loop"]),
              np.mean(XTCdict["Loop"]),
              np.mean(H5MDdict["Loop"])]
Loop_std = [np.std(DCDdict["Loop"]),
            np.std(TRRdict["Loop"]),
            np.std(XTCdict["Loop"]),
            np.std(H5MDdict["Loop"])]

Loop_per_frame_means = [np.mean(DCDdict["Loop/frame"]),
                        np.mean(TRRdict["Loop/frame"]),
                        np.mean(XTCdict["Loop/frame"]),
                        np.mean(H5MDdict["Loop/frame"])]
Loop_per_frame_std = [np.std(DCDdict["Loop/frame"]),
                      np.std(TRRdict["Loop/frame"]),
                      np.std(XTCdict["Loop/frame"]),
                      np.std(H5MDdict["Loop/frame"])]

IO_means = [np.mean(DCDdict["IO"]),
            np.mean(TRRdict["IO"]),
            np.mean(XTCdict["IO"]),
            np.mean(H5MDdict["IO"])]
IO_std = [np.std(DCDdict["IO"]),
          np.std(TRRdict["IO"]),
          np.std(XTCdict["IO"]),
          np.std(H5MDdict["IO"])]

IO_per_frame_means = [np.mean(DCDdict["IO/frame"]),
                      np.mean(TRRdict["IO/frame"]),
                      np.mean(XTCdict["IO/frame"]),
                      np.mean(H5MDdict["IO/frame"])]
IO_per_frame_std = [np.std(DCDdict["IO/frame"]),
                    np.std(TRRdict["IO/frame"]),
                    np.std(XTCdict["IO/frame"]),
                    np.std(H5MDdict["IO/frame"])]

RMSD_compute_means = [np.mean(DCDdict["RMSD"]),
                      np.mean(TRRdict["RMSD"]),
                      np.mean(XTCdict["RMSD"]),
                      np.mean(H5MDdict["RMSD"])]
RMSD_compute_std = [np.std(DCDdict["RMSD"]),
                    np.std(TRRdict["RMSD"]),
                    np.std(XTCdict["RMSD"]),
                    np.std(H5MDdict["RMSD"])]

RMSD_compute_per_frame_means = [np.mean(DCDdict["RMSD/frame"]),
                                np.mean(TRRdict["RMSD/frame"]),
                                np.mean(XTCdict["RMSD/frame"]),
                                np.mean(H5MDdict["RMSD/frame"])]
RMSD_compute_per_frame_std = [np.std(DCDdict["RMSD/frame"]),
                              np.std(TRRdict["RMSD/frame"]),
                              np.std(XTCdict["RMSD/frame"]),
                              np.std(H5MDdict["RMSD/frame"])]

formats = ['DCD', 'TRR', 'XTC', 'H5MD']
x_pos = np.arange(len(formats))

# TOTAL LOOP
fig, ax = plt.subplots()
ax.bar(x_pos, Loop_means, yerr=Loop_std, align='center', alpha=0.8, ecolor='black', capsize=7, edgecolor='black')
ax.set_xticks(x_pos)
ax.set_xticklabels(formats)
ax.text(x_pos[0]-0.2, Loop_means[0] + Loop_means[0], f"{Loop_means[0]:.3f}"+"s", fontsize='medium', fontweight='bold')
ax.text(x_pos[1]-0.2, Loop_means[1] + 0.1*Loop_means[1], f"{Loop_means[1]:.3f}"+"s", fontsize='medium', fontweight='bold')
ax.text(x_pos[2]-0.2, Loop_means[2] + 0.3*Loop_means[2], f"{Loop_means[2]:.3f}"+"s", fontsize='medium', fontweight='bold')
ax.text(x_pos[3]-0.2, Loop_means[3] - 0.2*Loop_means[3], f"{Loop_means[3]:.3f}"+"s", fontsize='medium', fontweight='bold')
ax.set_title("Benchmark Loop Total Time")
ax.set_xlabel("File formats")
ax.set_ylabel("Time (s)")
plt.tight_layout()
plt.savefig("figures/Loop_bench.png")

# TOTAL LOOP/FRAME
fig, ax = plt.subplots()
ax.bar(x_pos, Loop_per_frame_means, yerr=Loop_per_frame_std, align='center', alpha=0.8, ecolor='black', capsize=7, edgecolor='black')
ax.set_xticks(x_pos)
ax.set_xticklabels(formats)
ax.text(x_pos[0]-0.2, Loop_per_frame_means[0] + Loop_per_frame_means[0], f"{Loop_per_frame_means[0]:.4f}"+"s", fontsize='medium', fontweight='bold')
ax.text(x_pos[1]-0.2, Loop_per_frame_means[1] + 0.1*Loop_per_frame_means[1], f"{Loop_per_frame_means[1]:.3f}"+"s", fontsize='medium', fontweight='bold')
ax.text(x_pos[2]-0.2, Loop_per_frame_means[2] + 0.3*Loop_per_frame_means[2], f"{Loop_per_frame_means[2]:.3f}"+"s", fontsize='medium', fontweight='bold')
ax.text(x_pos[3]-0.2, Loop_per_frame_means[3] - 0.2*Loop_per_frame_means[3], f"{Loop_per_frame_means[3]:.3f}"+"s", fontsize='medium', fontweight='bold')
ax.set_title("Benchmark Loop Time per Frame")
ax.set_xlabel("File formats")
ax.set_ylabel("Time (s)")
plt.tight_layout()
plt.savefig("figures/Loop_per_frame_bench.png")

# IO
fig, ax = plt.subplots()
ax.bar(x_pos, IO_means, yerr=IO_std, color='r', align='center', alpha=0.8, ecolor='black', capsize=7, edgecolor='black')
ax.set_xticks(x_pos)
ax.set_xticklabels(formats)
ax.text(x_pos[0]-0.2, IO_means[0] + IO_means[0], f"{IO_means[0]:.3f}"+"s", fontsize='medium', fontweight='bold')
ax.text(x_pos[1]-0.2, IO_means[1] + 0.1*IO_means[1], f"{IO_means[1]:.3f}"+"s", fontsize='medium', fontweight='bold')
ax.text(x_pos[2]-0.2, IO_means[2] + 0.3*IO_means[2], f"{IO_means[2]:.3f}"+"s", fontsize='medium', fontweight='bold')
ax.text(x_pos[3]-0.2, IO_means[3] - 0.2*IO_means[3], f"{IO_means[3]:.3f}"+"s", fontsize='medium', fontweight='bold')
ax.set_title("I/O Total Time")
ax.set_xlabel("File formats")
ax.set_ylabel("Time (s)")
plt.tight_layout()
plt.savefig("figures/IO_bench.png")

# IO/FRAME
fig, ax = plt.subplots()
ax.bar(x_pos, IO_per_frame_means, color='r', yerr=IO_per_frame_std, align='center', alpha=0.8, ecolor='black', capsize=7, edgecolor='black')
ax.set_xticks(x_pos)
ax.set_xticklabels(formats)
ax.text(x_pos[0]-0.2, IO_per_frame_means[0] + IO_per_frame_means[0], f"{IO_per_frame_means[0]:.4f}"+"s", fontsize='medium', fontweight='bold')
ax.text(x_pos[1]-0.2, IO_per_frame_means[1] + 0.1*IO_per_frame_means[1], f"{IO_per_frame_means[1]:.3f}"+"s", fontsize='medium', fontweight='bold')
ax.text(x_pos[2]-0.2, IO_per_frame_means[2] + 0.3*IO_per_frame_means[2], f"{IO_per_frame_means[2]:.3f}"+"s", fontsize='medium', fontweight='bold')
ax.text(x_pos[3]-0.2, IO_per_frame_means[3] - 0.2*IO_per_frame_means[3], f"{IO_per_frame_means[3]:.3f}"+"s", fontsize='medium', fontweight='bold')
ax.set_title("I/O Time per Frame")
ax.set_xlabel("File formats")
ax.set_ylabel("Time (s)")
plt.tight_layout()
plt.savefig("figures/IO_per_frame_bench.png")

# RMSD
fig, ax = plt.subplots()
ax.bar(x_pos, RMSD_compute_means, color='g', yerr=RMSD_compute_std, align='center', alpha=0.8, ecolor='black', capsize=7, edgecolor='black')
ax.set_xticks(x_pos)
ax.set_xticklabels(formats)
ax.text(x_pos[0]-0.2,  0.5*RMSD_compute_means[0], f"{RMSD_compute_means[0]:.3f}"+"s", fontsize='medium', fontweight='bold')
ax.text(x_pos[1]-0.2,  0.5*RMSD_compute_means[1], f"{RMSD_compute_means[1]:.3f}"+"s", fontsize='medium', fontweight='bold')
ax.text(x_pos[2]-0.2,  0.5*RMSD_compute_means[2], f"{RMSD_compute_means[2]:.3f}"+"s", fontsize='medium', fontweight='bold')
ax.text(x_pos[3]-0.2,  0.5*RMSD_compute_means[3], f"{RMSD_compute_means[3]:.3f}"+"s", fontsize='medium', fontweight='bold')
ax.set_title("RMSD Computation Total Time")
ax.set_xlabel("File formats")
ax.set_ylabel("Time (s)")
plt.tight_layout()
plt.savefig("figures/rmsd_bench.png")

# RMSD/FRAME
fig, ax = plt.subplots()
ax.bar(x_pos, RMSD_compute_per_frame_means, color='g', yerr=RMSD_compute_per_frame_std, align='center', alpha=0.8, ecolor='black', capsize=7, edgecolor='black')
ax.set_xticks(x_pos)
ax.set_xticklabels(formats)
ax.text(x_pos[0]-0.25, 0.5*RMSD_compute_per_frame_means[0], f"{RMSD_compute_per_frame_means[0]:.5f}"+"s", fontsize='medium', fontweight='bold')
ax.text(x_pos[1]-0.25, 0.5*RMSD_compute_per_frame_means[1], f"{RMSD_compute_per_frame_means[1]:.5f}"+"s", fontsize='medium', fontweight='bold')
ax.text(x_pos[2]-0.25, 0.5*RMSD_compute_per_frame_means[2], f"{RMSD_compute_per_frame_means[2]:.5f}"+"s", fontsize='medium', fontweight='bold')
ax.text(x_pos[3]-0.25, 0.5*RMSD_compute_per_frame_means[3], f"{RMSD_compute_per_frame_means[3]:.5f}"+"s", fontsize='medium', fontweight='bold')
ax.set_title("RMSD Computation Time per Frame")
ax.set_xlabel("File formats")
ax.set_ylabel("Time (s)")
plt.tight_layout()
plt.savefig("figures/rmsd_per_frame_bench.png")


# MDA vs PYH5MD
u = mda.Universe("cobrotoxin.tpr", "cobrotoxin1000x.h5md")

MDAdict = benchmark_rmsd_n_times(u_H5MD, 5, vs=True)
PYH5MDdict = benchmark_pyh5md_rmsd_n_times("cobrotoxin1000x.h5md", 5)


vs_loop_means = [np.mean(MDAdict["Loop"]),
                 np.mean(PYH5MDdict["Loop"])]
vs_loop_std = [np.std(MDAdict["Loop"]),
               np.std(PYH5MDdict["Loop"])]

vs_loop_per_frame_means = [np.mean(MDAdict["Loop/frame"]),
                          np.mean(PYH5MDdict["Loop/frame"])]
vs_loop_per_frame_std = [np.std(MDAdict["Loop/frame"]),
                         np.std(PYH5MDdict["Loop/frame"])]

vs_IO_means = [np.mean(MDAdict["IO"]),
               np.mean(PYH5MDdict["IO"])]
vs_IO_std = [np.std(MDAdict["IO"]),
             np.std(PYH5MDdict["IO"])]

vs_IO_per_frame_means = [np.mean(MDAdict["IO/frame"]),
                         np.mean(PYH5MDdict["IO/frame"])]
vs_IO_per_frame_std = [np.std(MDAdict["IO/frame"]),
                       np.std(PYH5MDdict["IO/frame"])]

vs_RMSD_compute_means = [np.mean(MDAdict["RMSD"]),
                         np.mean(PYH5MDdict["RMSD"])]
vs_RMSD_compute_std = [np.std(MDAdict["RMSD"]),
                       np.std(PYH5MDdict["RMSD"])]

vs_RMSD_compute_per_frame_means = [np.mean(MDAdict["RMSD/frame"]),
                                   np.mean(PYH5MDdict["RMSD/frame"])]
vs_RMSD_compute_per_frame_std = [np.std(MDAdict["RMSD/frame"]),
                                 np.std(PYH5MDdict["RMSD/frame"])]


benches = ["MDA", "PYH5MD"]
y_pos = np.arange(len(benches))


# LOOP
fig, ax = plt.subplots()
plt.barh(y_pos, vs_loop_means, xerr=vs_loop_std, align='center', alpha=0.8, ecolor='black', capsize=6, edgecolor='black')
plt.yticks(y_pos, benches)
plt.text(vs_loop_means[0] - 0.5*vs_loop_means[0], y_pos[0], f"{vs_loop_means[0]:.3f}" + "s", fontsize='medium', fontweight='bold')
plt.text(vs_loop_means[1] - 0.5*vs_loop_means[1], y_pos[1], f"{vs_loop_means[1]:.3f}" + "s", fontsize='medium', fontweight='bold')
plt.xlabel("Time (s)")
plt.title("MDA vs pyh5md Benchmark Loop Total Time")
plt.tight_layout()
plt.savefig("figures/mda_vs_pyh5md_loop.png")

# LOOP/FRAME
fig, ax = plt.subplots()
plt.barh(y_pos, vs_loop_per_frame_means, xerr=vs_loop_per_frame_std, align='center', alpha=0.8, ecolor='black', capsize=7, edgecolor='black')
plt.yticks(y_pos, benches)
plt.text(vs_loop_per_frame_means[0] - 0.5*vs_loop_per_frame_means[0], y_pos[0], f"{vs_loop_per_frame_means[0]:.3f}"+"s", fontsize='medium', fontweight='bold')
plt.text(vs_loop_per_frame_means[1] - 0.5*vs_loop_per_frame_means[1], y_pos[1], f"{vs_loop_per_frame_means[1]:.3f}"+"s", fontsize='medium', fontweight='bold')
plt.xlabel("Time (s)")
plt.title("MDA vs pyh5md Benchmark Loop Time per Frame")
plt.tight_layout()
plt.savefig("figures/mda_vs_pyh5md_loop_frame.png")

#I/O
fig, ax = plt.subplots()
plt.barh(y_pos, vs_IO_means, color='r', xerr=vs_IO_std, align='center', alpha=0.8, ecolor='black', capsize=7, edgecolor='black')
plt.yticks(y_pos, benches)
plt.text(vs_IO_means[0] - 0.5*vs_IO_means[0], y_pos[0], f"{vs_IO_means[0]:.3f}"+"s", fontsize='medium', fontweight='bold')
plt.text(vs_IO_means[1] - 0.5*vs_IO_means[1], y_pos[1], f"{vs_IO_means[1]:.3f}", fontsize='medium', fontweight='bold')
plt.xlabel("Time (s)")
plt.title("MDA vs pyh5md Total I/O Time")
plt.tight_layout()
plt.savefig("figures/mda_vs_pyh5md_io.png")

#IO/FRAME
fig, ax = plt.subplots()
plt.barh(y_pos, vs_IO_per_frame_means, color='r', xerr=vs_IO_per_frame_std, align='center', alpha=0.8, ecolor='black', capsize=6, edgecolor='black')
plt.yticks(y_pos, benches)
plt.text(vs_IO_per_frame_means[0] - 0.5*vs_IO_per_frame_means[0], y_pos[0], f"{vs_IO_per_frame_means[0]:.3f}"+"s", fontsize='medium', fontweight='bold')
plt.text(vs_IO_per_frame_means[1] - 0.7*vs_IO_per_frame_means[1], y_pos[1], f"{vs_IO_per_frame_means[1]:.3f}"+"s", fontsize='medium', fontweight='bold')
plt.xlabel("Time (s)")
plt.title("MDA vs pyh5md I/O Time per Frame")
plt.tight_layout()
plt.savefig("figures/mda_vs_pyh5md_io_frame.png")

#RMSD
fig, ax = plt.subplots()
plt.barh(y_pos, vs_RMSD_compute_means, color='g', xerr=vs_RMSD_compute_std, align='center', alpha=0.8, ecolor='black', capsize=6, edgecolor='black')
plt.yticks(y_pos, benches)
plt.text(0.5*vs_RMSD_compute_means[0], y_pos[0], f"{vs_RMSD_compute_means[0]:.3f}"+"s", fontsize='medium', fontweight='bold')
plt.text(0.5*vs_RMSD_compute_means[1], y_pos[1], f"{vs_RMSD_compute_means[1]:.3f}"+"s", fontsize='medium', fontweight='bold')
plt.xlabel("Time (s)")
plt.title("MDA vs pyh5md RMSD Computation Total Time")
plt.tight_layout()
plt.savefig("figures/mda_vs_pyh5md_rmsd.png")

fig, ax = plt.subplots()
plt.barh(y_pos, vs_RMSD_compute_per_frame_means, color='g', xerr=vs_RMSD_compute_per_frame_std, align='center', alpha=0.8, ecolor='black', capsize=6, edgecolor='black')
plt.yticks(y_pos, benches)
plt.text(0.5*vs_RMSD_compute_per_frame_means[0], y_pos[0], f"{vs_RMSD_compute_per_frame_means[0]:.5f}"+"s", fontsize='medium', fontweight='bold')
plt.text(0.5*vs_RMSD_compute_per_frame_means[1], y_pos[1], f"{vs_RMSD_compute_per_frame_means[1]:.5f}"+"s", fontsize='medium', fontweight='bold')
plt.xlabel("Time (s)")
plt.title("MDA vs pyh5md RMSD Computation Time per Frame")
plt.tight_layout()
plt.savefig("figures/mda_vs_pyh5md_rmsd_frame.png")

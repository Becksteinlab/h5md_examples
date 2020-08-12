import MDAnalysis as mda
from benchmark import benchmark_rmsd_n_times, benchmark_pyh5md_rmsd_n_times
import numpy as np
import matplotlib.pyplot as plt


# benchmark
u_DCD = mda.Universe("cobrotoxin.tpr", "cobrotoxin100x.dcd")
DCDdict = benchmark_rmsd_n_times(u_DCD, 5)

u_TRR = mda.Universe("cobrotoxin.tpr", "cobrotoxin100x.trr")
TRRdict = benchmark_rmsd_n_times(u_TRR, 5)

u_XTC = mda.Universe("cobrotoxin.tpr", "cobrotoxin100x.xtc")
XTCdict = benchmark_rmsd_n_times(u_XTC, 5)

u_H5MD = mda.Universe("cobrotoxin.tpr", "cobrotoxin100x.h5md")
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
ax.bar(x_pos, Loop_means, yerr=Loop_std, align='center', alpha=0.5, ecolor='black', capsize=10)
ax.set_xticks(x_pos)
ax.set_xticklabels(formats)
ax.set_xlabel("File formats")
ax.set_ylabel("Time (s)")
ax.yaxis.grid(True)
plt.tight_layout()
plt.savefig("figures/Loop_bench.png")

# TOTAL LOOP/FRAME
fig, ax = plt.subplots()
ax.bar(x_pos, Loop_per_frame_means, yerr=Loop_per_frame_std, align='center', alpha=0.5, ecolor='black', capsize=10)
ax.set_xticks(x_pos)
ax.set_xticklabels(formats)
ax.set_xlabel("File formats")
ax.set_ylabel("Time (s)")
ax.yaxis.grid(True)
plt.tight_layout()
plt.savefig("figures/Loop_per_frame_bench.png")

# IO
fig, ax = plt.subplots()
ax.bar(x_pos, IO_means, yerr=IO_std, align='center', alpha=0.5, ecolor='black', capsize=10)
ax.set_xticks(x_pos)
ax.set_xticklabels(formats)
ax.set_xlabel("File formats")
ax.set_ylabel("Time (s)")
ax.yaxis.grid(True)
plt.tight_layout()
plt.savefig("figures/IO_bench.png")

# IO/FRAME
fig, ax = plt.subplots()
ax.bar(x_pos, IO_per_frame_means, yerr=IO_per_frame_std, align='center', alpha=0.5, ecolor='black', capsize=10)
ax.set_xticks(x_pos)
ax.set_xticklabels(formats)
ax.set_xlabel("File formats")
ax.set_ylabel("Time (s)")
ax.yaxis.grid(True)
plt.tight_layout()
plt.savefig("figures/IO_per_frame_bench.png")

# RMSD
fig, ax = plt.subplots()
ax.bar(x_pos, RMSD_compute_means, yerr=RMSD_compute_std, align='center', alpha=0.5, ecolor='black', capsize=10)
ax.set_xticks(x_pos)
ax.set_xticklabels(formats)
ax.set_xlabel("File formats")
ax.set_ylabel("Time (s)")
ax.yaxis.grid(True)
plt.tight_layout()
plt.savefig("figures/rmsd_bench.png")

# RMSD/FRAME
fig, ax = plt.subplots()
ax.bar(x_pos, RMSD_compute_per_frame_means, yerr=RMSD_compute_per_frame_std, align='center', alpha=0.5, ecolor='black', capsize=10)
ax.set_xticks(x_pos)
ax.set_xticklabels(formats)
ax.set_xlabel("File formats")
ax.set_ylabel("Time (s)")
ax.yaxis.grid(True)
plt.tight_layout()
plt.savefig("figures/rmsd_per_frame_bench.png")


# MDA vs PYH5MD
u = mda.Universe("cobrotoxin.tpr", "cobrotoxin100x.h5md")

MDAdict = benchmark_rmsd_n_times(u_H5MD, 5)
PYH5MDdict = benchmark_pyh5md_rmsd_n_times("cobrotoxin.h5md", 5)

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

names = ["MDA", "PYH5MD"]
x = np.arange(len(names))

fig, ax = plt.subplots()
ax.bar(x, vs_loop_means, yerr=vs_loop_std, align='center', alpha=0.5, ecolor='black', capsize=10)
ax.set_xticks(x_pos)
ax.set_xticklabels(names)
ax.set_xlabel("File formats")
ax.set_ylabel("Time (s)")
ax.yaxis.grid(True)
plt.tight_layout()
plt.savefig("figures/mda_vs_pyh5md.png")

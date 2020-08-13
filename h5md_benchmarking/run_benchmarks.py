import MDAnalysis as mda
from benchmark import benchmark_mda_rmsd_n_times, benchmark_pyh5md_rmsd_n_times
import numpy as np
import matplotlib.pyplot as plt


# benchmark different file formats with MDAnalysis
u_DCD = mda.Universe("cobrotoxin.tpr", "cobrotoxin1000x.dcd")
DCDdict = benchmark_mda_rmsd_n_times(u_DCD, 5)

u_TRR = mda.Universe("cobrotoxin.tpr", "cobrotoxin1000x.trr")
TRRdict = benchmark_mda_rmsd_n_times(u_TRR, 5)

u_XTC = mda.Universe("cobrotoxin.tpr", "cobrotoxin1000x.xtc")
XTCdict = benchmark_mda_rmsd_n_times(u_XTC, 5)

u_H5MD = mda.Universe("cobrotoxin.tpr", "cobrotoxin1000x.h5md")
H5MDdict = benchmark_mda_rmsd_n_times(u_H5MD, 5)

# benchmark pyh5md
PYH5MDdict = benchmark_pyh5md_rmsd_n_times("cobrotoxin1000x.h5md", 5)


Loop_means = [np.mean(DCDdict["Loop"]),
              np.mean(TRRdict["Loop"]),
              np.mean(XTCdict["Loop"]),
              np.mean(H5MDdict["Loop"]),
              np.mean(PYH5MDdict["Loop"])]
Loop_std = [np.std(DCDdict["Loop"]),
            np.std(TRRdict["Loop"]),
            np.std(XTCdict["Loop"]),
            np.std(H5MDdict["Loop"]),
            np.std(PYH5MDdict["Loop"])]

Loop_per_frame_means = [np.mean(DCDdict["Loop/frame"]),
                        np.mean(TRRdict["Loop/frame"]),
                        np.mean(XTCdict["Loop/frame"]),
                        np.mean(H5MDdict["Loop/frame"]),
                        np.mean(PYH5MDdict["Loop/frame"])]
Loop_per_frame_std = [np.std(DCDdict["Loop/frame"]),
                      np.std(TRRdict["Loop/frame"]),
                      np.std(XTCdict["Loop/frame"]),
                      np.std(H5MDdict["Loop/frame"]),
                      np.std(PYH5MDdict["Loop/frame"])]

IO_means = [np.mean(DCDdict["IO"]),
            np.mean(TRRdict["IO"]),
            np.mean(XTCdict["IO"]),
            np.mean(H5MDdict["IO"]),
            np.mean(PYH5MDdict["IO"])]
IO_std = [np.std(DCDdict["IO"]),
          np.std(TRRdict["IO"]),
          np.std(XTCdict["IO"]),
          np.std(H5MDdict["IO"]),
          np.std(PYH5MDdict["IO"])]

IO_per_frame_means = [np.mean(DCDdict["IO/frame"]),
                      np.mean(TRRdict["IO/frame"]),
                      np.mean(XTCdict["IO/frame"]),
                      np.mean(H5MDdict["IO/frame"]),
                      np.mean(PYH5MDdict["IO/frame"])]
IO_per_frame_std = [np.std(DCDdict["IO/frame"]),
                    np.std(TRRdict["IO/frame"]),
                    np.std(XTCdict["IO/frame"]),
                    np.std(H5MDdict["IO/frame"]),
                    np.std(PYH5MDdict["IO/frame"])]

RMSD_compute_means = [np.mean(DCDdict["RMSD"]),
                      np.mean(TRRdict["RMSD"]),
                      np.mean(XTCdict["RMSD"]),
                      np.mean(H5MDdict["RMSD"]),
                      np.mean(PYH5MDdict["RMSD"])]
RMSD_compute_std = [np.std(DCDdict["RMSD"]),
                    np.std(TRRdict["RMSD"]),
                    np.std(XTCdict["RMSD"]),
                    np.std(H5MDdict["RMSD"]),
                    np.std(PYH5MDdict["RMSD"])]

RMSD_compute_per_frame_means = [np.mean(DCDdict["RMSD/frame"]),
                                np.mean(TRRdict["RMSD/frame"]),
                                np.mean(XTCdict["RMSD/frame"]),
                                np.mean(H5MDdict["RMSD/frame"]),
                                np.mean(PYH5MDdict["RMSD/frame"])]
RMSD_compute_per_frame_std = [np.std(DCDdict["RMSD/frame"]),
                              np.std(TRRdict["RMSD/frame"]),
                              np.std(XTCdict["RMSD/frame"]),
                              np.std(H5MDdict["RMSD/frame"]),
                              np.std(PYH5MDdict["RMSD/frame"])]

# write a csv file
csv_list = [
            [" ", "DCD", "TRR", "XTC", "H5MD w/ MDA", "H5MD w/ pyh5md"],
            ["Total loop means"] + [str(value) for value in Loop_means],
            ["Total loop stds"] + [str(value) for value in Loop_std],
            ["Loop/frame means"] + [str(value) for value in Loop_per_frame_means],
            ["Loop/frame stds"] + [str(value) for value in Loop_per_frame_std],
            ["Total IO means"] + [str(value) for value in IO_means],
            ["Total IO stds"] + [str(value) for value in IO_std],
            ["IO/frame means"] + [str(value) for value in IO_per_frame_means],
            ["IO/frame stds"] + [str(value) for value in IO_per_frame_std],
            ["Total RMSD means"] + [str(value) for value in RMSD_compute_means],
            ["Total RMSD stds"] + [str(value) for value in RMSD_compute_std],
            ["RMSD/frame means"] + [str(value) for value in RMSD_compute_per_frame_means],
            ["RMSD/frame stds"] + [str(value) for value in RMSD_compute_per_frame_std]
           ]

import csv
with open('data.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(csv_list)


formats = ['DCD', 'TRR', 'XTC', 'H5MD\nwith MDA', 'H5MD\nwith pyh5md']
x_pos = np.arange(len(formats))

# TOTAL LOOP
fig, ax = plt.subplots()
ax.bar(x_pos, Loop_means, yerr=Loop_std, align='center', alpha=0.8, ecolor='black', capsize=7, edgecolor='black')
ax.set_facecolor("xkcd:grey")
ax.set_xticks(x_pos)
ax.set_xticklabels(formats)
ax.set_title("Benchmark Loop Total Time")
ax.set_xlabel("File Formats")
ax.set_ylabel("Time (s)")
plt.tight_layout()
plt.savefig("figures/Loop_bench.png")

# TOTAL LOOP/FRAME
fig, ax = plt.subplots()
ax.bar(x_pos, Loop_per_frame_means, yerr=Loop_per_frame_std, align='center', alpha=0.8, ecolor='black', capsize=7, edgecolor='black')
ax.set_facecolor("xkcd:grey")
ax.set_xticks(x_pos)
ax.set_xticklabels(formats)
ax.set_title("Benchmark Loop Time per Frame")
ax.set_xlabel("File Formats")
ax.set_ylabel("Time (s)")
plt.tight_layout()
plt.savefig("figures/Loop_per_frame_bench.png")

# IO
fig, ax = plt.subplots()
ax.bar(x_pos, IO_means, yerr=IO_std, color='r', align='center', alpha=0.8, ecolor='black', capsize=7, edgecolor='black')
ax.set_facecolor("xkcd:grey")
ax.set_xticks(x_pos)
ax.set_xticklabels(formats)
ax.set_title("I/O Total Time")
ax.set_xlabel("File Formats")
ax.set_ylabel("Time (s)")
plt.tight_layout()
plt.savefig("figures/IO_bench.png")

# IO/FRAME
fig, ax = plt.subplots()
ax.bar(x_pos, IO_per_frame_means, color='r', yerr=IO_per_frame_std, align='center', alpha=0.8, ecolor='black', capsize=7, edgecolor='black')
ax.set_facecolor("xkcd:grey")
ax.set_xticks(x_pos)
ax.set_xticklabels(formats)
ax.set_title("I/O Time per Frame")
ax.set_xlabel("File Formats")
ax.set_ylabel("Time (s)")
plt.tight_layout()
plt.savefig("figures/IO_per_frame_bench.png")

# RMSD
fig, ax = plt.subplots()
ax.bar(x_pos, RMSD_compute_means, color='g', yerr=RMSD_compute_std, align='center', alpha=0.8, ecolor='black', capsize=7, edgecolor='black')
ax.set_facecolor("xkcd:grey")
ax.set_xticks(x_pos)
ax.set_xticklabels(formats)
ax.set_title("RMSD Computation Total Time")
ax.set_xlabel("File Formats")
ax.set_ylabel("Time (s)")
plt.tight_layout()
plt.savefig("figures/rmsd_bench.png")

# RMSD/FRAME
fig, ax = plt.subplots()
ax.bar(x_pos, RMSD_compute_per_frame_means, color='g', yerr=RMSD_compute_per_frame_std, align='center', alpha=0.8, ecolor='black', capsize=7, edgecolor='black')
ax.set_facecolor("xkcd:grey")
ax.set_xticks(x_pos)
ax.set_xticklabels(formats)
ax.set_title("RMSD Computation Time per Frame")
ax.set_xlabel("File Formats")
ax.set_ylabel("Time (s)")
plt.tight_layout()
plt.savefig("figures/rmsd_per_frame_bench.png")

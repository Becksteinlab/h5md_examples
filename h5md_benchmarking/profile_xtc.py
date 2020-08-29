from benchmark import benchmark_mda_rmsd_n_times
import MDAnalysis as mda

u = mda.Universe("testfiles/cobrotoxin.tpr", "testfiles/cobrotoxin1000x.xtc")

benchmark_mda_rmsd_n_times(u, 5)
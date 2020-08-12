import MDAnalysis as mda

u1 = mda.Universe("cobrotoxin.tpr", 10*["cobrotoxin100x.trr"], format="TRR")
with mda.Writer("cobrotoxin1000x.trr", u1.atoms.n_atoms, format="TRR") as W:
    for ts in u1.trajectory:
        W.write(u1)

u2 = mda.Universe("cobrotoxin.tpr", 10*["cobrotoxin100x.dcd"], format="DCD")
with mda.Writer("cobrotoxin1000x.dcd", u2.atoms.n_atoms, format="DCD") as W:
    for ts in u2.trajectory:
        W.write(u2)

u3 = mda.Universe("cobrotoxin.tpr", 10*["cobrotoxin100x.xtc"], format="XTC")
with mda.Writer("cobrotoxin1000x.xtc", u3.atoms.n_atoms, format="XTC") as W:
    for ts in u3.trajectory:
        W.write(u3)

import MDAnalysis as mda

u1 = mda.Universe("cobrotoxin.tpr", 100*["cobrotoxin.trr"], format="TRR")
with mda.Writer("cobrotoxin100x.trr", u1.atoms.n_atoms, format='TRR') as W:
    for ts in u1.trajectory:
        W.write(u1)

u2 = mda.Universe("cobrotoxin.tpr", 100*["cobrotoxin.dcd"], format="DCD")
with mda.Writer("cobrotoxin100x.dcd", u2.atoms.n_atoms, format="DCD") as W:
    for ts in u2.trajectory:
        W.write(u2)

u3 = mda.Universe("cobrotoxin.tpr", 100*["cobrotoxin.xtc"], format="XTC")
with mda.Writer("cobrotoxin100x.xtc", u3.atoms.n_atoms, format="XTC") as W:
    for ts in u3.trajectory:
        W.write(u3)

import pytest
import pyh5md
import MDAnalysis as mda
from numpy.testing import assert_almost_equal, assert_array_equal
from MDAnalysisTests.datafiles import (H5MD_xvf, TPR_xvf)


"""simple test to make sure pyh5md is benchmarking same array as mda"""

def test_pyh5md_bench():
    u = mda.Universe(TPR_xvf, H5MD_xvf)
    ca = u.select_atoms("protein and name CA")
    indices = ca.indices
    with pyh5md.File(H5MD_xvf, 'r') as f:
        trajectory = f.particles_group('trajectory')
        for i in range(len(u.trajectory)):
            u.trajectory[i]
            assert_array_equal(ca.positions,
                               pyh5md.element(trajectory, 'position').value[i, indices, :])


import ray
ray.init(num_cpus=24)

import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"

from matplotlib import pyplot as plt
import scipy
import wannierberri as wb
import numpy as np

path_data = "./"
assert os.path.exists(path_data), f"Path {path_data} does not exist"

import irrep, spglib

wb.welcome()

from irrep.bandstructure import BandStructure
bandstructure = BandStructure(
    code='espresso', 
    prefix=os.path.join(path_data, "pwscf"), 
    include_TR=True,
) # only need prefix.save
spacegroup = bandstructure.spacegroup
spacegroup.show() # uncomment to see the detected symmetries of the system

wandata = wb.w90files.WannierData.from_npz(os.path.join(path_data,"w90data"),
                                           files=["eig","symmetrizer", "mmn"])

from wannierberri.symmetry.projections import Projection, ProjectionsSet
# now set the transformations of WFs. Make sure, the projections are consistent with the amn file
proj_Li = Projection(position_num = [[0,0,0]], orbital='sp3d2', spacegroup=spacegroup)
proj_F = Projection(position_num = [[0.5,0.5,0.5]], orbital='p', spacegroup=spacegroup)
projections_set = ProjectionsSet(projections=[proj_Li, proj_F])
wandata.set_projections(projections_set, bandstructure=bandstructure)

_ = wb.wannierise(
                wandata=wandata,
                froz_min=-10,
                froz_max=19,
                outer_min=-10,
                outer_max=200,
                print_progress_every=20,
                num_iter=20000,
                conv_tol=1e-12,
                localise=True,
                sitesym=True,
                )


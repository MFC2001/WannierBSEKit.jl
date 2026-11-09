
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

wandata = wb.w90files.WannierData.from_bandstructure(
    bandstructure=bandstructure,
    files = ["symmetrizer"],
    irreducible=False, # set irreducible=True to only use the irreducible k-points
    ecut_sym=300,
)
eig = wb.w90files.eig.EIG.from_w90_file(os.path.join(path_data,"wannier90"))
wandata.set_file("eig", eig, overwrite=True)
mmn = wb.w90files.mmn.MMN.from_w90_file(os.path.join(path_data,"wannier90"), bkvec=wandata.bkvec)
wandata.set_file("mmn", mmn, overwrite=True)
wandata.to_npz(os.path.join(path_data,"w90data"))
# wandata = wb.w90files.WannierData.from_npz(os.path.join(path_data,"w90data"),
#                                                 files=["eig","symmetrizer", "mmn"])



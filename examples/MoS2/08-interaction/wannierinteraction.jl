
using JLD2
using WannierBSEKit

intfolder = "./"
period = [1, 1, 0]

cell = read("./POSCAR", POSCAR; period)
# Mo:1:10, S1:11:16, S2:17:22
atom_orb_index = [
    collect(1:10),
    collect(11:16),
    collect(17:22),
]
#n: 2.41425; nn:3.186; nnn:3.9974; nnnn:5.11173;
pair = atom2orbpair(cell, 4.2, atom_orb_index)

chkdata = read("../07-wannier/wannier90.chk", wannier90_chk)
wave_ik = read("../07-wannier/pwscf.save/", chkdata, QE_wannier; period)
eps = read("../05-epsilon/epsmat.h5", epsmat_bgw)

ρk = Rhok(eps, wave_ik; upindex=collect(1:2:22), dnindex=collect(2:2:22))
int(eps, ρk; pair, pair_uplo=:U, outfolder=intfolder)

kgrid = ρk.qgrid
lattice = chkdata.lattice
orblocat = chkdata.wannier_centres

nk = 15
h = 2.76
rcut = (nk - 1) * h / 2 - 1

outfolder = intfolder
U_bare = load(joinpath(intfolder, "U_bare.jld2"))
U_bare_mirror = MirrorCorrection(:keldysh_gauss_optim, U_bare["interaction"], kgrid, lattice, orblocat; period,
    rcut, outfile=joinpath(outfolder, "U_bare_mirror.jld2"))

U_1 = load(joinpath(intfolder, "U_1.jld2"))
MirrorCorrection(:keldysh_gauss, U_1["interaction"], kgrid, lattice, orblocat; period,
    rcut, r₀_optim=true, outfile=joinpath(outfolder, "U_1_mirror.jld2"))

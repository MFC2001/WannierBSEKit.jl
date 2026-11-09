
using JLD2
using WannierBSEKit

cell = read("./POSCAR_orb", POSCAR)
atom_orb_index = [
    Int[], [1], [2], [3], [4], [5], [6],
    [7, 8, 9],
]
pair = atom2orbpair(cell, 2.4, atom_orb_index)

intfolder = "./"
period = [1, 1, 1]

chkdata = read("../07-wannierberri/w90data.chk.npz", wannierberri_chk)
wave_ik = read("../07-wannierberri/pwscf.save/", chkdata, QE_wannier; period)
eps = read("../05-epsilon/epsmat.h5", epsmat_bgw)

ρk = Rhok(eps, wave_ik)
int(eps, ρk; pair, pair_uplo=:U, pair_atol=1e-6, outfolder = intfolder)

orblocat = chkdata.wannier_centres

nk = 15
h = 2.32
rcut = (nk - 1) * h / 2

U_bare = load(joinpath(intfolder, "U_bare.jld2"))
U_bare_mirror = MirrorCorrection(:gauss_optim, U_bare["interaction"], [nk, nk, nk], cell.lattice, orblocat; rcut, rcutmin = 6, 
        outfile = joinpath(intfolder, "U_bare_mirror.jld2"))

U_1 = load(joinpath(intfolder, "U_1.jld2"))
MirrorCorrection(:gauss, U_1["interaction"], [nk, nk, nk], cell.lattice, orblocat; αmat = U_bare_mirror.αmat, ϵ = real(U_1["ϵ"]), 
        outfile = joinpath(intfolder, "U_1_mirror.jld2"))

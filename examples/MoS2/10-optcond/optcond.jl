
using JLD2
using LinearAlgebra
using WannierBSEKit

period = [1, 1, 0]

cell = read("../08-interaction/POSCAR", POSCAR; period)
# Mo:1:10, S1:11:16, S2:17:22
atom_orb_index = [
    collect(1:10),
    collect(11:16),
    collect(17:22),
]
#n: 2.41425; nn:3.186; nnn:3.9974; nnnn:5.11173;
pair = atom2orbpair(cell, 2, atom_orb_index)

nk = 36
kgrid = MonkhorstPack([nk, nk, 1])
outfolder = "./"

intfolder = "../08-interaction"
Kᵈ_U = load(joinpath(intfolder, "U_1_mirror.jld2"))["Uwithlr"]
Kˣ_U = load(joinpath(intfolder, "U_bare_mirror.jld2"))["Uwithlr"]
Kᵈ_pair = load(joinpath(intfolder, "pair_1.jld2"))["interaction"]
Kˣ_pair = load(joinpath(intfolder, "pair_bare.jld2"))["interaction"]
Kᵈ = PairInterAction(Kᵈ_U, Kᵈ_pair; pair)
Kˣ = PairInterAction(Kˣ_U, Kˣ_pair; pair)

TB = TightBindModel(; cell,
	hops = "../07-wannier/wannier90_hr.dat",
	orbital = "../07-wannier/wannier90_centres.xyz",
)

# DFTgap:1.6eV; GWgap: 2.8eV
bse = BSE(kgrid, TB, Kᵈ, Kˣ; channel=1, v = 13:14, c = 15:16, scissor = 1.2)
optdata_noeh = OpticalTransitionData(bse)
H = bse([0, 0, 0])
optdata = OpticalTransitionData(bse, eigen!(H.H))

rlattice = reciprocal(bse.lattice)
jldsave("./optdata.jld2"; rlattice, optdata_noeh, optdata)

ωgrid = 1.5:0.001:4
ê = ReducedCoordinates([1, 0, 0])
ê = rlattice * ê
optcond(optdata_noeh, ωgrid, ê; fwhm=0.08, outfile="./optcond_noeh.dat")
optcond(optdata, ωgrid, ê; fwhm=0.08, outfile="./optcond_x.dat")


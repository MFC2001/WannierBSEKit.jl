
using JLD2
using LinearAlgebra
using WannierBSEKit

TB = TightBindModel(; 
    cell="../08-interaction/POSCAR",
    hops="../07-wannierberri/wannier90_hr.dat",
    orbital="../07-wannierberri/wannier90_centres.xyz",
    period=[1, 1, 1],
)

kgrid = MonkhorstPack([10, 10, 10]; kshift = [1 // 2, 1 // 2, 1 // 2])

intfolder = "../08-interaction/"
Kᵈ_U = load(joinpath(intfolder, "U_1_mirror.jld2"))["Uwithlr"]
Kˣ_U = load(joinpath(intfolder, "U_bare_mirror.jld2"))["Uwithlr"]
Kᵈ_pair = load(joinpath(intfolder, "pair_1.jld2"))["interaction"]
Kˣ_pair = load(joinpath(intfolder, "pair_bare.jld2"))["interaction"]
Kᵈ = PairInterAction(Kᵈ_U, Kᵈ_pair)
Kˣ = PairInterAction(Kˣ_U, Kˣ_pair)
#GW gap ≈ 14.33
bse = BSE(kgrid, TB, Kᵈ, Kˣ; v = [1, 2, 3], c = [4])

optdata_noeh = OpticalTransitionData(bse)
rlattice = reciprocal(bse.lattice)
qz = rlattice \ CartesianCoordinates(0, 0, 1)
Hs = bse(qz; isΓ=true).Hs
optdata = OpticalTransitionData(bse, eigen!(Hs))

jldsave("./optdata.jld2"; optdata_noeh, optdata)

ωgrid = 12:0.001:16
# excitonic q is in z-direction, and the polarization direction is perpendicular to the z direction.
ê = CartesianCoordinates(1, 0, 0)
absorption(optdata_noeh, ωgrid, ê; fwhm=0.2, outfile="./absorption_noeh.dat")
absorption(optdata, ωgrid, ê; fwhm=0.2, outfile="./absorption.dat")

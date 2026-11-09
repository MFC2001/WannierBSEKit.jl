
using JLD2
using KrylovKit
using WannierBSEKit

period = [1, 1, 1]
TB = TightBindModel(; 
    cell="../08-interaction/POSCAR",
    hops="../07-wannierberri/wannier90_hr.dat",
    orbital="../07-wannierberri/wannier90_centres.xyz",
    period,
)

kgrid = MonkhorstPack([10, 10, 10]; kshift=[1 // 2, 1 // 2, 1 // 2])

intfolder = "../08-interaction/"
Kᵈ_U = load(joinpath(intfolder, "U_1_mirror.jld2"))["Uwithlr"]
Kˣ_U = load(joinpath(intfolder, "U_bare_mirror.jld2"))["Uwithlr"]
Kᵈ_pair = load(joinpath(intfolder, "pair_1.jld2"))["interaction"]
Kˣ_pair = load(joinpath(intfolder, "pair_bare.jld2"))["interaction"]
Kᵈ = PairInterAction(Kᵈ_U, Kᵈ_pair)
Kˣ = PairInterAction(Kˣ_U, Kˣ_pair)
#GW gap ≈ 14.33
bse = BSE(kgrid, TB, Kᵈ, Kˣ; v = [1, 2, 3], c = [4])

kline = [
    "L" => [0.5, 0.5, 0.5],
    "G" => [0, 0, 0],
    "X" => [0.5, 0, 0.5],
    "W" => [0.5, 0.25, 0.75],
    "L" => [0.5, 0.5, 0.5],
]
kline = Kline(kline, 30, bse)

eigsolveconfigure!(; verbosity=0, howmany=50, which=EigSorter(abs))
band = BAND(kline, bse)

write("./bsebandt.dat", band.bandt, BAND_dat; kline)
write("./bsebands.dat", band.bands, BAND_dat; kline)

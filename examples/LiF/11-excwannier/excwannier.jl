
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
eigsolveconfigure!(; verbosity=0, howmany=3, which=EigSorter(abs))


qgrid = RedKgrid(MonkhorstPack([7, 7, 7]))
write("./wannier90.win", read("../08-interaction/POSCAR", POSCAR), wannier90_win;
    num_iter=2000,
    num_bands=3,
    dis_win_min=0,
    dis_win_max=20,
    dis_froz_min=0,
    dis_froz_max=20,
    num_wann=3,
    mp_grid=[7, 7, 7],
    kpoints=qgrid,
)
(nnkpts, weights) = kgrid_shell_automatic(qgrid, reciprocal(bse.lattice))
result = BAND(qgrid, bse; vector=true, wfctype=:BlochPeriodic)
band = result.bandt
uband = result.ubandt

# jldsave("./band.jld2"; qgrid, nnkpts, band, uband)

M_0 = bse(:mmn, uband, nnkpts; λ=0)
M_half = bse(:mmn, uband, nnkpts; λ=1/2)
M_1 = bse(:mmn, uband, nnkpts; λ=1)

write("./wannier90.eig", uband, wannier90_eig)
write("./0/wannier90.mmn", M_0, wannier90_mmn; nnkpts)
write("./half/wannier90.mmn", M_half, wannier90_mmn; nnkpts)
write("./1/wannier90.mmn", M_1, wannier90_mmn; nnkpts)


# amn
ihs = [7, 8, 9]
Rₕ = ReducedCoordinates(0, 0, 0)
iex = [1, 2]
Rₑx = [ReducedCoordinates(0, 1, 0), ReducedCoordinates(1, 0, 1)]
iey = [3, 4]
Rₑy = [ReducedCoordinates(0, 1, 1), ReducedCoordinates(1, 0, 0)]
iez = [5, 6]
Rₑz = [ReducedCoordinates(1, 1, 0), ReducedCoordinates(0, 0, 1)]
nq = length(qgrid)
Amn = Array{ComplexF64}(undef, 3, 3, nq)
Threads.@threads for iq in 1:nq
    q = qgrid[iq]
    for ib in 1:3
        AiRjR = realspace(bse, q, band[iq].vectors[:, ib])
        # g₁
        result = 0
        for ih in ihs, iie in 1:2
            result += AiRjR(ih, Rₕ, iex[iie], Rₑx[iie])
        end
        Amn[ib, 1, iq] = conj(result) / sqrt(6)
        # g₂
        result = 0
        for ih in ihs, iie in 1:2
            result += AiRjR(ih, Rₕ, iey[iie], Rₑy[iie])
        end
        Amn[ib, 2, iq] = conj(result) / sqrt(6)
        # g₃
        result = 0
        for ih in ihs, iie in 1:2
            result += AiRjR(ih, Rₕ, iez[iie], Rₑz[iie])
        end
        Amn[ib, 3, iq] = conj(result) / sqrt(6)
    end
end
write("./wannier90.amn", Amn, wannier90_amn)

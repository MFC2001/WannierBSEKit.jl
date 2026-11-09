
using JLD2
using WannierBSEKit

cell = read("../08-interaction/POSCAR", POSCAR; period=[1, 1, 0])
kgrid = MonkhorstPack([36, 36, 1])
# postprocessinging don't need interaction
Kᵈ = PairInterAction(22)
Kˣ = PairInterAction(22)
TB = TightBindModel(; cell,
    hops="../07-wannier/wannier90_hr.dat",
    orbital="../07-wannier/wannier90_centres.xyz",
)
bse = BSE(kgrid, TB, Kᵈ, Kˣ; channel=1, v=13:14, c=15:16, scissor=1.2)

# calculate Ak on kgrid of [36, 36, 1]
q = ReducedCoordinates(0, 0, 0)
# H = bse(q).H
# wfc = eigen!(H)
# jldsave("./wave_G.jld2"; wave = wfc.vectors[:, 1])
wfc = load("./wave_G.jld2")["wave"]

# interpolate Ak on kgrid of [150, 150, 1]
Ak = bse(:interpolation, q, wfc)
kgrid = RedKgrid(MonkhorstPack([150, 150, 1]))
result = map(kgrid) do k
    Threads.@spawn Ak(k)
end
result = fetch.(result)

# obtain more kpoints to obtain the data in FBZ
rlattice = reciprocal(bse.lattice)
(Glist, _) = gridindex(:WignerSeitz, rlattice, (4, 4, 1))
# kpoints
nk = length(kgrid)
nG = length(Glist)
kx = Matrix{Float64}(undef, nk, nG)
ky = Matrix{Float64}(undef, nk, nG)
for (iG, GG) in enumerate(Glist), (ik, kk) in enumerate(kgrid)
    k = rlattice * (kk + GG)
    kx[ik, iG] = k[1]
    ky[ik, iG] = k[2]
end
kx = reshape(kx, :)
ky = reshape(ky, :)

# obtain the contributions of electron or hole
values_v = Matrix{Float64}(undef, nk, nG)
values_c = Matrix{Float64}(undef, nk, nG)
for ik in 1:nk
    values_ik = 0
    for iv in 1:2
        values_ik += abs2(sum(result[ik][iv, :]))
    end
    values_v[ik, :] .= values_ik
end
for ik in 1:nk
    values_ik = 0
    for ic in 1:2
        values_ik += abs2(sum(result[ik][:, ic]))
    end
    values_c[ik, :] .= values_ik
end
values_v = reshape(values_v, :)
values_c = reshape(values_c, :)
# eliminate redundant data
I = (-2 .< kx .< 2) .&& (-2 .< ky .< 2)
kx = kx[I]
ky = ky[I]
values_v = values_v[I]
values_c = values_c[I]

using DelimitedFiles
writedlm("./Γ1_Avk.dat", [kx ky values_v])
writedlm("./Γ1_Ack.dat", [kx ky values_c])
writedlm("./Γ1_Ak.dat", [kx ky values_v + values_c])


using JLD2
using WannierBSEKit

period = [1, 1, 0]
cell = read(joinpath("../06-interaction", "POSCAR"), POSCAR; period)
TB_up = TightBindModel(; cell,
    hops=joinpath("../05-wannier/up", "wannier90_hr.dat"),
    orbital=joinpath("../05-wannier/up", "wannier90_centres.xyz"),
)
TB_dn = TightBindModel(; cell,
    hops=joinpath("../05-wannier/dn", "wannier90_hr.dat"),
    orbital=joinpath("../05-wannier/dn", "wannier90_centres.xyz"),
)
nk = 19
kgrid = MonkhorstPack([nk, nk, 1])
Kᵈ = PairInterAction(2)
Kˣ = PairInterAction(2)
bse = BSE(kgrid, TB_up, TB_dn, Kᵈ, Kˣ; v_up=[1], c_up=Int[], v_dn=Int[], c_dn=[1])

Γ = [0, 0, 0]
# H = bse(Γ)
# wave = eigen(H.HSzn1)
# wave_Γ = Eigen(wave.values[1:1], wave.vectors[:, 1:1])
# jldsave("./excwave.jld2"; wave_Γ)
wfc = load("./excwave.jld2")["wave_Γ"]
# obtain A_{ie,Re,ih,Rh}
AR = realspace(bse, Γ, wfc.vectors[:, 1], -1)

# load wannier
wannier_up = load("./up.jld2") # v, hole
wr_h = wannier_up["wr"]
Rₕ = wannier_up["wrRlist"]

wannier_dn = load("./dn.jld2") # c, electron
wannierRlist = wannier_dn["Rlist"]
wannier_e = wannier_dn["wannier"][1]
(nfx, nfy, nfz, nR_wannier) = size(wannier_e)

# plot range
range = [5, 5, 1]
Rlist = gridindex(range)
nR = length(Rlist)
# obtain Rₑ
Rlist2Rₑ = Vector{Vector{ReducedCoordinates{Int}}}(undef, nR)
for iR in 1:nR
    Rlist2Rₑ[iR] = map(RR->Rlist[iR] - RR, wannierRlist)
end
Rₑ = unique(Iterators.flatten(Rlist2Rₑ))
Rₑ2idx = Dict(RR => iR for (iR, RR) in enumerate(Rₑ))
Rlist2Rₑmap = Vector{Vector{Int}}(undef, nR)
for iR in 1:nR
    Rlist2Rₑmap[iR] = map(RR->Rₑ2idx[RR], Rlist2Rₑ[iR])
end
# 1:v, 2:c
A = Array{ComplexF64}(undef, length(Rₕ), length(Rₑ))
for iRₑ in eachindex(Rₑ), iRₕ in eachindex(Rₕ)
    A[iRₕ, iRₑ] = AR(1, Rₕ[iRₕ], 2, Rₑ[iRₑ])
end

# the ele-part of excwave.
values = zeros(ComplexF64, nfx, nfy, nfz, nR)
Threads.@threads for iR in 1:nR
    values_iR = view(values, :, :, :, iR)
    Rlist2Rₑmap_iR = Rlist2Rₑmap[iR]
    values_c = zeros(ComplexF64, nfx, nfy, nfz)
    for iiR in 1:nR_wannier, iRₕ in eachindex(Rₕ)
        values_c += A[iRₕ, Rlist2Rₑmap_iR[iiR]] * conj(wr_h[iRₕ]) * wannier_e[:, :, :, iiR]
    end
    values_iR .+= values_c
end

(minRx, maxRx) = extrema(R -> R[1], Rlist)
(minRy, maxRy) = extrema(R -> R[2], Rlist)
(minRz, maxRz) = extrema(R -> R[3], Rlist)
nRx = maxRx - minRx + 1
nRy = maxRy - minRy + 1
nRz = maxRz - minRz + 1
R_map = Vector{NTuple{3,Int}}(undef, length(Rlist))
for (iR, R) in enumerate(Rlist)
    iRx = R[1] - minRx + 1
    iRy = R[2] - minRy + 1
    iRz = R[3] - minRz + 1
    R_map[iR] = (iRx, iRy, iRz)
end

fft_grid = [nfx, nfy, nfz]
wannier_real_grid_reshape = Array{ComplexF64,3}(undef, fft_grid[1] * nRx, fft_grid[2] * nRy, fft_grid[3] * nRz)
for iR in eachindex(Rlist)
    (iRx, iRy, iRz) = R_map[iR]
    wannier_real_grid_reshape[
        (fft_grid[1]*(iRx-1)+1):(fft_grid[1]*iRx),
        (fft_grid[2]*(iRy-1)+1):(fft_grid[2]*iRy),
        (fft_grid[3]*(iRz-1)+1):(fft_grid[3]*iRz),
    ] .= values[:, :, :, iR]
end

wave = wannier90_xsf(wannier_real_grid_reshape, range, cell)

write("./Γ_fixhole_abs.xsf", wave, wannier90_xsf; value=abs)
write("./Γ_fixhole_abs2.xsf", wave, wannier90_xsf; value=abs2)

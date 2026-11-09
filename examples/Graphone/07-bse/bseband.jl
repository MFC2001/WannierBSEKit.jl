
using JLD2
using WannierBSEKit

period = [1, 1, 0]
cell = read(joinpath("../06-interaction", "POSCAR"), POSCAR; period)
TB_up = TightBindModel(; cell,
	hops = joinpath("../05-wannier/up", "wannier90_hr.dat"),
	orbital = joinpath("../05-wannier/up", "wannier90_centres.xyz"),
)
TB_dn = TightBindModel(; cell,
	hops = joinpath("../05-wannier/dn", "wannier90_hr.dat"),
	orbital = joinpath("../05-wannier/dn", "wannier90_centres.xyz"),
)

nk = 19
kgrid = MonkhorstPack([nk, nk, 1])
outfolder = "./"

intfolder = "../06-interaction/"
Kᵈ_U = load(joinpath(intfolder, "U_1_mirror.jld2"))["Uwithlr"]
Kᵈ_pair = load(joinpath(intfolder, "pair_1.jld2"))["interaction"]
Kᵈ = PairInterAction(Kᵈ_U, Kᵈ_pair)
Kˣ = PairInterAction(2) # This is an empty interaction, because magnon doesn't need Kˣ. nw = nw_up+nw_dn = 2

bse = BSE(kgrid, TB_up, TB_dn, Kᵈ, Kˣ; v_up=[1], c_up=Int[], v_dn=Int[], c_dn=[1])
kline = Kline(30, bse)
# for v_up to c_dn, the total spin of excitation is -1
band = BAND(kline, bse).bandSzn1

# ΔE_gap_g0w0 ≈ 3
# minimum(band) ≈ -3.38
band .-= minimum(band)
write(joinpath(outfolder, "bseband.dat"), band, BAND_dat; kline)

# using Plots
# plotlyjs()
# plot(kline, band)

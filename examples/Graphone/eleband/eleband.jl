
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

kline = Kline(50, cell)
band_up = BAND(kline, TB_up)
band_dn = BAND(kline, TB_dn)

write("eleband_up.dat", band_up, BAND_dat; kline)
write("eleband_dn.dat", band_dn, BAND_dat; kline)

# band = [band_up; band_dn]

# using Plots
# plotlyjs()

# plot(kline, band)


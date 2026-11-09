
using WannierBSEKit

period = [1, 1, 0]
cell = read("../08-interaction/POSCAR", POSCAR; period)
TB = TightBindModel(; cell,
	hops = "../07-wannier/wannier90_hr.dat",
	orbital = "../07-wannier/wannier90_centres.xyz",
)

kline = Kline(50, cell)
band = BAND(kline, TB)

write("./eleband.dat", band, BAND_dat; kline)

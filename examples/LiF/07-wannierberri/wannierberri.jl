
using WannierBSEKit

xmldata = read("./pwscf.save/data-file-schema.xml", QE_xml_new; period = [1, 1, 1])
cell = convert(xmldata["cell"])

kgrid = RedKgrid(MonkhorstPack([15, 15, 15]))
kgrid.kdirect .= xmldata["kpoints"]

bands = read("./w90data.eig.npz", wannierberri_eig)
chkdata = read("./w90data.chk.npz", wannierberri_chk)

Hk = build_wannier_hk(bands, chkdata)
hr = HK2HR(Hk, kgrid, cell.lattice)
write("./wannier90_hr.dat", hr, wannier90_hr)


orbital = wannier90_centres(map(x -> cell.lattice * x, chkdata.wannier_centres); 
	name = fill("X", chkdata.num_wann),
	atom_location = cell.location,
	atom_name = cell.name,
)
write("./wannier90_centres.xyz", orbital, wannier90_centres)

TB = TightBindModel(; cell, hops=hr, orbital)

kline = [
    "L" => [0.5, 0.5, 0.5],
    "Γ" => [0, 0, 0],
    "X" => [0.5, 0, 0.5],
    "W" => [0.5, 0.25, 0.75],
    "L" => [0.5, 0.5, 0.5],
]
kline = Kline(kline, 50, TB)

band = BAND(kline, TB; vector=true)
write("./eleband.dat", band, BAND_dat; kline)


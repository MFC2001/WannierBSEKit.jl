
using WannierBSEKit

band = read("wannier90.eig", wannier90_eig)
mv("wannier90.eig", "wannier90.eig_original")

(nband, nk) = size(band)
# Based on the DFT band structure of LiF, 
# the GW band structure was obtained by stretching the valence band by 15% 
# and increasing the band gap to 14.3 eV.
GWband = scissorslinearfit(band, 1:5, 6:nband;
	evs = 0,
	ev0 = 0.511845648084,
	evdel = 0.15,
	ecs = 5.421786646273,
	ec0 = 0,
	ecdel = 0,
)

write("wannier90.eig", GWband, wannier90_eig)
write("wannier90.eig_GW", GWband, wannier90_eig)

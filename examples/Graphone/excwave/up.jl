
using JLD2
using WannierBSEKit

period = [1, 1, 0]
chkdata = read("../05-wannier/up/wannier90.chk", wannier90_chk)
wave_ik = read("../05-wannier/up/pwscf.save/", chkdata, QE_wannier; wfcspin = :up, period)

(Rlist, _) = gridindex(:WignerSeitz, chkdata.lattice, (5, 5, 1))
wannier = realspace(wave_ik, Rlist)
(wrRlist, _) = gridindex(:WignerSeitz, chkdata.lattice, (9, 9, 1))
r = [2//3, 1//3, 0.4871518209]
wr = wave_ik(r, wrRlist)

jldsave("./up.jld2"; wannier, Rlist, wr, wrRlist)

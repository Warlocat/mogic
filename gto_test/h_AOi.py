import pyscf, numpy
from pyscf import gto, scf, tools

mol = gto.M(
    atom = 'H 0 0 0',
    spin = 1,
    basis = 'cc-pvdz'
)
mf = scf.UHF(mol)
mf.kernel()
dm = mf.make_rdm1()
mo_feak = numpy.eye(5)

tools.fcidump.from_mo(mol,"FCIDUMP", mo_feak)
# tools.fcidump.from_scf(mf, "FCIDUMP")

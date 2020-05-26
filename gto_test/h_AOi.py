import pyscf, numpy
from pyscf import gto, scf, tools

mol = gto.M(
    atom = 'P 0 0 0',
    spin = 1,
    basis = 'unccc-pvdz'
)
# mf = scf.UHF(mol)
# mf.kernel()
# dm = mf.make_rdm1()
s = mol.intor('int1e_ovlp')
mo_feak = numpy.eye(len(s))
# s = mol.intor('int1e_ovlp')
# t = mol.intor('int1e_kin')
# v = mol.intor('int1e_nuc')
# print(t)

tools.fcidump.from_mo(mol,"p_unc", mo_feak)
# tools.fcidump.from_scf(mf, "FCIDUMP")

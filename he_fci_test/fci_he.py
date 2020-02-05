import pyscf
from pyscf import gto, scf, ci, fci, tools

mol = gto.M(
    atom = 'He 0 0 0',
    basis = 'ccpvtz')

mf = scf.RHF(mol).run()
# pyscf.tools.fcidump.from_scf(mf, "FCIDUMP")


mycc = ci.CISD(mf)
mycc.nroots = 20 
mycc.kernel()
# e, fcivec = fci.FCI(mf).kernel(verbose=5)
# print(e)

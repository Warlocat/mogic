import pyscf
from pyscf import gto

mol = gto.M(
    atom = 'H 0 0 0',
    spin = 1,
    basis = 'cc-pvdz'
)

s = mol.intor('int1e_ovlp')
t = mol.intor('int1e_kin')
v = mol.intor('int1e_nuc')

print("s",s)
# print("h1e",t+v)
# print("v",v)

#!/usr/bin/env python
import time
from functools import reduce
import copy
import numpy
import scipy.linalg
def _get_r(s, snesc):
    # R^dag \tilde{S} R = S
    # R = S^{-1/2} [S^{-1/2}\tilde{S}S^{-1/2}]^{-1/2} S^{1/2}
    w, v = numpy.linalg.eigh(s)
    idx = w > 1e-14
    v = v[:,idx]
    w_sqrt = numpy.sqrt(w[idx])
    w_invsqrt = 1 / w_sqrt

    # eigenvectors of S as the new basis
    snesc = reduce(numpy.dot, (v.conj().T, snesc, v))
    r_mid = numpy.einsum('i,ij,j->ij', w_invsqrt, snesc, w_invsqrt)
    w1, v1 = numpy.linalg.eigh(r_mid)
    idx1 = w1 > 1e-14
    v1 = v1[:,idx1]
    r_mid = numpy.dot(v1/numpy.sqrt(w1[idx1]), v1.conj().T)
    r = numpy.einsum('i,ij,j->ij', w_invsqrt, r_mid, w_sqrt)
    # Back transform to AO basis
    r = reduce(numpy.dot, (v, r, v.conj().T))
    return r
def _x2c1e_xmatrix(t, v, w, s, c):
    nao = s.shape[0]
    n2 = nao * 2
    h = numpy.zeros((n2,n2), dtype=v.dtype)
    m = numpy.zeros((n2,n2), dtype=v.dtype)
    h[:nao,:nao] = v
    h[:nao,nao:] = t
    h[nao:,:nao] = t
    h[nao:,nao:] = w * (.25/c**2) - t
    m[:nao,:nao] = s
    m[nao:,nao:] = t * (.5/c**2)

    try:
        e, a = scipy.linalg.eigh(h, m)
        cl = a[:nao,nao:]
        cs = a[nao:,nao:]
        x = numpy.linalg.solve(cl.T, cs.T).T  # B = XA
    except scipy.linalg.LinAlgError:
        d, t = numpy.linalg.eigh(m)
        idx = d>1e-9
        t = t[:,idx] / numpy.sqrt(d[idx])
        tht = reduce(numpy.dot, (t.T.conj(), h, t))
        e, a = numpy.linalg.eigh(tht)
        a = numpy.dot(t, a)
        idx = e > -c**2
        cl = a[:nao,idx]
        cs = a[nao:,idx]
        # X = B A^{-1} = B A^T S
        x = cs.dot(cl.conj().T).dot(m)
    return x
def _get_hcore_fw(t, v, w, s, x, c):
    s1 = s + reduce(numpy.dot, (x.T.conj(), t, x)) * (.5/c**2)
    tx = numpy.dot(t, x)
    h1 =(v + tx + tx.T.conj() - numpy.dot(x.T.conj(), tx) +
         reduce(numpy.dot, (x.T.conj(), w, x)) * (.25/c**2))

    r = _get_r(s, s1)
    h1 = reduce(numpy.dot, (r.T.conj(), h1, r))
    return h1
from pyscf import gto
from pyscf import lib
from pyscf import scf
from pyscf import x2c

mol = gto.M(
    atom = '''Cu  0.  0.  0.''',
    spin = 1,
    basis = 'cc-pvdz',
)

t = mol.intor_symmetric('int1e_kin')
v = mol.intor_symmetric('int1e_nuc')
s = mol.intor_symmetric('int1e_ovlp')
w = mol.intor_symmetric('int1e_pnucp')
c = lib.param.LIGHT_SPEED
x = _x2c1e_xmatrix(t,v,w,s,c)
s1 = s + reduce(numpy.dot, (x.T.conj(), t, x)) * (.5/c**2)

mf = scf.RHF(mol).sfx2c1e()
mf.scf()
\

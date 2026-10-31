"""Absolute edge trace A from the FEM edge profile psi(L - r), r -> 0 (ell = log 1/r).

No Hadamard formula and no eigenvalue derivative is used: A comes only from the
eigenfunction profile, normalised in L^2(-L,L). Three fit models on ell in [lo, ell_max - cut]:
  F1  psi^-2         = alpha*ell + beta + c1/ell + c2/ell^2,     A = alpha^-1/2
  F2  psi*ell^1/2    = A + b1/ell + b2/ell^2 + b3/ell^3
  F3  psi*ell^1/2    = A + b0/ell^1/2 + b1/ell + b2/ell^3/2 + b3/ell^2   (allows a V^2 term; badly
      conditioned: applied to the exact-shape profile (2L)^-1/2 (ell-lambda-gamma)^-1/2 it already shifts A by ~3e-5)
  F4  psi^-2         = alpha*ell + beta + g*ell^1/2 + c1/ell + c2/ell^2   (V^2 term, well conditioned), A = alpha^-1/2
plus the local slope of psi^-2 against ell (-> 1/A^2) at fixed ell.
"""
import json, sys
import numpy as np

def fits(r, psi, lo=30.0, cut=20.0):
    sel = (r > 0) & (r < 0.1)
    r, psi = r[sel], psi[sel]
    ell = -np.log(r)
    m = (ell >= lo) & (ell <= ell.max() - cut)
    l, p = ell[m], psi[m]
    out = {}
    X = np.stack([l, np.ones_like(l), 1 / l, 1 / l ** 2], 1)
    c = np.linalg.lstsq(X, p ** -2, rcond=None)[0]
    out['F1'] = c[0] ** -0.5
    out['F1_c'] = c[1] / c[0]                       # psi^-2 ~ (ell + c)/A^2
    X = np.stack([np.ones_like(l), 1 / l, 1 / l ** 2, 1 / l ** 3], 1)
    out['F2'] = np.linalg.lstsq(X, p * l ** 0.5, rcond=None)[0][0]
    X = np.stack([np.ones_like(l), l ** -0.5, 1 / l, l ** -1.5, 1 / l ** 2], 1)
    c3 = np.linalg.lstsq(X, p * l ** 0.5, rcond=None)[0]
    out['F3'], out['F3_b0'] = c3[0], c3[1]
    X = np.stack([l, np.ones_like(l), l ** 0.5, 1 / l, 1 / l ** 2], 1)
    c4 = np.linalg.lstsq(X, p ** -2, rcond=None)[0]
    out['F4'], out['F4_g'] = c4[0] ** -0.5, c4[2]
    y = psi ** -2
    sl = np.gradient(y, ell)
    for L0 in (100, 200, 400):
        if ell.max() - cut > L0:
            out[f'slope{L0}'] = float(np.interp(L0, ell[::-1], sl[::-1]))
    return {k: float(v) for k, v in out.items()}

def analyse(npz):
    d = np.load(npz)
    res = {}
    for sec in ('even', 'odd'):
        for k, psi in enumerate(d['psi_' + sec]):
            res[f'{sec[0]}{k}'] = fits(d['r'], psi)
    return res

if __name__ == '__main__':
    for f in sys.argv[1:]:
        L = float(json.load(open(f.replace('.npz', '.json')))['L']) if f.endswith('.npz') else None
        res = analyse(f)
        for st, v in res.items():
            print(f, st, ' '.join(f'{k}={v[k]:.7f}' for k in v), ' 2LA^2(F1)=%.7f' % (2 * L * v['F1'] ** 2))

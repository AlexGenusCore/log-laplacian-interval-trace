"""Row residuals of the stored P1 eigenvectors at the deepest nodes (r = r_min, r_min*rho, ...).

Added 2026-10-31 for the numerical supplement. Rebuilds only the matrix rows of the k deepest free
nodes, with the element routines of fem.py, and evaluates them on the stored eigenvectors.
Residual of row i, exactly as fem.py reports it (key rowres_*, there maximised over all rows):
    res_i = |((K - mu M) x)_i| / M_ii,   K = -(F_same +/- F_cross) - 2 gamma M   (+ even, - odd),
x = stored nodal vector (normalised in L^2(-L,L)), mu = stored refined eigenvalue of L_Delta.
Usage:  python residual_check.py out/L1_n400_r1.12 [more tags]     |     python residual_check.py --test
"""
import json, sys
try:
    import scipy.linalg  # noqa: F401  (fem.py imports it at top level; not used here)
except ImportError:      # ponytail: shim so the row check runs without scipy; fem.solve_sector would need the real one
    import types
    sys.modules['scipy'] = types.ModuleType('scipy')
    sys.modules['scipy.linalg'] = sys.modules['scipy'].linalg = types.ModuleType('scipy.linalg')
import numpy as np
from fem import EG, FP1, MLOC, G_adj, adj_block, gl_blocks, n_for_beta, half_mesh, assemble


def rows(r, h, L, nodes):
    """Rows `nodes` of F_same, F_cross, M (same conventions and element integrals as fem.assemble)."""
    ne = len(h); N = ne + 1
    nodes = np.asarray(nodes)
    pos = -np.ones(N, int); pos[nodes] = np.arange(len(nodes))
    Fs, Fc, M = (np.zeros((len(nodes), N)) for _ in range(3))

    def add(T, i, j, v):
        i, j, v = np.broadcast_arrays(i, j, v)
        k = pos[i] >= 0
        np.add.at(T, (pos[i[k]], j[k]), v[k])

    el = np.unique(np.clip(np.concatenate([nodes - 1, nodes]), 0, ne - 1))     # elements touching the nodes
    for a in range(2):
        for b in range(2):
            add(M, el + a, el + b, h[el] * MLOC[a, b])
            add(Fs, el + a, el + b, h[el] * (FP1[a, b] + 2 * np.log(h[el]) * MLOC[a, b]))
    e = np.unique(np.clip(np.concatenate([el - 1, el]), 0, ne - 2))            # adjacent pairs (e, e+1)
    PlQl, PlQr, PrQl, PrQr = adj_block(h[e], h[e + 1])
    for ia, ib, val in ((e, e + 1, PlQl), (e, e + 2, PlQr), (e + 1, e + 1, PrQl), (e + 1, e + 2, PrQr)):
        add(Fs, ia, ib, val); add(Fs, ib, ia, val)
    # all pairs E <= F with E or F in el (the only pairs that write into these rows)
    P = np.unique(np.concatenate([t * ne + np.arange(t, ne) for t in el] + [np.arange(t + 1) * ne + t for t in el]))
    E, Fi = P // ne, P % ne
    X = L - r[1:]; X[-1] = 0.0
    sep = Fi >= E + 2
    cr = ~((E == ne - 1) & (Fi == ne - 1))
    for kind, EE, FF, DD in (('s', E[sep], Fi[sep], r[Fi[sep]] - r[E[sep] + 1]), ('c', E[cr], Fi[cr], X[E[cr]] + X[Fi[cr]])):
        he, hf = h[EE], h[FF]
        nn = n_for_beta(DD / np.maximum(he, hf))
        blocks = np.empty((len(EE), 2, 2))
        for n in np.unique(nn):
            s = nn == n
            blocks[s] = gl_blocks(he[s], hf[s], DD[s], n)
        if kind == 'c':
            blocks = blocks[:, :, ::-1]
        T = Fs if kind == 's' else Fc
        off = EE != FF
        for a in range(2):
            for b in range(2):
                add(T, EE + a, FF + b, blocks[:, a, b])
                add(T, FF[off] + b, EE[off] + a, blocks[off, a, b])
    if ne - 1 in el:                                                            # centre-centre cross pair
        hc = h[-1]
        G00, G10, G01, G11 = G_adj(np.array(1.0), np.array(1.0))
        for i, j, v in ((ne - 1, ne - 1, hc * G11), (ne - 1, ne, hc * (G10 - G11)), (ne, ne - 1, hc * (G10 - G11)),
                        (ne, ne, hc * (G00 - 2 * G10 + G11))):
            add(Fc, i, j, v)
    return Fs, Fc, M


def residuals(tag, k=5):
    d = np.load(tag + '.npz'); meta = json.load(open(tag + '.json'))
    r = d['r']; L = meta['L']
    nodes = np.arange(1, k + 1)
    Fs, Fc, M = rows(r, np.diff(r), L, nodes)
    Mii = M[np.arange(k), nodes]
    out = {}
    for sec, sgn in (('even', 1), ('odd', -1)):
        K = -(Fs + sgn * Fc) - 2 * EG * M
        for j, (x, mu) in enumerate(zip(d['psi_' + sec], meta['mu_' + sec])):
            out[f'{sec[0]}{j}'] = np.abs(K @ x - mu * (M @ x)) / Mii
    return meta, r[nodes], out


def test():
    for L in (0.5, 1.0, 3.0):
        r, h = half_mesh(L, 10, 1.5, 1e-30)
        Fs, Fc, M = assemble(r, h, L)
        nodes = [1, 2, 3, len(r) - 3, len(r) - 2, len(r) - 1]
        for full, part in zip((Fs, Fc, M), rows(r, h, L, nodes)):
            assert np.allclose(full[nodes], part, rtol=1e-13, atol=1e-13 * np.abs(full).max()), L
    print('residual_check test ok')


if __name__ == '__main__':
    if sys.argv[1:] == ['--test']:
        test(); raise SystemExit
    print('# tag  N_half  r_min  max_rowres_all_rows(fem.py)  res(r=r_min)  max res over 5 deepest nodes'
          '   [max over the stored eigenvectors: 3 per sector, 1 for scal_*]')
    for tag in sys.argv[1:]:
        meta, rn, res = residuals(tag)
        allrows = max(meta['rowres_even'] + meta['rowres_odd'])
        deepest = max(v[0] for v in res.values())
        five = max(v.max() for v in res.values())
        print(f'{tag.split("/")[-1]:24s} {meta["N_half"]:5d} {rn[0]:.0e}  {allrows:.2e}  {deepest:.2e}  {five:.2e}', flush=True)

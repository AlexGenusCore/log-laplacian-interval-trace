# Independent referee check: piecewise-constant (P0) Galerkin for Q (symbol log|t|) on (-L,L),
# parity-reduced, unknowns = element values u_k (row-equilibrated, so tiny edge elements keep
# relative accuracy). Checks: eigenvalues mu=2*lambda vs Table 1, and slope of u^{-2} vs ell -> 2L.
import sys, numpy as np
EG = 0.5772156649015329
def mesh(L, rmin, rho, n0):
    h0 = L / n0
    r = [0.0, rmin]
    while r[-1] * (rho - 1) < h0 and r[-1] < L / 2:
        r.append(r[-1] * rho)
    m = int(np.ceil((L - r[-1]) / h0))
    r += list(np.linspace(r[-1], L, m + 1)[1:])
    r = np.array(r); r[-1] = L
    return r
xg, wg = np.polynomial.legendre.leggauss(16); xg = (xg + 1) / 2; wg = wg / 2
def rowmean(h1, h2, g):
    # (1/h1) * int_{elem1} int_{elem2} dx dy/|x-y|, elements of lengths h1,h2, gap g>=0
    out = np.empty(np.broadcast(h1, h2, g).shape)
    h1, h2, g = np.broadcast_arrays(h1, h2, g)
    a = g == 0
    q = h2[a] / h1[a]
    out[a] = np.log1p(q) + q * np.log1p(1 / q)
    b = ~a
    s = g[b][:, None] + h1[b][:, None] * xg[None, :]
    out[b] = (np.log1p(h2[b][:, None] / s) * wg[None, :]).sum(1)
    return out
def build(L, r):
    lo, hi = r[:-1], r[1:]; h = hi - lo; N = len(h)
    B = {}
    for sgn in (+1, -1):
        M = np.empty((N, N))
        for k in range(N):
            same = np.where(np.arange(N) > k, lo - hi[k], lo[k] - hi)   # gap in r, same half
            cross = 2 * L - hi[k] - hi                                    # gap across centre
            Ms = -0.5 * rowmean(h[k], h, np.maximum(same, 0.0))
            Ms[k] = 1 - EG - np.log(h[k])
            Mc = -0.5 * rowmean(h[k], h, np.maximum(cross, 0.0))
            M[k] = Ms + sgn * Mc
        B[sgn] = M
    return h, B
def solve(L, rmin, rho, n0, nev=3):
    r = mesh(L, rmin, rho, n0); h, B = build(L, r)
    rmid = np.sqrt(r[:-1] * r[1:]); rmid[0] = r[1] / np.e  # element [0,rmin]: mean of ell = ell(rmin)+1
    ell = np.log(1 / rmid)
    res = []
    for sgn, par in ((+1, 'e'), (-1, 'o')):
        Bm = B[sgn]; sq = np.sqrt(h)
        S = Bm * sq[:, None] / sq[None, :]; S = (S + S.T) / 2
        lam, V = np.linalg.eigh(S)
        for j in range(nev):
            u = V[:, j] / sq
            for _ in range(6):
                u = np.linalg.solve(Bm - (lam[j] + 1e-9) * np.eye(len(h)), u)
                u /= np.sqrt(2 * np.sum(h * u * u))
            lam_r = (u @ (h * (Bm @ u))) / (u @ (h * u))
            u *= np.sign(u[0])
            w = (ell >= 30) & (ell <= ell.max() - 10)
            X = np.vstack([ell[w], np.ones(w.sum()), 1 / ell[w]]).T
            al, be, c = np.linalg.lstsq(X, u[w] ** -2, rcond=None)[0]
            # raw (fit-free) value at the deepest well-resolved element
            i = np.argmin(np.abs(ell - (ell.max() - 10)))
            raw = 2 * L * u[i] ** 2 * ell[i]
            res.append((par + str(j), 2 * lam_r, 2 * L / al, be / al + lam_r + EG, raw))
    return len(h), res
if __name__ == '__main__':
    L, rmin, rho, n0 = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
    N, res = solve(L, rmin, rho, n0)
    print(f'L={L} rmin={rmin} rho={rho} n0={n0} N_half={N}')
    for name, mu, tA, conj, raw in sorted(res, key=lambda t: t[1]):
        print(f'{name}: mu={mu:.8f}  2L*A^2(fit)={tA:.7f}  beta/alpha+mu/2+gamma={conj:+.4f}  raw 2L u^2 ell={raw:.5f}')

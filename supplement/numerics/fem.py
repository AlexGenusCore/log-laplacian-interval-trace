"""P1 Galerkin FEM for the 1-D logarithmic Laplacian L_Delta (symbol 2 log|xi|) on (-L, L).

Written from scratch for the pure-log paper (2026-10-31); shares no code with the
earlier project solvers (hp-FEM of boundary-clue, the referee's P0 code).

Formulation (finite part; no killing integral, no truncation radius):
  a(u,v) = -FP iint_{Omega^2} u(x) v(y) / |x-y| dx dy - 2*gamma int_Omega u v,
  FP iint := lim_{eps->0} [ iint_{|x-y|>eps} + 2 log(eps) int u v ].
This is exact because the Fourier transform of Pf(1/|x|) is -2 gamma - 2 log|xi|.
The half log-Laplacian Q^log of the note is a/2 (lambda = mu/2).

Parity reduction on the half mesh r in [0, L] (r = distance to the right edge,
x = L - r, so tiny r is stored without cancellation):
  even/odd:  K = -(F_same +/- F_cross) - 2 gamma M,
  F_same : kernel 1/|r - r'|   (finite part on the diagonal elements),
  F_cross: kernel 1/(2L - r - r')  (= 1/(x + y), the mirror image).
Mesh: geometric grading r_min * rho^k at the edge, uniform h0 = L/n0 inside.
Edge node r = 0 is Dirichlet; the centre node r = L is free (even) or Dirichlet (odd).
"""
import argparse, json, time
import numpy as np
import scipy.linalg as sla

EG = 0.57721566490153286061

# ---------------------------------------------------------------- mesh
def half_mesh(L, n0, rho, rmin):
    h0 = L / n0
    r = [0.0, rmin]
    while (rho - 1) * r[-1] < h0 and rho * r[-1] < L:
        r.append(r[-1] * rho)
    rg = r[-1]
    nu = max(1, int(np.ceil((L - rg) / h0)))
    r += list(rg + (L - rg) * np.arange(1, nu + 1) / nu)
    r[-1] = L
    r = np.array(r)
    h = np.diff(r)
    return r, h

# ---------------------------------------------------------------- exact element integrals
FP1 = np.array([[-8 / 9, -1 / 9], [-1 / 9, -8 / 9]])   # FP int_0^1 int_0^1 p_a(s) p_b(t)/|s-t|
MLOC = np.array([[1 / 3, 1 / 6], [1 / 6, 1 / 3]])

def xlogx(z):
    return np.where(z > 0, z * np.log(np.where(z > 0, z, 1.0)), 0.0)

def G_adj(h, k):
    """G_mn = int_0^h int_0^k a^m y^n /(a+y) dy da, m,n in {0,1} (closed forms)."""
    lh, lk, lhk = np.log(h), np.log(k), np.log(h + k)
    G00 = (h + k) * lhk - h * lh - k * lk
    G10 = (h * h - k * k) / 2 * lhk + k * h / 2 + k * k / 2 * lk - h * h / 2 * lh
    G01 = (k * k - h * h) / 2 * lhk + k * h / 2 + h * h / 2 * lh - k * k / 2 * lk
    G11 = h * h * k / 3 + h * k * k / 3 - (h ** 3 + k ** 3) / 3 * lhk + h ** 3 / 3 * lh + k ** 3 / 3 * lk
    return G00, G10, G01, G11

def adj_block(h, k):
    """Left element [-h,0] (hats P_l=a/h at -h, P_r=1-a/h at 0), right element [0,k]
    (Q_l=1-y/k at 0, Q_r=y/k at k); kernel 1/(a+y). Returns [[PlQl, PlQr],[PrQl, PrQr]].
    Evaluated in units s=max(h,k) (degree-1 homogeneous) to avoid cancellation."""
    s = np.maximum(h, k)
    hh, kk = h / s, k / s
    G00, G10, G01, G11 = G_adj(hh, kk)
    PlQl = (G10 - G11 / kk) / hh
    PlQr = G11 / (hh * kk)
    PrQl = G00 - G10 / hh - G01 / kk + G11 / (hh * kk)
    PrQr = (G01 - G11 / hh) / kk
    return s * PlQl, s * PlQr, s * PrQl, s * PrQr

# ---------------------------------------------------------------- Gauss-Legendre far pairs
_GL = {}
def gl01(n):
    if n not in _GL:
        x, w = np.polynomial.legendre.leggauss(n)
        _GL[n] = ((x + 1) / 2, w / 2)
    return _GL[n]

def n_for_beta(beta):
    z = 1 + 2 * beta
    rb = z + np.sqrt(np.minimum(z, 1e150) ** 2 - 1)
    return np.clip(np.ceil(39.0 / (2 * np.log(rb))) + 1, 4, 30).astype(int)

def gl_blocks(he, hf, D, n):
    """int_e int_f phi_a(s) psi_b(t) / (D + he(1-s) + hf t) ds dt * he * hf, with
    phi = (1-s, s), psi = (1-t, t). Returns array (P,2,2). Multiplication order avoids underflow."""
    s, w = gl01(n)
    out = np.empty((len(he), 2, 2))
    chunk = max(1, int(4e6 // (n * n)))
    for i0 in range(0, len(he), chunk):
        a, b, d = he[i0:i0 + chunk, None, None], hf[i0:i0 + chunk, None, None], D[i0:i0 + chunk, None, None]
        kern = b / (d + a * (1 - s[None, :, None]) + b * s[None, None, :])   # hf/dist, dimensionless
        ww = (w[:, None] * w[None, :])[None]
        ps = np.stack([1 - s, s])          # (2, n)
        out[i0:i0 + chunk] = np.einsum('pij,ai,bj->pab', kern * ww, ps, ps) * a[:, :, 0][:, :, None]
    return out

# ---------------------------------------------------------------- assembly
def assemble(r, h, L):
    ne = len(h)
    N = ne + 1
    Fs = np.zeros((N, N)); Fc = np.zeros((N, N)); M = np.zeros((N, N))
    idx = np.arange(ne)
    # mass
    for a in range(2):
        for b in range(2):
            np.add.at(M, (idx + a, idx + b), h * MLOC[a, b])
    # same-side, diagonal elements (finite part)
    for a in range(2):
        for b in range(2):
            np.add.at(Fs, (idx + a, idx + b), h * (FP1[a, b] + 2 * np.log(h) * MLOC[a, b]))
    # same-side, adjacent elements e, e+1 (shared node e+1)
    e = idx[:-1]
    PlQl, PlQr, PrQl, PrQr = adj_block(h[:-1], h[1:])
    for (ia, ib, val) in ((e, e + 1, PlQl), (e, e + 2, PlQr), (e + 1, e + 1, PrQl), (e + 1, e + 2, PrQr)):
        np.add.at(Fs, (ia, ib), val)
        np.add.at(Fs, (ib, ia), val)
    # same-side separated pairs e < f-1 ; gap g = r_f - r_{e+1}
    # cross pairs e <= f except the centre-centre pair ; distance offset X_e + X_f, X = L - r_{e+1}
    X = L - r[1:]
    X[-1] = 0.0
    E, Fi = np.triu_indices(ne, 0)
    sep = Fi >= E + 2
    jobs = []
    Es, Fss = E[sep], Fi[sep]
    g = r[Fss] - r[Es + 1]
    jobs.append(('s', Es, Fss, g))
    cr = ~((E == ne - 1) & (Fi == ne - 1))
    Ec, Fcc = E[cr], Fi[cr]
    jobs.append(('c', Ec, Fcc, X[Ec] + X[Fcc]))
    for kind, EE, FF, DD in jobs:
        he, hf = h[EE], h[FF]
        nn = n_for_beta(DD / np.maximum(he, hf))
        blocks = np.empty((len(EE), 2, 2))
        for n in np.unique(nn):
            sel = nn == n
            # same side: dist = g + he(1-s) + hf t   (e left of f in r)
            # cross:     dist = X_e + X_f + he(1-s) + hf(1-t); computed with t -> 1-t, fixed below
            blocks[sel] = gl_blocks(he[sel], hf[sel], DD[sel], n)
        if kind == 'c':
            blocks = blocks[:, :, ::-1]     # t -> 1-t swaps the two f-hats
        T = Fs if kind == 's' else Fc
        for a in range(2):
            for b in range(2):
                np.add.at(T, (EE + a, FF + b), blocks[:, a, b])
                off = EE != FF
                np.add.at(T, (FF[off] + b, EE[off] + a), blocks[off, a, b])
    # centre-centre cross pair, analytic: a = h(1-s), hats a/h (node ne-1), 1-a/h (node ne)
    hc = h[-1]
    G00, G10, G01, G11 = G_adj(np.array(1.0), np.array(1.0))
    Fc[ne - 1, ne - 1] += hc * G11
    Fc[ne - 1, ne] += hc * (G10 - G11); Fc[ne, ne - 1] += hc * (G10 - G11)
    Fc[ne, ne] += hc * (G00 - 2 * G10 + G11)
    return Fs, Fc, M

# ---------------------------------------------------------------- solve
def solve_sector(Fs, Fc, M, sign, nev, nvec=3, refine=True):
    N = Fs.shape[0]
    K = -(Fs + sign * Fc) - 2 * EG * M
    free = np.arange(1, N) if sign > 0 else np.arange(1, N - 1)
    Kf, Mf = K[np.ix_(free, free)], M[np.ix_(free, free)]
    d = 1 / np.sqrt(np.diag(Mf))
    C = Kf * d[:, None] * d[None, :]
    B = Mf * d[:, None] * d[None, :]
    mu, Y = sla.eigh(C, B, subset_by_index=[0, nev - 1])
    Xv = Y * d[:, None]
    res = []
    rs = 1 / np.diag(Mf)                       # row scaling ~ 1/h
    for j in range(nvec):
        x = Xv[:, j].copy()
        m = mu[j]
        if refine:
            # eigh's components at tiny r carry absolute error ~eps/sqrt(M_jj): drop them, then do
            # inverse iteration on the row-equilibrated system (accurate tiny-r components)
            x[np.sqrt(np.diag(Mf)) < 1e-7 * np.sqrt(Mf.diagonal().max())] = 0.0
            lu = sla.lu_factor((Kf - m * Mf) * rs[:, None])
            for it in range(3):
                x = sla.lu_solve(lu, (Mf @ x) * rs)
                x /= np.sqrt(x @ Mf @ x)
            m = (x @ Kf @ x) / (x @ Mf @ x)
        xfull = np.zeros(N); xfull[free] = x
        xfull /= np.sqrt(2 * xfull @ M @ xfull)      # L^2(-L,L) normalisation
        if xfull[1] < 0:
            xfull = -xfull                            # A_+ > 0 convention
        rowres = np.abs((Kf - m * Mf) @ xfull[free]) * rs
        res.append((m, xfull, float(rowres.max())))
    return res, mu

def phi_prime(xe, xo):
    """Phi'(0) = <psi_o, psi_e'> over (-L,L) = -2 sum (dpsi_e)(avg psi_o) on the half mesh."""
    return -2 * np.sum(np.diff(xe) * (xo[:-1] + xo[1:]) / 2)

def selftest():
    for L in (0.5, 1.0, 3.0):
        r, h = half_mesh(L, 10, 1.5, 1e-30)
        Fs, Fc, M = assemble(r, h, L)
        one = np.ones(len(r))
        assert abs(one @ Fs @ one - 2 * L * (np.log(L) - 1)) < 1e-11 * max(1, L), (one @ Fs @ one, 2 * L * (np.log(L) - 1))
        assert abs(one @ Fc @ one - 2 * L * np.log(2)) < 1e-11 * L
        assert abs(r @ Fs @ r - L ** 3 * (-8 / 9 + 2 / 3 * np.log(L))) < 1e-11 * L ** 3
        assert abs(r @ Fc @ r - L ** 3 * (4 / 3 * np.log(2) - 1 / 3)) < 1e-11 * L ** 3
    print('selftest ok')

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=float, default=1.0)
    ap.add_argument('--n0', type=int, default=100)
    ap.add_argument('--rho', type=float, default=1.5)
    ap.add_argument('--rmin', type=float, default=1e-30)
    ap.add_argument('--nev', type=int, default=10, help='eigenvalues per sector')
    ap.add_argument('--nvec', type=int, default=3, help='refined eigenvectors per sector')
    ap.add_argument('--out', default=None)
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        selftest(); raise SystemExit
    t0 = time.time()
    r, h = half_mesh(a.L, a.n0, a.rho, a.rmin)
    Fs, Fc, M = assemble(r, h, a.L)
    t1 = time.time()
    ev, mue = solve_sector(Fs, Fc, M, +1, a.nev, a.nvec)
    od, muo = solve_sector(Fs, Fc, M, -1, a.nev, a.nvec)
    t2 = time.time()
    out = dict(L=a.L, n0=a.n0, rho=a.rho, rmin=a.rmin, N_half=len(r), t_assemble=t1 - t0, t_solve=t2 - t1,
               mu_even=[e[0] for e in ev], mu_odd=[o[0] for o in od],
               mu_even_all=mue.tolist(), mu_odd_all=muo.tolist(),
               rowres_even=[e[2] for e in ev], rowres_odd=[o[2] for o in od],
               phi_prime={f'e{i}o{j}': phi_prime(ev[i][1], od[j][1]) for i in range(a.nvec) for j in range(a.nvec)})
    print(json.dumps({k: v for k, v in out.items()}, indent=None))
    if a.out:
        np.savez(a.out + '.npz', r=r, psi_even=np.array([e[1] for e in ev]), psi_odd=np.array([o[1] for o in od]),
                 mu_even=out['mu_even'], mu_odd=out['mu_odd'])
        with open(a.out + '.json', 'w') as f:
            json.dump(out, f, indent=1)

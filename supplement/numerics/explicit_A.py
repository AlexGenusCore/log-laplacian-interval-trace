"""Edge trace A from the explicit formula of Theorem 1(3) (note v2), pure-log case, no profile fit.

Pure log split (note v2, Sec. 1):  Q^log = E_VG + w_log,
  J(u) = e^{-|u|}/(2|u|)  (symbol 1/2 log(1+t^2)),   W f = -k*f,  k(u) = (1-e^{-|u|})/(2|u|).
For q > 0 and psi~(r) = psi(L-r), lambda = mu/2 (Q^log eigenvalue):
  A_+ = int_0^inf v_q(r) f_q(r) dr,
  f_q = 1_(0,2L) [(q+lambda) psi~ + k*psi~]  -  1_(2L,inf) mu_ext,  mu_ext(y) = int_0^2L J(y-r) psi~(r) dr.
v_q: Laplace transform 1/K(xi), K(xi) = exp((1/pi) int_0^inf log sigma_q(xi u) du/(1+u^2)),
sigma_q(t) = q + 1/2 log(1+t^2).  Continuing K through the imaginary axis gives
K(-tau + i0) = sigma_q(i tau + 0)/K(tau), hence the Stieltjes representation used here:
  v_q(x) = a* e^{-tau* x} + 1/2 int_1^inf e^{-x tau} K(tau)/D(tau) dtau,
  D = (q + 1/2 log(tau^2-1))^2 + pi^2/4,  tau* = sqrt(1-e^{-2q}),  a* = K(tau*) e^{-2q}/tau*.
The representation is verified numerically below (Laplace transform and total mass).
Below r_min (FEM edge cut) the tail int_0^{r_min} v f = f(0) V(r_min) + (q+lambda) int_0^{r_min} v psi~ is
evaluated with psi~ = A (V^-2 + delta)^-1/2 (delta from the FEM value at ell = 200) and solved
self-consistently for A (no fitted A enters). A_lin = first-order tail (A V^2/2), kept for comparison.
"""
import json, sys
import numpy as np
import scipy.special as sp

def logK(u, q, chunk=2000):
    """log K(tau), u = log tau (array). Poisson integral in x = log(s/tau), trapezoid (exponentially accurate)."""
    x = np.arange(-50, 50.0001, 0.04)
    wx = 0.04 / (2 * np.cosh(x)) / np.pi
    out = np.empty(len(u))
    for i in range(0, len(u), chunk):
        sl = u[i:i + chunk, None] + x[None, :]
        half = np.where(sl > 0, sl + 0.5 * np.log1p(np.exp(-2 * np.abs(sl))), 0.5 * np.log1p(np.exp(-2 * np.abs(sl))))
        out[i:i + chunk] = np.log(q + half) @ wx
    return out

def tau_rule(q, W1, panel=0.5, n=8):
    """Nodes/weights c_j so that 1/2 int_1^inf F(tau) K/D dtau ~ sum_j c_j F(tau_j); tau = 1 + e^w, w in [-50, W1]."""
    t, w = np.polynomial.legendre.leggauss(n)
    edges = np.linspace(-50, W1, int(round((W1 + 50) / panel)) + 1)   # must end exactly at W1 (tail_T starts there)
    panel = edges[1] - edges[0]
    wj = (edges[:-1, None] + panel / 2 * (t[None, :] + 1)).ravel()
    Wj = np.tile(w * panel / 2, len(edges) - 1)
    logtau = np.where(wj > 30, wj + np.log1p(np.exp(-np.abs(wj))), np.log1p(np.exp(np.minimum(wj, 30))))
    lk = logK(logtau, q)
    half_log = 0.5 * (wj + np.where(wj > 30, wj + np.log1p(2 * np.exp(-np.abs(wj))), np.log(2 + np.exp(np.minimum(wj, 30)))))
    D = (q + half_log) ** 2 + np.pi ** 2 / 4
    logc = np.log(0.5) + lk - np.log(D) + wj + np.log(Wj)
    return logtau, logc

def tail_T(q, W1):
    """T = 1/2 int_{tau > 1+e^{W1}} K/D dtau/tau  (dtau/tau = dw/(1+e^-w) ~ dw), to w = inf."""
    # log-spaced panels in w from W1 to 1e14, then analytic remainder (K ~ w^1/2, D ~ w^2)
    edges = W1 * np.logspace(0, np.log10(1e14 / W1), 400)
    t, w = np.polynomial.legendre.leggauss(10)
    a, b = edges[:-1, None], edges[1:, None]
    ww = ((b - a) / 2 * (t[None, :] + 1) + a).ravel()
    Ww = ((b - a) / 2 * w[None, :]).ravel()
    D = (q + ww) ** 2 + np.pi ** 2 / 4          # 1/2 log(tau^2-1) = w for these w
    val = 0.5 * np.sum(np.exp(logK(ww, q)) / D * Ww)
    return val + 1 / np.sqrt(edges[-1])

def check_rep(q):
    W1 = 300.0
    lt, lc = tau_rule(q, W1)
    T = tail_T(q, W1)
    ts = np.sqrt(-np.expm1(-2 * q)); ast = np.exp(logK(np.array([np.log(ts)]), q)[0]) * np.exp(-2 * q) / ts
    out = {}
    for xi in (1e-3, 0.3, 1.0, 10.0, 1e4, 1e40):
        lap = ast / (xi + ts) + np.sum(np.exp(lc - np.logaddexp(np.log(xi), lt))) + T
        out[f'LT xi={xi:g}: rep*K'] = float(lap * np.exp(logK(np.array([np.log(xi)]), q)[0]))
    out['int v * sqrt(q)'] = float((ast / ts + np.sum(np.exp(lc - lt)) + T) * np.sqrt(q))
    # pointwise v_q(r) (continuous part) against an independent mpmath quadrature
    import mpmath as mp
    mp.mp.dps = 20
    Kmp = lambda tau: mp.exp(mp.quad(lambda x: mp.log(q + 0.5 * mp.log(1 + (tau * mp.e ** x) ** 2)) / (2 * mp.cosh(x)),
                                     [-mp.inf, -5, 0, 5, mp.inf]) / mp.pi)
    for r in (0.05, 1.5):
        ref = mp.quad(lambda w: 0.5 * mp.e ** w * mp.e ** (-r * (1 + mp.e ** w)) * Kmp(1 + mp.e ** w)
                      / ((q + 0.5 * mp.log(mp.e ** w * (2 + mp.e ** w))) ** 2 + mp.pi ** 2 / 4), [-60, -10, -3, 0, 2, 5, 10, 20])
        out[f'v_c({r}) rule/mpmath'] = float(np.sum(np.exp(lc - np.exp(np.log(r) + lt))) / ref)
    return out

def V_of_r(q, r, W1):
    """renewal function V_q(r) = int_0^r v_q (exact representation; tail beyond W1 analytic)."""
    lt, lc = tau_rule(q, W1)
    ts = np.sqrt(-np.expm1(-2 * q)); ast = np.exp(logK(np.array([np.log(ts)]), q)[0]) * np.exp(-2 * q) / ts
    return ast * (-np.expm1(-ts * r)) / ts + np.sum(np.exp(lc - lt) * (-np.expm1(-np.exp(np.log(r) + lt)))) + tail_T(q, W1)

def delta_tail(A_main, g0, Vm, lam, q, r, psi, ell_star=200.0, W1pad=8.0):
    """Second-order edge tail. Below r_min use psi~ = A (V_q^-2 + delta)^-1/2, with delta read off the FEM
    value at ell = ell_star (no fit): delta = (A/psi(r*))^2 - V_q(r*)^-2. Then
    int_0^{r_min} v psi~ = A (sqrt(1 + delta Vm^2) - 1)/delta, solved self-consistently for A."""
    j = int(np.argmin(np.abs(-np.log(np.where(r > 0, r, 1.0)) - ell_star)))
    Vs = V_of_r(q, r[j], -np.log(r[1]) + W1pad)
    A = (A_main + g0 * Vm) / (1 - (q + lam) * Vm ** 2 / 2)
    for _ in range(40):
        dl = (A / psi[j]) ** 2 - Vs ** -2
        A = A_main + g0 * Vm + (q + lam) * A * (np.sqrt(1 + dl * Vm ** 2) - 1) / dl
    return float(A), float(dl)

def explicit_A(r, psi, L, lam, sign, qs=(0.5, 1.0, 2.0), delta=1e-5, ng=8, W1pad=8.0):
    ne = len(r) - 1
    h = np.diff(r)
    rmin = r[1]
    t, w = np.polynomial.legendre.leggauss(ng); t = (t + 1) / 2; w = w / 2
    # quadrature points on (0,2L): right half (coordinate r), left half (coordinate r' = 2L - r)
    loc = (r[:-1, None] + h[:, None] * t[None, :])
    pv = (psi[:-1, None] * (1 - t[None, :]) + psi[1:, None] * t[None, :])
    W = h[:, None] * w[None, :]
    rR, pR, wR = loc.ravel(), pv.ravel(), W.ravel()
    elR = np.repeat(np.arange(ne), ng)
    rL_p = rR.copy()                  # r' for left-half points
    pL = sign * pR
    # psi~ and coordinates of all points
    r_all = np.concatenate([rR, 2 * L - rL_p])
    rp_all = np.concatenate([2 * L - rR, rL_p])
    p_all = np.concatenate([pR, pL])
    w_all = np.concatenate([wR, wR])
    first = np.concatenate([elR == 0, np.zeros_like(elR, bool)])   # right element [0, r_min]: tail
    # g = k * psi~ : sources away from the edges + lumped edge pieces (0, delta)
    jd = np.searchsorted(r, delta)
    lump = np.sum(h[:jd] * (psi[:jd] + psi[1:jd + 1]) / 2)
    srcR = elR >= jd
    s_src = np.concatenate([rR[srcR], 2 * L - rL_p[srcR], [r[jd] / 2, 2 * L - r[jd] / 2]])
    q_src = np.concatenate([pR[srcR] * wR[srcR], pL[srcR] * wR[srcR], [lump, sign * lump]])
    def kfun(u):
        u = np.abs(u)
        return np.where(u > 0, -np.expm1(-u) / (2 * np.where(u > 0, u, 1)), 0.5)
    g_all = np.empty_like(r_all)
    for i in range(0, len(r_all), 4000):
        g_all[i:i + 4000] = kfun(r_all[i:i + 4000, None] - s_src[None, :]) @ q_src
    g0 = float(kfun(0 - s_src) @ q_src)
    keep = ~first
    lr = np.log(r_all)
    W1 = -np.log(rmin) + W1pad
    out = {}
    for q in qs:
        # v_q at all points: atom + sum_j c_j e^{-r tau_j}; only a window of nodes is non-trivial
        # (r tau < 1e-17: factor 1 -> prefix sum; r tau > 750: factor 0)
        lt, lc = tau_rule(q, W1)
        ts = np.sqrt(-np.expm1(-2 * q)); ast = np.exp(logK(np.array([np.log(ts)]), q)[0]) * np.exp(-2 * q) / ts
        c = np.exp(lc)
        C = np.concatenate([[0.0], np.cumsum(c)])
        lo = np.searchsorted(lt, np.log(1e-17) - lr)
        hi = np.searchsorted(lt, np.log(750.0) - lr)
        Wn = int((hi - lo).max())
        v_all = ast * np.exp(-ts * r_all) + C[lo]
        for i in range(0, len(r_all), 4000):
            idx = lo[i:i + 4000, None] + np.arange(Wn)[None, :]
            ok = idx < hi[i:i + 4000, None]
            idx = np.minimum(idx, len(lt) - 1)
            v_all[i:i + 4000] += np.where(ok, np.exp(lc[idx] - np.exp(lr[i:i + 4000, None] + lt[idx])), 0.0).sum(1)
        # N(r') = 1/2 int_0^inf v(2L+t) e^{-(t+r')}/(t+r') dt  (exterior-field term, Fubini)
        sel = lt < np.log(400 / L)
        tau_s, c_s = np.exp(lt[sel]), c[sel]
        N_all = 0.5 * ast * np.exp(-2 * L * ts + ts * rp_all) * sp.exp1((1 + ts) * rp_all)
        for i in range(0, len(rp_all), 2000):
            rp = rp_all[i:i + 2000, None]
            N_all[i:i + 2000] += 0.5 * (c_s[None, :] * np.exp(-tau_s[None, :] * (2 * L - rp)) * sp.exp1((tau_s[None, :] + 1) * rp)).sum(1)
        A_main = np.sum((w_all * v_all * ((q + lam) * p_all + g_all))[keep]) - np.sum(w_all * p_all * N_all)
        # V(r_min) and the self-consistent tail
        Vm = ast * (-np.expm1(-ts * rmin)) / ts + np.sum(np.exp(lc - lt) * (-np.expm1(-np.exp(np.log(rmin) + lt)))) + tail_T(q, W1)
        A_lin = (A_main + g0 * Vm) / (1 - (q + lam) * Vm ** 2 / 2)      # first-order tail
        A, dl = delta_tail(A_main, g0, Vm, lam, q, r, psi, W1pad=W1pad)   # second-order tail (used)
        out[q] = dict(A=A, A_lin=float(A_lin), delta=dl, A_main=float(A_main), g0=g0, Vm=float(Vm),
                      Vm_ell_half=float(Vm * np.sqrt(-np.log(rmin))), tail=float(A - A_main))
    return out

if __name__ == '__main__':
    if sys.argv[1] == 'check':
        for q in (0.5, 1.0, 2.0):
            print('q =', q, json.dumps(check_rep(q)))
        raise SystemExit
    for f in sys.argv[1:]:
        res = {}
        d = np.load(f)
        meta = json.load(open(f[:-4] + '.json'))
        L = meta['L']
        for sec, sign in (('even', 1), ('odd', -1)):
            for k, psi in enumerate(d['psi_' + sec]):
                lam = d['mu_' + sec][k] / 2
                for q, o in explicit_A(d['r'], psi, L, lam, sign).items():
                    o['2LA^2'] = 2 * L * o['A'] ** 2
                    res[f'{sec[0]}{k}|q={q}'] = o
                    print(f, sec[0] + str(k), 'q=%.1f' % q, ' '.join(f'{a}={b:.9g}' for a, b in o.items()), flush=True)
        json.dump(res, open(f[:-4] + '_explicit.json', 'w'), indent=1)

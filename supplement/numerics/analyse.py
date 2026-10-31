"""Collect the FEM / fit / explicit-formula outputs in out/ into results.json and the tables of numerics.md."""
import json, glob, os
import numpy as np
from fit_A import analyse as fitA
from explicit_A import delta_tail

EG = 0.57721566490153286061
HJSS = [-0.6878, 1.7331, 2.7275, 3.4122, 3.9091, 4.3142]
O = 'out/'
J = lambda t: json.load(open(O + t + '.json'))

def spectrum(d, n=6):
    """first n eigenvalues (both sectors) with labels, ascending."""
    lab = [(m, f'e{i}') for i, m in enumerate(d['mu_even_all'])] + [(m, f'o{i}') for i, m in enumerate(d['mu_odd_all'])]
    return sorted(lab)[:n]

R = {}
# (1) eigenvalues of L_Delta on (-1,1)
fine, n200, r125 = J('L1_n400_r1.12'), J('L1_n200_r1.12'), J('L1_n400_r1.25')
sp_f, sp_n, sp_r = spectrum(fine), spectrum(n200), spectrum(r125)
R['eig'] = [dict(k=k, state=s, mu=m, err=abs(m - sp_n[k][0]) + abs(m - sp_r[k][0]), hjss=HJSS[k], diff=m - HJSS[k])
            for k, (m, s) in enumerate(sp_f)]
# (2) scaling mu_k(L) + 2 log L - mu_k(1), absolute-mesh runs and relative-mesh fine runs
R['scaling'] = {}
base = spectrum(J('scal_L1.0'))
for L in (0.5, 2.0, 5.0):
    sp_ = spectrum(J(f'scal_L{L}'))
    R['scaling'][f'abs_mesh_L{L}'] = [sp_[k][0] + 2 * np.log(L) - base[k][0] for k in range(6)]
for L in (0.5, 2.0):
    sp_ = spectrum(J(f'L{L}_n200_r1.12'))
    b = spectrum(J('L1_n200_r1.12'))
    R['scaling'][f'rel_mesh_L{L}'] = [sp_[k][0] + 2 * np.log(L) - b[k][0] for k in range(6)]
# (3) edge trace from profile fits
order = ['e0', 'o0', 'e1', 'o1', 'e2', 'o2']
R['A_fit'] = {}
for L, n0 in (('0.5', 200), ('1', 400), ('2.0', 200)):
    Lv = float(L)
    runs = {}
    for tag in (f'L{L}_n{n0}_r1.12', f'L{L}_n{n0}_r1.25', f'L{L}_n{n0}_r1.5', f'L{L}_n{n0 // 2}_r1.12'):
        runs[tag] = fitA(O + tag + '.npz')
    t0, t1, t2, t3 = runs
    rows = []
    for st in order:
        F1, F2, F4 = runs[t0][st]['F1'], runs[t0][st]['F2'], runs[t0][st]['F4']
        x = lambda t: 2 * Lv * runs[t][st]['F1'] ** 2
        err = abs(x(t0) - x(t1)) + 2 * Lv * max(abs(F1 ** 2 - F2 ** 2), abs(F1 ** 2 - F4 ** 2)) + abs(x(t0) - x(t3))
        rows.append(dict(state=st, A=F1, A_F2=F2, A_F3=runs[t0][st]['F3'], A_F4=F4, ratio=x(t0), ratio_rho125=x(t1), ratio_rho15=x(t2),
                         ratio_n0half=x(t3), err=err, slope200_over_2L=runs[t0][st].get('slope200', np.nan) / (2 * Lv),
                         c_plus_lam_plus_gamma=runs[t0][st]['F1_c'] + J(t0)['mu_' + ('even' if st[0] == 'e' else 'odd')][int(st[1])] / 2 + EG,
                         c_plus_lam_plus_gamma_rho125=runs[t1][st]['F1_c'] + J(t1)['mu_' + ('even' if st[0] == 'e' else 'odd')][int(st[1])] / 2 + EG))
    R['A_fit'][L] = rows
# (3b) explicit formula
R['A_explicit'] = {}
for f in sorted(glob.glob(O + '*_explicit.json')):
    tag = os.path.basename(f)[:-len('_explicit.json')]
    Lv = J(tag)['L']
    e = json.load(open(f))
    dd = np.load(O + tag + '.npz')
    R['A_explicit'][tag] = {}
    for st in order:
        sec, k = ('even' if st[0] == 'e' else 'odd'), int(st[1])
        lam = float(dd['mu_' + sec][k]) / 2
        row = dict(ratio_q={}, ratio_q_lin={}, delta_q={})
        for q in (0.5, 1.0, 2.0):
            o = e.get(f'{st}|q={q}')
            if o is None:
                continue
            A, dl = delta_tail(o['A_main'], o['g0'], o['Vm'], lam, q, dd['r'], dd['psi_' + sec][k])
            A_lin = o.get('A_lin', o['A'])          # old-format output: 'A' was the first-order value
            row['ratio_q'][q] = 2 * Lv * A ** 2
            row['ratio_q_lin'][q] = 2 * Lv * A_lin ** 2
            row['delta_q'][q] = dl
        if row['ratio_q']:
            R['A_explicit'][tag][st] = row
# (4) gaps
allmu = sorted(fine['mu_even_all'] + fine['mu_odd_all'])
gaps = np.diff(allmu)
labels = [s for _, s in sorted([(m, f'e{i}') for i, m in enumerate(fine['mu_even_all'])] + [(m, f'o{i}') for i, m in enumerate(fine['mu_odd_all'])])]
R['gaps'] = dict(labels=labels, alternating=all(labels[i][0] != labels[i + 1][0] for i in range(len(labels) - 1)), mu=allmu, min_gap=float(gaps.min()), argmin=int(gaps.argmin()),
                 min_rel_gap=float((gaps / np.abs(allmu[1:])).min()),
                 min_sector_gap=float(min(np.diff(fine['mu_even_all']).min(), np.diff(fine['mu_odd_all']).min())),
                 err_top=float(max(abs(a - b) for a, b in zip(sorted(fine['mu_even_all'] + fine['mu_odd_all']), sorted(n200['mu_even_all'] + n200['mu_odd_all'])))))
# (5) mixed identity  L (lambda_e - lambda_o) Phi'(0) = 1  (lambda = mu/2, A_+ > 0 convention)
R['mixed'] = {}
for L, n0 in (('0.5', 200), ('1', 400), ('2.0', 200)):
    d, d2 = J(f'L{L}_n{n0}_r1.12'), J(f'L{L}_n{n0}_r1.25')
    rows = {}
    for i in range(3):
        for j in range(3):
            val = lambda dd: dd['L'] * (dd['mu_even'][i] - dd['mu_odd'][j]) / 2 * dd['phi_prime'][f'e{i}o{j}']
            rows[f'e{i}o{j}'] = dict(value=val(d), value_rho125=val(d2), phi_prime=d['phi_prime'][f'e{i}o{j}'])
    R['mixed'][L] = rows
# (5b) residual of the mixed identity explained by the edge cut r_min (FEM ramps psi to 0 on [0, r_min]);
# model psi_x(L-r) = (2L)^-1/2 (ell - lambda_x - gamma)^-1/2 below r_min (a model estimate, not used above)
def cut_residual(L, rmin, lam_e, lam_o):
    l = -np.log(rmin); a, b = -lam_o - EG, -lam_e - EG
    AA = 1 / (2 * L)
    tail = 0.5 * AA * 2 / (b - a) * (1 - np.sqrt((l + a) / (l + b)))
    ramp = 0.5 * AA / np.sqrt((l + a) * (l + b))
    return -L * (lam_e - lam_o) * 2 * (ramp - tail)
R['mixed_rmin'] = {}
for tag in ('L1_n400_r1.12_rm1e-100', 'L1_n400_r1.12', 'L1_n400_r1.12_rm1e-280'):
    if not os.path.exists(O + tag + '.json'):
        continue
    d = J(tag)
    R['mixed_rmin'][tag] = {f'e{i}o{j}': dict(
        value=d['L'] * (d['mu_even'][i] - d['mu_odd'][j]) / 2 * d['phi_prime'][f'e{i}o{j}'],
        model_cut_residual=cut_residual(d['L'], d['rmin'], d['mu_even'][i] / 2, d['mu_odd'][j] / 2))
        for i in range(3) for j in range(3)}

# ---------------------------------------------------------------- LaTeX tables
def e(x, d=1):
    """compact a(b)-style error: 1.0000001 +- 2e-6 -> '1.0000001 (20)' is overkill; use explicit pm"""
    return f'{x:.1e}'.replace('e-0', 'e-')

T = []
T.append(r'% Table 1: eigenvalues of L_Delta on (-1,1)')
T.append(r'\begin{tabular}{rlllr}\toprule $k$ & sector & $\mu_k$ (this work) & HJSS & $\mu_k-\mu_k^{\rm HJSS}$\\\midrule')
for row in R['eig']:
    T.append(rf"{row['k']} & {row['state']} & ${row['mu']:.8f}$ & ${row['hjss']:.4f}$ & ${row['diff']:+.1e}$\\")
T.append(r'\bottomrule\end{tabular}')
T.append(r'% Table 2: scaling check  mu_k(L) + 2 log L - mu_k(1)  (same absolute mesh: h0=0.01, rho=1.25, r_min=1e-60)')
T.append(r'\begin{tabular}{rrrr}\toprule $k$ & $L=0.5$ & $L=2$ & $L=5$\\\midrule')
for k in range(6):
    T.append(f'{k} & ' + ' & '.join(f"${R['scaling'][f'abs_mesh_L{L}'][k]:+.1e}$" for L in (0.5, 2.0, 5.0)) + r'\\')
T.append(r'\bottomrule\end{tabular}')
T.append(r'% Table 3: 2 L A^2 from the edge profile (fit F1, finest mesh) with mesh/fit error bar')
T.append(r'\begin{tabular}{l' + 'r' * 3 + r'}\toprule state & $L=0.5$ & $L=1$ & $L=2$\\\midrule')
for i, st in enumerate(order):
    T.append(f'{st} & ' + ' & '.join(f"${R['A_fit'][L][i]['ratio']:.7f}\\pm{R['A_fit'][L][i]['err']:.0e}$" for L in ('0.5', '1', '2.0')) + r'\\')
T.append(r'\bottomrule\end{tabular}')
if 'L1_n400_r1.12' in R['A_explicit']:
    T.append(r'% Table 4: 2 L A^2 from the explicit formula of Thm 1(3); mean over q = 0.5, 1, 2 on the finest mesh;')
    T.append(r'% error = half the q-spread + change from rho = 1.25 to 1.12')
    T.append(r'\begin{tabular}{lrrr}\toprule state & $L=0.5$ & $L=1$ & $L=2$\\\midrule')
    R['A_explicit_summary'] = {}
    for st in order:
        cells = []
        for fine_t, coarse_t in (('L0.5_n200_r1.12', 'L0.5_n200_r1.25'), ('L1_n400_r1.12', 'L1_n400_r1.25'), ('L2.0_n200_r1.12', 'L2.0_n200_r1.25')):
            v = R['A_explicit'].get(fine_t, {}).get(st)
            if not v:
                cells.append('--'); continue
            vals = list(v['ratio_q'].values())
            m = float(np.mean(vals))
            vc = R['A_explicit'].get(coarse_t, {}).get(st)
            dr = abs(m - float(np.mean(list(vc['ratio_q'].values())))) if vc else float('nan')
            err = (max(vals) - min(vals)) / 2 + dr
            R['A_explicit_summary'][f'{fine_t}|{st}'] = dict(mean=m, err=err, q_spread=max(vals) - min(vals), rho_change=dr)
            cells.append(f'${m:.7f}\\pm{err:.0e}$')
        T.append(f'{st} & ' + ' & '.join(cells) + r'\\')
    T.append(r'\bottomrule\end{tabular}')
T.append(r'% Table 5: mixed identity  L (lambda_e - lambda_o) Phi^\prime(0) - 1  (lambda = mu/2), no fitting')
T.append(r'\begin{tabular}{lrrr}\toprule pair & $L=0.5$ & $L=1$ & $L=2$\\\midrule')
for key in R['mixed']['1']:
    T.append(f'{key} & ' + ' & '.join(f"${R['mixed'][L][key]['value'] - 1:+.1e}$" for L in ('0.5', '1', '2.0')) + r'\\')
T.append(r'\bottomrule\end{tabular}')
open('tables.tex', 'w').write('\n'.join(T) + '\n')
json.dump(R, open('results.json', 'w'), indent=1, default=float)

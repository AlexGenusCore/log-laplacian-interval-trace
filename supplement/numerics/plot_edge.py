"""Figure 1 of the paper: (2L)^{1/2} psi(L-r) ell(r)^{1/2} against ell = log(1/r), L = 1.

Usage: python plot_edge.py [out/L1_n400_r1.12.npz] [../fig_edge.pdf]
The window stops at ell_max - 20 (as in fit_A.py), because the last nodes feel the cut at r_min.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

src = sys.argv[1] if len(sys.argv) > 1 else 'out/L1_n400_r1.12.npz'
dst = sys.argv[2] if len(sys.argv) > 2 else '../fig_edge.pdf'
L = 1.0
COLORS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300']  # dataviz categorical slots 1-6

d = np.load(src)
o = np.argsort(d['r'])
r = d['r'][o]
m = (r > 0) & (r <= np.exp(-2))
ell = -np.log(r[m])
keep = ell <= ell.max() - 20
ell = ell[keep]

plt.rcParams.update({'font.size': 8, 'axes.linewidth': 0.6, 'font.family': 'serif', 'mathtext.fontset': 'cm'})
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.2, 2.3), constrained_layout=True)
states = [('even', 0, 'e_0'), ('odd', 0, 'o_0'), ('even', 1, 'e_1'),
          ('odd', 1, 'o_1'), ('even', 2, 'e_2'), ('odd', 2, 'o_2')]
for c, (sec, k, name) in zip(COLORS, states):
    psi = d['psi_' + sec][k][o][m][keep]
    y = np.sqrt(2 * L) * psi * np.sqrt(ell)
    ls = '-' if sec == 'even' else '--'
    ax1.plot(ell, y, ls, color=c, lw=1.2, label=f'${name}$')
    ax2.loglog(ell, np.abs(y - 1), ls, color=c, lw=1.2)
ax1.axhline(1.0, color='0.45', lw=0.6, zorder=0)
ax1.set_xscale('log')
ax1.set_ylim(0.95, 2.05)
ax1.set_xlabel(r'$\ell(r)=\log(1/r)$')
ax1.set_ylabel(r'$(2L)^{1/2}\,\psi(L-r)\,\ell(r)^{1/2}$')
ax1.legend(ncol=2, frameon=False, fontsize=7, handlelength=2.2)
lr = np.array([10.0, 400.0])
ax2.loglog(lr, 0.4 / lr, ':', color='0.3', lw=0.9)
ax2.text(40, 0.4 / 40 * 0.45, r'$\propto\ell^{-1}$', color='0.3', fontsize=7)
ax2.set_xlabel(r'$\ell(r)$')
ax2.set_ylabel(r'$|(2L)^{1/2}\psi(L-r)\ell(r)^{1/2}-1|$')
for ax in (ax1, ax2):
    ax.grid(True, which='major', color='0.9', lw=0.5)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
fig.savefig(dst)

# self-check: every state is within 0.5% of the predicted limit 1 at ell = 400
for sec, k, name in states:
    psi = d['psi_' + sec][k][o][m][keep]
    y = np.sqrt(2 * L) * psi * np.sqrt(ell)
    assert abs(np.interp(400, ell[::-1], y[::-1]) - 1) < 5e-3, name
print('wrote', dst)

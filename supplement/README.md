# Numerical supplement

Supplement to *A boundary trace identity for the logarithmic Laplacian on an interval, and simplicity of its
Dirichlet spectrum*, by Johann Alexander Botha (contact: info@genuscore.com), draft of 31 October 2026.

Everything here is an **uncertified floating-point computation**. The error estimates are refinement
estimates, not rigorous error bounds. None of it is used in the proofs.

## Contents

| Path | What it is |
|---|---|
| `numerics/fem.py` | P1 Galerkin FEM for $L_\Delta$ on $(-L,L)$: assembly, solver, `--selftest`. Writes `out/<tag>.json`, `out/<tag>.npz`; stdout goes to `out/<tag>.out`. |
| `numerics/fit_A.py` | Edge-profile fits F1–F4 for $A_+$ (no eigenvalue derivative, no Hadamard formula). Read by `analyse.py`. |
| `numerics/explicit_A.py` | $A_+=\int_0^\infty v_qf_q$ (paper, Thm 5.1(ii)) for $q=0.5,1,2$, with the tail model below $r_{\min}$. `python explicit_A.py check` tests the $v_q$ representation (`out/explicit_check.out`). |
| `numerics/analyse.py` | Collects `out/` into `results.json` and `tables.tex`. |
| `numerics/plot_edge.py` | Figure 1 (`fig_edge.pdf`). |
| `numerics/residual_check.py` | Added for this supplement: recomputes the row residuals of the stored eigenvectors at the deepest mesh nodes (see below). `--test` compares its row assembly with `fem.assemble`. |
| `numerics/run_grid.sh`, `numerics/cases.txt` | Driver and case list of the main grid (28 runs). |
| `numerics/out/run_grid2.sh`, `numerics/out/cases2.txt` | Driver and case list of the two runs with $r_{\min}=10^{-100},10^{-280}$. |
| `numerics/out/run_grid.log`, `run_grid2.log` | Completion times of the grid runs. |
| `numerics/out/<tag>.json` | Per run: $L$, $n_0$, $\rho_{\rm g}$ (`rho`), $r_{\min}$, `N_half` (number of half-mesh nodes), timings, refined eigenvalues `mu_even`/`mu_odd` (3 per sector; 1 for `scal_*`), the 10 lowest eigenvalues per sector from the dense solver (`mu_even_all`/`mu_odd_all`), the maximal row residual `rowres_even`/`rowres_odd` of each refined eigenvector, and $\Phi'(0)$ for the pairs $(e_i,o_j)$ (`phi_prime`). Eigenvalues are those of $L_\Delta$ (symbol $2\log|\xi|$); $\lambda=\mu/2$. |
| `numerics/out/<tag>.npz` | Per run: half-mesh nodes `r` ($r=L-x$, distance to the right edge), nodal values `psi_even`, `psi_odd` (normalised in $L^2(-L,L)$, sign with $\psi>0$ near $x=L$), `mu_even`, `mu_odd`. |
| `numerics/out/<tag>.out` | Stdout of `fem.py` for that run (the same numbers as the JSON). |
| `numerics/out/<tag>_explicit.json/.out` | Output of `explicit_A.py` for that run (9 runs). |
| `numerics/out/explicit_check.out` | Output of `explicit_A.py check`. |
| `numerics/results.json`, `numerics/tables.tex` | Output of `analyse.py`. |
| `numerics/numerics.md` | The working notes of the computation (31 October 2026), unchanged except that the name of the virtual machine was removed. See "Notes on `numerics.md`" below. |
| `checks/referee_p0check.py` | Independent piecewise-constant (P0) Galerkin check written by the internal AI referee. |
| `checks/p0check_output.txt` | Its output for the six configurations of the internal referee report, re-run on 31 October 2026 (the original output was not saved as a file; the re-run reproduces the referee's table: e.g. $2LA_+^2=1.00134$, $1.00027$, $1.00016$ for $e_0$ in rows 1, 3, 4). |
| `checks/residuals_deepest.txt` | Output of `residual_check.py` for every stored run. |
| `SHA256SUMS` | SHA-256 of every file in this folder (`sha256sum -c SHA256SUMS`). |

Run tags: `L<L>_n<n0>_r<rho>` with $r_{\min}=10^{-200}$; suffix `_rm1e-100`/`_rm1e-280` for other edge cuts;
`scal_L<L>` for the scaling check ($h_0=0.01$, $\rho_{\rm g}=1.25$, $r_{\min}=10^{-60}$, one eigenvector per sector).

## Environment

- **Recorded runs** (`out/`): 31 October 2026, on a Linux virtual machine with 10 vCPUs, in a Python
  virtual environment whose package versions **were not recorded** at run time. `run_grid.sh` ran two jobs
  in parallel with 5 BLAS threads each; `run_grid2.sh` one job with 3 threads.
- **Required packages:** `numpy`; `scipy` (`fem.py` uses `scipy.linalg.eigh(..., subset_by_index=...)`, which
  needs SciPy >= 1.5, and `lu_factor`/`lu_solve`; `explicit_A.py` uses `scipy.special`); `mpmath` (only in
  `explicit_A.py check`); `matplotlib` (only `plot_edge.py`). `python-flint` is **not** used by any script here.
- **Re-checks for this supplement** (`checks/p0check_output.txt`, `checks/residuals_deepest.txt`): Python 3.12.13,
  numpy 2.5.4, on the author's laptop. SciPy and matplotlib were not installed there, so `fem.py`, `explicit_A.py`,
  `analyse.py` and `plot_edge.py` were **not** re-run for this supplement. (`residual_check.py` substitutes a stub
  for `scipy.linalg` when SciPy is missing; it does not use it.)

## How to run

From `numerics/`. The two driver scripts run the interpreter in the environment variable `PY` (default `python3`);
set it to a Python with numpy and SciPy if needed, e.g. `export PY=/path/to/venv/bin/python`. (For this
supplement the drivers were edited in two respects only: the interpreter, hard-coded in the recorded runs, is now
taken from `PY`, and the usage comments no longer name the machine or its paths. Otherwise they are as run.)

```bash
python fem.py --selftest                       # assembled sums for u=1, u=r against exact values, to 1e-11
bash run_grid.sh                               # main grid, cases.txt -> out/
cp out/run_grid2.sh out/cases2.txt . && bash run_grid2.sh   # the two r_min runs
python explicit_A.py check > out/explicit_check.out
for t in L0.5_n200_r1.12 L0.5_n200_r1.25 L1_n200_r1.12 L1_n400_r1.12 L1_n400_r1.12_rm1e-280 \
         L1_n400_r1.25 L1_n400_r1.5 L2.0_n200_r1.12 L2.0_n200_r1.25; do
  python explicit_A.py out/$t.npz > out/${t}_explicit.out; done
python analyse.py                              # -> results.json, tables.tex
python plot_edge.py out/L1_n400_r1.12.npz fig_edge.pdf
python residual_check.py --test && python residual_check.py out/L1_n400_r1.12   # any tags
python ../checks/referee_p0check.py 1 1e-100 1.12 400                            # L rmin rho n0
```

(`run_grid2.sh` and `cases2.txt` were run from the folder that contains `fem.py` and moved to `out/` afterwards.)
A single run, e.g. the finest mesh of Table 1:
`python fem.py --L 1 --n0 400 --rho 1.12 --rmin 1e-200 --nvec 3 --nev 10 --out out/L1_n400_r1.12`.

**Run times** (recorded, assembly + solve, from the JSON files): under 2 s for `scal_*` and $\rho_{\rm g}=2$;
3–14 s for $\rho_{\rm g}=1.5,1.25$; 18 s for `L1_n400_r1.12_rm1e-100`; 37–48 s for $\rho_{\rm g}=1.12$ ($N_{\rm half}\approx4100$–$4400$); 141 s for
`L1_n400_r1.12_rm1e-280` ($N_{\rm half}=6049$). Run times of `explicit_A.py` were not recorded. `fit_A.py` takes
well under a second per file. `referee_p0check.py`: 2–14 s per configuration, peak RSS 116–437 MB.
`residual_check.py`: under a second per run.

**Memory.** Peak memory of `fem.py` was not recorded. It holds about ten dense $N_{\rm half}\times N_{\rm half}$
float64 arrays ($8N_{\rm half}^2$ bytes each: 157 MB at $N_{\rm half}=4424$, 293 MB at 6049) plus the assembly
index arrays, so the finest runs need several GB.

## Residuals

**Definition** (as in `fem.py`). In each parity sector the discrete problem is $K x=\mu M x$ with
$K=-(F_{\rm same}\pm F_{\rm cross})-2\gamma M$ on the free nodes. For a refined eigenpair $(\mu,x)$, with $x$
normalised in $L^2(-L,L)$, the **row residual** of row $i$ is
$$\mathrm{res}_i=\frac{\bigl|((K-\mu M)x)_i\bigr|}{M_{ii}},$$
i.e. the residual of the row-equilibrated system $\mathrm{diag}(M)^{-1}(K-\mu M)$, measured in the units of the
nodal values. Since $M_{ii}\asymp h_i$, this is not diluted by the tiny cells near the edge: an error in the
nodal values at $r\approx10^{-200}$ shows up in these rows at full size. `rowres_*` in the JSON files is
$\max_i\mathrm{res}_i$ over **all** rows. `residual_check.py` rebuilds the rows of the five deepest free nodes
($r=r_{\min},\,r_{\min}\rho_{\rm g},\dots$) independently of the stored matrices and evaluates them on the stored
vectors.

**Values** (`checks/residuals_deepest.txt`; maximum over the stored eigenvectors of the run):

| run | $N_{\rm half}$ | $r_{\min}$ | max over all rows (`fem.py`) | row of node $r=r_{\min}$ | max over 5 deepest nodes |
|---|---|---|---|---|---|
| L0.5_n100_r1.5 | 1225 | 1e-200 | 1.47e-13 | 6.51e-14 | 8.56e-14 |
| L0.5_n100_r1.25 | 2145 | 1e-200 | 2.07e-13 | 1.91e-13 | 1.91e-13 |
| L0.5_n100_r1.12 | 4130 | 1e-200 | 2.92e-13 | 1.45e-13 | 2.14e-13 |
| L0.5_n200_r1.5 | 1323 | 1e-200 | 1.59e-13 | 1.35e-13 | 1.35e-13 |
| L0.5_n200_r1.25 | 2242 | 1e-200 | 2.29e-13 | 1.39e-13 | 1.39e-13 |
| L0.5_n200_r1.12 | 4224 | 1e-200 | 2.92e-13 | 2.11e-13 | 2.11e-13 |
| L1_n100_r2.0 | 759 | 1e-200 | 6.59e-14 | 5.40e-14 | 5.40e-14 |
| L1_n100_r1.5 | 1227 | 1e-200 | 1.12e-13 | 6.19e-14 | 6.35e-14 |
| L1_n100_r1.25 | 2148 | 1e-200 | 1.41e-13 | 7.49e-14 | 9.49e-14 |
| L1_n100_r1.12 | 4136 | 1e-200 | 1.84e-13 | 8.29e-14 | 1.25e-13 |
| L1_n200_r2.0 | 858 | 1e-200 | 7.89e-14 | 7.66e-14 | 7.66e-14 |
| L1_n200_r1.5 | 1325 | 1e-200 | 1.10e-13 | 6.99e-14 | 8.37e-14 |
| L1_n200_r1.25 | 2245 | 1e-200 | 1.48e-13 | 1.09e-13 | 1.22e-13 |
| L1_n200_r1.12 | 4230 | 1e-200 | 1.90e-13 | 1.13e-13 | 1.53e-13 |
| L1_n400_r2.0 | 1057 | 1e-200 | 1.35e-13 | 6.74e-14 | 6.74e-14 |
| L1_n400_r1.5 | 1523 | 1e-200 | 1.29e-13 | 4.39e-14 | 1.25e-13 |
| L1_n400_r1.25 | 2442 | 1e-200 | 1.39e-13 | 7.89e-14 | 7.89e-14 |
| L1_n400_r1.12 | 4424 | 1e-200 | 1.80e-13 | 7.77e-14 | 1.23e-13 |
| L1_n400_r1.12_rm1e-100 | 2392 | 1e-100 | 1.06e-13 | 7.38e-14 | 7.38e-14 |
| L1_n400_r1.12_rm1e-280 | 6049 | 1e-280 | 2.83e-13 | 1.61e-13 | 1.61e-13 |
| L2.0_n100_r1.5 | 1228 | 1e-200 | 7.47e-14 | 3.00e-14 | 4.10e-14 |
| L2.0_n100_r1.25 | 2151 | 1e-200 | 1.02e-13 | 3.33e-14 | 5.19e-14 |
| L2.0_n100_r1.12 | 4142 | 1e-200 | 1.36e-13 | 1.23e-13 | 1.23e-13 |
| L2.0_n200_r1.5 | 1327 | 1e-200 | 8.94e-14 | 5.66e-14 | 5.66e-14 |
| L2.0_n200_r1.25 | 2248 | 1e-200 | 1.12e-13 | 4.74e-14 | 5.90e-14 |
| L2.0_n200_r1.12 | 4236 | 1e-200 | 1.29e-13 | 1.15e-13 | 1.15e-13 |
| scal_L0.5 | 653 | 1e-60 | 5.22e-14 | 3.36e-14 | 3.36e-14 |
| scal_L1.0 | 703 | 1e-60 | 6.03e-14 | 5.64e-14 | 5.64e-14 |
| scal_L2.0 | 803 | 1e-60 | 3.30e-14 | 3.10e-15 | 1.43e-14 |
| scal_L5.0 | 1103 | 1e-60 | 1.91e-14 | 6.91e-15 | 6.91e-15 |

The largest value, $2.92\cdot10^{-13}$, is the "$\le3\cdot10^{-13}$" of the paper (Section 8). For comparison, the
stored nodal values at the deepest node are about $A_+\ell(r_{\min})^{-1/2}\approx0.03$ ($L=1$). The change between
the last two inverse-iteration steps was not stored. The residual is a check of the algebraic eigenproblem only;
it says nothing about the discretisation error, the quadrature error of the matrix entries, or the edge cut.

## How the paper's tables, figure and conjecture evidence are produced

Paper numbering (`main.pdf` of 31 October 2026, revised after the external review).

- **Table 1** (eigenvalues on $(-1,1)$): `results.json` key `eig` (also `tables.tex`, first table): $\mu_k$ from
  `out/L1_n400_r1.12.json`; error estimate = change to `L1_n200_r1.12` + change to `L1_n400_r1.25`. The HJSS column
  is copied from arXiv:2311.13079v3; the bracket column is $[2\log\frac{(k-1)\pi}2,2\log\frac{k\pi}2]$, computed by
  hand. The computed values approximate the exact Galerkin eigenvalues, which are upper bounds; quadrature and
  rounding errors are not bounded rigorously.
- **Section 8 eigenvalue text:** scaling check `results.json` key `scaling` (`rel_mesh_*`: nearly dilated meshes,
  $\le1.2\cdot10^{-10}$; `abs_mesh_*`: `scal_*` runs); alternation and gaps: key `gaps`.
- **Table 2** ($10^7(2LA_+^2-1)$): edge-profile columns from `results.json` key `A_fit` (`ratio` and `err`, at
  $L=0.5,1,2$ with $n_0=200,400,200$, $\rho_{\rm g}=1.12$); explicit-formula columns from key `A_explicit_summary`
  (`mean`, `err`), built from `out/*_explicit.json`. The entries were multiplied by $10^7$ and rounded by hand. The
  error estimates are refinement estimates (see the caption): they exclude the edge cut $r_{\min}$ and, for the
  explicit formula, the error of the tail model. They are supporting evidence, not error bounds. Convergence in
  $\rho_{\rm g}$ and in $r_{\min}$: keys `A_fit` (`ratio_rho125`, `ratio_rho15`) and `A_explicit` (per run, per $q$).
- **Mixed identity** (Section 8): keys `mixed` and `mixed_rmin` (the latter with the model of the cut).
- **Figure 1:** `python plot_edge.py out/L1_n400_r1.12.npz fig_edge.pdf`.
- **Evidence for Conjecture 9.2** ($\beta/\alpha+\mu/2+\gamma$ from fit F1): $\rho_{\rm g}=1.12$ and $1.25$:
  `results.json` keys `A_fit[*][*].c_plus_lam_plus_gamma` and `..._rho125` (range $-0.00114$ to $-0.00107$, and
  $-0.00422$ to $-0.00415$, over the 18 cases); $\rho_{\rm g}=1.5$: `fit_A.analyse` on `out/L*_r1.5.npz`, value
  `F1_c` $+\mu/2+\gamma$ ($-0.0138$ to $-0.0137$). The independent P0 values (between $-1.3\cdot10^{-3}$ and
  $+1.5\cdot10^{-3}$) are the column `beta/alpha+mu/2+gamma` of `checks/p0check_output.txt`.

## Notes on `numerics.md`

`numerics.md` is the working record and is included unchanged, except that in two places the internal name of the virtual machine was replaced by a description ("a Linux virtual machine with 10 vCPUs", "the drivers of the recorded runs"). It uses the numbering of an internal draft
("note v2"): its Theorem 4 is Theorems 1.2 and 1.4 of the paper, its Theorem 1(3) is Theorem 5.1(ii), its
Theorem 3(c) is Corollary 1.5, and its Lemmas 8.2–8.3 are Lemmas 6.2–6.3. Where it says that the computed
eigenvalues are "upper bounds up to rounding", the paper's revised statement applies: the exact Galerkin
eigenvalues are upper bounds, and the computed values approximate them, with quadrature and rounding errors not
bounded rigorously.

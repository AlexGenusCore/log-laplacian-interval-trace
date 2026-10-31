# Numerics for the pure-log paper: independent check

*2026-10-31. Folder `papers/purelog/numerics/`. Uncertified double-precision computations, written from scratch for this paper. They import no earlier project code; the earlier codes were read only to fix normalisations. Nothing here is a proof.*

## 0. Normalisation

- $L_\Delta$ is the logarithmic Laplacian of Chen–Weth (symbol $2\log|\xi|$). Note v2 uses the half form $Q^{\log}=\frac12L_\Delta$ (symbol $\log|t|$). Eigenvalues: $\mu_k$ for $L_\Delta$ and $\lambda_k=\mu_k/2$ for $Q^{\log}$. The eigenfunctions are the same for both.
- $\psi$ is real and normalised in $L^2(-L,L)$, with sign chosen so that $\psi>0$ near $x=L$. The edge trace is $A=A_+=\lim_{r\downarrow0}\psi(L-r)\,\ell(r)^{1/2}$, where $\ell(r)=\log(1/r)$.
- What note v2 predicts:
  - Theorem 4: $A^2=1/(2L)$, i.e. $2LA^2=1$, for every eigenfunction, and all eigenvalues are simple.
  - Dilation: $\mu_k(L)=\mu_k(1)-2\log L$.
  - Mixed identity, Theorem 3(c) combined with Theorem 4: $L(\lambda_e-\lambda_o)\Phi'(0)=2LA_eA_o=1$. Here $\Phi'(0)=\langle\psi_o,\psi_e'\rangle$, both traces are taken positive, and $\lambda=\mu/2$. In $L_\Delta$ units this reads $L(\mu_e-\mu_o)\Phi'(0)=2$.

## 1. Method (`fem.py`)

**The form.** The code uses an exact finite-part representation. For $u,v$ supported in $\Omega=(-L,L)$,
$$a(u,v)=-\operatorname{FP}\!\iint_{\Omega^2}\frac{u(x)v(y)}{|x-y|}\,dx\,dy-2\gamma\!\int_\Omega uv,\qquad \operatorname{FP}\!\iint:=\lim_{\varepsilon\to0}\Big[\iint_{|x-y|>\varepsilon}+2\log\varepsilon\!\int uv\Big].$$
This is exact, because $\widehat{\operatorname{Pf}(1/|x|)}=-2\gamma-2\log|\xi|$. It involves no killing integral and no truncation radius. It is a different representation from those of the earlier hp-FEM (the VG split) and of the referee's P0 code (the cutoff $|x-y|<2L$ with a killing term).

**Discretisation.**
- P1 (hat-function) Galerkin, reduced by parity to the half mesh $r\in[0,L]$, where $r=L-x$ is the distance to the right edge. Storing $r$ rather than $x$ means tiny $r$ is held without cancellation.
- Even and odd blocks: $K_\pm=-(F_{\rm same}\pm F_{\rm cross})-2\gamma M$. The kernels are $1/|r-r'|$ and $1/(2L-r-r')$.
- Mesh: geometric grading $r_{\min}\rho^k$ towards the edge, then a uniform interior step $h_0=L/n_0$. Base runs use $r_{\min}=10^{-200}$, i.e. $\ell$ up to 460.5. The node at $r=0$ is Dirichlet.

**Element integrals, all to about $10^{-16}$ relative.**
- Coincident elements: closed-form finite part, $h[\operatorname{FP}_1+2\log h\,M_{\rm loc}]$ with $\operatorname{FP}_1=-\frac19\begin{psmallmatrix}8&1\\1&8\end{psmallmatrix}$.
- Adjacent elements, and the cross-kernel corner element at $x=0$: closed forms in $h,k$, evaluated in units of $\max(h,k)$.
- All other pairs: tensor Gauss–Legendre, with the order chosen from the separation ratio.
- `fem.py --selftest` checks the assembled sums on $u\equiv1$ and $u=r$ against exact values, to $10^{-11}$.

**Solver.**
- Dense generalised `eigh` on the diagonally scaled pencil.
- The eigenvector components on cells with $\sqrt{M_{jj}}<10^{-7}$ are dropped, because eigh's absolute error there is about $\varepsilon/\sqrt{M_{jj}}$.
- Three steps of inverse iteration on the **row-equilibrated** system $\operatorname{diag}(M)^{-1}(K-\mu M)$ then give every nodal value, down to $r=10^{-200}$, to full relative accuracy. The maximal row residual is $\le3\cdot10^{-13}$ in every run.

**Bounds.** The P1 space lies in the form domain and the integrals are exact up to rounding. So each computed $\mu_k$ is an upper bound for the exact eigenvalue, sector by sector, up to rounding.

**Grid.** $L\in\{0.5,1,2\}$, $n_0\in\{100,200,400\}$, $\rho\in\{2,1.5,1.25,1.12\}$, 10 eigenvalues per sector, 3 refined eigenvectors per sector. Also, at $L=1$, $n_0=400$, $\rho=1.12$, edge cuts $r_{\min}=10^{-100}$ and $10^{-280}$. Runs are on a Linux virtual machine with 10 vCPUs: $N_{\rm half}\le6049$, at most 2.5 min each. Raw data: `out/*.json`, `out/*.npz`.

## 2. Results

### (1) Eigenvalues of $L_\Delta$ on $(-1,1)$

Finest mesh: $n_0=400$, $\rho=1.12$, $N_{\rm half}=4424$. The error is the sum of the changes when $n_0$ is halved and when $\rho$ goes from 1.12 to 1.25. It is $\le6\cdot10^{-8}$ for $k\le5$.

| $k$ | state | $\mu_k$ (this work) | est. error | HJSS v3 | difference |
|---|---|---|---|---|---|
| 0 | e0 | $-0.68800259$ | $1\cdot10^{-9}$ | $-0.6878$ | $-2.0\cdot10^{-4}$ |
| 1 | o0 | $1.73292027$ | $3\cdot10^{-9}$ | $1.7331$ | $-1.8\cdot10^{-4}$ |
| 2 | e1 | $2.72741038$ | $6\cdot10^{-9}$ | $2.7275$ | $-0.9\cdot10^{-4}$ |
| 3 | o1 | $3.41209699$ | $1\cdot10^{-8}$ | $3.4122$ | $-1.0\cdot10^{-4}$ |
| 4 | e2 | $3.90890346$ | $3\cdot10^{-8}$ | $3.9091$ | $-2.0\cdot10^{-4}$ |
| 5 | o2 | $4.31398631$ | $5\cdot10^{-8}$ | $4.3142$ | $-2.1\cdot10^{-4}$ |

- The values agree with arXiv 2311.13079 v3, which gives 4 decimals, to within $2.2\cdot10^{-4}$.
- All six of our values lie below the published ones. They are also below the final referee's P0 values ($-0.68799,\dots,4.31449$, note v2 §13.1).
- Both facts are consistent with our values being finer upper bounds.

### (2) Scaling $\mu_k(L)=\mu_k(1)-2\log L$

The table shows $\mu_k(L)+2\log L-\mu_k(1)$. All runs use the **same absolute mesh parameters** ($h_0=0.01$, $\rho=1.25$, $r_{\min}=10^{-60}$), so the meshes are not dilates of each other.

| $k$ | $L=0.5$ | $L=2$ | $L=5$ |
|---|---|---|---|
| 0 | $+1.8\cdot10^{-8}$ | $-5.2\cdot10^{-9}$ | $-7.3\cdot10^{-9}$ |
| 1 | $+1.3\cdot10^{-7}$ | $-1.8\cdot10^{-8}$ | $-2.2\cdot10^{-8}$ |
| 2 | $+5.5\cdot10^{-7}$ | $-5.6\cdot10^{-8}$ | $-6.5\cdot10^{-8}$ |
| 3 | $+1.7\cdot10^{-6}$ | $-1.5\cdot10^{-7}$ | $-1.7\cdot10^{-7}$ |
| 4 | $+4.0\cdot10^{-6}$ | $-3.5\cdot10^{-7}$ | $-3.8\cdot10^{-7}$ |
| 5 | $+8.1\cdot10^{-6}$ | $-6.9\cdot10^{-7}$ | $-7.6\cdot10^{-7}$ |

- The residuals are discretisation error. They are largest at $L=0.5$, where $h_0=0.01=L/50$ is the coarsest relative mesh.
- With meshes that are (nearly) dilates ($h_0=L/200$, $\rho=1.12$), the residual is $\le1.2\cdot10^{-10}$.
- *What this tests:* the law is exact by dilation covariance, so this is a consistency check of the code, not of new mathematics.

### (3) Edge trace $A$, absolute determination, $2LA^2$ (prediction 1)

**(a) Profile fit (`fit_A.py`).** No eigenvalue derivative and no Hadamard formula is used. $A$ comes only from the $L^2$-normalised nodal values $\psi(L-r)$, on the window $\ell\in[30,\ell_{\max}-20]$, about 3700 nodes at $\rho=1.12$. Fit models:
- **F1**: $\psi^{-2}=\alpha\ell+\beta+c_1/\ell+c_2/\ell^2$, with $A=\alpha^{-1/2}$;
- **F2**: $\psi\ell^{1/2}=A+\sum_{k=1}^3b_k\ell^{-k}$;
- **F4**: F1 plus an $\ell^{1/2}$ term. This admits a $V^2$ (i.e. $\ell^{-1/2}$) correction to $\psi\ell^{1/2}$;
- the local slope $d\psi^{-2}/d\ell$ at $\ell=200$, which should equal $2L$.

The table gives F1 on the finest mesh ($\rho=1.12$; $n_0=400$ at $L=1$, $n_0=200$ at $L=0.5,2$). The error bar is the sum of three terms:
- the $\rho$-change $1.25\to1.12$;
- $\max(|{\rm F1}-{\rm F2}|,|{\rm F1}-{\rm F4}|)$;
- the $n_0$-halving change.

| state | $L=0.5$ | $L=1$ | $L=2$ |
|---|---|---|---|
| e0 | $1.0000000\pm1\cdot10^{-6}$ | $1.0000001\pm8\cdot10^{-7}$ | $1.0000001\pm9\cdot10^{-7}$ |
| o0 | $1.0000001\pm2\cdot10^{-6}$ | $1.0000000\pm1\cdot10^{-6}$ | $1.0000001\pm1\cdot10^{-6}$ |
| e1 | $1.0000001\pm2\cdot10^{-6}$ | $1.0000001\pm2\cdot10^{-6}$ | $1.0000001\pm2\cdot10^{-6}$ |
| o1 | $1.0000001\pm3\cdot10^{-6}$ | $1.0000001\pm2\cdot10^{-6}$ | $1.0000002\pm2\cdot10^{-6}$ |
| e2 | $1.0000002\pm4\cdot10^{-6}$ | $1.0000001\pm2\cdot10^{-6}$ | $1.0000002\pm3\cdot10^{-6}$ |
| o2 | $1.0000003\pm6\cdot10^{-6}$ | $1.0000001\pm3\cdot10^{-6}$ | $1.0000003\pm4\cdot10^{-6}$ |

Equivalently, $A=1,\ 0.7071068,\ 0.5$ at $L=0.5,1,2$, each to within $3.1\cdot10^{-6}$ relative (value plus error bar).

**Convergence in $\rho$** (F1, $L=1$, $n_0=400$; e0 / o2):

| $\rho$ | e0 | o2 |
|---|---|---|
| 2 | $1.0000520$ | $1.0000861$ |
| 1.5 | $1.0000072$ | $1.0000127$ |
| 1.25 | $1.0000008$ | $1.0000015$ |
| 1.12 | $1.00000005$ | $1.00000008$ |

**Other checks.**
- The local slope at $\ell=200$, divided by $2L$, is $1.0000000$–$1.0000003$ in all 18 cases.
- Changing the edge cut to $r_{\min}=10^{-100}$ or $10^{-280}$ changes $2LA^2$ by $\le2\cdot10^{-7}$.
- The fit F3 ($\psi\ell^{1/2}=A+b_0\ell^{-1/2}+\dots$, polynomial in $\ell^{-1/2}$ up to $\ell^{-2}$) shifts $A$ by up to $4\cdot10^{-5}$. That basis is badly conditioned on this data: applied to the exact-shape profile $(2L)^{-1/2}(\ell-\lambda-\gamma)^{-1/2}$ it produces most of the same shift. The well-conditioned F4, which admits the same $V^2$ term, moves $2LA^2$ by $\le1.2\cdot10^{-6}$. So no $\ell^{-1/2}$ term is resolved.

**(b) Explicit formula of Theorem 1(3) (`explicit_A.py`).**

*What is computed.* The pure-log split (J_VG, $\tfrac12\log(1+t^2)$, $w_{\log}$) is used, for $q=0.5,1,2$:
$$A=\int_0^\infty v_q f_q,\qquad f_q=\mathbf 1_{(0,2L)}\big[(q+\lambda)\tilde\psi+k*\tilde\psi\big]-\mathbf 1_{(2L,\infty)}\mu_{\rm ext}.$$
No profile fit enters.

*The kernel $v_q$.* It is computed from a Stieltjes representation. From [Kw11, Lemma 3.8, (3.4)], $\psi^\dagger(\xi)\psi^\dagger(-\xi)=\psi(-\xi^2)$; taking boundary values gives $K(-\tau+i0)=\big(q+\tfrac12\log(1-\tau^2+i0)\big)/K(\tau)$, and therefore
$$v_q(x)=a_*e^{-\tau_*x}+\tfrac12\int_1^\infty e^{-x\tau}\,\frac{K(\tau)}{(q+\frac12\log(\tau^2-1))^2+\pi^2/4}\,d\tau,$$
with $\tau_*=\sqrt{1-e^{-2q}}$ (the continuation of $1/K$ has a simple pole at $-\tau_*$, which gives the atom) and $a_*=K(\tau_*)e^{-2q}/\tau_*$.
- This is the $q>0$ analogue of Kwaśnicki–Małecki–Ryznar, *Suprema of Lévy processes*, Ann. Probab. 41 (2013), Prop. 4.5, eq. (4.6), which is stated there for $q=0$. The Stieltjes-inversion derivation above is sketched in Section 8 of the paper; it justifies the numerical method and is not used in the proofs. The representation is also checked numerically (`out/explicit_check.out`): $\int e^{-\xi x}v_q\,dx=1/K(\xi)$ at six sampled values spanning $[10^{-3},10^{40}]$ ($\xi=10^{-3},0.3,1,10,10^4,10^{40}$), with relative discrepancies below $1.2\cdot10^{-13}$ (no uniform check over the interval); $\int v_q=q^{-1/2}$ to $10^{-13}$; and the continuous part pointwise at $r=0.05$ and $1.5$ against an independent mpmath quadrature.
- As a by-product, $V_q(r)=(\ell(r)+q-\gamma)^{-1/2}(1+O(\ell^{-2}))$ numerically, in line with note v2 N8.

*The edge tail.* The FEM stops at $r_{\min}=10^{-200}$, where $V_q\approx0.047$. Below it,
$$\int_0^{r_{\min}}v_qf_q=f_q(0)V_q(r_{\min})+(q+\lambda)\int_0^{r_{\min}}v_q\tilde\psi .$$
- The last integral uses the tail model $\tilde\psi=A(V_q^{-2}+\delta)^{-1/2}$. This gives $A\,(\sqrt{1+\delta V_m^2}-1)/\delta$.
- $\delta$ is read off the FEM value at $\ell=200$, as $\delta=(A/\psi)^2-V_q^{-2}$. This is one data point, not a fit; $\ell=100$ or $300$ gives the same result.
- The equation is solved for $A$ self-consistently, so no fitted $A$ enters.
- We find $\delta=-(q+\lambda)-0.001$, consistent with the observation in §3.
- The first-order tail $(q+\lambda)AV_m^2/2$ alone misses $(q+\lambda)^2AV_m^4/8$. That is up to $10^{-5}$ in $2LA^2$, and it shows up as a $q$-dependence. Both values are in `results.json` (`ratio_q`, `ratio_q_lin`).

*The exterior term.* It is computed through Fubini, as $\int\tilde\psi(r)N(2L-r)\,dr$ with an $E_1$ kernel.

*Three defects found and fixed while validating this code, all by its own checks:*
- the τ-rule overlapped the analytic tail, so $V(r_{\min})$ was $5\cdot10^{-4}$ too large;
- a panel-width mismatch made $v_q$ wrong by $3\cdot10^{-5}$ (caught by the mpmath pointwise check);
- the first-order tail was not enough (caught by the $q$-spread).

The values below are from the fixed code.

| state | $L=0.5$ | $L=1$ | $L=2$ |
|---|---|---|---|
| e0 | $1.0000001\pm3\cdot10^{-7}$ | $1.0000001\pm2\cdot10^{-7}$ | $1.0000001\pm5\cdot10^{-8}$ |
| o0 | $1.0000003\pm6\cdot10^{-7}$ | $1.0000002\pm4\cdot10^{-7}$ | $1.0000002\pm3\cdot10^{-7}$ |
| e1 | $1.0000003\pm8\cdot10^{-7}$ | $1.0000002\pm6\cdot10^{-7}$ | $1.0000002\pm5\cdot10^{-7}$ |
| o1 | $1.0000004\pm1\cdot10^{-6}$ | $1.0000003\pm7\cdot10^{-7}$ | $1.0000003\pm7\cdot10^{-7}$ |
| e2 | $1.0000005\pm1\cdot10^{-6}$ | $1.0000003\pm8\cdot10^{-7}$ | $1.0000004\pm9\cdot10^{-7}$ |
| o2 | $1.0000006\pm1\cdot10^{-6}$ | $1.0000003\pm9\cdot10^{-7}$ | $1.0000005\pm1\cdot10^{-6}$ |

- Each entry is the mean over $q\in\{0.5,1,2\}$ on the finest mesh. The error is half the $q$-spread plus the change from $\rho=1.25$ to $1.12$.
- The largest $q$-spread is $1.4\cdot10^{-7}$.
- Convergence at $L=1$, $n_0=400$, all states and $q$: $2LA^2-1$ lies in $[1.8\cdot10^{-7},6.7\cdot10^{-6}]$ at $\rho=1.5$, in $[4.8\cdot10^{-8},1.4\cdot10^{-6}]$ at $\rho=1.25$, and in $[2.2\cdot10^{-8},4.0\cdot10^{-7}]$ at $\rho=1.12$. With $r_{\min}=10^{-280}$ it lies in $[1.6\cdot10^{-8},2.5\cdot10^{-7}]$.
- With only the first-order tail, $2LA^2-1\in[-2.7\cdot10^{-5},6.7\cdot10^{-8}]$, with a visible $q$-dependence. This is why the second-order tail is used.
- *What this tests.* This is a second absolute route to $A$. It is fit-free in $A$ and gives $2LA^2=1$ to $\le1.9\cdot10^{-6}$ in all 18 cases. It uses the same FEM eigenfunctions as (3a), but samples them globally, weighted by $v_q$, rather than only at the edge.
- *Relation to κ.* Theorem 1(3) rests on the half-line inversion of note v2 §5 and does not use the layer Lemmas 8.2–8.3. So (3b), like (3a), tests $\kappa=1$ against the exact $\lambda'=-1/L$.
- *Caveat.* (3b) is also a consistency check of the eigenfunctions and of Theorem 1(3). It is not an independent proof of anything.

### (4) Simplicity: eigenvalue gaps ($L=1$, $n_0=400$, $\rho=1.12$; 10 eigenvalues per sector)

- The 20 lowest eigenvalues alternate e0, o0, e1, o1, …, e9, o9 (interlacing).
- The smallest gap is $\mu_{19}-\mu_{18}=0.10550$, between $6.73859$ and $6.84409$. The smallest relative gap is $0.0154$.
- The smallest gap within one sector is $0.2164$.
- For these eigenvalues the discretisation error is $\le5.1\cdot10^{-6}$, so every computed gap exceeds it by a factor $>2\cdot10^4$.
- So the computed spectrum is simple, with a wide margin. This is a numerical observation consistent with Theorem 4; it is not a proof of it.

### (5) Mixed identity without fitting: $L(\lambda_e-\lambda_o)\Phi'(0)=1$

$\Phi'(0)=\langle\psi_o,\psi_e'\rangle$ is computed exactly for the P1 functions, as $-2\sum(\Delta\psi_e)\,\overline{\psi_o}$. The table shows the value minus 1 ($\rho=1.12$, $r_{\min}=10^{-200}$):

| pair | $L=0.5$ | $L=1$ | $L=2$ |
|---|---|---|---|
| e0o0 | $+8.7\cdot10^{-7}$ | $+8.6\cdot10^{-7}$ | $+8.6\cdot10^{-7}$ |
| e0o1 | $+2.5\cdot10^{-6}$ | $+2.5\cdot10^{-6}$ | $+2.5\cdot10^{-6}$ |
| e0o2 | $+3.8\cdot10^{-6}$ | $+3.7\cdot10^{-6}$ | $+3.8\cdot10^{-6}$ |
| e1o0 | $+1.5\cdot10^{-7}$ | $+1.5\cdot10^{-7}$ | $+1.5\cdot10^{-7}$ |
| e1o1 | $+8.4\cdot10^{-8}$ | $+7.1\cdot10^{-8}$ | $+8.4\cdot10^{-8}$ |
| e1o2 | $+4.5\cdot10^{-7}$ | $+3.8\cdot10^{-7}$ | $+4.4\cdot10^{-7}$ |
| e2o0 | $+7.5\cdot10^{-7}$ | $+7.1\cdot10^{-7}$ | $+7.4\cdot10^{-7}$ |
| e2o1 | $+6.9\cdot10^{-8}$ | $+4.0\cdot10^{-8}$ | $+6.9\cdot10^{-8}$ |
| e2o2 | $+8.9\cdot10^{-8}$ | $+3.0\cdot10^{-8}$ | $+8.9\cdot10^{-8}$ |

- The residual comes entirely from the edge cut. The FEM ramps $\psi$ to 0 on $[0,r_{\min}]$, which misses a piece of $\int\psi_o\psi_e'$ of order $(\lambda_o-\lambda_e)/\ell(r_{\min})^2$.
- The residual scales as $\ell(r_{\min})^{-2}$. At $L=1$, for e0o2: $+1.47\cdot10^{-5}$, $+3.70\cdot10^{-6}$ and $+1.90\cdot10^{-6}$ at $r_{\min}=10^{-100},10^{-200},10^{-280}$.
- A model of the cut, which uses the asymptotic profile $(2L)^{-1/2}(\ell-\lambda-\gamma)^{-1/2}$ (a model estimate, not part of the fit-free number), predicts $+1.49\cdot10^{-5}$, $+3.71\cdot10^{-6}$ and $+1.89\cdot10^{-6}$. After subtracting it, all 9 pairs agree with 1 to $\le1.1\cdot10^{-8}$ at $r_{\min}=10^{-200}$ and $10^{-280}$. At $10^{-100}$ the agreement is $\le2\cdot10^{-7}$, because the model's own higher-order error is larger there (`results.json`, key `mixed_rmin`).
- *What this tests* (note v2 §13.1): the mixed identity and Theorem 4 jointly. It is **blind to the layer constant κ**, since $\kappa$ cancels. Only (3) tests $\kappa=1$.

## 3. What these numbers do and do not show

- **The constant.** (3a) and (3b) are absolute determinations of $A$, made without any Hadamard formula, compared with the exact $\lambda'=-1/L$. They give $2LA^2=1$:
  - to $\le7\cdot10^{-6}$, profile fit, 18 cases;
  - to $\le2\cdot10^{-6}$, explicit formula, 18 cases, 3 values of $q$.

  This confirms $\kappa=1$ more than two orders of magnitude more tightly than the referee's P0 check, which reached the $10^{-3}$ level.
- **Independence.** This code is independent of the hp-FEM and of the referee's P0 code: different representation of the form, different elements (P1), different solver.
- **Shared assumptions.** It shares the definition of the operator and the meaning of $A$ with them.
- **Status.** All results are uncertified floating-point computations. The error bars are refinement estimates, not bounds. The eigenvalues are Galerkin upper bounds up to rounding.
- **Observation, not proved.** On the finest meshes, $\psi(L-r)^{-2}=2L\,(\ell(r)+c)+O(\ell^{-1})$ with $c+\lambda+\gamma=-0.0042,\,-0.0011$ at $\rho=1.25,1.12$ (it was $-0.0137$ at $\rho=1.5$). That is, $c\to-\lambda-\gamma$ under refinement, for all states and all $L$. This matches the heuristic $V_q\approx(\ell+q-\gamma)^{-1/2}$ with $q=-\lambda$.

## 4. Files

- `fem.py`: assembly and solver (`--selftest`).
- `fit_A.py`: profile fits.
- `explicit_A.py`: Theorem 1(3) (`check` mode verifies the $v_q$ representation).
- `analyse.py`: writes `results.json` and `tables.tex`.
- `run_grid.sh` with `cases.txt`, and `out/run_grid2.sh` with `out/cases2.txt`: the drivers of the recorded runs.
- `out/`: raw `.json`, `.npz` and `.out` for every run, plus `*_explicit.json`/`.out` and `explicit_check.out`.
- **Reproduce:** `python fem.py --L 1 --n0 400 --rho 1.12 --rmin 1e-200 --nvec 3 --out out/L1_n400_r1.12` (about 1 min on 5 cores), then `python analyse.py`.

## 5. LaTeX tables (`tables.tex`, generated)

```latex
% Table 1: eigenvalues of L_Delta on (-1,1)
\begin{tabular}{rlllr}\toprule $k$ & sector & $\mu_k$ (this work) & HJSS & $\mu_k-\mu_k^{\rm HJSS}$\\\midrule
0 & e0 & $-0.68800259$ & $-0.6878$ & $-2.0e-04$\\
1 & o0 & $1.73292027$ & $1.7331$ & $-1.8e-04$\\
2 & e1 & $2.72741038$ & $2.7275$ & $-9.0e-05$\\
3 & o1 & $3.41209699$ & $3.4122$ & $-1.0e-04$\\
4 & e2 & $3.90890346$ & $3.9091$ & $-2.0e-04$\\
5 & o2 & $4.31398631$ & $4.3142$ & $-2.1e-04$\\
\bottomrule\end{tabular}
% Table 2: scaling check  mu_k(L) + 2 log L - mu_k(1)  (same absolute mesh: h0=0.01, rho=1.25, r_min=1e-60)
\begin{tabular}{rrrr}\toprule $k$ & $L=0.5$ & $L=2$ & $L=5$\\\midrule
0 & $+1.8e-08$ & $-5.2e-09$ & $-7.3e-09$\\
1 & $+1.3e-07$ & $-1.8e-08$ & $-2.2e-08$\\
2 & $+5.5e-07$ & $-5.6e-08$ & $-6.5e-08$\\
3 & $+1.7e-06$ & $-1.5e-07$ & $-1.7e-07$\\
4 & $+4.0e-06$ & $-3.5e-07$ & $-3.8e-07$\\
5 & $+8.1e-06$ & $-6.9e-07$ & $-7.6e-07$\\
\bottomrule\end{tabular}
% Table 3: 2 L A^2 from the edge profile (fit F1, finest mesh) with mesh/fit error bar
\begin{tabular}{lrrr}\toprule state & $L=0.5$ & $L=1$ & $L=2$\\\midrule
e0 & $1.0000000\pm1e-06$ & $1.0000001\pm8e-07$ & $1.0000001\pm9e-07$\\
o0 & $1.0000001\pm2e-06$ & $1.0000000\pm1e-06$ & $1.0000001\pm1e-06$\\
e1 & $1.0000001\pm2e-06$ & $1.0000001\pm2e-06$ & $1.0000001\pm2e-06$\\
o1 & $1.0000001\pm3e-06$ & $1.0000001\pm2e-06$ & $1.0000002\pm2e-06$\\
e2 & $1.0000002\pm4e-06$ & $1.0000001\pm2e-06$ & $1.0000002\pm3e-06$\\
o2 & $1.0000003\pm6e-06$ & $1.0000001\pm3e-06$ & $1.0000003\pm4e-06$\\
\bottomrule\end{tabular}
% Table 4: 2 L A^2 from the explicit formula of Thm 1(3); mean over q = 0.5, 1, 2 on the finest mesh;
% error = half the q-spread + change from rho = 1.25 to 1.12
\begin{tabular}{lrrr}\toprule state & $L=0.5$ & $L=1$ & $L=2$\\\midrule
e0 & $1.0000001\pm3e-07$ & $1.0000001\pm2e-07$ & $1.0000001\pm5e-08$\\
o0 & $1.0000003\pm6e-07$ & $1.0000002\pm4e-07$ & $1.0000002\pm3e-07$\\
e1 & $1.0000003\pm8e-07$ & $1.0000002\pm6e-07$ & $1.0000002\pm5e-07$\\
o1 & $1.0000004\pm1e-06$ & $1.0000003\pm7e-07$ & $1.0000003\pm7e-07$\\
e2 & $1.0000005\pm1e-06$ & $1.0000003\pm8e-07$ & $1.0000004\pm9e-07$\\
o2 & $1.0000006\pm1e-06$ & $1.0000003\pm9e-07$ & $1.0000005\pm1e-06$\\
\bottomrule\end{tabular}
% Table 5: mixed identity  L (lambda_e - lambda_o) Phi^\prime(0) - 1  (lambda = mu/2), no fitting
\begin{tabular}{lrrr}\toprule pair & $L=0.5$ & $L=1$ & $L=2$\\\midrule
e0o0 & $+8.7e-07$ & $+8.6e-07$ & $+8.6e-07$\\
e0o1 & $+2.5e-06$ & $+2.5e-06$ & $+2.5e-06$\\
e0o2 & $+3.8e-06$ & $+3.7e-06$ & $+3.8e-06$\\
e1o0 & $+1.5e-07$ & $+1.5e-07$ & $+1.5e-07$\\
e1o1 & $+8.4e-08$ & $+7.1e-08$ & $+8.4e-08$\\
e1o2 & $+4.5e-07$ & $+3.8e-07$ & $+4.4e-07$\\
e2o0 & $+7.5e-07$ & $+7.1e-07$ & $+7.4e-07$\\
e2o1 & $+6.9e-08$ & $+4.0e-08$ & $+6.9e-08$\\
e2o2 & $+8.9e-08$ & $+3.0e-08$ & $+8.9e-08$\\
\bottomrule\end{tabular}
```

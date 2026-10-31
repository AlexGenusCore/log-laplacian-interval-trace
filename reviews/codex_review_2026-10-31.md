I recommend **public release as a preprint after minor corrections**. I read all 745 lines of [main.tex](../paper/main.tex), including every proof and the appendix. I found no fatal or major mathematical gap: the Wiener–Hopf inversion, boundary estimate, dilation identity, trace normalisation and simplicity argument withstand this review. The numerical section generally distinguishes illustration from proof and discloses its limitations well. The principal attribution correction is that the withdrawn Pohozaev formula formally implies the separate endpoint identities, not merely their sum. For journal consideration, the contribution should remain centred on the different proof, quantitative estimate and additional identities, since the principal conclusions already appear in acknowledged public work. I found no claim of progress on RH supported or suggested by the results themselves. All checks were read-only.

1. **Minor — priority accuracy: Remark 6.4, line 544; related wording in the abstract and §§1.1–1.5.**  
   The assertion that the separate endpoint identities “are not contained” in AGV24 is incorrect. Its announced formula, with \(f(x,u)=\mu u\), gives on \((0,d)\)
   \[
   2\|u\|_2^2=2d\,A_+(u)^2.
   \]
   Reflection gives the other endpoint. Equivalently, translating the origin in the formula on \((-L,L)\) recovers both identities. This is a consequence of the *announced statement*, without accepting its withdrawn proof. [AGV24, Theorem 2.3 and Corollary 3.8](https://arxiv.org/html/2411.15985v1).

   **Fix:** Say that the withdrawn formula formally contains the full eigenfunction trace identity after translation, while its proof and the required trace existence were not established. The exact finite-dilation identity (6.1) remains a distinct contribution.

2. **Minor — numerical validation: §8, line 627; `explicit_A.py`, function `check_rep`.**  
   The claim that the Laplace transform matches \(1/K_q\) “to \(10^{-13}\) on \([10^{-3},10^{40}]\)” overstates the recorded check. The code tests six points:
   \[
   10^{-3},\ 0.3,\ 1,\ 10,\ 10^4,\ 10^{40}.
   \]
   The saved relative discrepancies reach approximately \(1.18\times10^{-13}\).

   **Fix:** Write “at six sampled values spanning this range, with relative discrepancies below \(1.2\times10^{-13}\).” No uniform validation over the interval was performed.

3. **Minor — justification of the numerical kernel: §8, line 627.**  
   The paper correctly discloses that KMR13 Proposition 4.5 concerns the unkilled process and does not literally justify the \(q>0\) implementation. The cited proposition and equation number (4.6) are correct. [KMR13, Proposition 4.5](https://arxiv.org/pdf/1103.0935).

   **Fix:** Include the short Stieltjes-inversion derivation of the formula actually used:
   \[
   v_q(x)=a_*e^{-s_*x}
   +\frac12\int_1^\infty
   \frac{e^{-sx}K_q(s)}
   {(q+\frac12\log(s^2-1))^2+\pi^2/4}\,ds,
   \]
   where
   \[
   s_*=\sqrt{1-e^{-2q}},\qquad
   a_*=\frac{e^{-2q}K_q(s_*)}{s_*}.
   \]
   The atom comes from the pole at \(-s_*\); the continuous density comes from the boundary value across \((-\infty,-1)\). This follows from the factorisation identity in Kw11 Lemma 3.8. I checked the residue and density against the implementation. This omission affects the presentation of the numerical method, not the main proofs.

4. **Minor — domain of an expression: §6 opening, line 464.**  
   The definition
   \[
   (J_0*\psi)(y)=\int\frac{\psi(x)}{2|x-y|}\,dx
   \]
   is written without restricting \(y\). It generally diverges inside the support, and can diverge at an endpoint even for the eigenfunctions considered here.

   **Fix:** Define this expression for \(y\notin\overline I\). State that endpoint values are irrelevant to the exterior integrals. All subsequent uses are valid exterior uses.

5. **Minor — explicit constant bookkeeping: proof of Theorem 1.1, line 454.**  
   Theorem 5.1 supplies a bound \(C_E\ell^{-1}\) at each endpoint, whereas Theorem 1.1 bounds the **sum** of the two errors.

   **Fix:** Explicitly choose \(C=2C_E\|u\|_2\), after making the stated parameter-only replacements. The existential theorem is unaffected; this makes its advertised constant unambiguous.

6. **Minor — exposition: abstract, §1.5 and repeated priority discussions.**  
   The mathematical development is readable, but the abstract and introduction are overloaded with repeated descriptions of research notes, review history and qualifications. This makes the actual contribution harder for an outside specialist to identify.

   **Fix:** Retain a concise, prominent acknowledgement of BT26 and Mathia26, then state precisely what differs here. Move the detailed note-by-note comparison to a separate appendix or supplementary discussion. Preserve the qualifications about unverified proofs and the withdrawn formula.

7. **Typo — §2.5 proof, line 256.**  
   Simultaneously translating \(f\) and \(g\) leaves \(\hat f\,\overline{\hat g}\) unchanged; it does not multiply that product by a nontrivial unimodular factor.

   **Fix:** Say that the two Fourier phases cancel. Also specify \(z>0\) for the inequality \(\Ein(z)\le z\) in line 165.

**Claims independently verified, and how**

- **Proposition 2.1 and Lemmas 2.2–2.5:** I rederived the Fourier/position-space identity, the variance-gamma decomposition, truncation, density and dilation formulas. The constants
  \[
  c_1=1,\qquad \rho_1=-2\gamma,\qquad \mathcal E_{\log}=2\mathcal Q
  \]
  are consistent. The disjoint-support interaction kernel is exactly \(1/(2|x-y|)\) for \(\mathcal Q\).

- **Lemma 2.6 and Appendix A:** I checked the contraction inequalities, resolvent-to-semigroup domination, Duhamel identity and \(L^p\) bootstrap. The chosen \(\tau\) ensures \(\|T_\tau\|_{r\to r}\le1/2\), including \(r=\infty\). The Gamma-function bound and number of bootstrap steps yield the stated \(B_\infty\).

- **Lemmas 3.1–3.3:** I checked the complete Bernstein/Stieltjes representations, boundary factorisation, total mass of \(v_q\), asymptotic estimates and renewal-function inequalities. In particular,
  \[
  0\le\frac1{1+s}-e^{-s}\le\frac{2s}{(1+s)^2}
  \]
  gives the asserted renewal estimate. Independent NumPy quadrature at 30 parameter pairs also satisfied the bounds for \(\log K_q\) and \(\beta\). The relevant Kw11 numbering and hypotheses match [arXiv:1006.0524v4](https://arxiv.org/html/1006.0524v4).

- **Proposition 4.1:** I checked every step of the inversion, including Fourier signs and the Hardy-space support argument. The identity
  \[
  \sigma_q|\widehat{Gf}|^2
  =|\widehat{P_+h}|^2
  \]
  supplies the claimed form-domain membership. The defect pairing vanishes with the stated supports; density then extends the equation to all \(D_+\).

- **Lemma 4.2 and Theorem 5.1:** The convolution continuity argument, tail estimates and \(q\)-independence of the trace are sound. I independently checked
  \[
  \int_0^\infty E_1(x)^2\,dx=2\log2,
  \qquad
  \sup_{\ell\ge1}\ell e^{-\ell/2}=2/e,
  \]
  which produce the displayed constants. Theorem 1.1 follows, with issue 5’s bookkeeping clarification.

- **Lemmas 6.1–6.3 and Theorem 1.2:** The truncated functions are admissible tests, and the cross-boundary integrals converge absolutely. The exterior-field coefficient is exactly \(1\): the kernel’s \(1/2\) cancels the factor \(2\) from integrating \(1/(r\sqrt{\ell(r)})\). Dividing (6.1) by \(2L(\theta-1)\) therefore gives
  \[
  A_\pm^2=\frac{\|\psi\|_2^2}{2L}.
  \]
  There is no missing factor of two. The layer estimate is sufficient for the limits actually taken.

- **Theorem 1.4, Proposition 7.1 and Corollary 1.5:** Simplicity follows from injectivity of the linear scalar trace on an eigenspace. Parity is subsequently deduced, so there is no circularity. I checked the translation signs and the conversion \(\mu=2\lambda\); the mixed identity has the stated factor \(4\) and magnitude \(2/L\). Its formulation through \(\Phi'(0)\) avoids assuming unproved derivative integrability.

- **Corollary 9.3 and the background comparisons:** The parity limit argument uses the appropriate fractional eigenfunction convergence. FJW22 Theorem 1.1(v), Corollary 1.4, and Zhang’s Lemma 6.1 have the stated content. The discussion correctly distinguishes weak bracketing from simplicity and correctly excludes direct application of complete Bernstein results to \(\log\xi\). [FJW22](https://arxiv.org/html/2010.10448), [Zhang](https://arxiv.org/html/2608.23457v2), [Kwaśnicki–Wszoła](https://arxiv.org/html/2609.24892v1).

- **Numerical tables and diagnostics:** I reproduced all 18 profile-fit central values and their stated refinement-estimate calculations from the stored nodal data. I independently integrated the P1 functions to confirm normalisation and recomputed all 27 mixed-identity residuals, obtaining
  \[
  3.0172\times10^{-8}\le\text{residual}\le3.8069\times10^{-6}.
  \]
  Rebuilding the five deepest rows in all 30 stored runs gave a maximum residual \(2.1402\times10^{-13}\), consistent with the text. The assembly self-test passed.

- **Eigenvalues:** The six entries in Table 1 match the raw data and lie inside the stated brackets. The maximum reported refinement estimate is \(5.315\times10^{-8}\), and the minimum gap among the first 20 values is \(0.1054953\). A separately assembled uniform P0 discretisation with 160, 320 and 640 cells converged towards the reported eigenvalues; its first value was respectively
  \[
  -0.68759663,\quad -0.68785350,\quad -0.68794665.
  \]
  This is a coarse independent check, not confirmation of eight-digit accuracy. The cited HJSS25 comparison values are transcribed correctly. [HJSS25](https://arxiv.org/html/2311.13079).

- **Numerical honesty and supplement integrity:** All 128 supplied SHA-256 checksums passed. The manuscript discloses the stronger fitting assumptions, excluded tail-model errors, one-sided discrepancies and unsuccessful P0/F3 behaviour. The two trace-estimation routes share FEM data and modelling assumptions; they should not be counted as independent validations at the quoted precision.

- **Priority:** The accessible rStation manuscript dated 27 October 2026 explicitly contains trace existence, the endpoint squares and interval simplicity in §3.3 and Corollary 2.3. Those acknowledgements are accurate. The cached Mathia notes contain the earlier trace/remainder claims described here. I found no additional uncited overlapping result in the searches completed, but this is not an exhaustive novelty finding. [BT26](https://rstation.ai/problems/riemann_hypothesis/papers/67). AGV24’s withdrawal date and reason are correctly reported. [Withdrawal notice](https://arxiv.org/abs/2411.15985v2).

**What I could not check**

- Direct requests to the specified rStation APIs failed, and pagination was partly inaccessible. I inspected the accessible manuscript and listings, plus identified cached API responses; I could not complete a fresh platform-wide API search.
- Zenodo access and GitHub content/search access were limited. Some comparisons therefore used cached primary PDFs and repository files. I did not independently authenticate every historical commit timestamp or the account of independent discovery.
- SciPy and mpmath were unavailable. I did not rerun the full finest-mesh eigensolver or the explicit-formula quadrature. My numerical checks used NumPy, stored eigenvectors, rebuilt residual rows and independent P0 assembly.
- I did not verify every peripheral citation or original book theorem number against its printed edition. The principal accessible theorem references used in the proof chain were checked.
- I did not certify the proofs of BT26 or Mathia26, rigorous numerical error bounds, Conjecture 9.2, derivative asymptotics, or extensions to higher dimensions and the full Weil operator. Those are not established by this review.
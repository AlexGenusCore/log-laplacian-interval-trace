# Response to the external review of 31 October 2026

Paper: `papers/purelog/main.tex` / `main.pdf` (now 24 pp.). The versions before these changes are
`review_ext1/main_before_ext1.tex` and `review_ext1/refs_before_ext1.bib`. Statement numbers below are those of the
compiled PDF. They are unchanged from the reviewed version, except that Conjecture 9.3 is now Corollary 9.3. (The
review's "Thm 1.3" is our Thm 1.4 and its "Cor 1.4" is our Cor 1.5.) The PDF compiles with tectonic with no errors,
no undefined references or citations, and no overfull or underfull boxes.

We thank the reviewer for the careful reading, for the checks of the normalisation, the Hardy-space step and the
simplicity deduction, and above all for pointing us to the overlapping manuscript.

## Main finding: the earlier manuscript by bramble_trace

**Found.** The manuscript is "Boundary traces and support variation for localized Weil operators" by the agent
handle bramble_trace: rStation paper riemann_hypothesis_67, version 1, submitted 2026-10-27 19:01:48 UTC and
published 19:07:56 UTC, https://rstation.ai/problems/riemann_hypothesis/papers/67/v1. rStation is a platform for AI
research agents, and the paper was accepted there by the platform's AI reviewer agents. Its Corollary 2.3 is as the
reviewer describes. With $a$ the half-length, $|\tau_\pm|^2=1/(2a)$ is exactly our $A_\pm^2=\|u\|^2/|I|$, because the
trace weight $\ell^{1/2}$ is the same; there is no factor 2 to reconcile. The search also found earlier public
research notes, S. Garcia Murillo's "mathia" repository (notes PL-319 to PL-323, 13–15 September 2026). These derive
the existence of the trace, with a remainder $o(\ell^{-1})$, for eigenstates of the localized Weil operator. Full
statement-by-statement comparison: `review_ext1/comparison.md`; search record: `hunt_log.md`, `bramble_summary.md`,
`overlap_sweep.md`.

**Changes.**
- **Abstract** (p. 1): "with an explicit rate" became "with a remainder of order $\log(1/\dist(x,\partial I))^{-1}$".
  New sentences state that existence of the trace, the trace identity and simplicity were obtained earlier in the
  rStation manuscript (the last two from a Hadamard-type support-variation formula), and that existence of the trace
  had also been derived in public research notes. They also say that we give a different, Pohozaev-type proof, and
  what is new to our knowledge: the quantitative remainder, the one-endpoint dilation identity, and the even–odd
  identity.
- **§1.1** (pp. 1–2): "three questions were open" is now restricted to the journal and arXiv literature, with a
  pointer to the recent unrefereed answers. In item (1), "the existence … was not established" is replaced by
  "derived in [Mathia26, PL-319–PL-323] and proved in [BT26, §3.3]". Item (2) adds that the 1-D eigenfunction case
  of the AGV identity follows from [BT26, Cor. 2.3]. Item (3), "Only $\mu_1$ was known to be simple", is now
  qualified and cites [BT26, Cor. 2.3].
- **§1.2** (pp. 2–3): new paragraph describing [BT26] (§3.3, Thms 2.1–2.2, Cor. 2.3) and [Mathia26] (PL-297,
  PL-306, PL-319–PL-323), plus [Ve26]. Also new citations of Zhang [Zh26, Lemma 6.1] and Chandru [Ch26].
- **§1.4** (Method (1)): notes that [BT26, (3.19)] is a limit formula of the same kind as ours, without a rate.
- **Rem 1.3**: the boundary-trace form of the Hadamard formula, for clusters, one-sided variations and the Weil
  operator, is credited to [BT26, Thms 2.1–2.2], and Mathia's conditional zero-order Hadamard law to PL-306.
- **§1.5**, retitled "Relation to earlier work, novelty and status" (pp. 4–5), is rewritten. It now:
  - describes the overlap first;
  - removes "we know of no proof" of the summed AGV identity;
  - presents the paper as an independent Pohozaev-type proof with a quantitative edge law, not a first proof;
  - limits "to our knowledge new" to the identity (6.1), Prop. 7.1 for different eigenvalues, and Cor. 1.5;
  - states that for equal eigenvalues Prop. 7.1 reduces to [BT26, (6.3)];
  - says the rate is quantitative, not "stronger": Mathia claims $o(\ell^{-1})$ without constants;
  - says we knew parts of [Mathia26] while working on the paper, so no independence is claimed from it.
- **Rem 5.2(4)**, the remark on [Kw23] in §7, **Rem 7.3** and **Rem 9.4(2)**: cite [BT26 §3.3], [Mathia26 PL-297,
  PL-319–PL-323] and [Ch26] where they bear on the statement.
- **Conjecture 9.3 is now Corollary 9.3** (parity alternation), with a three-line proof:
  - [Zh26, Lemma 6.1] gives the parity $(-1)^{k+1}$ of the fractional eigenfunctions;
  - [FJW22, Thm 1.1(v)] gives $L^2$ convergence as $s\to0$;
  - our Thm 1.4 gives simplicity.

  Both cited statements were checked against the papers. The same deduction, by a different route, was posted
  earlier as a forum comment by bramble_trace [BT26c], which is credited.
- **§8**: the doubled Galerkin values of [BT26, §10] agree with Table 1 to within $4\cdot10^{-8}$, and those of
  [Ve26] to within $3.4\cdot10^{-6}$.
- **refs.bib**: new entries `BT26`, `BT26c`, `Mathia26`, `Zh26`, `Ch26` and `Ve26`, built only from checked
  metadata; `KW26` now names its version (v1, 21 September 2026).

**Open points for the author before any public version** (marked as `% AUTHOR TO CONFIRM` / `% AUTHOR TO DECIDE`
comments in `main.tex`):
1. **Independence.** §1.5 states that our proofs were found independently of [BT26]. The project record supports
   this: the proofs date from 30 October, and no project file mentions rStation before 31 October. Only the author
   can confirm it. If it is not confirmed, delete that sentence. The abstract makes no independence claim.
2. **Reviewer identity.** The Verification section uses a wording that does not depend on whether the reviewer was a
   human or an AI. If a human expert, and if they agree, the acknowledgement can be strengthened (the alternative
   sentence is in the comment).
3. **BT26 v2.** Re-check [BT26] for a v2 before submission; its statement numbers may change. Open the rStation
   forum URL in `BT26c` once in a browser, since it was built from the API and never opened.
4. **Title.** It was kept. A retitle, e.g. "A dilation proof of the boundary trace identity …", is optional.

## Requested amendments

**1. Lemma 6.3: uniformity in $\theta$.** *Changed* (pp. 13–14).
- The error is now written $\varepsilon(\delta,\theta)$.
- The lemma states that $\varepsilon(\delta,\theta)\to0$ as $\delta\to0$ for each fixed $\theta$, uniformly for
  bounded $\theta$, and whenever $\log\theta=o(\log(1/\delta))$. Because of the term $|AB|\log\theta/\Lambda$, the
  bound does not tend to 0 uniformly in unrestricted $\theta\ge1$.
- A sentence after the proof records the only two uses, and both have bounded $\theta$:
  - the proof of Thm 1.2, with $\theta=1+\delta/(2L)\le1+\rho_0/(2L)$, which the proof now says explicitly;
  - the proof of Prop. 7.1, with $\theta=1$, where the error is now written $h\,\varepsilon(h,1)$.

**2. Galerkin upper bounds.** *Changed* (§8 Method, p. 16; Table 1 caption and text, p. 17).
- The text now says that the exact Galerkin eigenvalues of the conforming P1 space are upper bounds, sector by
  sector, by min-max.
- It also says that the computed values only approximate them: the far-field element integrals use quadrature, and
  neither quadrature nor rounding error is bounded rigorously, so the computed values are not certified upper
  bounds.
- "consistent with being finer Galerkin upper bounds" became "consistent with approximating Galerkin upper bounds
  on finer meshes".
- No other place in the paper calls the computed values upper bounds. The working notes `numerics.md` do, and are
  kept unchanged as a record; the supplement README points out the correction.

**3. Table 2 error estimates.** *Changed* (§8 introduction, p. 16; trace-constant text and Table 2 caption, p. 18).
- The text and the caption now state that the error estimates are refinement estimates, and that they exclude both
  the edge cut $r_{\min}$ and, for the explicit formula, the error of the tail model.
- They also state that Table 2 is supporting evidence, not a rigorous error bound.

**4. Hopf lemma plus scaling.** *Changed* (§1.1 item (1), p. 2; the paragraph after the proof of Thm 1.4, p. 15).
- **Deduction verified.** Take $0<\theta\le e^{\mu_1(I)/2}$. Then $u_1(\cdot/\theta)$ is a positive first
  eigenfunction on $\theta I$ with eigenvalue $\mu_1(I)-2\log\theta\ge0$ (Lemma 2.5), so it is a nonnegative
  supersolution and [HLS25, Thm 1.4] applies. Since $\ell(\theta d)=\ell(d)+\log(1/\theta)$, we have
  $\ell(\theta d)/\ell(d)\to1$, which gives $\liminf u_1\ell^{1/2}>0$ on every interval.
- **The old restriction is removed.** The earlier text limited the citation to $|I|\le1.418$. It now gives the
  rescaling argument in one line, credited as an observation from HLS and Lemma 2.5 that the reviewer pointed out.
- **New content, stated plainly.** Relative to these bounds, what is new is the existence of the exact limit, with
  rate and value, for every eigenfunction, including the sign-changing ones, to which the Hopf lemma does not apply.
  The text adds that existence and value are also in [BT26], without a rate.

**5. Numerical supplement.** *Done.* It is in `papers/purelog/supplement/` and `papers/purelog/supplement.zip` (3.5
MB; 128 files with SHA256SUMS). Nothing was uploaded.
- **Code:** `fem.py`, `fit_A.py`, `explicit_A.py`, `analyse.py`, `plot_edge.py`, the drivers and the case lists.
- **Outputs and notes:** all raw outputs in `out/` (`.json`, `.out` and `.npz`, 4.1 MB in total, so the `.npz` files
  are included), plus `results.json`, `tables.tex` and `numerics.md`.
- **Checks:** the internal referee's P0 code with a fresh re-run of its six configurations, which reproduces the
  referee's table.
- **New residual check:** `residual_check.py` rebuilds the matrix rows of the five deepest nodes independently,
  tested against the full assembly, and evaluates them on the stored vectors.
- **README:** covers how to run, the environment, run times, what each file is, the residual definition, a table
  of deepest-node residuals for all 30 runs, and how Tables 1–2, Figure 1 and the Conjecture 9.2 evidence are
  produced.
- **Residual results:** all-row residuals are at most $2.92\cdot10^{-13}$. At the five deepest nodes the residuals
  are at most $2.14\cdot10^{-13}$, and at $r=r_{\min}$ at most $2.11\cdot10^{-13}$.
- **In the paper:** §8 now points to the supplement for these deepest-node residuals, and **Data availability**
  describes the supplement (p. 21).
- **Limits we state plainly:**
  - the package versions of the VM environment used for the recorded runs were not recorded;
  - SciPy and matplotlib are not installed on the laptop, so the FEM, explicit-formula and plotting scripts were
    not re-run for the supplement; only the P0 check and the residual check were;
  - the change between the last two inverse-iteration steps was not stored.

## Verification statement

Rewritten (p. 21).
- It now describes the external review: what it covered and what it confirmed, that it is not a certification, and
  that it did not reproduce the numerics.
- It lists the changes the review requested, all of which have been made.
- It states that no one has certified the proofs and that the code has not been checked outside the AI-assisted
  process. It names the Hardy-space step and Appendix A as the parts most in need of expert reading.
- It says that [BT26] theorem numbers refer to v1.

## Final pass after the two check reports (31 October 2026)

Sources: `check_novelty.md` (verdict MAJOR) and `check_amendments.md` (verdict minor). Backups made before editing:
`main_before_final2.tex`, `refs_before_final2.bib`, `supplement_before_final2.zip`, `numerics_md_before_final2.md`,
`responses_ext1_before_final2.md`. Page numbers below refer to the recompiled `main.pdf` (24 pages; tectonic run
with no errors, no undefined references or citations, no overfull or underfull boxes, clean BibTeX log). The
`% AUTHOR TO CONFIRM` and `% AUTHOR TO DECIDE` comment blocks and the author block are unchanged.

**Correction to the response above.** The answer above says that Mathia gives o(ℓ⁻¹) "without constants". That
is wrong. PL-329, PL-336 and PL-337 state a uniform O(ℓ⁻¹) remainder, as described in item 1.

1. **The O(ℓ⁻¹) rate is no longer claimed as new (MAJOR).** We read the notes in the local clone
   `sources/mathia_fresh/repo`.
   - **PL-329** (ac5835e30, 16 Sep 2026). It is a conditional statement about endpoint profiles g that satisfy the
     exact endpoint identity of PL-319, with the envelope |g| ≤ AQ^{-1/2} and a source bound S. The conclusion
     is |g − CQ^{-1/2}| ≤ K(1+S)/Q + Ke^{−η(Q−R)}, with K, η and R independent of the family parameter. The
     quantitative maximum-principle step (PL329-20) is asserted ("The same argument is quantitative"), not
     written out.
   - **PL-336** (2a3652475, 17 Sep). It applies PL-329 to the odd eigenbranch of Suzuki's
     A_a = ½L_Δ + γI + K_a at the first prime activation, conditionally on a uniform L^∞ bound. The results are
     (PL336-4) and u₀(1−2e^{−Q}) = C₀Q^{−1/2} + O(Q^{−1}) (PL336-6).
   - **PL-337** (8b910524f and ce507a6a8, 17 Sep). It proves ‖u‖_∞ ≤ C(J,M) for every normalised eigenpair with
     a in a compact set J and |λ| ≤ M.
   - **None of these is refereed.** Q = ℓ + log 2, so O(Q⁻¹) is O(ℓ⁻¹).

   Changes made:
   - **Abstract (p. 1).** "The quantitative remainder bound" is removed from the "new" list. The abstract now says
     that the research notes "also state a remainder of the same order O(ℓ⁻¹)", and that our existence proof
     follows essentially the same Wiener–Hopf route as BT26 and gives a self-contained proof of the remainder
     bound for the pure logarithmic operator.
   - **§1.5 (p. 5).** The notes are described as above: hypotheses, the asserted step, conditional versus proved,
     and the operator family. The sentence "Neither [BT26] nor [Mathia26] gives a quantitative bound of this kind;
     … qualitative … without a constant" is removed. It is replaced by: our proof is self-contained, the constant
     is explicit, BT26 gives no rate, and "a remainder of this order is stated in [PL-329, PL-336, PL-337] for the
     localized Weil operator". PL-323's o(ℓ⁻¹) is kept as a description of PL-323.
   - **Further citations:**
     - §1.1(1) (p. 2);
     - §1.2 (p. 3);
     - Rem. 5.2(4) (p. 13), for the bounded-forcing extension;
     - the paragraph after the proof of Thm 1.4 (p. 15);
     - Rem. 9.4(2) (p. 20).
   - **§8, model F1 (p. 17), optional item 7.** A parenthesis notes that PL-323's o(ℓ⁻¹), if correct and
     applicable, would exclude an ℓ^{1/2} term, consistently with Conj. 9.2.
   - **`refs.bib`, Mathia26.** The title now lists PL-329, PL-336 and PL-337, and the note gives "committed 13–17
     September 2026" with the three commit hashes.
2. **Abstract rate wording (p. 1).** It now reads: with ℓ(x) = log(1/dist(x,∂I)), the limit of uℓ^{1/2} exists, and
   u(x) = A_±(u)ℓ(x)^{−1/2} + O(ℓ(x)^{−1}), with an explicit constant (C_E of eq. (5.1), Thm 5.1(iii)).
3. **PS25 lemma numbers.** Checked against the text of arXiv:2407.11718v1 (`sources/checknov_ps25_2407.11718v1`).
   - **Lemma 3.11** is the Hopf lemma for L_Δu ≥ Vu, u ≥ 0, with ‖V‖ ≤ λ₁(B) on a ball B ⊂ Ω touching ∂Ω. It
     gives liminf u₁ℓ^{1/2} > 0 directly. We added this to §1.1(1) (p. 2), keeping the scaling argument credited to
     the reviewer, and to the paragraph after Thm 1.4 (p. 15).
   - **Lemma 3.4** is the antisymmetric Hopf lemma. We use it for the eigenfunction of the lowest odd eigenvalue.
     A two-line proof is given on p. 15 that this eigenfunction has one sign on (0, L). By (2.1) with R ≥ |I|, the
     odd-function form is ½∬_{x,y>0}[(f(x)−f(y))²K(x−y) + (f(x)+f(y))²K(x+y)], with K decreasing. So
     sgn(x)|u(x)| is also a minimiser, and Thm 1.4 gives u = ±sgn·|u|.
   - **Numerical check.** The stored first odd eigenvectors (L = 0.5, 1, 2) have no sign change on the half-mesh.
   - **Sign-changing claim, qualified.** It now reads: "apart from the eigenfunction of the lowest odd eigenvalue,
     no lower bound follows from Hopf-type lemmas".
   - **§1.2 (p. 2)** cites Lemmas 3.4 and 3.11.
   - **Verification section.** It says that PS25 lemma numbers refer to arXiv v1.
4. **Mathia scope (pp. 3, 5).**
   - The notes are stated for Suzuki's A_a = ½L_Δ + γI + K_a, with arithmetic K_a. That the argument covers the
     pure log case (it only uses the source L_Δu = μu, which is bounded and tends to 0) is marked as our inference.
   - PL-297 is now described exactly. It derives τ_s/√s → R^{−1/2} for the principal eigenfunction on intervals
     (−R, R) with μ₁ > 0, from DFW's Hadamard formula and FJW. This "formally corresponds to" A² = 1/(2R), and the
     value of the log trace is not stated there.
   - The words "formal limit" were removed in §1.2, §1.5 and §6.
5. **Same route as BT26.**
   - Method item (1) (p. 4): "Essentially the same route … is taken in [BT26, §3.3, (3.19)]".
   - §1.5 (p. 5): "Our proof of the existence of the trace takes essentially the same Wiener–Hopf route as
     [BT26, §3.3]; only the proof of the trace value, through the dilation identity (6.1), is different."
   - Abstract: the "different, Pohozaev-type" claim is restricted to the trace identity and simplicity.
6. **§1.1 opening (p. 2).** "All three have been answered" is replaced by "Statements and proofs addressing them have
   been posted in recent unrefereed work", followed by: for question (2), that work, like this paper, covers only
   the one-dimensional eigenfunction case of the AGV24 identity.
7. **Ve26 (p. 17).** The text now says that the values of Ve26, "which are given there to six decimals", agree "to
   within 3.5·10⁻⁶". The actual difference is 3.41·10⁻⁶. (Corrected by the final check: an intermediate wording "to
   the six decimals given there" was wrong, since the sixth decimals differ.)
8. **Lemma 6.3 and §1.1(1).**
   - After Lemma 6.3 (p. 14): when AB ≠ 0, the non-uniformity is genuine. For fixed δ, U((δ−t)/θ) → 0 uniformly,
     so the integral tends to 0 and the error tends to |AB|δ.
   - §1.1(1): now says "distances to the boundary scale by θ", and states the conclusion liminf u₁ℓ^{1/2} > 0 at
     both endpoints of every interval.
9. **Verification section (p. 22).** "The release of the numerical code and raw output" now reads "the preparation
   of a numerical supplement with the code and raw output (see Data availability)".
10. **Supplement.**
    - **Drivers.** `run_grid.sh` and `out/run_grid2.sh`, in both `numerics/` and `supplement/numerics/`, now
      `export PY=${PY:-python3}` and run `"$PY" fem.py …`. This was tested: `PY=echo` is honoured, and the default
      is `python3`. The usage comments no longer name the VM or a home-directory path, and the run_grid2 comment now names
      `run_grid2.sh`.
    - **`numerics.md`** (both copies). The host name of the VM is replaced by "a Linux virtual machine with 10 vCPUs" and by
      "the drivers of the recorded runs".
    - **README.** The how-to-run section now explains `PY` and the two edits to the drivers. The `numerics.md` notes
      record the host-name edit.
    - **Checksums and zip.** `SHA256SUMS` was regenerated (128 files, all OK), and `supplement.zip` was rebuilt
      (133 entries; testzip passes; byte-identical to the folder).
    - **Leak check.** Grepping the zip for host names, user names, home-directory paths, private IP prefixes and private e-mail addresses gives no hits, including
      inside the `.npz` members. No `~/` or `.venv` remains in any text file. The only path left is the
      project-relative `papers/purelog/numerics/` in the `numerics.md` header.

**Not done:**
- the two `% AUTHOR` blocks are left to the author;
- nothing was pushed, published or sent.

# Responses to the Codex review of 31 October 2026 (codex_review_2026-10-31.md)

Reviewer: Codex (GPT-6.1 Sol, high setting), an AI system. Verdict: release as a preprint after minor corrections; no fatal or major gap.
Backups taken before editing: `main_before_codex1.tex`, `refs_before_codex1.bib`, `supplement_before_codex1.zip`, `SHA256SUMS_before_codex1` (all in this folder). `refs.bib` is unchanged.

## 1. Priority accuracy: AGV24 contains the separate endpoint identities (Remark 6.4, abstract, §§1.1-1.5)

**Accepted.** Checked against the arXiv v1 text (arXiv:2411.15985v1, 24 Nov 2024, downloaded fresh). Thm 2.3: for a bounded C^{1,1} domain, assuming the boundary limit of u^2 ln(delta^{-2}) exists and L_Delta u is pointwise bounded,
N∫(2F - uf) + 2∫x·F_x + 2∫u^2 = ∫_{∂Ω} u^2 ln(δ^{-2}) (x·ν). Cor 3.8 states the same for bounded weak solutions; its proof is "follows from [25, Thm 1.1] and Theorem 2.3", so it does not show that the limit exists.
With N=1 and f = μu, F = μu²/2, F_x = 0: the interior side is 2∫u², and the boundary value is 2A². On (0,d) (x·ν = 0 at 0, = d at d): 2∫u² = 2d A_+², i.e. A_+² = ||u||²/|I|. On (-d,0): the same for A_-². No star-shapedness or origin condition is imposed in Thm 2.3, so the translates are allowed. The reviewer is right.

Changes:
- Remark 6.4 rewritten: it shows the computation on (-L,L), (0,d) and (-d,0); says the announced, withdrawn formula formally contains each endpoint identity and hence all of Theorem 1.2; says its proof was withdrawn and that the trace existence was assumed (Thm 2.3) and not established (Cor 3.8); says the finite-dilation identity (6.1) is not contained there.
- The sentence after Theorem 1.2 now says the formula, applied with the origin at an endpoint, formally contains Theorem 1.2, with the same two caveats.
- §1.1 item (2) adds that the announced identity on translates of I formally gives each endpoint identity (Remark 6.4).
- Abstract and §1.5: the "summed identity only" wording is replaced by "formally contains the trace identity; its proof and the existence of the trace were not established there".
- The novelty claim for (6.1) is unchanged.

## 2. §8: overstated v_q check

**Accepted.** Verified in `numerics/out/explicit_check.out` (output of `check_rep`): six ξ values (1e-3, 0.3, 1, 10, 1e4, 1e40) for q = 0.5, 1, 2. The largest relative discrepancy is 1.183e-13 (q=2, ξ=1e40).
The text now says: "at six sampled values spanning [1e-3,1e40] (ξ = 1e-3, 0.3, 1, 10, 1e4, 1e40), ... relative discrepancies below 1.2e-13; no uniform check over the interval was made."
The same overstatement in `numerics/numerics.md` (also in the supplement copy) was corrected in the same way. The `numerics/numerics.md` checksum in `supplement/SHA256SUMS` was updated and `supplement.zip` refreshed. All 128 checksums pass, and the unpacked zip matches the directory.

## 3. §8: derivation of the v_q formula actually used

**Accepted.** A new display (8.1) gives the representation, followed by a short Stieltjes-inversion derivation. By (K4), 1/K_q is Stieltjes with measure ν = π^{-1} m̃_0, and v_q = ∫e^{-sx}ν(ds). Kw11 Lemma 3.8, (3.4) (checked against the arXiv text) gives K_q(ξ)K_q(-ξ) = q + ½log(1-ξ²) on C\R. Taking the boundary value at -s+i0:
- for s > 1 the denominator is q + ½log(s²-1) + iπ/2, which gives the density ½K_q(s)/D(s);
- on (0,1) the boundary value is real, with one simple pole at -s_*, where s_* = √(1-e^{-2q}). The derivative of the denominator there is s_* e^{2q}, so the atom is a_* = e^{-2q}K_q(s_*)/s_*;
- there is no mass at s = 1.
The text says this justifies the numerical method and is not used in the proofs. It also keeps the KMR13 (q = 0) attribution.
Verification:
- The formula matches `explicit_A.py` exactly: docstring, `tau_rule` density, and `ts`/`ast` in `check_rep`.
- An independent mpmath check (`codex1_vq_check.py`, output `codex1_vq_check.out`, both in this folder and not in the supplement) computed K_q directly from (3.1) and the Laplace transform of (8.1) by separate mpmath quadrature. For q = 0.5 and 1 at ξ = 0.3, 1, 10, |LT·K - 1| ≤ 1.1e-16 (dps 15).
- Discrimination check: the atom alone accounts for 11% of LT·K at ξ = 1 (q = 1), so the test would catch an error in it.

## 4. §6 opening: domain of J_0*ψ

**Accepted.** J_0*ψ is now defined only for y ∈ R \ closure(I). Absolute convergence follows from ψ ∈ L¹(I) and |x-y| ≥ dist(y,I) > 0. The text says the function is used only in integrals over R\I and at the points L+t, -L-t with t > 0, and that endpoint values (a null set, where the integral may diverge) play no role. All later uses (Lemmas 6.1, 6.2, Prop. 7.1) are exterior uses.

## 5. Proof of Theorem 1.1: explicit constant

**Accepted.** After the parameter-only replacements, the constant is called C_E^*(L,λ). Theorem 5.1(iii) is applied to u/||u|| and to its reflection, which has the same L, λ and sup norm. This bounds each endpoint error by C_E^*||u||ℓ^{-1} on 0 < r ≤ ρ_0. So C = 2C_E^*||u||_2, i.e. C(b-a,μ) = 2C_E^*((b-a)/2, μ/2).

## 6. Exposition: concise abstract and introduction

**Accepted.**
- **Abstract.** The earlier-work part is now two sentences: (i) one acknowledgement of BT26 (rStation, 27 Oct 2026, unrefereed), Mathia26 (trace existence derived, remainder of the same order stated, localized Weil operator) and the withdrawn AGV24 formula with its caveats; (ii) a precise statement of what differs: the self-contained remainder bound with explicit constant, the different Pohozaev-type proof of the trace identity and simplicity, and the new dilation and even-odd identities.
- **§1.2.** The long BT26/Mathia26 paragraph is replaced by a pointer to §1.5 and Appendix B. The Veltri sentence is kept.
- **§1.5** is rewritten:
  - one paragraph acknowledging BT26 and Mathia26 and their status (AI-agent acceptance, not checked by a human, not established results; the same holds for this paper; we have not verified all their proofs);
  - a list (a)-(c) of what differs: rate and explicit constant; dilation-identity route versus collar/Schur; new identities, with the BT26 (6.3) reduction;
  - AGV24, BT26c parity and "none of the tools is new";
  - independent discovery, limited searches, no certification.
- **New Appendix B, "Detailed comparison with earlier public work".** It holds the note-by-note material formerly in §§1.2 and 1.5, with every qualification kept:
  - BT26: handle, AI-reviewer acceptance, not human-checked, Cor 2.3 in its notation, GR12 route, (3.19) without rate, Thms 2.1-2.2 collar/Schur, Hadamard form, (6.3), §10 Galerkin values;
  - Mathia26: A_a definition, PL-319-323 with o(ℓ^{-1}), "our inference" for the pure-log case, PL-329/336/337 details including "asserted rather than written out", conditional L∞ bound, pointwise realisation not written out, "not checked line by line", PL-297 and PL-306;
  - AGV24: withdrawn, formally contains Thm 1.2, trace existence not established, (6.1) not contained.
- **Organisation paragraph.** Now mentions Appendix B.

## 7. Typos

**Accepted.**
- Lemma 2.5 proof: a simultaneous translation multiplies f̂ and ĝ by the same phase. In f̂ times the conjugate of ĝ the phases cancel, so the product is unchanged. Reflection replaces t by -t, and log|t| is even.
- §2.1: "Ein(z) ≤ z for z > 0".

## Verification section

A new paragraph covers the second external review: by Codex, an AI system of the same family (6.1 Sol, high setting), 31 October 2026, on the revised version. It records:
- that it is not a human expert and not a certification;
- what it checked and what it did not rerun;
- the verdict: release as a preprint after minor corrections;
- the list of corrections made in this version.

The first-review paragraph now says its changes "were made in the revised version" and points to §1.5 and Appendix B.

## Build

`tectonic main.tex`: no errors, no undefined references or citations, and no overfull or underfull boxes; 26 pages. An overfull \vbox on page 1 from the first abstract draft was removed by tightening the abstract. No forbidden strings (private e-mail addresses, host names, local paths) appear in main.tex, refs.bib or the PDF text.

## Release check (31 October 2026)

A final check of these responses (verdict: minor) confirmed points 1-5 and 7 against the sources and the stored output, including an independent check of the v_q representation (q = 0.5, 1, 2; ξ from 0.01 to 1e20; |LT·K - 1| ≤ 1.8e-9). It found one regression from point 6 and one cosmetic issue, both fixed before release:
- **Abstract.** The tightened abstract had dropped the credit that our existence proof follows essentially the same Wiener–Hopf route as BT26. Restored: "What differs here is a self-contained proof, along essentially the same Wiener–Hopf route as the rStation manuscript, of the remainder bound with an explicit constant, ..."
- **Appendix B, BT26 paragraph.** Cor. 2.3 was stated twice in a row, the first time with a used before it was defined. The first statement was removed; the one that defines the notation is kept.
- **Lemma 2.5 proof.** The phase sentence now reads: the Fourier phase e^{ict} of f̂ cancels the phase e^{-ict} of conj(ĝ).

Recompiled with tectonic: no errors, no undefined references, no overfull or underfull boxes; 26 pages.

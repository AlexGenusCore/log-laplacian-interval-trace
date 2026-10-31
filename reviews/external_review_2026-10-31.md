# External review of papers/purelog/main.pdf (21 pp), received 2026-10-31 (pasted by the author)

Assessment: central argument coherent; no fatal gap found in this review; novelty claim needs revision; several smaller
statements need tightening. A mathematical review, not a certification of the proof. All 21 pages and key references
checked; numerics not reproduced (code/raw outputs not attached).

Main finding: earlier overlapping manuscript "Boundary traces and support variation for localized Weil operators",
posted 27 October 2026 under bramble_trace; its Corollary 2.3 states every interval eigenvalue of the compressed log|D|
is simple and normalized eigenfunctions have endpoint trace squares 1/(2a). Reviewer did not validate its proof.
"If this is your earlier version, cite it and explain that the present paper supplies a more focused proof and an
explicit boundary error estimate. If it is independent work, compare the arguments and narrow the priority claim."

Checks passed: AGV withdrawal (date, reason); HLS optimal boundary bound + Hopf lower bound (not the exact limit);
Kwasnicki arXiv v4 Lemma 3.8, Example 6.4; fractional one-endpoint Pohozaev simplicity argument credited correctly;
Table 1 vs HJSS max discrepancy 2.1369e-4 (< 2.2e-4); KW26 exists (abstract only). Normalisation: exterior kernel
factor 1/2 consistent with Q = E_log/2; (J0*u)(L+t) ~ A_+ sqrt(log(1/t)); layer gives A_+^2 delta, other side
||u||^2 delta/(2L); hence A_+^2 = A_-^2 = ||u||^2/|I|, no missing factor 2. Prop 4.1 Hardy-step support signs
consistent. Simplicity deduction valid given the trace identity.

Requested changes:
1. Lemma 6.3 (p.12): error has log(theta)/log(1/delta); limit not uniform in unrestricted theta >= 1. State for fixed
   theta, uniformly for bounded theta, or under log theta = o(log(1/delta)); main proof uses theta = 1 + delta/(2L).
2. p.15: exact conforming Galerkin eigenvalues are upper bounds; computed (quadrature) values approximate them.
3. Table 2 (pp.15-16): error estimates exclude the cutoff; explicit formula uses a tail model; evidence, not rigorous.
4. p.1: Hopf applies directly when mu_1 >= 0; with the scaling law it extends to the first eigenfunction on any
   interval (reviewer's deduction). Genuine additional claim: existence of the exact trace, esp. sign-changing.
5. Provide the numerical supplement.
Conclusion: proceed to specialist review after these amendments; the overlapping manuscript must be addressed first.

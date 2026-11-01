# A boundary trace identity for the logarithmic Laplacian on an interval, and simplicity of its Dirichlet spectrum

Johann Alexander Botha (GenusCore, info@genuscore.com). Version 1.0.0, 31 October 2026.

The paper is [`paper/main.pdf`](paper/main.pdf) (source `paper/main.tex`, `paper/refs.bib`, `paper/fig_edge.pdf`).
The numerical supplement is in [`supplement/`](supplement/README.md). The review records are in [`reviews/`](reviews/).

> **Status. Please read this first.**
> - The paper is **unrefereed**. No human expert has checked the proofs or the code.
> - It was produced **with substantial AI assistance** (Anthropic's Claude, used through Claude Code), directed by
>   the author, who takes full responsibility for the content.
> - The only reviews so far are by AI systems: several internal referee passes by Claude-based referees, and two
>   external reviews on 31 October 2026 by a separate AI system (6.1 Sol, high setting; the second one run through
>   Codex). None of them is a certification. The second external review recommended release as a preprint after
>   minor corrections, which have been made (see [`reviews/`](reviews/)).
> - The parts that most need an expert reading are Step 4 of Proposition 4.1 (the Hardy-space step) and Appendix A.
> - The main statements are **not new as statements**; see "Earlier work" below.
> - The paper is about the logarithmic Laplacian. It makes no claim about the Riemann hypothesis.

## Progress

Public releases across the author's three related repositories, as of 1 November 2026. All are **unrefereed** and
**AI-assisted** (Anthropic's Claude, used through Claude Code, directed by the author). Every review so far was run by AI
systems: internal Claude-based referees, and for the later items a second AI system, OpenAI's Codex. **No human expert has
checked any of them.** Independent review is welcome: info@genuscore.com.

| Date | Release | Result | Review status |
|---|---|---|---|
| 2026-10-25 | [weil-ground-state-simple-even](https://github.com/AlexGenusCore/weil-ground-state-simple-even), first release (commit `26fd9e6`; no tag, no DOI) | Certified simple, even ground state of the localized Weil form for `L` in `(0, 1.09]` | Internal AI review rounds only; the windows 1.03–1.09 had no review round |
| 2026-10-26 | weil-ground-state-simple-even, version of 2026-10-26 (commit `8adfbe5`; no tag, no DOI) | Range extended to `(0, 1.19]` | Internal AI review only; the windows 1.03–1.19 had no review round |
| 2026-10-31 | weil-ground-state-simple-even v1.1.0, archived as v1.1.1 (credits corrected): [10.5281/zenodo.23832221](https://doi.org/10.5281/zenodo.23832221) | Range extended to `(0, 1.4035]` (Piece III, with Elston's certifier) | Internal AI review only; Piece III not audited beyond `L = 1.19` |
| 2026-10-31 | [weil-edge-laws](https://github.com/AlexGenusCore/weil-edge-laws) v0.1.2, research notes: [10.5281/zenodo.23832560](https://doi.org/10.5281/zenodo.23832560) | Edge law, Hadamard formula, kink (echo) and tail laws | Internal AI referee rounds only |
| 2026-10-31 | [log-laplacian-interval-trace](https://github.com/AlexGenusCore/log-laplacian-interval-trace) v1.0.0, paper: [10.5281/zenodo.23832744](https://doi.org/10.5281/zenodo.23832744) | Boundary trace identity and simple Dirichlet spectrum for the logarithmic Laplacian on an interval | Internal AI referees; two reviews by GPT-6.1 Sol (the second through Codex): release as a preprint after minor corrections, which were made |
| 2026-10-31 | weil-edge-laws v0.2.0, paper: [10.5281/zenodo.23834265](https://doi.org/10.5281/zenodo.23834265) | The notes consolidated into a paper with complete proofs, pending external review | Internal AI referees; Codex (GPT-6.1 Sol): release as a preprint after minor fixes, which were made; the changes were not reviewed again |
| 2026-11-01 | weil-ground-state-simple-even v2.0.0: [10.5281/zenodo.23852617](https://doi.org/10.5281/zenodo.23852617) | Range extended to `(0, 2]` (Piece IV) | Fresh-agent AI audits of Piece IV; two Codex (GPT-6.1 Sol) rounds, the first "not yet" (three major points), the second "release as a preprint after minor corrections"; all points addressed; the round-two fixes were not reviewed again. Neither round re-ran the large certificates |

Several of these results, or parts of them, were obtained earlier by others; each repository's README and paper give
the credits (among them Zhu, Kramarenko-Byrd, Elston, the rStation authors bramble_trace and lihongmin, and the mathia
notes). None of the results is a result about the Riemann hypothesis. Work on extending the range and on computer-light
proofs continues; nothing beyond the releases above is claimed.

## Results in plain terms

The logarithmic Laplacian $L_\Delta$ is the operator with Fourier symbol $2\log|\xi|$ (the derivative at $s=0$ of
the fractional Laplacian $(-\Delta)^s$). Consider its Dirichlet eigenfunctions on a bounded interval $I$. They
vanish at the endpoints, but only logarithmically slowly. Write $\ell(x)=\log(1/\mathrm{dist}(x,\partial I))$.

1. **Boundary trace (Theorem 1.1).** Every eigenfunction $u$ has an exact edge law: at each endpoint the limit
   $A_\pm(u)$ of $u(x)\,\ell(x)^{1/2}$ exists, and $u(x)=A_\pm(u)\,\ell(x)^{-1/2}+O(\ell(x)^{-1})$ with an explicit
   constant that depends only on $|I|$, the eigenvalue and $\|u\|_{L^2}$.
2. **Trace identity (Theorem 1.2).** $A_+(u)^2=A_-(u)^2=\|u\|_{L^2}^2/|I|$ for every eigenfunction.
3. **Simplicity (Theorem 1.4).** Every Dirichlet eigenvalue of $L_\Delta$ on an interval is simple, so every
   eigenfunction is even or odd. As in Kwaśnicki's short proof for the fractional Laplacian, this follows from the
   trace identity.
4. **Further identities.** An exact one-endpoint dilation identity (6.1), which compares an eigenfunction with its
   dilation about one endpoint and from which Theorem 1.2 follows; a translation identity for pairs of
   eigenfunctions (Proposition 7.1); and a mixed identity for even–odd pairs, $(\mu_e-\mu_o)\Phi'(0)=\pm2/L$
   on $(-L,L)$ (Corollary 1.5). Parity alternation (Corollary 9.3) follows from known results and Theorem 1.4.
5. **Numerics (Section 8).** Uncertified P1 finite-element computations illustrate the results: eigenvalues on
   $(-1,1)$, and $2LA_+^2-1$ for the six lowest eigenfunctions at $L=0.5,1,2$. They also support a conjectured
   second-order term (Conjecture 9.2). None of the numerics is used in the proofs.

The proofs combine the Wiener–Hopf factorisation of the variance-gamma symbol on a half-line with the fact that
the interaction kernel of $L_\Delta$ between disjoint sets is exactly $|x-y|^{-1}$.

## Earlier work, and what is ours

**Existence of the trace, the trace identity and simplicity were obtained earlier by others.** Full details are in
§1.5 and Appendix B of the paper.

- **bramble_trace, "Boundary traces and support variation for localized Weil operators"**, rStation paper
  riemann_hypothesis 67, v1, published 27 October 2026,
  <https://rstation.ai/problems/riemann_hypothesis/papers/67/v1> (cited as [BT26]). rStation is a research
  platform for AI agents, and the paper was accepted there by the platform's AI reviewer agents. As far as we
  know, no human has checked its proofs, so it is an overlapping public claim whose proof still needs checking, as
  is this paper. It proves the existence of the boundary traces (§3.3, without a rate), and, for $\log|D|=\frac12
  L_\Delta$ on an interval, that every eigenvalue is simple and every normalised eigenfunction has squared traces
  $1/|I|$ at both endpoints (Corollary 2.3). That is our Theorem 1.2 and the simplicity part of Theorem 1.4. Its
  (6.3) is the equal-eigenvalue case of our Proposition 7.1, and its forum comment [BT26c] posted the parity
  alternation earlier.
- **Sergio Garcia Murillo, Mathia research notes**, prime_lattice line, <https://github.com/murillo128/mathia>
  (cited as [Mathia26]). Notes PL-319 to PL-323 (September 2026) derived the existence of the trace for
  eigenstates of the localized Weil operator, and notes PL-329, PL-336 and PL-337 state a remainder of the same
  order $O(\ell^{-1})$ as ours. These notes are not refereed either.
- A Pohozaev identity announced in arXiv:2411.15985 (since withdrawn by its authors, with a proof that was withdrawn
  and without establishing that the trace exists) formally contains the trace identity.

**What differs here:**
- a self-contained proof of the remainder bound with an explicit constant, for the pure logarithmic operator,
  along essentially the same Wiener–Hopf route as [BT26, §3.3];
- a different, Pohozaev-type proof of the trace identity and of simplicity, through the exact dilation identity
  (6.1), instead of the support-variation (Hadamard-type) route of [BT26];
- to our knowledge new: the dilation identity (6.1), Proposition 7.1 for eigenfunctions with different eigenvalues,
  and the even–odd identity of Corollary 1.5.

None of the tools is new. Our literature searches were limited, and not finding a result is not proof that it is
new. The proofs here were found independently of [BT26], which we learned of from an external review; we knew
parts of [Mathia26] while working on the paper.

## Reproducing the numerics

Everything is in `supplement/`; its [README](supplement/README.md) lists every file, the environment, run times
and how Tables 1–2, Figure 1 and the evidence for Conjecture 9.2 are produced. In short:

```bash
cd supplement && sha256sum -c SHA256SUMS      # 128 files
cd numerics
python fem.py --selftest                      # needs numpy and SciPy >= 1.5
bash run_grid.sh                              # main grid (28 runs) -> out/
python explicit_A.py check > out/explicit_check.out   # needs mpmath
python analyse.py                             # -> results.json, tables.tex
python plot_edge.py out/L1_n400_r1.12.npz fig_edge.pdf  # needs matplotlib
python residual_check.py --test
python ../checks/referee_p0check.py 1 1e-100 1.12 400   # independent P0 check
```

The raw output of every recorded run is included in `supplement/numerics/out/`, so the tables can be rebuilt with
`analyse.py` without rerunning the solver. Package versions of the recorded runs were not recorded. The
computations are floating-point and uncertified; their error estimates are refinement estimates, not rigorous
bounds.

To rebuild the paper: `cd paper && tectonic main.tex` (or `pdflatex main && bibtex main && pdflatex main && pdflatex main`).

## Reviews

[`reviews/`](reviews/) holds the summary of the first external review, the second external review (run through
Codex) in full, and our responses to both. The responses mention working files (search logs, comparison notes,
backups) that are not part of this repository.

## How to cite

See [`CITATION.cff`](CITATION.cff), or:

```bibtex
@techreport{Botha2026LogLaplacianTrace,
  author      = {Botha, Johann Alexander},
  title       = {A boundary trace identity for the logarithmic {L}aplacian on an interval, and simplicity of its {D}irichlet spectrum},
  institution = {GenusCore},
  type        = {Preprint},
  year        = {2026},
  month       = oct,
  note        = {Version 1.0.0, 31 October 2026. Unrefereed},
  url         = {https://github.com/AlexGenusCore/log-laplacian-interval-trace}
}
```

Please also cite [BT26] and [Mathia26] for the results they obtained first.

## Related repositories

- <https://github.com/AlexGenusCore/weil-ground-state-simple-even>: a computer-assisted proof that the ground state
  of the localized Weil quadratic form is simple and even on a range of window sizes.
- <https://github.com/AlexGenusCore/weil-edge-laws>: research notes on the edge law, the Hadamard formula, and
  echo and tail laws for eigenfunctions of the localized Weil quadratic form, whose archimedean part is a bounded
  perturbation of $\frac12L_\Delta$.

## License

The paper and all non-code files are under CC BY 4.0 ([`paper/LICENSE`](paper/LICENSE)); the code in
`supplement/` is under the MIT License ([`LICENSE`](LICENSE)).

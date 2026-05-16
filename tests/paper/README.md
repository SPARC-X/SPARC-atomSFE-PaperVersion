# Paper results

LaTeX fragments, figures, and matching tests for the SPARC-atomSFE paper. Each topic has its own subdirectory under `paper/`; shared figures are in [`figures/`](figures/). Full text, equations, tables, and captions are in the linked `.tex` files.

---

## Radial Schrödinger

Verifies the spectral finite-element framework for the all-electron radial Schrödinger equation by comparing occupied eigenvalues for neutral uranium to hydrogenic reference energies, with errors below 0.1 nano-Hartree.

| | |
|--|--|
| TeX | [`schrodinger/z92_hydrogenic.tex`](schrodinger/z92_hydrogenic.tex) |
| Standalone PDF preview | [`schrodinger/z92_hydrogenic_standalone.tex`](schrodinger/z92_hydrogenic_standalone.tex) |
| Test | [`../test_z92_schrodinger.py`](../test_z92_schrodinger.py) |
| Figure | [`figures/test_z92_schrodinger_summary.pdf`](figures/test_z92_schrodinger_summary.pdf) |

Refresh the figure with `REGENERATE_SUMMARY_PDF = True`, then `python ../test_z92_schrodinger.py`. Preview TeX from `schrodinger/`: `pdflatex z92_hydrogenic_standalone.tex` (requires the PDF in `figures/`).

---

## All-electron discretization convergence

Narrative and figures for Kohn--Sham mesh convergence ($R_{max}$, $N_{fe}$) in the all-electron setting. The pseudopotential companion fragment is [`convergence/pseudopotential_convergence.tex`](convergence/pseudopotential_convergence.tex).

| | |
|--|--|
| TeX | [`convergence/all_electron_convergence.tex`](convergence/all_electron_convergence.tex) |
| Standalone PDF preview | [`convergence/all_electron_convergence_standalone.tex`](convergence/all_electron_convergence_standalone.tex) |
| Figures | [`figures/gga_pbe_convergence_test_summary.pdf`](figures/gga_pbe_convergence_test_summary.pdf), [`figures/boundary_subset_rule_ae.pdf`](figures/boundary_subset_rule_ae.pdf), [`figures/lda_svwn_convergence_test_featom_summary.pdf`](figures/lda_svwn_convergence_test_featom_summary.pdf), [`figures/lda_svwn_convergence_test_summary.pdf`](figures/lda_svwn_convergence_test_summary.pdf), [`figures/rscan_convergence_test_summary.pdf`](figures/rscan_convergence_test_summary.pdf) |

From `tests/data/compare/`, run `python gga_pbe_convergence_test.py`, `python lda_svwn_convergence_test_featom.py`, `python lda_svwn_convergence_test.py`, and `python rscan_convergence_test.py`; copy the four matching `*_summary.pdf` files into `paper/figures/`. Place `boundary_subset_rule_ae.pdf` in `paper/figures/` (schematic of nested AE mesh subsets). Preview: `cd convergence` then `pdflatex all_electron_convergence_standalone.tex`.

---

## Pseudopotential discretization convergence

Same layout as the all-electron subsection: two paragraphs of main text, PBE two-panel figure, supplementary nested-mesh schematic, then LDA-SVWN and rSCAN panels.

| | |
|--|--|
| TeX | [`convergence/pseudopotential_convergence.tex`](convergence/pseudopotential_convergence.tex) |
| Standalone PDF preview | [`convergence/pseudopotential_convergence_standalone.tex`](convergence/pseudopotential_convergence_standalone.tex) |
| Figures | [`figures/pseudo_gga_pbe_convergence_test_summary.pdf`](figures/pseudo_gga_pbe_convergence_test_summary.pdf), [`figures/boundary_subset_rule_psp.pdf`](figures/boundary_subset_rule_psp.pdf), [`figures/pseudo_lda_svwn_convergence_test_summary.pdf`](figures/pseudo_lda_svwn_convergence_test_summary.pdf), [`figures/pseudo_rscan_convergence_test_summary.pdf`](figures/pseudo_rscan_convergence_test_summary.pdf) |

From `tests/data/compare/`, run `python pseudo_gga_pbe_convergence_test.py`, `python pseudo_lda_svwn_convergence_test.py`, and `python pseudo_rscan_convergence_test.py`; copy the three `*_summary.pdf` files into `paper/figures/`. Place `boundary_subset_rule_psp.pdf` in `paper/figures/`. Preview: `cd convergence` then `pdflatex pseudopotential_convergence_standalone.tex`.

---

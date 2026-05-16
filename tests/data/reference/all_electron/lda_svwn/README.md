# LDA_SVWN reference dataset

## `featom_atoms_lda.json`

- Citation: Čertík, O.; Pask, J. E.; Fernando, I.; Goswami, R.; Sukumar, N.; Collins, L. A.; Manzini, G.; Vackář, J. "High-order finite element method for atomic structure calculations." *Comput. Phys. Comm.* **297**, 109051 (2024).
- Source: generated locally with `featom` (`solve_schroed`) for `Z=1..92`.

At the top level, the JSON stores metadata (`generator`, `xc_label`), run settings (`settings`), summary counts (`n_success`, `n_failed`), failed cases (`failures`), and the per-atom data list (`records`).

The `settings` block records the run parameters used to generate this dataset: `z_min`/`z_max` are the atomic-number range; `eps` is the SCF convergence tolerance; `rmin` and `rmax` are the radial-domain bounds (bohr); `a` is the exponential-mesh concentration parameter; `Ne` is the number of finite elements; `Nq` is the number of quadrature points per element; and `p` is the polynomial order.

### Per-species fields mapping (`records` entries)

| Physical meaning | Mathematical expression | Variable name (JSON key) |
| --- | --- | --- |
| Atomic number | $Z$ | `Z` |
| Number of occupied states | $N_{\mathrm{occ}}$ | `n_occupied_states` |
| Total energy | $E_{\mathrm{tot}}$ (a.u.) | `total_energy_au` |
| Occupied eigenvalues | $\epsilon_i$ for occupied $i$ (a.u.) | `occupied_eigenvalues_au` |

**Notes**

- Failed atoms (non-converged SCF in the selected settings) are recorded under `failures` with the corresponding `Z` and tail error text.
- The dataset is not restricted to closed-shell species; it follows a full `Z=1..92` sweep and keeps whichever cases converge under the run settings.

## `featom_z92_fe16_R040_reference.json`

- Source: FEATOM (`solve_schroed`) **single-Z** run for **uranium (Z = 92)** at the mesh used in `atom/tests/data/lda_svwn/finite_element_sweep/fe16_R040` (see `settings.dataset_dir`).
- Purpose: occupied eigenvalues (and total energy) as the reference vector in `tests/data/compare/lda_svwn_convergence_test_featom.py` for `configuration_092` in LDA_SVWN summaries.

The `settings` block uses the **same field names** as sweep summary `input_parameters` where applicable (`atomic_number`, `scf_tolerance`, `domain_size`, `finite_element_number`, `polynomial_order`, `quadrature_point_number`, `mesh_concentration`), plus `dataset_dir`, `runner`, and `generated_at_utc`.

Schema matches the multi-species file (`records` with one entry, `n_success` / `failures` at top level).

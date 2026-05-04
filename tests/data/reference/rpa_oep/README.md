# RPA-OEP reference dataset

## `shubhang_z_1_to_18_rpa_oep.json`

- Citation: Trivedi, S. K.; Suryanarayana, P. "Spectral finite-element formulation of the optimized effective potential method for atomic structure in the random phase approximation." *arXiv e-prints* (2025), arXiv:2512.
- Source: compiled from RPA-OEP reference summaries for `Z=1..18`.

At the top level, the JSON stores dataset metadata (`dataset_dir`), run settings (`input_parameters`), summary counts (`n_configurations`, `n_summarized`), and the per-configuration data list (`config_summaries`).

The `input_parameters` block records the run parameters used to generate this dataset: `all_electron_flag` and `spin_polarized_flag` are model flags; `xc_functional` is the XC label; `use_oep` indicates OEP usage; `domain_size`, `finite_element_number`, `polynomial_order`, `quadrature_point_number`, and `oep_basis_number` define the radial/OEP discretization; `mesh_type`, `mesh_concentration`, and `mesh_spacing` define the radial mesh; `scf_tolerance` is the SCF convergence tolerance; `frequency_quadrature_point_number` and `angular_momentum_cutoff` control the RPA integration/cutoff settings; and `double_hybrid_flag` and `enable_parallelization` are workflow toggles.

### Per-configuration fields mapping (`config_summaries` entries)


| Physical meaning                | Mathematical expression                     | Variable name (JSON key)      |
| ------------------------------- | ------------------------------------------- | ----------------------------- |
| Configuration label             | —                                           | `configuration`               |
| Atomic number                   | $Z$                                         | `atomic_number`               |
| Electron count                  | $N_e$                                       | `n_electrons`                 |
| Convergence status              | —                                           | `converged`                   |
| Total energy                    | $E_{\mathrm{tot}}$ (Ha)                     | `total_energy_ha`             |
| Number of occupied states       | $N_{\mathrm{occ}}$                          | `occupied_state_count`        |
| Occupied eigenvalues            | $\epsilon_i$ for occupied $i$ (Ha)          | `occupied_eigenvalues_ha`     |
| Lowest 5 unoccupied eigenvalues | $\epsilon_a$ for lowest unoccupied $a$ (Ha) | `lowest_5_unoccupied_eigs_ha` |


### Energy-component fields mapping (`energies_ha`)


| Physical meaning          | Mathematical expression              | Variable name (JSON key) |
| ------------------------- | ------------------------------------ | ------------------------ |
| Radial kinetic energy     | $T_r$ (Ha)                           | `kinetic_radial`         |
| Angular kinetic energy    | $T_l$ (Ha)                           | `kinetic_angular`        |
| Total kinetic energy      | $T=T_r+T_l$ (Ha)                     | `total_kinetic`          |
| External potential energy | $E_{\mathrm{ext}}$ (Ha)              | `external_potential`     |
| Hartree energy            | $E_{\mathrm{H}}$ (Ha)                | `hartree`                |
| Exchange energy           | $E_{\mathrm{x}}$ (Ha)                | `exchange`               |
| Correlation energy        | $E_{\mathrm{c}}$ (Ha)                | `correlation`            |
| Exact exchange (HF-style) | $E_{\mathrm{x}}^{\mathrm{HF}}$ (Ha)  | `exact_exchange_hf`      |
| RPA correlation           | $E_{\mathrm{c}}^{\mathrm{RPA}}$ (Ha) | `rpa_correlation`        |
| Total potential energy    | $V_{\mathrm{tot}}$ (Ha)              | `total_potential`        |
| Total energy              | $E_{\mathrm{tot}}$ (Ha)              | `total_energy`           |


**Notes**

- Some optional terms (for example `exact_exchange_hf` or `rpa_correlation`) can be `null` depending on how the source summaries were exported.
- `config_summaries` contains one entry per configuration in the selected `Z=1..18` range.


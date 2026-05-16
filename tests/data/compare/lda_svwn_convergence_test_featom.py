r"""
FEATOM reference vs LDA_SVWN sweep: per-state absolute eigenvalue error.

Left: domain-radius sweep — one curve per committed case ($N_{\mathrm{fe}}=16$ omitted; odd/even retained).
Right: fixed $R$, vary $N_{\mathrm{fe}}$ (even $N_{\mathrm{fe}}$ only; $N_{\mathrm{fe}}=16$ omitted).

Reference occupied eigenvalues (Hartree, 18 states) are read from
``tests/data/reference/all_electron/lda_svwn/featom_z92_fe16_R040_reference.json``
(single-Z FEATOM run; ``settings`` keys match sweep ``input_parameters`` naming).
Summaries are read from committed files under
``tests/data/summary/all_electron/lda_svwn/<sweep>/`` as flat
``fe*_R*__*.json`` files (see ``build_summary_from_out.py`` / ``summary_naming``).
For each sweep case file, the row matching ``configuration`` (default U,
``configuration_092``) supplies ``occupied_eigenvalues_ha`` for comparison to
the reference.

Outputs PDF only (same font defaults as ``pseudo_gga_pbe_convergence_test.py``).

Run from anywhere (defaults read ``atomSFE/tests/data/summary/all_electron/lda_svwn/``,
write ``atomSFE/tests/data/compare/lda_svwn_convergence_test_featom_summary.pdf``)::

    python atomSFE/tests/data/compare/lda_svwn_convergence_test_featom.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Iterator

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (
    LogFormatterMathtext,
    LogLocator,
    NullLocator,
)

_DATA_DIR = Path(__file__).resolve().parent.parent
if str(_DATA_DIR) not in sys.path:
    sys.path.insert(0, str(_DATA_DIR))
from summary_naming import glob_sweep_summaries, mesh_tag_from_summary_path

_COMPARE_DIR = Path(__file__).resolve().parent
_SUMMARY_DIR = _DATA_DIR / "summary"
_DEFAULT_LDA_SVWN_ROOT = _SUMMARY_DIR / "all_electron" / "lda_svwn"
_DEFAULT_FEATOM_REF_JSON = (
    _DATA_DIR / "reference" / "all_electron" / "lda_svwn" / "featom_z92_fe16_R040_reference.json"
)
_DEFAULT_CONFIGURATION = "configuration_092"
_DEFAULT_FIXED_R = 40
_DEFAULT_OUT_PDF = _COMPARE_DIR / "lda_svwn_convergence_test_featom_summary.pdf"

plt.rcParams.update(
    {
        "font.family": "serif",
        "mathtext.fontset": "stix",
        "font.size": 15,
    }
)

AXIS_LABEL_FONTSIZE = 18
X_AXIS_LABEL_FONTSIZE = 18
TICK_LABEL_FONTSIZE = 15
LEGEND_FONTSIZE = 13  # single-column legend: smaller than axis titles so many entries fit

# Figure size scaled up vs ``(12, 4.8)`` in ``pseudo_gga_pbe_convergence_test.py`` for the larger PDF fonts.
_FIG_WIDTH_IN = 16.0
_FIG_HEIGHT_IN = 6.4

# Distinct curve colors (matplotlib tab10); no extra repo deps beyond ``tests/data``.
_COLORS: tuple[str, ...] = (
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
)


def _finalize_paper_semilogy_axes(ax: plt.Axes, *, x_mode: str = "fe") -> None:
    """Log-y grid/ticks in a publication-friendly style; no minor ticks on x (``x_mode='fe'``)."""
    _ = x_mode  # reserved for parity with older callers
    ax.grid(True, which="major", linestyle="-", linewidth=0.75, alpha=0.32)
    ax.grid(True, which="minor", axis="y", linestyle=":", linewidth=0.75, alpha=0.45)
    ax.yaxis.set_major_locator(LogLocator(base=10.0))
    ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", labelsize=TICK_LABEL_FONTSIZE)
    ax.tick_params(axis="both", which="minor", labelsize=TICK_LABEL_FONTSIZE - 1)


def load_z92_featom_reference_eigenvalues(
    path: Path = _DEFAULT_FEATOM_REF_JSON,
    *,
    z: int = 92,
) -> np.ndarray:
    """Occupied eigenvalues (Ha) from a FEATOM single-Z JSON (``records``), fixed list order."""
    data = json.loads(path.read_text(encoding="utf-8"))
    for rec in data.get("records", []):
        if int(rec.get("Z", 0)) != int(z):
            continue
        occ = rec.get("occupied_eigenvalues_au")
        if occ is None:
            continue
        return np.asarray(occ, dtype=float)
    raise ValueError(f"No occupied_eigenvalues_au for Z={z} in {path}")


ref_eigenvalues = load_z92_featom_reference_eigenvalues()
N_REF = int(ref_eigenvalues.size)


def parse_mesh_tag(mesh: str) -> tuple[int, int]:
    """Parse ``fe12_R040``-style mesh token into (finite_element_number, domain_radius_bohr)."""
    parts = mesh.split("_", 1)
    if len(parts) != 2:
        raise ValueError(f"Unrecognized mesh tag: {mesh!r}")
    fe_txt, r_txt = parts
    return int(fe_txt.replace("fe", "")), int(round(float(r_txt.replace("R", ""))))


def read_occupied_eigenvalues_from_summary_json(
    summary_path: Path, configuration: str
) -> np.ndarray | None:
    if not summary_path.is_file():
        return None
    data = json.loads(summary_path.read_text(encoding="utf-8"))
    for row in data.get("config_summaries", []):
        if row.get("configuration") != configuration:
            continue
        if not row.get("converged", True):
            return None
        occ = row.get("occupied_eigenvalues_ha")
        if occ is None:
            return None
        return np.asarray(occ, dtype=float)
    return None


def iter_sweep_summary_paths(sweep_root: Path) -> Iterator[Path]:
    if not sweep_root.is_dir():
        return
    yield from glob_sweep_summaries(sweep_root)


def _include_fe_for_domain_radius_sweep(fe: int) -> bool:
    """Keep all ``N_fe`` except 16 (same mesh as the FEATOM reference JSON)."""
    return int(fe) != 16


def _include_fe_for_finite_element_sweep(fe: int) -> bool:
    """Even ``N_fe`` only; omit 16 (avoids the trivial FEATOM-reference mesh)."""
    n = int(fe)
    return n % 2 == 0 and n != 16


def collect_radius_curves(
    sweep_dir: Path,
    configuration: str,
    ref: np.ndarray,
) -> list[tuple[int, float, np.ndarray]]:
    """Return list of (N_fe, R_bohr, abs_error[0:n_ref]); omit only ``N_fe=16``.

    Sorted by R then N_fe. (Domain-radius summaries often use odd ``N_fe`` per radius; do not require even here.)
    """
    rows: list[tuple[int, float, np.ndarray]] = []
    for summary_path in iter_sweep_summary_paths(sweep_dir):
        try:
            fe, r_bohr = parse_mesh_tag(mesh_tag_from_summary_path(summary_path))
        except ValueError:
            continue
        if not _include_fe_for_domain_radius_sweep(fe):
            continue
        ev = read_occupied_eigenvalues_from_summary_json(summary_path, configuration)
        if ev is None or ev.size < ref.size:
            continue
        err = np.abs(ev[: ref.size] - ref)
        rows.append((int(fe), float(r_bohr), err))
    rows.sort(key=lambda t: (t[1], t[0]))
    return rows


def collect_fe_curves(
    sweep_dir: Path,
    configuration: str,
    fixed_r: int,
    ref: np.ndarray,
) -> list[tuple[float, np.ndarray]]:
    """Return list of (N_fe, abs_error[0:n_ref]) sorted by N_fe; even N_fe only, excluding N_fe=16."""
    rows: list[tuple[float, np.ndarray]] = []
    for summary_path in iter_sweep_summary_paths(sweep_dir):
        try:
            fe, r_bohr = parse_mesh_tag(mesh_tag_from_summary_path(summary_path))
        except ValueError:
            continue
        if not _include_fe_for_finite_element_sweep(fe):
            continue
        if int(r_bohr) != int(fixed_r):
            continue
        ev = read_occupied_eigenvalues_from_summary_json(summary_path, configuration)
        if ev is None or ev.size < ref.size:
            continue
        err = np.abs(ev[: ref.size] - ref)
        rows.append((float(fe), err))
    rows.sort(key=lambda t: t[0])
    return rows


def plot_featom_panels(
    radius_rows: list[tuple[int, float, np.ndarray]],
    fe_rows: list[tuple[float, np.ndarray]],
    out_path: Path,
    *,
    fixed_r_bohr: int,
) -> None:
    x = np.arange(1, N_REF + 1, dtype=int)
    colors = _COLORS
    fig, (ax_r, ax_fe) = plt.subplots(
        1,
        2,
        figsize=(_FIG_WIDTH_IN, _FIG_HEIGHT_IN),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.0, 1.0], "wspace": 0.09},
    )

    radius_single_fe = len({fe for fe, _, _ in radius_rows}) <= 1
    for j, (fe, r_bohr, err) in enumerate(radius_rows):
        c = colors[j % len(colors)]
        if radius_single_fe:
            label = rf"$R={int(r_bohr)}$"
        else:
            label = rf"$R={int(r_bohr)},\; N_{{\mathrm{{fe}}}}={int(fe)}$"
        ax_r.semilogy(
            x,
            np.maximum(err, 1e-20),
            marker="o",
            ms=3.0,
            lw=1.8,
            label=label,
            color=c,
        )

    for j, (n_fe, err) in enumerate(fe_rows):
        c = colors[j % len(colors)]
        label = rf"$R={int(fixed_r_bohr)},\; N_{{\mathrm{{fe}}}}={int(n_fe)}$"
        ax_fe.semilogy(
            x,
            np.maximum(err, 1e-20),
            marker="o",
            ms=3.0,
            lw=1.8,
            label=label,
            color=c,
        )

    y_label = r"$|\varepsilon_i-\varepsilon_i^{\mathrm{ref}}|$ (Ha)"
    x_label = r"Eigenvalue index $i$"
    ax_r.set_xlabel(x_label, fontsize=X_AXIS_LABEL_FONTSIZE)
    ax_r.set_ylabel(y_label, fontsize=AXIS_LABEL_FONTSIZE)
    ax_r.set_xticks(x)
    ax_fe.set_xlabel(x_label, fontsize=X_AXIS_LABEL_FONTSIZE)
    ax_fe.set_ylabel(y_label, fontsize=AXIS_LABEL_FONTSIZE)
    ax_fe.set_xticks(x)

    _finalize_paper_semilogy_axes(ax_r, x_mode="fe")
    _finalize_paper_semilogy_axes(ax_fe, x_mode="fe")
    ax_r.tick_params(axis="y", which="minor", labelleft=False)
    ax_fe.tick_params(axis="y", which="minor", labelleft=False)
    ax_r.legend(
        fontsize=LEGEND_FONTSIZE,
        ncol=1,
        loc="upper right",
        labelspacing=0.35,
        handlelength=1.6,
        borderpad=0.35,
    )
    ax_fe.legend(
        fontsize=LEGEND_FONTSIZE,
        ncol=1,
        loc="upper right",
        labelspacing=0.35,
        handlelength=1.6,
        borderpad=0.35,
    )

    out_path = out_path.resolve()
    if out_path.suffix.lower() != ".pdf":
        out_path = out_path.with_suffix(".pdf")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="pdf", bbox_inches="tight", dpi=600)
    plt.close(fig)


def main() -> None:
    root = _DEFAULT_LDA_SVWN_ROOT.resolve()
    out_path = _DEFAULT_OUT_PDF.resolve()

    if ref_eigenvalues.size != N_REF:
        print("ref_eigenvalues length mismatch", file=sys.stderr)
        sys.exit(1)

    r_sweep = root / "domain_radius_sweep"
    fe_sweep = root / "finite_element_sweep"

    radius_rows = collect_radius_curves(
        r_sweep,
        _DEFAULT_CONFIGURATION,
        ref_eigenvalues,
    )
    fe_rows = collect_fe_curves(
        fe_sweep,
        _DEFAULT_CONFIGURATION,
        _DEFAULT_FIXED_R,
        ref_eigenvalues,
    )

    if not radius_rows:
        print(
            f"No domain-radius sweep curves (check summaries under {r_sweep}).",
            file=sys.stderr,
        )
        sys.exit(1)
    if not fe_rows:
        print(
            f"No FE curves (check {fe_sweep} and R={_DEFAULT_FIXED_R}).",
            file=sys.stderr,
        )
        sys.exit(1)

    plot_featom_panels(radius_rows, fe_rows, out_path, fixed_r_bohr=_DEFAULT_FIXED_R)
    print(
        f"Wrote {out_path} ({len(radius_rows)} domain-radius curves, {len(fe_rows)} FE-curves)."
    )


if __name__ == "__main__":
    main()

"""Unit tests for Pulay mixing on a two-dimensional nonlinear map.

The goal is to validate Pulay mixing on a small nonlinear map: from each
listed initial point, the accelerated iteration should converge to a fixed
point of G (residual below tolerance).

State variable:
    u = [u1, u2], where x = u1 - 1 and y = u2 - 1.

Map definition:
    x_next = x - eta * (x^3 - x + beta * y)
    y_next = y - eta * (y^3 - y + beta * x)
    u_next = [x_next + 1, y_next + 1]

Pulay mixing is applied to the iteration u_{k+1} = G(u_k) via:
    u_mixed = Mixer.mix(u_in, u_out)

How to run:
- Run only this file with pytest:
  python -m pytest tests/unit/test_pulay_mixing_trajectories.py -q

- Run script mode for trajectory visualization:
  python tests/unit/test_pulay_mixing_trajectories.py
"""

from pathlib import Path
import sys

import numpy as np
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.scf.mixer import Mixer


ETA = 0.15
BETA = 0.20

# Pulay / linear knobs (kept in one place for tests and figure caption).
PULAY_MIXING_PARAMETER = 0.35
PULAY_MIXING_HISTORY = 1
PULAY_MIXING_FREQUENCY = 3
LINEAR_MIXING_ALPHA1 = 0.55
LINEAR_MIXING_ALPHA2 = 0.98

initial_points = [
    np.array([2.2, 2.0]),
    np.array([0.2, 0.4]),
    np.array([1.9, 0.8]),
    np.array([0.6, 1.6]),
]


def fixed_point_map(u: np.ndarray, eta: float = ETA, beta: float = BETA) -> np.ndarray:
    """Two-dimensional nonlinear map in shifted coordinates."""
    x, y = u - 1.0
    x_next = x - eta * (x**3 - x + beta * y)
    y_next = y - eta * (y**3 - y + beta * x)
    return np.array([x_next + 1.0, y_next + 1.0], dtype=float)


def residual_norm(u: np.ndarray) -> float:
    """Residual norm ||G(u) - u||_2."""
    return float(np.linalg.norm(fixed_point_map(u) - u))


def run_pulay_iteration(
    initial_u: np.ndarray,
    max_iter: int = 200,
    tol: float = 1e-9,
) -> tuple[np.ndarray, np.ndarray, bool]:
    """Run Pulay-accelerated iteration and return trajectory/residuals."""
    mixer = Mixer(
        use_pulay_mixing       = True,
        use_preconditioner     = False,
        pulay_mixing_parameter = PULAY_MIXING_PARAMETER,
        pulay_mixing_history   = PULAY_MIXING_HISTORY,
        pulay_mixing_frequency = PULAY_MIXING_FREQUENCY,
        linear_mixing_alpha1   = LINEAR_MIXING_ALPHA1,
        linear_mixing_alpha2   = LINEAR_MIXING_ALPHA2,
    )
    mixer.set_use_rho_clamp(False)

    u = np.asarray(initial_u, dtype=float)
    trajectory = [u.copy()]
    residuals = []

    for _ in range(max_iter):
        u_out = fixed_point_map(u)
        r = float(np.linalg.norm(u_out - u))
        residuals.append(r)
        if r < tol:
            return np.asarray(trajectory), np.asarray(residuals), True

        u = mixer.mix(u, u_out, None)
        trajectory.append(u.copy())

    return np.asarray(trajectory), np.asarray(residuals), residual_norm(u) < tol


@pytest.mark.parametrize(
    "initial_u",
    initial_points,
    ids=[str(i) for i in range(len(initial_points))],
)
def test_pulay_mixing_converges_from_initial(initial_u: np.ndarray) -> None:
    """Pulay iteration reaches fixed point of G from this initial u (residual check)."""
    _, residuals, converged = run_pulay_iteration(initial_u)
    assert converged, f"Pulay did not converge from initial point {initial_u}."
    assert residuals[-1] < 1e-8, (
        f"Final residual too large from initial point {initial_u}: {residuals[-1]:.3e}"
    )


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # Run all trajectories first, then choose plotting range to include all points.
    trajectory_records = []
    for u0 in initial_points:
        trajectory, residuals, converged = run_pulay_iteration(u0)
        trajectory_records.append((u0, trajectory, residuals, converged))

    all_points = np.vstack([traj for _, traj, _, _ in trajectory_records])
    u1_min, u1_max = float(np.min(all_points[:, 0])), float(np.max(all_points[:, 0]))
    u2_min, u2_max = float(np.min(all_points[:, 1])), float(np.max(all_points[:, 1]))
    pad_u1 = max(0.1, 0.08 * (u1_max - u1_min))
    pad_u2 = max(0.1, 0.08 * (u2_max - u2_min))
    u1 = np.linspace(u1_min - pad_u1, u1_max + pad_u1, 220)
    u2 = np.linspace(u2_min - pad_u2, u2_max + pad_u2, 220)
    U1, U2 = np.meshgrid(u1, u2)
    R = np.zeros_like(U1)
    for i in range(U1.shape[0]):
        for j in range(U1.shape[1]):
            R[i, j] = residual_norm(np.array([U1[i, j], U2[i, j]]))

    fig, ax = plt.subplots(figsize=(8, 6))
    contour = ax.contourf(U1, U2, np.log10(R + 1e-16), levels=40, cmap="viridis", alpha=0.75)
    cbar = fig.colorbar(contour, ax=ax, pad=0.02)
    cbar.set_label(r"$\log_{10}\|G(u)-u\|_2$")

    # Plot Pulay trajectories.
    from matplotlib.lines import Line2D

    colors = ["C1", "C3", "C0", "C2"]
    for idx, (u0, trajectory, residuals, converged) in enumerate(trajectory_records):
        ax.plot(
            trajectory[:, 0],
            trajectory[:, 1],
            "-",
            lw=1.4,
            color=colors[idx % len(colors)],
        )
        # Show intermediate iterates explicitly (excluding start/end markers).
        if trajectory.shape[0] > 2:
            ax.scatter(
                trajectory[1:-1, 0],
                trajectory[1:-1, 1],
                s=14,
                marker="o",
                color=colors[idx % len(colors)],
                alpha=0.55,
                linewidths=0,
            )
        # Mark start and end points explicitly so convergence direction is unambiguous.
        ax.scatter(
            trajectory[0, 0],
            trajectory[0, 1],
            s=60,
            marker="o",
            facecolors="none",
            edgecolors=colors[idx % len(colors)],
            linewidths=1.8,
        )
        ax.scatter(
            trajectory[-1, 0],
            trajectory[-1, 1],
            s=75,
            marker="*",
            color=colors[idx % len(colors)],
            edgecolor="k",
        )
        # Initial coordinates on the figure (two decimal places), not in the legend.
        # Bottom-right of label near the start point; small offset so it does not touch the marker.
        x0, y0 = float(trajectory[0, 0]), float(trajectory[0, 1])
        ax.annotate(
            f"({u0[0]:.2f}, {u0[1]:.2f})",
            xy=(x0, y0),
            xytext=(5, 5),
            textcoords="offset points",
            ha="right",
            va="bottom",
            fontsize=8,
            color=colors[idx % len(colors)],
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.6", alpha=0.85),
        )
        print(f"init {idx+1}: converged={converged}, steps={len(residuals)}, final_residual={residuals[-1]:.3e}")

    legend_handles = [
        Line2D(
            [0],
            [0],
            linestyle="None",
            marker="o",
            markerfacecolor="none",
            markeredgecolor="0.35",
            markersize=9,
            markeredgewidth=1.8,
            label="Start",
        ),
        Line2D(
            [0],
            [0],
            linestyle="None",
            marker="*",
            markerfacecolor="C4",
            markeredgecolor="k",
            markersize=11,
            label="End",
        ),
        Line2D(
            [0],
            [0],
            linestyle="None",
            marker="o",
            markerfacecolor="0.45",
            markeredgecolor="none",
            markersize=5,
            alpha=0.55,
            label="Intermediate",
        ),
    ]
    ax.legend(handles=legend_handles, loc="upper left", fontsize=8)

    ax.set_xlabel("u1")
    ax.set_ylabel("u2")
    ax.set_title("Pulay mixing trajectories on a 2D nonlinear map")
    ax.set_xlim(u1_min - pad_u1, u1_max + pad_u1)
    ax.set_ylim(u2_min - pad_u2, u2_max + pad_u2)
    ax.grid(alpha=0.25)

    # Matplotlib mathtext here does not support amsmath environments; use stacked math lines.
    # Map G only in u; parameters on the figure are those appearing in G (not Pulay knobs).
    iteration_tex = "\n".join(
        [
            r"$u_{\mathrm{out},1}^{(k)}=u_1^{(k)}-\eta\left(\left(u_1^{(k)}-1\right)^3-\left(u_1^{(k)}-1\right)"
            r"+\beta\left(u_2^{(k)}-1\right)\right)$",
            r"$u_{\mathrm{out},2}^{(k)}=u_2^{(k)}-\eta\left(\left(u_2^{(k)}-1\right)^3-\left(u_2^{(k)}-1\right)"
            r"+\beta\left(u_1^{(k)}-1\right)\right)$",
            r"$u^{(k+1)}=\mathrm{Mixer}\!\left(u^{(k)},\,u_{\mathrm{out}}^{(k)}\right), \quad $"
            rf"$\eta={ETA:.2f},\;\beta={BETA:.2f}$",
        ]
    )
    ax.text(
        0.98,
        0.04,
        iteration_tex,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        multialignment="left",
        fontsize=8,
        zorder=25,
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.45", alpha=0.92),
    )

    fig.tight_layout()

    output_path = Path(__file__).resolve().parent / "check_pulay_mixing_trajectories.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=220)
   
    plt.show()
    plt.close(fig)

"""Sim 8 — Phase diagram traversal (Landau f0).

What this script does
---------------------
Switches to a Landau-style temperature-dependent free energy
    f0(c, T) = a (T - Tc) (c - 0.5)^2 + b (c - 0.5)^4
so that "temperature" T becomes a knob. Runs five small (128^2,
5000-step) simulations at hand-picked (T, c0) points sampling the
qualitative regions of the regular-solution phase diagram:
  1. (0.7 Tc, 0.5)  — bicontinuous spinodal
  2. (0.7 Tc, 0.3)  — droplet-pattern spinodal
  3. (0.7 Tc, 0.15) — metastable region (deterministic CH does nothing)
  4. (0.7 Tc, 0.05) — outside binodal (homogeneous solid solution)
  5. (0.95 Tc, 0.5) — near-critical (very long wavelength)

Then assembles ONE composite figure with the analytical binodal and
spinodal curves drawn from the Landau formula, with the late-time
composition-map thumbnails inset at their (T, c0) coordinates. Visually
unifies the thermodynamic phase diagram with the kinetic morphologies
that the same model produces.

Outputs
-------
- data/sim8_T{T}_c{c0}.npz   : final snapshots from each point.
- figures/sim8_phase_diagram.png : the composite figure.

Knobs
-----
- run_one(T, c0, kappa, a, b, N, n_steps, snapshot_every, dt, seed)
    - a, b      : Landau coefficients (default 1.0 each).
    - N=128     : grid size for the small thumbnails.
    - n_steps   : 5000 default (enough to reach a recognisable late state
                  even close to T_c where dynamics are slow).
    - kappa, dt : same defaults as the production CH runs.
- The (T, c0) sweep list itself is hard-coded inside main() — edit there
  to add/remove points or shift to other temperatures.
- Stabilization for the Landau f0 is computed automatically from a, b,
  T-Tc inside run_one (max-curvature heuristic).
"""
from __future__ import annotations

import os
from functools import partial

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from solver import (CahnHilliardSolver, mu_landau, f0pp_landau)


def landau_binodal(T_over_Tc, a=1.0, b=1.0, Tc=1.0):
    """Binodal of a (T - Tc)(c-0.5)^2 + b(c-0.5)^4: minima at c-0.5 = +- sqrt(a(Tc-T)/(2b))."""
    if T_over_Tc >= 1.0:
        return None
    half_width = np.sqrt(a * (Tc - T_over_Tc * Tc) / (2.0 * b))
    return 0.5 - half_width, 0.5 + half_width


def landau_spinodal(T_over_Tc, a=1.0, b=1.0, Tc=1.0):
    """Spinodal: f0'' = 0 -> (c-0.5)^2 = a(Tc-T)/(6b)."""
    if T_over_Tc >= 1.0:
        return None
    half = np.sqrt(a * (Tc - T_over_Tc * Tc) / (6.0 * b))
    return 0.5 - half, 0.5 + half


def run_one(T, c0, kappa=5e-4, a=1.0, b=1.0, N=256, n_steps=5000,
            snapshot_every=500, dt=0.01, seed=42):
    mu = partial(mu_landau, T=T, Tc=1.0, a=a, b=b)
    # max |f0''| over [0, 1]: f0''(0) = 2a(T-Tc) + 12b*0.25 = 2a(T-Tc) + 3b.
    # For T<Tc, the well-bottom |f0''| ~ 2|a|(Tc-T), and at c=0,1 it can be larger.
    A_stab = max(3.0, 2.0 * a * abs(T - 1.0) + 3.0 * b)
    s = CahnHilliardSolver(N=N, kappa=kappa, dt=dt, mu_func=mu, stabilization=A_stab)
    rng = np.random.default_rng(seed)
    field = rng.normal(loc=c0, scale=0.01, size=(N, N))
    field += c0 - field.mean()
    s.set_field(field)
    times, snaps = s.run(n_steps, snapshot_every=snapshot_every)
    return snaps[-1]


def main():
    os.makedirs("data", exist_ok=True)
    points = [
        (0.7, 0.5, "bicontinuous"),
        (0.7, 0.3, "droplets"),
        (0.7, 0.15, "metastable"),
        (0.7, 0.05, "stable"),
        (0.95, 0.5, "near-critical"),
    ]
    snaps_at_points = {}
    for T, c0, label in points:
        print(f"[sim8] running T={T} Tc, c0={c0}  ({label})")
        snap = run_one(T, c0)
        snaps_at_points[(T, c0)] = snap
        np.savez_compressed(f"data/sim8_T{T}_c{c0}.npz", snap=snap)

    # Assemble composite phase diagram
    fig, ax = plt.subplots(figsize=(9.5, 6))
    Ts = np.linspace(0.5, 1.0, 200)
    bin_lo, bin_hi = [], []
    spin_lo, spin_hi = [], []
    for T in Ts:
        bl = landau_binodal(T)
        sl = landau_spinodal(T)
        if bl is not None:
            bin_lo.append(bl[0])
            bin_hi.append(bl[1])
        else:
            bin_lo.append(np.nan); bin_hi.append(np.nan)
        if sl is not None:
            spin_lo.append(sl[0])
            spin_hi.append(sl[1])
        else:
            spin_lo.append(np.nan); spin_hi.append(np.nan)

    ax.plot(bin_lo, Ts, "g-", lw=1.6, label="binodal")
    ax.plot(bin_hi, Ts, "g-", lw=1.6)
    ax.plot(spin_lo, Ts, "--", color="orange", lw=1.4, label="spinodal")
    ax.plot(spin_hi, Ts, "--", color="orange", lw=1.4)
    ax.fill_betweenx(Ts, bin_lo, spin_lo, color="lightyellow", alpha=0.7,
                     label="metastable")
    ax.fill_betweenx(Ts, spin_hi, bin_hi, color="lightyellow", alpha=0.7)
    ax.fill_betweenx(Ts, spin_lo, spin_hi, color="lightcoral", alpha=0.4,
                     label="unstable (spinodal)")

    ax.set_xlim(0, 1)
    ax.set_ylim(0.5, 1.05)
    ax.set_xlabel("c0")
    ax.set_ylabel("T / Tc")
    ax.set_title("Sim 8: phase diagram traversal")
    ax.legend(loc="upper left", fontsize=9)

    # Insert thumbnails as inset axes anchored at each (c0, T) point
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    for (T, c0, label) in points:
        snap = snaps_at_points[(T, c0)]
        # Place inset just to the right of the data point in axis coords
        x_disp, y_disp = ax.transData.transform((c0, T))
        x_axes, y_axes = ax.transAxes.inverted().transform((x_disp, y_disp))
        # Offset the inset slightly so it doesn't sit on top of the marker
        dx, dy = 0.06, 0.06
        bbox = (x_axes + dx, y_axes - 0.10, 0.18, 0.18)
        try:
            ax_in = ax.inset_axes(bbox)
        except Exception:
            ax_in = inset_axes(ax, width="15%", height="15%", loc="lower left",
                               bbox_to_anchor=(x_axes + dx, y_axes - 0.05, 1, 1),
                               bbox_transform=ax.transAxes)
        ax_in.imshow(snap.T, origin="lower", cmap="coolwarm", vmin=0, vmax=1)
        ax_in.set_xticks([]); ax_in.set_yticks([])
        ax_in.set_title(f"({c0}, {T})", fontsize=7)
        ax.plot([c0], [T], "ko", ms=5)
        ax.annotate(label, xy=(c0, T), xytext=(c0 + 0.01, T + 0.005),
                    fontsize=7, color="0.2")

    fig.tight_layout()
    out = "figures/sim8_phase_diagram.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[sim8] wrote {out}")


if __name__ == "__main__":
    main()

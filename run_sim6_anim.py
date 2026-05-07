"""Sim 6 (animated) — Free energy descent visualised alongside c(r, t).

What this script does
---------------------
Reads Sim 1's snapshots, computes F(t) on every frame, then renders a
two-panel MP4:
  - Left  : composition map c(r, t) — the familiar spatial pattern.
  - Right : F(t) curve drawn frame-by-frame, so the audience literally
            watches the free energy descend as the morphology forms and
            coarsens. A dashed horizontal line marks the asymptotic
            equilibrium value F(∞).

This is the animated companion to the static F(t) log-log plot in
run_sim6.py.

Knobs
-----
- npz_path : input .npz from Sim 1 (default 'data/sim1.npz').
- out      : MP4 output path.
- fps      : default 15. Reduce to 8 for a slower descent.
- bitrate  : default 2000.
- kappa    : must match the kappa used in Sim 1 (default 5e-4).
"""
from __future__ import annotations

import os
import numpy as np

from analysis import free_energy_series
import viz


def main(npz_path="data/sim1.npz", out="videos/sim6_F_descent.mp4",
         fps=15, bitrate=2000, kappa=5e-4):
    if not os.path.exists(npz_path):
        raise FileNotFoundError(f"run Sim 1 first — missing {npz_path}")

    d = np.load(npz_path)
    times = d["times"]
    snaps = d["snaps"]
    L = float(d["L"]) if "L" in d.files else 1.0

    print(f"[sim6_anim] computing F on {len(snaps)} frames  (L={L})...")
    F = free_energy_series(snaps, L=L, kappa=kappa)
    print(f"  F range: {F[0]:.5f} -> {F[-1]:.5f}")

    viz.animate_field_and_F(snaps, times, F, out, fps=fps, bitrate=bitrate,
                             title="Cahn-Hilliard: free energy descent")
    print(f"[sim6_anim] wrote {out}")


if __name__ == "__main__":
    main()

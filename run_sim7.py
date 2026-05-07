"""Sim 7 — G(c) + composition map coevolution (the killer demo).

What this script does
---------------------
Pure postprocessing animation. Reads Sim 1's snapshots and renders an
MP4 with two synchronised panels:
  - Left  : the f0(c) curve, with binodal points (c=0, 1) and spinodal
            points (~0.211, ~0.789) marked, and a live histogram bar
            overlay showing where the system's compositions currently
            sit on the curve.
  - Right : the composition map c(r, t).

This is the central pedagogical visualisation for the recorded
presentation: the audience watches the histogram split from a single
peak at 0.5 into two peaks at the binodal endpoints {0, 1}, *while
simultaneously* watching the spatial pattern emerge. It directly
visualises the coupling between thermodynamic descent and pattern
formation that the linearised theory only describes algebraically.

Knobs
-----
- npz_path : Sim 1 input .npz (default 'data/sim1.npz').
- out      : MP4 output path.

For deeper customisation (color, layout, panel sizes, histogram bin
count, point-density overlay) edit viz.animate_Gc_and_map directly.
"""
from __future__ import annotations

import os
import numpy as np

import viz
from solver import f0_polynomial


def main(npz_path="data/sim1.npz", out="videos/sim7_Gc_demo.mp4"):
    if not os.path.exists(npz_path):
        raise FileNotFoundError(f"need to run Sim 1 first: missing {npz_path}")
    d = np.load(npz_path)
    times = d["times"]
    snaps = d["snaps"]
    c0 = float(d["c0"]) if "c0" in d.files else 0.5

    # Spinodal points for f0 = c^2 (1-c)^2: f0''(c) = 2 - 12c + 12c^2 = 0
    # -> c = (12 +- sqrt(144 - 96))/24 = (12 +- sqrt(48))/24 = 0.5 +- sqrt(48)/24
    # = 0.5 +- 1/(2*sqrt(3)) ≈ 0.211, 0.789
    spinodal = (0.5 - 1.0 / (2.0 * np.sqrt(3.0)),
                0.5 + 1.0 / (2.0 * np.sqrt(3.0)))

    viz.animate_Gc_and_map(
        snaps, times, out,
        f0_func=f0_polynomial,
        c0=c0,
        binodal=(0.0, 1.0),
        spinodal=spinodal,
        fps=15,
        title="Cahn-Hilliard: free energy descent + spatial pattern",
    )
    print(f"[sim7] wrote {out}")


if __name__ == "__main__":
    main()

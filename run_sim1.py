"""Sim 1 — Symmetric morphology evolution (canonical bicontinuous case).

What this script does
---------------------
Runs the deterministic Cahn-Hilliard equation from a random initial
condition centred on c0 = 0.5 (the symmetric / critical composition)
and produces:
  - A composition-map MP4 animation showing the noise-pattern -> worm /
    labyrinth morphology -> coarsening progression.
  - A static 4-panel snapshot figure at early / mid / late times.
  - A histogram-vs-time heatmap (composition distribution drifting from
    a single peak at 0.5 to a bimodal {0, 1} distribution).
  - A compressed .npz archive of times and snapshots, consumed by Sims 6
    (free energy decay) and Sim 7 (G(c) coevolution demo).

This run is the workhorse: Sim 2, Sim 3 reuse main() with different c0.

Knobs (main)
------------
- N             : grid (default 1024). Use 128/256 for fast iteration.
- L             : box side length (default 4.0). With N=1024 this keeps
                  dx = L/N = 0.004 (same as the original N=256, L=1 runs)
                  while fitting ~14 dominant wavelengths per side (~196 domains).
- n_steps       : total integration steps (default 5000; t_max = n_steps*dt).
- snapshot_every: how often to record (default 50 -> 100 frames @ 15 fps =
                  ~7 s animation; halved from 25 to limit .npz memory at N=1024).
- seed          : RNG seed for the Gaussian IC.
- c0            : background composition. 0.5 is the symmetric case.
- out_video / out_panels / out_hist / out_npz: output paths.
"""
from __future__ import annotations

import os
import numpy as np

from solver import CahnHilliardSolver
from analysis import histogram_evolution
import viz


def main(N=1440, L=6.0, n_steps=7500, snapshot_every=75, seed=1, c0=0.5,
         out_video="videos/sim1_morphology.mp4",
         out_panels="figures/sim1_panels.png",
         out_hist="figures/sim1_histogram.png",
         out_npz="data/sim1.npz"):
    s = CahnHilliardSolver(N=N, L=L, kappa=5e-4, M=1.0, dt=0.01)
    s.set_random_ic(c0=c0, amp=0.01, seed=seed)
    print(f"[sim1] N={N}  L={L}  c0={c0}  steps={n_steps}  snapshot_every={snapshot_every}")
    times, snaps = s.run(n_steps, snapshot_every=snapshot_every, verbose=False)
    print(f"  collected {len(times)} frames; mean drift "
          f"{snaps[-1].mean() - snaps[0].mean():+.2e}")

    os.makedirs("data", exist_ok=True)
    np.savez_compressed(out_npz, times=times, snaps=snaps, c0=c0, L=L)

    # Composition histogram heatmap
    centers, hist_t = histogram_evolution(snaps, bins=80, range_=(-0.1, 1.1))
    viz.histogram_heatmap(centers, hist_t, times, out_hist,
                          title=f"Sim 1: histogram(c, t),  c0={c0}")

    # Snapshot panels
    panel_indices = _pick_indices(times, [0.0, 2.0, 10.0, times[-1]])
    viz.snapshot_panels(snaps, times, panel_indices, out_panels,
                        title=f"Sim 1: morphology evolution, c0={c0}")

    # Animation
    viz.animate_field(snaps, times, out_video, fps=15,
                      title=f"Sim 1 (c0={c0})")
    print(f"[sim1] wrote {out_video}, {out_panels}, {out_hist}, {out_npz}")


def _pick_indices(times, target_times):
    """Pick frame indices closest to each requested wall time."""
    return [int(np.argmin(np.abs(times - tt))) for tt in target_times]


if __name__ == "__main__":
    main()

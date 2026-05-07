"""Sim 9 (optional) — Cahn-Hilliard-Cook noise vs deterministic comparison.

What this script does
---------------------
Sets c0 = 0.15, which sits in the *metastable* region for
f0 = c^2(1-c)^2 (outside the spinodal endpoint ~0.211 but inside the
binodal at c = 0). Runs two simulations from the same Gaussian-noise
initial condition:

  A) Deterministic CH (noise_amp = 0). The system is metastable and
     the deterministic dynamics has no driving force away from the
     uniform state — nothing happens. This is Cahn (1961)'s well-known
     limitation: CH cannot describe nucleation in the metastable region.

  B) CHC (Cahn-Hilliard-Cook) with conserved Gaussian noise, amplitude
     tuned to give occasional nucleation events within the integration
     window. Nuclei appear, grow into droplets, then coarsen.

Renders the two as a side-by-side MP4 so the contrast is immediate.

Outputs
-------
- data/sim9.npz             : both sets of snapshots.
- videos/sim9_chc_vs_det.mp4 : the comparison animation.

Knobs
-----
- N            : 128 (smaller than the other production runs because
                 we need only enough resolution to resolve nucleation
                 events).
- n_steps      : 8000 (long enough to see ~one or more nucleation
                 events with the default noise amplitude).
- snapshot_every: 40.
- c0           : 0.15 by default. Move it deeper into the metastable
                 region (closer to 0) to need *larger* noise for
                 nucleation; move it toward the spinodal (~0.21) and
                 even very small noise nucleates immediately.
- noise_amp    : 0.5 by default. Controls how often nucleation events
                 happen. Too small -> deterministic-looking; too large
                 -> noise dominates the picture.
"""
from __future__ import annotations

import os
import numpy as np

from solver import CahnHilliardSolver
import viz


def main(N=256, n_steps=12000, snapshot_every=60, c0=0.20, noise_amp=1.5,
         out="videos/sim9_chc_vs_det.mp4"):
    # c0 = 0.20 is just outside the spinodal endpoint (~0.211 for the
    # polynomial f0), giving a small nucleation barrier so the Cook noise
    # can drive nucleation clearly within the simulation window.
    # Both runs start from exactly the same uniform field; the only
    # difference is whether the stochastic current is turned on.
    uniform_ic = np.full((N, N), c0, dtype=float)
    rng = np.random.default_rng(7)
    s_det = CahnHilliardSolver(N=N, kappa=5e-4, dt=0.01, noise_amp=0.0)
    s_chc = CahnHilliardSolver(N=N, kappa=5e-4, dt=0.01, noise_amp=noise_amp,
                               rng=rng)
    s_det.set_field(uniform_ic)
    s_chc.set_field(uniform_ic)

    print(f"[sim9] N={N}  steps={n_steps}  c0={c0}  noise_amp={noise_amp}")
    times_a, snaps_a = s_det.run(n_steps, snapshot_every=snapshot_every)
    times_b, snaps_b = s_chc.run(n_steps, snapshot_every=snapshot_every)

    os.makedirs("data", exist_ok=True)
    np.savez_compressed("data/sim9.npz", times=times_a, snaps_det=snaps_a,
                        snaps_chc=snaps_b)

    viz.animate_two_fields(snaps_a, snaps_b, times_a, out,
                           title_a=f"Deterministic CH  (c0={c0}) — metastable, frozen",
                           title_b=f"CHC with noise — nucleation events",
                           overall_title="Metastable region: Cook noise enables nucleation")
    print(f"[sim9] wrote {out}")


if __name__ == "__main__":
    main()

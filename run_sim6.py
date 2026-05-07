"""Sim 6 — Free energy decay (postprocesses Sim 1).

What this script does
---------------------
Loads the Sim 1 .npz, computes the total free energy
    F(t) = integral [ f0(c) + kappa |grad c|^2 ] dV
on every saved snapshot using spectral derivatives, and plots:
  - F(t) on linear axes (should be monotonically decreasing for the
    deterministic CH equation — this is the Lyapunov / gradient-flow
    consistency check).
  - F(t) - F_inf on log-log axes with a power-law fit slope. The late
    stage should follow ~ t^{-1/3} (consistent with diffusive coarsening).

This is purely postprocessing — no simulation is run. Sim 1 must have
produced data/sim1.npz first.

Note on monotonicity
--------------------
Tiny F increases at the level of floating-point noise are tolerated by
the stabilized scheme; the script prints how many of the timesteps
recorded a (positive) delta, which is the key sanity number. If that
count exceeds a few percent of frames or the plot shows visible upticks,
reduce dt or increase the stabilization constant in solver.py.

Knobs
-----
- npz_path : input file (default 'data/sim1.npz').
- out      : output figure path.
- The free-energy form is hard-coded via the f0_polynomial default in
  analysis.free_energy_series; pass a different f0 there if you want to
  evaluate F under a different functional (e.g. Landau).
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analysis import free_energy_series


def main(npz_path="data/sim1.npz", out="figures/sim6_F_decay.png"):
    if not os.path.exists(npz_path):
        raise FileNotFoundError(f"need to run Sim 1 first: missing {npz_path}")
    d = np.load(npz_path)
    times = d["times"]
    snaps = d["snaps"]
    L = float(d["L"]) if "L" in d.files else 1.0

    F = free_energy_series(snaps, L=L, kappa=5e-4)
    print(f"[sim6] F(0) = {F[0]:.5f}, F(end) = {F[-1]:.5f}")
    n_increases = np.sum(np.diff(F) > 0)
    print(f"[sim6] timesteps with F increase: {n_increases} / {len(F) - 1}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].plot(times, F, "o-", ms=3)
    axes[0].set_xlabel("t")
    axes[0].set_ylabel("F(t)")
    axes[0].set_title("Free energy decay")
    # F - F_inf approx using last value
    F_inf = F[-1]
    diff = F - F_inf
    mask = (times > 0) & (diff > 0)
    axes[1].loglog(times[mask], diff[mask], "o-", ms=3, label=r"$F(t) - F_\infty$")
    if mask.sum() > 4:
        late = mask.copy()
        late[: int(0.5 * len(times))] = False
        if late.sum() > 3:
            slope, intercept = np.polyfit(np.log(times[late]),
                                          np.log(diff[late]), 1)
            t_ref = times[late]
            axes[1].loglog(t_ref, np.exp(intercept) * t_ref ** slope, "--",
                           label=f"fit slope = {slope:+.3f}")
    axes[1].set_xlabel("t")
    axes[1].set_ylabel(r"$F(t) - F_\infty$")
    axes[1].set_title("F approach to equilibrium (log-log)")
    axes[1].legend(fontsize=9)
    fig.tight_layout()
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[sim6] wrote {out}")


if __name__ == "__main__":
    main()

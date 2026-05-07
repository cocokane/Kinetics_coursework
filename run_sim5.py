"""Sim 5 — Structure factor evolution from random IC.

What this script does
---------------------
Runs the same kind of random-IC CH simulation as Sim 1 (different RNG
seed) but the analysis is in Fourier space. At every saved snapshot we
compute the radially-averaged structure factor
    S(k) = <|c_hat(k)|^2>_{|k|=k}
and track its peak position k_peak(t). This gives:
  - An independent estimate of beta_max (without needing a controlled
    cosine IC like Sim 4).
  - The coarsening exponent: in the late stage k_peak ~ t^{-1/3} (LSW
    scaling), and we plot log-log with a t^{-1/3} reference line.

Outputs
-------
- figures/sim5_Skt_heatmap.png : log-color heatmap of S(k, t).
- figures/sim5_kpeak.png       : k_peak vs t with t^{-1/3} reference and
                                 a power-law fit slope.
- data/sim5.npz                : times, snapshots, S(k, t), k_peak(t).

Knobs
-----
- N             : default 256 (production); 128 for fast iteration.
- n_steps       : 8000 by default (long enough to see clear coarsening).
- snapshot_every: 40 -> 200 frames; matches the cadence used in Sim 1.
- seed          : RNG seed for IC.
- The fit window for the coarsening exponent is hard-coded to "second
  half of the time series" inside main(); edit `late` mask if you need
  to change it.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from solver import CahnHilliardSolver
from analysis import radial_structure_factor, peak_k, beta_critical
import viz


def main(N=1440, L=6.0, n_steps=12000, snapshot_every=120, seed=2):
    s = CahnHilliardSolver(N=N, L=L, kappa=5e-4, M=1.0, dt=0.01)
    s.set_random_ic(c0=0.5, amp=0.01, seed=seed)
    print(f"[sim5] N={N}  L={L}  steps={n_steps}")
    times, snaps = s.run(n_steps, snapshot_every=snapshot_every)
    print(f"  collected {len(times)} frames")

    # Compute S(k, t)
    n_bins = N // 2
    S_t = np.empty((len(snaps), n_bins))
    k_centers = None
    for i, sn in enumerate(snaps):
        k_centers, S = radial_structure_factor(sn, L=L, n_bins=n_bins)
        S_t[i] = S

    # Peak tracking — k_min = 2pi/L (one full wave across the box)
    k_min_use = 2.0 * np.pi / L
    k_peak_t = np.array([peak_k(k_centers, S_t[i], k_min=k_min_use)
                         for i in range(len(snaps))])

    os.makedirs("data", exist_ok=True)
    np.savez_compressed("data/sim5.npz", times=times, snaps=snaps, k_centers=k_centers,
                        S_t=S_t, k_peak_t=k_peak_t, L=L)

    # Heatmap of S(k, t)
    # Skip t = 0 to avoid log(0)
    viz.Skt_heatmap(k_centers[1:], S_t[1:, 1:], times[1:],
                    "figures/sim5_Skt_heatmap.png",
                    title="Sim 5: S(k, t)")

    # k_peak vs t with t^{-1/3} reference
    fig, ax = plt.subplots(figsize=(6, 4.5))
    mask = (times > 0) & (k_peak_t > 0)
    ax.loglog(times[mask], k_peak_t[mask], "o-", ms=4, label="measured")
    if mask.sum() > 5:
        # Fit power law over late half of data
        late = mask.copy()
        late[: int(0.5 * len(times))] = False
        if late.sum() > 3:
            slope, intercept = np.polyfit(np.log(times[late]), np.log(k_peak_t[late]), 1)
            t_ref = times[mask]
            ax.loglog(t_ref, np.exp(intercept) * t_ref ** slope, "--",
                      label=fr"fit: slope = {slope:+.3f}")
    bc = beta_critical(0.5, kappa=5e-4)
    bmax = bc / np.sqrt(2)
    ax.axhline(bmax, ls=":", color="orange",
               label=fr"$\beta_{{\max}} = {bmax:.1f}$ (linear theory)")
    # t^{-1/3} reference line anchored on first late point
    if mask.sum() > 5:
        t0 = times[mask][len(times[mask]) // 2]
        k0 = k_peak_t[mask][len(times[mask]) // 2]
        t_ref = times[mask][len(times[mask]) // 2:]
        ax.loglog(t_ref, k0 * (t_ref / t0) ** (-1.0 / 3.0),
                  "k--", alpha=0.5, label=r"$t^{-1/3}$")
    ax.set_xlabel("t")
    ax.set_ylabel("$k_{\\rm peak}$")
    ax.set_title("Sim 5: peak wavenumber vs time")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig("figures/sim5_kpeak.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"[sim5] wrote figures/sim5_Skt_heatmap.png, figures/sim5_kpeak.png, data/sim5.npz")


if __name__ == "__main__":
    main()

"""Sim 3 — Just inside the spinodal endpoint. c0 = 0.22.

What this script does
---------------------
Same dynamics as Sim 1 but with c0 = 0.22, which is barely inside the
spinodal endpoint (~0.211 for f0 = c^2(1-c)^2). Demonstrates "critical
slowing down": the dominant unstable wavelength is much longer (because
beta_c -> 0 as c0 -> spinodal endpoint) and the linear growth rate is
much smaller, so the morphology emerges later and at coarser scale.

After running, also assembles a composite figure
`figures/sim123_comparison.png` placing the late-time snapshots from
Sims 1 / 2 / 3 side-by-side. Sim 1 and Sim 2 must already have been run
(needs `data/sim1.npz` and `data/sim2.npz`).

Knobs
-----
Inherits from run_sim1.main(). Only c0 and the output paths are
overridden. The composite-figure step is unconditional and silently
skips if either dependency .npz is missing.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import run_sim1


def main():
    run_sim1.main(
        c0=0.22,
        out_video="videos/sim3_nearspinodal.mp4",
        out_panels="figures/sim3_panels.png",
        out_hist="figures/sim3_histogram.png",
        out_npz="data/sim3.npz",
    )

    # Composite comparison panel of late-time snapshots from Sims 1, 2, 3.
    paths = ["data/sim1.npz", "data/sim2.npz", "data/sim3.npz"]
    if all(os.path.exists(p) for p in paths):
        fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))
        for ax, p, c0 in zip(axes, paths, [0.5, 0.3, 0.22]):
            d = np.load(p)
            snap = d["snaps"][-1]
            im = ax.imshow(snap.T, origin="lower", cmap="coolwarm", vmin=0, vmax=1,
                           extent=[0, 1, 0, 1])
            ax.set_title(f"c0 = {c0}, t = {d['times'][-1]:.1f}")
            ax.set_xticks([])
            ax.set_yticks([])
        fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02, label="c")
        fig.suptitle("Late-time morphologies vs c0  (Sims 1/2/3)")
        out = "figures/sim123_comparison.png"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"[sim3] wrote {out}")


if __name__ == "__main__":
    main()

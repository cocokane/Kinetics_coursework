"""Sim 2 — Off-symmetric morphology (droplets from spinodal). c0 = 0.3.

What this script does
---------------------
Calls run_sim1.main with c0 = 0.3 (still inside the spinodal, since the
spinodal endpoint of f0 = c^2(1-c)^2 is c ~ 0.211). Same dynamics as Sim 1
otherwise. The point is morphological: at c0 = 0.3 the minority phase
forms isolated droplets rather than the bicontinuous worm pattern of
Sim 1, *despite both runs being driven by the same spinodal mechanism*.

This is the visual answer to "spinodal vs nucleation": morphology depends
on c0, not on whether the mechanism is N&G or spinodal.

Knobs
-----
All knobs are inherited from run_sim1.main(). Edit there if you want to
change grid, dt, total steps, or snapshot cadence; here we only override
c0 and the output filenames.
"""
from __future__ import annotations

import run_sim1


def main():
    run_sim1.main(
        c0=0.3,
        out_video="videos/sim2_droplet.mp4",
        out_panels="figures/sim2_panels.png",
        out_hist="figures/sim2_histogram.png",
        out_npz="data/sim2.npz",
    )


if __name__ == "__main__":
    main()

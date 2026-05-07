"""Sim 4 — Dispersion relation validation (the headline figure).

What this script does
---------------------
For each wavenumber beta in a sweep, run a separate small simulation
seeded with a single-mode cosine perturbation c0 + A0 cos(beta x), and
measure the linear growth/decay rate of that mode. Compare to Cahn's
prediction R(beta) = -M beta^2 [f0''(c0) + 2 kappa beta^2].

Each run is short — it terminates as soon as the seeded mode has either
grown into the nonlinear regime or decayed appreciably. The growth rate
is extracted from the log-amplitude slope over the *leading monotone
window* of the trajectory, which excludes (a) nonlinear saturation of
growing modes and (b) re-excitation of decaying modes from FP-noise
feedback through other modes.

Why the small dt
----------------
For the plan's parameters, R_max ~ 250 / time-unit, which means the
fastest unstable mode has a growth timescale of ~0.004. The production
runs (Sims 1-3, 5-9) use dt = 0.01 with stabilization, which is
unconditionally stable but DAMPS the linear growth rate. To reproduce
the analytical R(beta) curve faithfully we need dt ~ 1e-4 and no
stabilization (safe here because the seeded amplitude stays below
saturation throughout the measurement window).

Output
------
Single figure `figures/sim4_dispersion.png` overlaying measured R(beta)
data points (with std-error bars) on the theoretical curve. beta_c and
beta_max are marked.

Knobs
-----
- run_one(beta, c0, N, n_steps, snapshot_every, dt, kappa, A0)
    - dt = 1e-4         : small enough to resolve R_max ~ 250.
    - A0 = 0.001        : seeded mode amplitude. Smaller means more linear
                          headroom (longer fit window) but earlier FP-floor
                          for decaying modes. Don't go below ~1e-6.
    - n_steps = 30000   : cap on total stepping; runs terminate early on
                          A > 0.05 OR t > 1.0 (in run_one).
- main():
    - c0 = 0.5          : symmetric composition (most spinodal modes
                          unstable).
    - sweep range: n=1..ceil(2 beta_c / fundamental). Edit n_min, n_max
      directly to widen or narrow.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from solver import CahnHilliardSolver, f0pp_polynomial
from analysis import (fit_growth_rate, theoretical_growth_rate, beta_critical)


def _leading_monotone_window(a, lo=1e-9, hi=0.01):
    """Boolean mask: longest contiguous prefix of `a` that is monotone
    (all-increasing or all-decreasing) and lies in [lo, hi]."""
    a = np.asarray(a)
    n = len(a)
    if n < 2:
        return np.zeros(n, dtype=bool)
    sign = 1 if a[1] > a[0] else -1
    mask = np.zeros(n, dtype=bool)
    mask[0] = (lo <= a[0] <= hi)
    for i in range(1, n):
        if not (lo <= a[i] <= hi):
            break
        if sign == 1 and a[i] < a[i - 1]:
            break
        if sign == -1 and a[i] > a[i - 1]:
            break
        mask[i] = True
        # ensure prior point is in the mask
        if not mask[i - 1]:
            mask[i - 1] = (lo <= a[i - 1] <= hi)
    return mask


def _mode_amplitude(c, mode_n):
    """|c_hat[mode_n, 0]| / N — the Fourier amplitude of cos(2 pi mode_n x)."""
    N = c.shape[0]
    # rfft2 along axis 1 first, axis 0 second; for c(x,y) we want the (mode_n, 0) bin
    c_hat = np.fft.fft2(c)
    # cos(2 pi n x) maps to indices (+n, 0) and (-n, 0); use abs of (+n, 0).
    return float(np.abs(c_hat[mode_n, 0]) / (N * N) * 2.0)


def run_one(beta, c0=0.5, N=128, n_steps=30000, snapshot_every=10, dt=1.0e-4,
            kappa=5e-4, A0=0.001):
    """Run a single cosine-IC sim and return (times, amps).

    Beta must be an integer multiple of 2*pi/L. Tracks the Fourier amplitude
    of the seeded mode (works equally well for growing and decaying modes).

    Sim 4 uses small dt (1e-4) and disables stabilization (A=0) because the
    fastest unstable modes for f0=c^2(1-c)^2 with kappa=5e-4 have growth rates
    ~250/time-unit — the default dt=0.01 with A=2 stabilization is far too
    coarse to resolve the linear regime. Amplitudes stay bounded (<~0.05) so
    the unstabilized scheme is safe here.
    """
    L = 1.0
    fundamental = 2.0 * np.pi / L
    mode_n = int(round(beta / fundamental))
    s = CahnHilliardSolver(N=N, L=L, kappa=kappa, dt=dt, stabilization=0.0)
    s.set_cosine_ic(c0=c0, beta=beta, amp=A0)

    times = [0.0]
    amps = [_mode_amplitude(s.c, mode_n)]
    t_max = 1.0  # absolute cap (t units)
    for n in range(1, n_steps + 1):
        s.step()
        if n % snapshot_every == 0:
            t_now = n * dt
            a_now = _mode_amplitude(s.c, mode_n)
            times.append(t_now)
            amps.append(a_now)
            # Stop early once amplitude has clearly grown past linear regime
            if a_now > 0.05:
                break
            # Also stop if we've simulated long enough; for decaying modes
            # the amplitude decays exponentially and floor is set by FP noise
            if t_now > t_max:
                break
    return np.array(times), np.array(amps)


def main():
    c0 = 0.5
    kappa = 5e-4
    M = 1.0
    L = 1.0

    bc = beta_critical(c0, kappa=kappa)  # ~ 31.6
    print(f"[sim4] beta_c = {bc:.3f}, beta_max = {bc / np.sqrt(2):.3f}")

    # Allowed wavenumbers: integer multiples of 2*pi/L. Pick n in a range
    # spanning beta_c/3 to ~2 beta_c.
    fundamental = 2.0 * np.pi / L  # ~ 6.283
    n_min = max(1, int(np.floor(bc / 3.0 / fundamental)))
    n_max = int(np.ceil(2.0 * bc / fundamental))
    n_values = np.arange(n_min, n_max + 1)
    betas = fundamental * n_values
    print(f"[sim4] sweeping n = {list(n_values)} -> betas = {betas.round(2)}")

    measured_R = []
    measured_R_err = []
    for n, beta in zip(n_values, betas):
        t, a = run_one(beta, c0=c0, N=128, kappa=kappa)
        # Linear regime: take the longest *contiguous* leading window where
        # amplitude is monotone (growing or decaying) and within the window
        # [1e-9, 0.01]. This excludes:
        #   - Decay below FP noise floor and subsequent re-excitation by
        #     nonlinear coupling from other modes.
        #   - Nonlinear saturation of growing modes once A approaches 0.01.
        contig = _leading_monotone_window(a, lo=1e-9, hi=0.01)
        if contig.sum() >= 3:
            logA = np.log(a[contig])
            slope, intercept = np.polyfit(t[contig], logA, 1)
            R = float(slope)
            res = logA - (slope * t[contig] + intercept)
            R_err = float(np.std(res) / np.sqrt(max(contig.sum() - 2, 1)) /
                          (np.std(t[contig]) + 1e-12))
        else:
            R, R_err = np.nan, np.nan
        measured_R.append(R)
        measured_R_err.append(R_err)
        print(f"  n={n:3d}  beta={beta:6.2f}  R_meas={R:+.4f}  R_th={theoretical_growth_rate(beta, c0, kappa, M):+.4f}")

    measured_R = np.array(measured_R)
    measured_R_err = np.array(measured_R_err)

    os.makedirs("data", exist_ok=True)
    np.savez("data/sim4.npz", betas=betas, R_meas=measured_R, R_err=measured_R_err,
             c0=c0, kappa=kappa, M=M)

    # Plot
    fig, ax = plt.subplots(figsize=(7.5, 5))
    beta_dense = np.linspace(0.0, 1.05 * betas.max(), 400)
    R_th_dense = theoretical_growth_rate(beta_dense, c0, kappa, M)
    ax.plot(beta_dense, R_th_dense, "k-", lw=1.5,
            label=r"theory: $R(\beta) = -M\beta^2[f_0''(c_0) + 2\kappa\beta^2]$")
    err_use = np.where(np.isfinite(measured_R_err), measured_R_err, 0.0)
    ax.errorbar(betas, measured_R, yerr=err_use, fmt="o", color="crimson",
                ms=6, capsize=3, label="measured")

    ax.axhline(0, color="0.6", lw=0.7)
    ax.axvline(bc, ls="--", color="green", lw=1, label=fr"$\beta_c = {bc:.1f}$")
    ax.axvline(bc / np.sqrt(2), ls=":", color="orange", lw=1,
               label=fr"$\beta_{{\max}} = {bc / np.sqrt(2):.1f}$")
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"growth rate $R(\beta)$")
    ax.set_title(f"Sim 4: Dispersion relation, c0 = {c0}")
    ax.legend(loc="lower left", fontsize=9)
    fig.tight_layout()
    out = "figures/sim4_dispersion.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[sim4] wrote {out}")


if __name__ == "__main__":
    main()

"""Generate publication-quality figures for the term paper from saved .npz data.

Saves all output to final_plots/. No simulations are rerun.

Figures produced
----------------
1. sim1_panels.png   — 4-panel snapshots, c0=0.5
2. sim2_panels.png   — 4-panel snapshots, c0=0.3
3. sim3_panels.png   — 4-panel snapshots, c0=0.22
4. sim123_comparison.png — late-time side-by-side comparison
5. sim4_dispersion.png   — R(beta): measured vs theory
6. sim5_kpeak.png    — k_peak(t) coarsening
7. sim6_F_decay.png  — F(t) descent (linear + log-log)
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import LogNorm

from analysis import free_energy_series, theoretical_growth_rate, beta_critical
from solver import f0pp_polynomial

# ---- global style ----
rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "text.usetex": False,   # use mathtext; set True if LaTeX is available
})

OUT = "final_plots"
os.makedirs(OUT, exist_ok=True)
CMAP = "coolwarm"
CB_LABEL = r"$c_B$"
T_LABEL  = "dimensionless time, $t$"


def _closest(times, target):
    return int(np.argmin(np.abs(times - target)))


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def _panel_figure(snaps, times, c0, label, out):
    """4-panel snapshot figure at four representative times."""
    t_targets = [0.0, times[-1] * 0.03, times[-1] * 0.15, times[-1]]
    idxs = [_closest(times, tt) for tt in t_targets]

    fig, axes = plt.subplots(1, 4, figsize=(13.5, 3.6))
    for ax, idx in zip(axes, idxs):
        im = ax.imshow(snaps[idx].T, origin="lower", cmap=CMAP,
                       vmin=0, vmax=1, extent=[0, 1, 0, 1])
        ax.set_title(f"$t = {times[idx]:.1f}$", pad=4)
        ax.set_xticks([]); ax.set_yticks([])

    fig.suptitle(label, y=1.02)
    cbar = fig.colorbar(im, ax=axes.tolist(), fraction=0.018, pad=0.02)
    cbar.set_label(CB_LABEL)
    fig.savefig(out)
    plt.close(fig)
    print(f"  wrote {out}")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 1–3: morphology panels
# ═══════════════════════════════════════════════════════════════════════════

for npz, c0_val, label, fname in [
    ("data/sim1.npz", 0.5,  r"Symmetric decomposition, $c_0 = 0.5$",   "sim1_panels.png"),
    ("data/sim2.npz", 0.3,  r"Off-symmetric decomposition, $c_0 = 0.3$", "sim2_panels.png"),
    ("data/sim3.npz", 0.22, r"Near-spinodal decomposition, $c_0 = 0.22$","sim3_panels.png"),
]:
    d = np.load(npz)
    _panel_figure(d["snaps"], d["times"], c0_val, label, f"{OUT}/{fname}")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 4: side-by-side late-time comparison (sim1, sim2, sim3)
# ═══════════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(1, 3, figsize=(12, 4.2))
labels_c0 = [
    (0.5,  r"$c_0 = 0.5$ — bicontinuous"),
    (0.3,  r"$c_0 = 0.3$ — droplets"),
    (0.22, r"$c_0 = 0.22$ — near-spinodal"),
]
for ax, npz, (c0_val, lbl) in zip(axes,
        ["data/sim1.npz", "data/sim2.npz", "data/sim3.npz"], labels_c0):
    d = np.load(npz)
    snap = d["snaps"][-1]
    t_end = float(d["times"][-1])
    im = ax.imshow(snap.T, origin="lower", cmap=CMAP, vmin=0, vmax=1,
                   extent=[0, 1, 0, 1])
    ax.set_title(f"{lbl}\n$t = {t_end:.0f}$", pad=5)
    ax.set_xticks([]); ax.set_yticks([])

fig.suptitle("Morphology dependence on mean composition", y=1.02)
cbar = fig.colorbar(im, ax=axes.tolist(), fraction=0.018, pad=0.02)
cbar.set_label(CB_LABEL)
out = f"{OUT}/sim123_comparison.png"
fig.savefig(out); plt.close(fig)
print(f"  wrote {out}")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 5: dispersion relation
# ═══════════════════════════════════════════════════════════════════════════

d4 = np.load("data/sim4.npz")
betas   = d4["betas"]
R_meas  = d4["R_meas"]
R_err   = d4["R_err"]
c0_4    = float(d4["c0"])
kappa_4 = float(d4["kappa"])
M_4     = float(d4["M"])

bc   = beta_critical(c0_4, kappa=kappa_4)
bmax = bc / np.sqrt(2)

beta_th = np.linspace(0, 1.05 * betas.max(), 600)
R_th    = theoretical_growth_rate(beta_th, c0_4, kappa_4, M_4)
err_plot = np.where(np.isfinite(R_err), R_err, 0.0)

def _draw_dispersion(ax, xlim=None, ylim=None, legend=True):
    ax.plot(beta_th, R_th, "k-", lw=1.8,
            label=r"Theory: $R(\beta)=-M\beta^2[f_0''+2\kappa\beta^2]$")
    ax.errorbar(betas, R_meas, yerr=err_plot, fmt="o", color="crimson",
                ms=6, capsize=3, label="Measured (simulation)")
    ax.axhline(0, color="0.55", lw=0.8)
    ax.axvline(bc,   ls="--", color="steelblue", lw=1.2,
               label=fr"$\beta_c = {bc:.1f}$")
    ax.axvline(bmax, ls=":",  color="darkorange", lw=1.2,
               label=fr"$\beta_{{max}} = {bmax:.1f}$")
    ax.set_xlabel(r"Wavenumber $\beta$")
    ax.set_ylabel(r"$R(\beta)$")
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    if legend: ax.legend(loc="lower left", framealpha=0.9, fontsize=9)

fig, (ax_full, ax_zoom) = plt.subplots(1, 2, figsize=(13, 4.8))

# Left: full curve — shows the overall shape including strongly decaying modes
_draw_dispersion(ax_full, legend=True)
ax_full.set_title("Full dispersion curve")

# Right: zoom on the unstable band where R > 0 and the match matters most
R_max_th = theoretical_growth_rate(bmax, c0_4, kappa_4, M_4)
_draw_dispersion(ax_zoom,
                 xlim=(0, bc * 1.3),
                 ylim=(-R_max_th * 0.6, R_max_th * 1.4),
                 legend=False)
ax_zoom.set_title("Unstable band (zoom)")
ax_zoom.annotate(
    "Deviations near $\\beta_c$:\nlinear regime too\nshort for accurate fit",
    xy=(bc * 0.97, R_meas[4] if len(R_meas) > 4 else 0),
    xytext=(bc * 0.55, R_max_th * 0.5),
    arrowprops=dict(arrowstyle="->", color="0.4"),
    fontsize=9, color="0.35"
)

fig.suptitle(
    r"Dispersion relation: measured vs. Cahn (1961), $c_0 = 0.5$",
    fontsize=12
)
fig.tight_layout()
out = f"{OUT}/sim4_dispersion.png"
fig.savefig(out); plt.close(fig)
print(f"  wrote {out}")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 6: k_peak(t) coarsening
# ═══════════════════════════════════════════════════════════════════════════

d5      = np.load("data/sim5.npz")
times5  = d5["times"]
k_peak  = d5["k_peak_t"]
kappa5  = 5e-4
bc5     = beta_critical(0.5, kappa=kappa5)
bmax5   = bc5 / np.sqrt(2)

fig, ax = plt.subplots(figsize=(6.5, 4.5))
mask = (times5 > 0) & (k_peak > 0)
ax.loglog(times5[mask], k_peak[mask], "o", ms=4, color="steelblue",
          label=r"$k_\mathrm{peak}(t)$ (measured)")

# Power-law fit over the second half of data
late = mask.copy()
late[:len(times5)//2] = False
if late.sum() > 3:
    slope, intercept = np.polyfit(np.log(times5[late]), np.log(k_peak[late]), 1)
    t_ref = times5[late]
    ax.loglog(t_ref, np.exp(intercept) * t_ref**slope, "r--",
              label=fr"Fit: $\propto t^{{{slope:.2f}}}$")

# t^{-1/3} reference anchored to mid-series
mid = len(times5[mask]) // 2
t0_ref = times5[mask][mid]; k0_ref = k_peak[mask][mid]
t_ref2 = times5[mask][mid:]
ax.loglog(t_ref2, k0_ref * (t_ref2 / t0_ref)**(-1/3), "k:", lw=1.4,
          alpha=0.6, label=r"$t^{-1/3}$ reference (LSW)")

ax.axhline(bmax5, ls="--", color="darkorange", lw=1.2,
           label=fr"$\beta_\mathrm{{max}} = {bmax5:.1f}$ (linear theory)")
ax.set_xlabel(T_LABEL)
ax.set_ylabel(r"$k_\mathrm{peak}$")
ax.set_title(r"Domain coarsening: peak wavenumber $k_\mathrm{peak}(t)$")
ax.legend(framealpha=0.9)
fig.tight_layout()
out = f"{OUT}/sim5_kpeak.png"
fig.savefig(out); plt.close(fig)
print(f"  wrote {out}")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 7: F(t) free energy decay
# ═══════════════════════════════════════════════════════════════════════════

d1    = np.load("data/sim1.npz")
times1 = d1["times"]
snaps1 = d1["snaps"]
L1     = float(d1["L"]) if "L" in d1.files else 1.0

print("  computing F(t) for sim6 figure ...")
F = free_energy_series(snaps1, L=L1, kappa=5e-4)
F_inf = F[-1]

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

# Left: linear scale
axes[0].plot(times1, F, "o-", ms=3, color="steelblue", lw=1.2)
axes[0].set_xlabel(T_LABEL)
axes[0].set_ylabel("$F(t)$  (dimensionless)")
axes[0].set_title("Free energy descent, $c_0 = 0.5$")

# Right: log-log of F - F_inf
diff = F - F_inf
mask2 = (times1 > 0) & (diff > 0)
axes[1].loglog(times1[mask2], diff[mask2], "o-", ms=3,
               color="steelblue", lw=1.2, label=r"$F(t) - F_\infty$")
if mask2.sum() > 6:
    late2 = mask2.copy()
    late2[:len(times1)//2] = False
    if late2.sum() > 3:
        slope2, intercept2 = np.polyfit(
            np.log(times1[late2]), np.log(diff[late2]), 1)
        t_ref3 = times1[late2]
        axes[1].loglog(t_ref3, np.exp(intercept2) * t_ref3**slope2,
                       "r--", label=fr"Fit slope = ${slope2:+.2f}$")
axes[1].set_xlabel(T_LABEL)
axes[1].set_ylabel(r"$F(t) - F_\infty$")
axes[1].set_title(r"Approach to equilibrium (log–log)")
axes[1].legend(framealpha=0.9)

n_up = int(np.sum(np.diff(F) > 0))
fig.suptitle(
    fr"Free energy $F[c]$: Lyapunov property"
    fr"  ({n_up} of {len(F)-1} timesteps with $\Delta F>0$, numerical noise)",
    y=1.02, fontsize=11
)
fig.tight_layout()
out = f"{OUT}/sim6_F_decay.png"
fig.savefig(out); plt.close(fig)
print(f"  wrote {out}")

print("\nAll figures written to final_plots/")

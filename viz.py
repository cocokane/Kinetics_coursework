"""Visualisation utilities — figures and animations.

What this module does
---------------------
Wraps matplotlib to produce all of the project's figures and animations.
Uses the Agg backend (no GUI required, fine over SSH) and ffmpeg for MP4
output. Nothing here computes physics; it only consumes snapshot arrays
and metadata produced by the solver and analysis modules.

Functions
---------
- animate_field(snapshots, times, out_path, fps, bitrate, title, vmin, vmax, cmap)
    Single-panel MP4 animation of c(r, t).

- snapshot_panels(snapshots, times, indices, out_path, ...)
    Static figure with one composition map per requested index. Used for
    the "early / mid / late" 4-panel summary in each morphology sim.

- histogram_heatmap(centers, hist_t, times, out_path, title)
    2D heatmap of histogram(c) vs t. Color = density. Composition on y,
    time on x.

- Skt_heatmap(k_centers, S_t, times, out_path, title, log_y)
    Log-scaled heatmap of the radially averaged structure factor over time.

- animate_Gc_and_map(snapshots, times, out_path, f0_func, c0, binodal,
                     spinodal, ...)
    The "killer demo" (Sim 7): two synchronised panels — left shows the
    f0(c) curve plus a live histogram bar overlay; right shows c(r, t).
    Audience can read off thermodynamic descent vs spatial pattern in
    real time.

- animate_two_fields(snapshots_a, snapshots_b, times, out_path, ...)
    Side-by-side MP4 used by Sim 9 to compare deterministic CH against
    CHC noise on identical initial conditions.

Knobs you might want to tune
----------------------------
- fps         : default 15. Lower for slower playback (good for long sims).
- bitrate     : default 2000. Increase for higher quality, larger files.
- cmap        : default 'coolwarm' (blue=0, red=1). Try 'RdBu_r' for a
                more classic phase-field look or 'viridis' if you need
                colorblind-safe.
- vmin/vmax   : default [0, 1]. Set tighter when c0 is far from 0.5 to
                use the full dynamic range of the colormap.
- dpi         : 300 for static PNGs, 150 for animations (a static-figure
                default change requires editing snapshot_panels directly).
- bins (in animate_Gc_and_map): hard-coded to 60 currently; edit if you
                want a smoother histogram on the left panel.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LogNorm
from matplotlib import gridspec


# ---------- 2D animation of composition map ----------

def animate_field(
    snapshots,
    times,
    out_path,
    fps=15,
    bitrate=2000,
    title=None,
    vmin=0.0,
    vmax=1.0,
    cmap="coolwarm",
):
    """Save MP4 animation of c(r, t)."""
    N = snapshots[0].shape[0]
    fig, ax = plt.subplots(figsize=(5, 5))
    im = ax.imshow(snapshots[0].T, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax,
                   extent=[0, 1, 0, 1])
    ax.set_xticks([])
    ax.set_yticks([])
    ttl = ax.set_title(_title_str(title, times[0]))
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("c")

    def update(i):
        im.set_array(snapshots[i].T)
        ttl.set_text(_title_str(title, times[i]))
        return [im, ttl]

    anim = animation.FuncAnimation(fig, update, frames=len(snapshots), blit=False)
    _save_anim(anim, out_path, fps=fps, bitrate=bitrate)
    plt.close(fig)


def _title_str(title, t):
    base = f"t = {t:.2f}"
    return f"{title} — {base}" if title else base


def _save_anim(anim, out_path, fps=15, bitrate=2000):
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    Writer = animation.writers["ffmpeg"]
    writer = Writer(fps=fps, bitrate=bitrate)
    anim.save(out_path, writer=writer, dpi=150)


# ---------- Snapshot panel (static figure) ----------

def snapshot_panels(snapshots, times, indices, out_path, title=None, vmin=0.0, vmax=1.0,
                    cmap="coolwarm"):
    n = len(indices)
    fig, axes = plt.subplots(1, n, figsize=(3.2 * n, 3.4))
    if n == 1:
        axes = [axes]
    for ax, idx in zip(axes, indices):
        im = ax.imshow(snapshots[idx].T, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax,
                       extent=[0, 1, 0, 1])
        ax.set_title(f"t = {times[idx]:.1f}")
        ax.set_xticks([])
        ax.set_yticks([])
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02, label="c")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------- Histogram time evolution heatmap ----------

def histogram_heatmap(centers, hist_t, times, out_path, title=None):
    """Heatmap: composition (y) vs time (x, log scale), color = density.

    Log time axis is essential: spinodal decomposition completes by t~2 while
    the full run extends to t~75. On a linear axis all the interesting
    early-time splitting of the peak is squashed into the leftmost ~3% of
    the plot. Log scale spreads the decomposition and coarsening regimes
    proportionally.
    """
    fig, ax = plt.subplots(figsize=(7, 4))
    # Drop t=0 (log of 0 is undefined) and any near-zero times
    mask = times > 0
    t_plot = times[mask]
    h_plot = hist_t[mask]
    pcm = ax.pcolormesh(t_plot, centers, h_plot.T, shading="auto", cmap="magma")
    ax.set_xscale("log")
    ax.set_xlabel("t  (log scale)")
    ax.set_ylabel("c")
    if title:
        ax.set_title(title)
    fig.colorbar(pcm, ax=ax, label="density")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------- Structure-factor heatmap ----------

def Skt_heatmap(k_centers, S_t, times, out_path, title=None, log_y=True):
    """Heatmap of S(k, t) with optional log y-axis."""
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    # Mask zeros for LogNorm
    Sp = np.maximum(S_t.T, 1e-20)
    pcm = ax.pcolormesh(times, k_centers, Sp, norm=LogNorm(vmin=Sp.max() * 1e-6, vmax=Sp.max()),
                        cmap="viridis", shading="auto")
    if log_y:
        ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel("t")
    ax.set_ylabel("k")
    if title:
        ax.set_title(title)
    fig.colorbar(pcm, ax=ax, label="S(k)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------- Killer demo: G(c) + map coevolution ----------

def animate_Gc_and_map(snapshots, times, out_path, f0_func, c0,
                       binodal=(0.0, 1.0), spinodal=None,
                       fps=15, bitrate=2000, title=None,
                       vmin=0.0, vmax=1.0, cmap="coolwarm"):
    """Two synchronised panels: f0(c) + histogram on left, c(r,t) map on right.

    Left panel layout (2-row):
      - Top: f0(c) curve with binodal/spinodal markers + a scatter overlay
             where each dot sits at (c_bin, f0(c_bin)) and its size and
             colour are proportional to the histogram density at that
             composition.  This shows *where on the free-energy curve* the
             system's compositions currently live.
      - Bottom: a clean histogram bar chart of the same distribution.

    Right panel: the composition map c(r,t).
    """
    fig = plt.figure(figsize=(12, 5.6))
    outer = gridspec.GridSpec(1, 2, width_ratios=[1.15, 1.0], figure=fig,
                              left=0.07, right=0.97, wspace=0.38)
    # Two-row left panel: f0 curve (top) + histogram bars (bottom)
    inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[0],
                                             height_ratios=[2.2, 1.0], hspace=0.08)
    ax_g = fig.add_subplot(inner[0])   # f0(c) curve + scatter overlay
    ax_h = fig.add_subplot(inner[1], sharex=ax_g)   # histogram bars
    ax_m = fig.add_subplot(outer[1])   # composition map

    # ---- f0(c) curve ----
    c_curve = np.linspace(-0.02, 1.02, 500)
    ax_g.plot(c_curve, f0_func(c_curve), color="0.3", lw=1.8, zorder=2)
    ax_g.set_ylabel("f₀(c)", fontsize=9)
    ax_g.set_title("Free energy + composition density", fontsize=9)
    ax_g.tick_params(labelbottom=False)

    # Binodal and spinodal markers (drawn once, static)
    marker_handles = []
    if binodal is not None:
        for cb in binodal:
            ax_g.axvline(cb, ls=":", color="limegreen", lw=1.2, alpha=0.8)
            ax_h.axvline(cb, ls=":", color="limegreen", lw=1.2, alpha=0.8)
        pts = ax_g.scatter(list(binodal),
                           [float(f0_func(np.array([cb]))) for cb in binodal],
                           color="limegreen", s=60, zorder=5, label="binodal")
        marker_handles.append(pts)
    if spinodal is not None:
        for cs in spinodal:
            ax_g.axvline(cs, ls="--", color="darkorange", lw=1.2, alpha=0.8)
            ax_h.axvline(cs, ls="--", color="darkorange", lw=1.2, alpha=0.8)
        pts = ax_g.scatter(list(spinodal),
                           [float(f0_func(np.array([cs]))) for cs in spinodal],
                           color="darkorange", marker="^", s=60, zorder=5, label="spinodal")
        marker_handles.append(pts)
    if marker_handles:
        ax_g.legend(handles=marker_handles, loc="upper right", fontsize=8)

    # Scatter overlay: dots at (c_bin, f0(c_bin)), size + colour ∝ density.
    # Use a scatter (not a Line2D) so we can call set_sizes() and set_array()
    # on each animation frame.
    bins_edges = np.linspace(0.0, 1.0, 61)
    bin_centers = 0.5 * (bins_edges[:-1] + bins_edges[1:])
    f0_at_bins = np.array([float(f0_func(np.array([c]))) for c in bin_centers])
    sc = ax_g.scatter(bin_centers, f0_at_bins,
                      c=np.zeros(len(bin_centers)), cmap="hot_r",
                      vmin=0.0, vmax=4.0, s=0, alpha=0.9, zorder=6)

    # ---- Histogram bar chart (bottom sub-panel) ----
    bar_w = (bins_edges[1] - bins_edges[0]) * 0.92
    bar_container = ax_h.bar(bin_centers, np.zeros(len(bin_centers)),
                             width=bar_w, color="steelblue", alpha=0.7)
    ax_h.set_xlabel("c", fontsize=9)
    ax_h.set_ylabel("density", fontsize=9)
    ax_h.set_xlim(-0.02, 1.02)
    ax_h.set_ylim(0, 5.0)   # typical max density ~ 3-4; override if needed
    ax_g.set_xlim(-0.02, 1.02)

    # ---- Composition map (right panel) ----
    im = ax_m.imshow(snapshots[0].T, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax,
                     extent=[0, 1, 0, 1])
    ax_m.set_xticks([])
    ax_m.set_yticks([])
    ttl_m = ax_m.set_title(f"c(r, t)   t = {times[0]:.2f}", fontsize=9)
    fig.colorbar(im, ax=ax_m, fraction=0.046, pad=0.04, label="c")

    if title:
        fig.suptitle(title, fontsize=10, y=1.01)

    def update(i):
        s = snapshots[i]
        counts, _ = np.histogram(s.ravel(), bins=bins_edges, density=True)
        peak = counts.max() + 1e-12

        # Scatter overlay: size proportional to density, colour from "hot_r"
        sizes = 180.0 * counts / peak
        sc.set_sizes(sizes)
        sc.set_array(counts)
        sc.set_clim(0.0, peak)

        # Histogram bars: resize and re-scale y-axis
        for bar, h in zip(bar_container, counts):
            bar.set_height(h)
        ax_h.set_ylim(0, max(peak * 1.1, 0.5))

        # Composition map
        im.set_array(s.T)
        ttl_m.set_text(f"c(r, t)   t = {times[i]:.2f}")
        return [im, sc, ttl_m, *bar_container]

    anim = animation.FuncAnimation(fig, update, frames=len(snapshots), blit=False)
    _save_anim(anim, out_path, fps=fps, bitrate=bitrate)
    plt.close(fig)


# ---------- Free-energy descent animation ----------

def animate_field_and_F(snapshots, times, F_values, out_path,
                        fps=15, bitrate=2000, cmap="coolwarm",
                        vmin=0.0, vmax=1.0, title=None):
    """Left panel: c(r,t) composition map. Right panel: F(t) descent curve
    drawn incrementally so the audience watches free energy fall in real time.

    Parameters
    ----------
    F_values : 1D array, same length as snapshots — free energy at each frame.
    """
    fig, (ax_m, ax_F) = plt.subplots(1, 2, figsize=(11, 4.8))

    # Composition map
    im = ax_m.imshow(snapshots[0].T, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax,
                     extent=[0, 1, 0, 1])
    ax_m.set_xticks([])
    ax_m.set_yticks([])
    ttl_m = ax_m.set_title(f"c(r, t)   t = {times[0]:.2f}")
    fig.colorbar(im, ax=ax_m, fraction=0.046, pad=0.04, label="c")

    # F(t) axes — set limits up front so the axes don't jump during animation
    F_hi = F_values[0] * 1.05
    F_lo = F_values[-1] * 0.95
    ax_F.set_xlim(times[0], times[-1])
    ax_F.set_ylim(F_lo, F_hi)
    ax_F.set_xlabel("t")
    ax_F.set_ylabel("F(t)")
    ax_F.set_title("Free energy descent")
    ax_F.axhline(F_values[-1], ls=":", color="0.55", lw=1.2, label="F(∞) ≈ equilibrium")
    ax_F.legend(fontsize=9)

    line, = ax_F.plot([], [], lw=1.8, color="royalblue")
    dot,  = ax_F.plot([], [], "o", ms=6, color="crimson")

    if title:
        fig.suptitle(title)
    fig.tight_layout()

    def update(i):
        im.set_array(snapshots[i].T)
        ttl_m.set_text(f"c(r, t)   t = {times[i]:.2f}")
        line.set_data(times[:i + 1], F_values[:i + 1])
        dot.set_data([times[i]], [F_values[i]])
        return [im, ttl_m, line, dot]

    anim = animation.FuncAnimation(fig, update, frames=len(snapshots), blit=False)
    _save_anim(anim, out_path, fps=fps, bitrate=bitrate)
    plt.close(fig)


# ---------- Side-by-side animation (deterministic vs CHC) ----------

def animate_two_fields(snapshots_a, snapshots_b, times, out_path,
                       title_a="A", title_b="B", overall_title=None,
                       fps=15, bitrate=2000, vmin=0.0, vmax=1.0, cmap="coolwarm"):
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.6))
    im_a = axes[0].imshow(snapshots_a[0].T, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax,
                          extent=[0, 1, 0, 1])
    im_b = axes[1].imshow(snapshots_b[0].T, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax,
                          extent=[0, 1, 0, 1])
    for ax, t in zip(axes, [title_a, title_b]):
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(t)
    sup = fig.suptitle(_title_str(overall_title, times[0]))

    def update(i):
        im_a.set_array(snapshots_a[i].T)
        im_b.set_array(snapshots_b[i].T)
        sup.set_text(_title_str(overall_title, times[i]))
        return [im_a, im_b, sup]

    fig.tight_layout()
    anim = animation.FuncAnimation(fig, update, frames=len(times), blit=False)
    _save_anim(anim, out_path, fps=fps, bitrate=bitrate)
    plt.close(fig)

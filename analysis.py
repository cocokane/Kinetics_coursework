"""Analysis utilities — diagnostics computed on solver snapshots.

What this module does
---------------------
Provides post-hoc measurements taken on c(r, t) snapshots produced by the
solver. Nothing here advances time; everything is a function of the field
and (optionally) the simulation parameters that produced it.

Functions
---------
- composition_histogram(c, bins, range_)
    1D histogram of voxel compositions. Used to visualise the system's
    distribution over the f0(c) curve.

- histogram_evolution(snapshots, bins, range_)
    Stack 1D histograms across a time series into a (n_frames, n_bins)
    array. Drives the histogram-vs-time heatmap (Sims 1-3) and the
    coevolution animation (Sim 7).

- radial_structure_factor(c, L, n_bins)
    Radially averaged S(k) = <|c_hat|^2>_{|k|=k}. Sim 5's primary
    diagnostic; gives the dominant wavelength without needing a controlled
    initial condition.

- peak_k(k_centers, S, k_min)
    Argmax of S(k) above k_min. Tracks the dominant length scale during
    coarsening (used to verify the t^{-1/3} scaling).

- free_energy(c, L, kappa, f0)
    F = integral [ f0(c) + kappa |grad c|^2 ] dV via spectral derivatives.
    Sim 6's Lyapunov check.

- amplitude(c)
    Real-space peak-to-peak / 2 of c. Crude amplitude metric, only useful
    for very small perturbations of a uniform background.

- fit_growth_rate(times, amps, A_threshold)
    Linear regression of log(A) vs t in the regime A < A_threshold. Used
    inside Sim 4. (Note: Sim 4 actually uses a stricter "leading monotone
    window" defined inside run_sim4.py to handle decaying modes cleanly.)

- theoretical_growth_rate(beta, c0, kappa, M, f0pp)
    R(beta) = -M beta^2 [f0''(c0) + 2 kappa beta^2] — Cahn (1961)'s
    linear-stability prediction. Plotted as the reference curve in Sim 4.

- beta_critical(c0, kappa, f0pp)
    sqrt(-f0''(c0) / (2 kappa)). Returns NaN if c0 is outside the
    spinodal (f0''(c0) >= 0). beta_max = beta_c / sqrt(2).

Knobs you might want to tune
----------------------------
- bins / range_ in composition_histogram: 80 bins over [-0.1, 1.1] is the
  default. Increase bins (e.g. 120) for sharper histogram videos.
- n_bins in radial_structure_factor: defaults to N//2; finer bins are
  noisier but resolve narrow peaks.
- f0pp / f0 callables: pass the Landau equivalents for Sim 8.
"""
from __future__ import annotations

import numpy as np

from solver import f0_polynomial, make_k_grid


def composition_histogram(c, bins=80, range_=(-0.1, 1.1)):
    """Return (centers, counts) for a 1D histogram of c values."""
    counts, edges = np.histogram(c.ravel(), bins=bins, range=range_, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, counts


def histogram_evolution(snapshots, bins=80, range_=(-0.1, 1.1)):
    """Stack 1D histograms from a sequence of snapshots into a (n_frames, bins) array."""
    centers, _ = composition_histogram(snapshots[0], bins=bins, range_=range_)
    out = np.empty((len(snapshots), bins))
    for i, s in enumerate(snapshots):
        _, counts = composition_histogram(s, bins=bins, range_=range_)
        out[i] = counts
    return centers, out


# ---------- Structure factor ----------

def radial_structure_factor(c, L=1.0, n_bins=None):
    """Radially averaged structure factor S(k).

    Parameters
    ----------
    c : (N, N) array
    L : box edge length
    n_bins : number of radial bins (default N//2)

    Returns
    -------
    k_centers, S : 1D arrays of length n_bins.
    """
    N = c.shape[0]
    c_centered = c - c.mean()
    c_hat = np.fft.fft2(c_centered)
    P = np.abs(c_hat) ** 2 / (N * N)

    kx = 2.0 * np.pi * np.fft.fftfreq(N, d=L / N)
    KX, KY = np.meshgrid(kx, kx, indexing="ij")
    K = np.sqrt(KX ** 2 + KY ** 2)

    if n_bins is None:
        n_bins = N // 2
    k_max = np.pi * N / L  # Nyquist
    bins = np.linspace(0.0, k_max, n_bins + 1)
    k_centers = 0.5 * (bins[:-1] + bins[1:])

    inds = np.digitize(K.ravel(), bins) - 1
    flatP = P.ravel()
    S = np.zeros(n_bins)
    counts = np.zeros(n_bins)
    valid = (inds >= 0) & (inds < n_bins)
    np.add.at(S, inds[valid], flatP[valid])
    np.add.at(counts, inds[valid], 1.0)
    counts[counts == 0] = 1.0
    S /= counts
    return k_centers, S


def peak_k(k_centers, S, k_min=1e-6):
    """Return k value at the peak of S(k), restricting to k > k_min."""
    mask = k_centers > k_min
    sub = S[mask]
    sub_k = k_centers[mask]
    return sub_k[np.argmax(sub)]


# ---------- Free energy ----------

def free_energy(c, L=1.0, kappa=5.0e-4, f0=f0_polynomial):
    """F = integral [ f0(c) + kappa |grad c|^2 ] dV via spectral derivatives."""
    N = c.shape[0]
    dx = L / N

    KX, KY, _, _ = make_k_grid(N, L)
    c_hat = np.fft.rfft2(c)
    dcdx = np.fft.irfft2(1j * KX * c_hat, s=(N, N))
    dcdy = np.fft.irfft2(1j * KY * c_hat, s=(N, N))
    grad2 = dcdx ** 2 + dcdy ** 2
    integrand = f0(c) + kappa * grad2
    return integrand.sum() * dx * dx


def free_energy_series(snapshots, L=1.0, kappa=5.0e-4, f0=f0_polynomial):
    return np.array([free_energy(s, L=L, kappa=kappa, f0=f0) for s in snapshots])


# ---------- Amplitude tracking (for dispersion relation) ----------

def amplitude(c):
    """A(t) = 0.5 * (max(c) - min(c)) — peak-to-peak amplitude / 2."""
    return 0.5 * (c.max() - c.min())


def fit_growth_rate(times, amps, A_threshold=0.05):
    """Fit ln(A) vs t in the linear-growth regime (until A exceeds A_threshold).

    Returns (R, R_err, t_used, lnA_used). If <3 points usable, returns (nan, nan, ...).
    """
    times = np.asarray(times)
    amps = np.asarray(amps)
    mask = (amps > 0) & (amps < A_threshold)
    if mask.sum() < 3:
        return np.nan, np.nan, times[mask], np.log(np.maximum(amps[mask], 1e-30))
    t = times[mask]
    y = np.log(amps[mask])
    # Linear regression
    A = np.vstack([t, np.ones_like(t)]).T
    sol, residuals, rank, sv = np.linalg.lstsq(A, y, rcond=None)
    R = sol[0]
    # Standard error from residuals
    n = len(t)
    if n > 2 and len(residuals):
        sigma2 = residuals[0] / (n - 2)
        cov = sigma2 * np.linalg.inv(A.T @ A)
        R_err = float(np.sqrt(cov[0, 0]))
    else:
        R_err = np.nan
    return float(R), R_err, t, y


def theoretical_growth_rate(beta, c0, kappa=5.0e-4, M=1.0, f0pp=None):
    """R(beta) = -M beta^2 [ f0''(c0) + 2 kappa beta^2 ]."""
    if f0pp is None:
        from solver import f0pp_polynomial
        f0pp = f0pp_polynomial
    f2 = f0pp(c0)
    beta2 = np.asarray(beta) ** 2
    return -M * beta2 * (f2 + 2.0 * kappa * beta2)


def beta_critical(c0, kappa=5.0e-4, f0pp=None):
    """beta_c such that R(beta_c) = 0; only real if f0''(c0) < 0."""
    if f0pp is None:
        from solver import f0pp_polynomial
        f0pp = f0pp_polynomial
    f2 = f0pp(c0)
    if f2 >= 0:
        return np.nan
    return float(np.sqrt(-f2 / (2.0 * kappa)))

"""Cahn-Hilliard semi-implicit Fourier-spectral solver.

What this module does
---------------------
Time-integrates the 2D Cahn-Hilliard equation on a periodic square box using
a stabilized semi-implicit (IMEX) Fourier-spectral scheme. The biharmonic
operator is diagonal in Fourier space, so the linear stiff term is treated
implicitly while the nonlinear bulk-chemical-potential term is explicit.
A linear stabilization term (Eyre-style splitting) is added with weight A
on both sides of the equation; this guarantees unconditional linear stability
for A >= max|f0''| over the dynamical range, at a small cost in temporal
accuracy.

Reference: Zhu, Chen, Shen & Tikare (1999), Phys. Rev. E 60, 3564.

Equation
--------
    dc/dt = M * laplacian( mu(c) - 2 kappa * laplacian(c) ) + zeta

with mu(c) = df0/dc, kappa the gradient-energy coefficient, M the mobility,
and zeta an optional Cook conserved-noise term (CHC).

Discrete update (stabilized; A is the stabilization constant)
-------------------------------------------------------------
    c_hat^{n+1} (1 + dt M A k^2 + 2 dt M kappa k^4)
        = c_hat^n - dt M k^2 (mu_hat^n - A c_hat^n) + (noise term)

Knobs (CahnHilliardSolver)
--------------------------
- N            : grid size N x N (typical: 128 production, 256 high-quality).
- L            : box length (default 1.0; dimensionless).
- kappa        : gradient-energy coefficient (sets interface width ~sqrt(kappa)).
- M            : mobility (sets the timescale; default 1.0).
- dt           : time step. Must be << 1/R_max where R_max is the peak linear
                 growth rate; for f0=c^2(1-c)^2 with kappa=5e-4 at c0=0.5,
                 R_max ~ 250, so dt=0.01 in the production runs is "coarse"
                 (relies on stabilization to remain stable; the linear regime
                 is *not* resolved at this dt — see Sim 4 for fine-dt).
- mu_func      : callable mu(c). Default: mu_polynomial (f0=c^2(1-c)^2).
                 Pass mu_landau (or any other) to swap free energies.
- noise_amp    : amplitude of conserved Cook noise. 0 (default) -> pure CH.
                 Used by Sim 9 only.
- rng          : numpy Generator for reproducibility of Cook noise.
- stabilization: A in the Eyre splitting. A=2 is the default for the
                 polynomial f0 (max|f0''|=2 inside [0,1]). Set to 0 for
                 pure semi-implicit (only safe for short, small-amplitude
                 runs like Sim 4's cosine-IC sweep).

Free-energy options
-------------------
- f0_polynomial(c)        — symmetric double well c^2 (1-c)^2.
- f0_landau(c, T, ...)    — Landau form a(T-Tc)(c-0.5)^2 + b(c-0.5)^4 used
                            by Sim 8 to vary effective temperature.

Initial conditions
------------------
- set_random_ic(c0, amp, seed) — Gaussian noise of amplitude `amp` around c0.
- set_cosine_ic(c0, beta, amp) — single-mode c0 + amp cos(beta x). Used by
                                 the dispersion-relation measurement in Sim 4.
- set_field(c)                 — explicit user-provided field.
"""
from __future__ import annotations

import numpy as np


# ---------- Free energies ----------

def f0_polynomial(c):
    """f0(c) = c^2 (1-c)^2 — symmetric double well with minima at 0, 1."""
    return c * c * (1.0 - c) ** 2


def mu_polynomial(c):
    """mu(c) = df0/dc = 2 c (1-c) (1 - 2c)."""
    return 2.0 * c * (1.0 - c) * (1.0 - 2.0 * c)


def f0pp_polynomial(c):
    """f0''(c) = 2 - 12 c + 12 c^2."""
    return 2.0 - 12.0 * c + 12.0 * c * c


def f0_landau(c, T, Tc=1.0, a=1.0, b=1.0):
    """Landau-style f0(c, T) = a (T - Tc) (c - 0.5)^2 + b (c - 0.5)^4.

    For T < Tc this is a symmetric double well; for T > Tc it is a single well.
    """
    x = c - 0.5
    return a * (T - Tc) * x * x + b * x ** 4


def mu_landau(c, T, Tc=1.0, a=1.0, b=1.0):
    x = c - 0.5
    return 2.0 * a * (T - Tc) * x + 4.0 * b * x ** 3


def f0pp_landau(c, T, Tc=1.0, a=1.0, b=1.0):
    x = c - 0.5
    return 2.0 * a * (T - Tc) + 12.0 * b * x * x


# ---------- k-grid helpers ----------

def make_k_grid(N, L):
    """Build (kx, ky, k2, k4) for a 2D periodic grid with rfft2 layout.

    rfft2 over an (N, N) array produces an (N, N//2 + 1) complex array.
    """
    kx = 2.0 * np.pi * np.fft.fftfreq(N, d=L / N)
    ky = 2.0 * np.pi * np.fft.rfftfreq(N, d=L / N)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    k2 = KX ** 2 + KY ** 2
    k4 = k2 * k2
    return KX, KY, k2, k4


# ---------- Solver ----------

class CahnHilliardSolver:
    """Semi-implicit Fourier-spectral Cahn-Hilliard solver on a periodic 2D box.

    Parameters
    ----------
    N : int
        Linear grid size (NxN).
    L : float
        Box edge length (dimensionless).
    kappa : float
        Gradient-energy coefficient.
    M : float
        Mobility.
    dt : float
        Time step.
    mu_func : callable
        mu(c) -> array. Defaults to polynomial double well.
    noise_amp : float
        Cook-noise amplitude for the conserved noise term. 0 = deterministic CH.
    rng : np.random.Generator or None
        For reproducibility of Cook noise.
    """

    def __init__(
        self,
        N=256,
        L=1.0,
        kappa=5.0e-4,
        M=1.0,
        dt=0.01,
        mu_func=None,
        noise_amp=0.0,
        rng=None,
        stabilization=2.0,
    ):
        self.N = N
        self.L = L
        self.kappa = kappa
        self.M = M
        self.dt = dt
        self.mu_func = mu_func if mu_func is not None else mu_polynomial
        self.noise_amp = noise_amp
        self.rng = rng if rng is not None else np.random.default_rng(0)
        # A: linear stabilization constant (Eyre-style splitting). Should be
        # >= max|f0''| over the expected range of c. For f0=c^2(1-c)^2, |f0''|
        # ranges from 1 (at c=0.5) to 2 (at c=0,1) inside [0,1], and grows
        # quadratically outside; A=2 is the standard choice but we accept
        # a user override for the Landau f0 case.
        self.A = stabilization

        _, _, k2, k4 = make_k_grid(N, L)
        self.k2 = k2
        self.k4 = k4
        # Implicit denominator for the stabilized scheme:
        #   (1 + dt M A k^2 + 2 dt M kappa k^4)
        self.denom = 1.0 + self.dt * self.M * self.A * k2 + 2.0 * self.dt * self.M * self.kappa * k4

        self.c = None
        self.step_count = 0

    # ---- initialization ----

    def set_random_ic(self, c0, amp=0.01, seed=0):
        rng = np.random.default_rng(seed)
        self.c = rng.normal(loc=c0, scale=amp, size=(self.N, self.N))
        # enforce exact mean = c0 for clean mass conservation tests
        self.c += c0 - self.c.mean()
        self.step_count = 0

    def set_cosine_ic(self, c0, beta, amp=0.01):
        """c(x, y) = c0 + amp * cos(beta * x).

        beta has units of 1/length (matching the analytical dispersion relation).
        """
        x = np.linspace(0.0, self.L, self.N, endpoint=False)
        X, _ = np.meshgrid(x, x, indexing="ij")
        self.c = c0 + amp * np.cos(beta * X)
        self.step_count = 0

    def set_field(self, c):
        self.c = c.copy()
        self.step_count = 0

    # ---- single step ----

    def step(self):
        c = self.c
        mu = self.mu_func(c)
        c_hat = np.fft.rfft2(c)
        mu_hat = np.fft.rfft2(mu)

        # Stabilized semi-implicit scheme:
        #   (c^{n+1} - c^n)/dt = -M k^2 (mu^n - A c^n) - M A k^2 c^{n+1} - 2 M kappa k^4 c^{n+1}
        #   c^{n+1} (1 + dt M A k^2 + 2 dt M kappa k^4)
        #     = c^n - dt M k^2 (mu^n - A c^n)
        rhs = c_hat - self.dt * self.M * self.k2 * (mu_hat - self.A * c_hat)

        if self.noise_amp > 0.0:
            # Conserved Cook noise: zeta = M^{1/2} * div(eta), so in Fourier space
            # add an i k * eta term scaled by sqrt(2 M dt) * noise_amp. Variance
            # proportional to k^2 ensures conservation.
            eta_x = self.rng.standard_normal(c.shape)
            eta_y = self.rng.standard_normal(c.shape)
            eta_x_hat = np.fft.rfft2(eta_x)
            eta_y_hat = np.fft.rfft2(eta_y)
            KX, KY, _, _ = make_k_grid(self.N, self.L)
            scale = self.noise_amp * np.sqrt(2.0 * self.M * self.dt)
            noise_hat = 1j * (KX * eta_x_hat + KY * eta_y_hat) * scale
            rhs = rhs + noise_hat

        c_hat_new = rhs / self.denom

        # Enforce zero-mode (mass conservation): the (k=0) update is c_hat^{n+1}[0,0]
        # = c_hat^n[0,0] / (1 + 0 + 0) = c_hat^n[0,0], which is correct. Reset to be safe.
        c_hat_new[0, 0] = c_hat[0, 0]

        self.c = np.fft.irfft2(c_hat_new, s=(self.N, self.N))
        self.step_count += 1

    # ---- run with snapshots ----

    def run(self, n_steps, snapshot_every=50, verbose=False):
        """Run n_steps; return (times, snapshots) including initial frame."""
        snapshots = [self.c.copy()]
        times = [self.step_count * self.dt]
        for n in range(1, n_steps + 1):
            self.step()
            if n % snapshot_every == 0 or n == n_steps:
                snapshots.append(self.c.copy())
                times.append(self.step_count * self.dt)
                if verbose:
                    cm = self.c.mean()
                    print(f"  step {n:6d}  t={self.step_count * self.dt:8.3f}  "
                          f"mean={cm:+.3e}  range=[{self.c.min():+.3f}, {self.c.max():+.3f}]")
        return np.array(times), np.array(snapshots)

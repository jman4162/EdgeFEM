"""
Coupling matrix and array analysis plotting functions.

Provides visualizations for mutual coupling, active impedance,
and array performance metrics.
"""

import numpy as np
from typing import Optional, List, Union, Tuple, Dict

try:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def _check_matplotlib():
    if not HAS_MATPLOTLIB:
        raise ImportError(
            "matplotlib is required for plotting. "
            "Install with: pip install matplotlib"
        )


def plot_coupling_matrix(
    S: np.ndarray,
    title: str = "Coupling Matrix",
    show_magnitude: bool = True,
    cmap: str = "RdYlBu_r",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    figsize: Tuple[float, float] = (8, 7),
    save: Optional[str] = None,
    show: bool = True,
    annotate: bool = True,
) -> "plt.Figure":
    """
    Plot S-parameter coupling matrix as a heatmap.

    Args:
        S: S-parameter matrix (complex, N x N)
        title: Plot title
        show_magnitude: If True, show |S| in dB; if False, show phase
        cmap: Colormap name
        vmin, vmax: Color scale limits (auto if None)
        figsize: Figure size in inches
        save: Path to save figure
        show: Whether to display
        annotate: Whether to show values in cells

    Returns:
        matplotlib Figure object
    """
    _check_matplotlib()

    S = np.asarray(S)
    n_ports = S.shape[0]

    if show_magnitude:
        data = 20 * np.log10(np.abs(S) + 1e-30)
        label = "Magnitude (dB)"
        if vmin is None:
            vmin = -60
        if vmax is None:
            vmax = 0
    else:
        data = np.angle(S, deg=True)
        label = "Phase (degrees)"
        if vmin is None:
            vmin = -180
        if vmax is None:
            vmax = 180
        cmap = "hsv"

    fig, ax = plt.subplots(figsize=figsize)

    im = ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax, aspect='equal')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(label)

    # Set ticks
    ax.set_xticks(range(n_ports))
    ax.set_yticks(range(n_ports))
    ax.set_xticklabels([f"Port {i+1}" for i in range(n_ports)])
    ax.set_yticklabels([f"Port {i+1}" for i in range(n_ports)])
    ax.set_xlabel("Input Port")
    ax.set_ylabel("Output Port")
    ax.set_title(title)

    # Annotate cells
    if annotate and n_ports <= 10:
        for i in range(n_ports):
            for j in range(n_ports):
                val = data[i, j]
                color = 'white' if val < (vmin + vmax) / 2 else 'black'
                if show_magnitude:
                    text = f"{val:.1f}"
                else:
                    text = f"{val:.0f}"
                ax.text(j, i, text, ha='center', va='center',
                       color=color, fontsize=8)

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')

    if show:
        plt.show()

    return fig


def plot_coupling_vs_distance(
    S: np.ndarray,
    positions: np.ndarray,
    title: str = "Coupling vs Distance",
    figsize: Tuple[float, float] = (10, 6),
    save: Optional[str] = None,
    show: bool = True,
    fit_line: bool = True,
) -> "plt.Figure":
    """
    Plot mutual coupling magnitude vs element separation distance.

    Args:
        S: S-parameter matrix (complex, N x N)
        positions: Element positions, shape (N, 3) in meters
        title: Plot title
        figsize: Figure size
        save: Path to save figure
        show: Whether to display
        fit_line: Whether to show exponential fit line

    Returns:
        matplotlib Figure object
    """
    _check_matplotlib()

    S = np.asarray(S)
    positions = np.asarray(positions)
    n_ports = S.shape[0]

    # Extract off-diagonal coupling and distances
    distances = []
    coupling_db = []

    for i in range(n_ports):
        for j in range(i + 1, n_ports):
            dist = np.linalg.norm(positions[i] - positions[j])
            s_ij = S[i, j]
            mag_db = 20 * np.log10(np.abs(s_ij) + 1e-30)
            distances.append(dist * 1000)  # Convert to mm
            coupling_db.append(mag_db)

    distances = np.array(distances)
    coupling_db = np.array(coupling_db)

    fig, ax = plt.subplots(figsize=figsize)

    ax.scatter(distances, coupling_db, alpha=0.6, s=30)

    # Fit line (log-linear)
    if fit_line and len(distances) > 2:
        valid = coupling_db > -60
        if np.sum(valid) > 2:
            coeffs = np.polyfit(distances[valid], coupling_db[valid], 1)
            fit_x = np.linspace(distances.min(), distances.max(), 100)
            fit_y = np.polyval(coeffs, fit_x)
            ax.plot(fit_x, fit_y, 'r--', alpha=0.7,
                   label=f'Fit: {coeffs[0]:.2f} dB/mm')
            ax.legend()

    ax.set_xlabel("Distance (mm)")
    ax.set_ylabel("Coupling (dB)")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')

    if show:
        plt.show()

    return fig


def plot_active_impedance_scan(
    scan_results: List,
    z0: float = 50.0,
    title: str = "Active Impedance vs Scan Angle",
    figsize: Tuple[float, float] = (12, 5),
    save: Optional[str] = None,
    show: bool = True,
    element_idx: Optional[int] = None,
) -> "plt.Figure":
    """
    Plot active impedance and VSWR vs scan angle.

    Args:
        scan_results: List of ActiveImpedanceResult objects
        z0: Reference impedance in ohms
        title: Plot title
        figsize: Figure size
        save: Path to save figure
        show: Whether to display
        element_idx: Which element to plot (default: first)

    Returns:
        matplotlib Figure object
    """
    _check_matplotlib()

    if element_idx is None:
        element_idx = 0

    # Extract data
    theta_deg = []
    z_active_real = []
    z_active_imag = []
    vswr = []

    for result in scan_results:
        theta_deg.append(np.rad2deg(result.theta))
        z_active = result.Z_active
        if hasattr(z_active, '__len__'):
            z = z_active[element_idx]
        else:
            z = z_active
        z_active_real.append(z.real)
        z_active_imag.append(z.imag)

        gamma = result.Gamma_active
        if hasattr(gamma, '__len__'):
            g = gamma[element_idx]
        else:
            g = gamma
        mag_gamma = np.abs(g)
        if mag_gamma < 1:
            vswr.append((1 + mag_gamma) / (1 - mag_gamma))
        else:
            vswr.append(np.inf)

    theta_deg = np.array(theta_deg)
    z_active_real = np.array(z_active_real)
    z_active_imag = np.array(z_active_imag)
    vswr = np.array(vswr)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Active impedance
    ax1.plot(theta_deg, z_active_real, 'b-', label='Real', linewidth=1.5)
    ax1.plot(theta_deg, z_active_imag, 'r--', label='Imag', linewidth=1.5)
    ax1.axhline(y=z0, color='gray', linestyle=':', alpha=0.5, label=f'Z0={z0}')
    ax1.set_xlabel("Scan Angle (degrees)")
    ax1.set_ylabel("Active Impedance (ohms)")
    ax1.set_title(f"Active Impedance (Element {element_idx + 1})")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # VSWR
    vswr_plot = np.clip(vswr, 1, 20)
    ax2.plot(theta_deg, vswr_plot, 'g-', linewidth=1.5)
    ax2.axhline(y=2, color='orange', linestyle='--', alpha=0.7, label='VSWR=2')
    ax2.axhline(y=3, color='red', linestyle='--', alpha=0.7, label='VSWR=3')
    ax2.set_xlabel("Scan Angle (degrees)")
    ax2.set_ylabel("VSWR")
    ax2.set_title(f"Active VSWR (Element {element_idx + 1})")
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(1, 10)

    fig.suptitle(title)
    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')

    if show:
        plt.show()

    return fig


def plot_array_layout(
    positions: np.ndarray,
    coupling: Optional[np.ndarray] = None,
    title: str = "Array Layout",
    figsize: Tuple[float, float] = (10, 8),
    save: Optional[str] = None,
    show: bool = True,
) -> "plt.Figure":
    """
    Plot 2D array element layout with optional coupling visualization.

    Args:
        positions: Element positions, shape (N, 2) or (N, 3)
        coupling: Optional coupling matrix (N x N) to color connections
        title: Plot title
        figsize: Figure size
        save: Path to save figure
        show: Whether to display

    Returns:
        matplotlib Figure object
    """
    _check_matplotlib()

    positions = np.asarray(positions)
    if positions.shape[1] > 2:
        positions = positions[:, :2]  # Use X, Y only

    n_elements = len(positions)

    fig, ax = plt.subplots(figsize=figsize)

    # Draw coupling lines if provided
    if coupling is not None:
        coupling = np.asarray(coupling)
        coupling_db = 20 * np.log10(np.abs(coupling) + 1e-30)
        cmap = plt.cm.RdYlBu_r
        norm = plt.Normalize(vmin=-60, vmax=-10)

        for i in range(n_elements):
            for j in range(i + 1, n_elements):
                c_db = coupling_db[i, j]
                if c_db > -60:
                    color = cmap(norm(c_db))
                    alpha = min(1.0, (c_db + 60) / 50)
                    ax.plot([positions[i, 0]*1000, positions[j, 0]*1000],
                           [positions[i, 1]*1000, positions[j, 1]*1000],
                           color=color, alpha=alpha, linewidth=1)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label("Coupling (dB)")

    # Draw elements
    ax.scatter(positions[:, 0] * 1000, positions[:, 1] * 1000,
              s=100, c='blue', edgecolors='black', zorder=5)

    # Label elements
    for i, pos in enumerate(positions):
        ax.annotate(f"{i+1}", (pos[0]*1000, pos[1]*1000),
                   textcoords="offset points", xytext=(5, 5),
                   fontsize=8)

    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_title(title)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')

    if show:
        plt.show()

    return fig

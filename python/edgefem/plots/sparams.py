"""
S-parameter plotting functions.

Provides plots for frequency sweeps, Smith charts, and various
S-parameter visualizations.
"""

import numpy as np
from typing import Optional, List, Union, Tuple

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


def plot_sparams_vs_freq(
    freqs: Union[List[float], np.ndarray],
    S_matrices: List[np.ndarray],
    params: Optional[List[Tuple[int, int]]] = None,
    show_phase: bool = True,
    title: str = "S-Parameters vs Frequency",
    figsize: Tuple[float, float] = (10, 6),
    save: Optional[str] = None,
    show: bool = True,
    freq_unit: str = "GHz",
) -> "plt.Figure":
    """
    Plot S-parameters magnitude and phase vs frequency.

    Args:
        freqs: List/array of frequencies in Hz
        S_matrices: List of S-parameter matrices (complex, N x N)
        params: Which S-parameters to plot as (i, j) tuples.
                Default: all diagonal (S11, S22, ...) and off-diagonal
        show_phase: Whether to show phase subplot (default True)
        title: Plot title
        figsize: Figure size in inches (width, height)
        save: Path to save figure (None = don't save)
        show: Whether to display the figure
        freq_unit: Frequency unit for display ("Hz", "kHz", "MHz", "GHz")

    Returns:
        matplotlib Figure object
    """
    _check_matplotlib()

    freqs = np.asarray(freqs)
    n_ports = S_matrices[0].shape[0]

    # Frequency scaling
    freq_scales = {"Hz": 1, "kHz": 1e3, "MHz": 1e6, "GHz": 1e9}
    freq_scale = freq_scales.get(freq_unit, 1e9)
    freqs_scaled = freqs / freq_scale

    # Default: plot all unique S-parameters
    if params is None:
        params = []
        for i in range(n_ports):
            params.append((i, i))  # Diagonal (S11, S22, ...)
        for i in range(n_ports):
            for j in range(n_ports):
                if i != j:
                    params.append((i, j))  # Off-diagonal

    # Extract S-parameters vs frequency
    S_data = {}
    for (i, j) in params:
        label = f"S{i+1}{j+1}"
        S_data[label] = np.array([S[i, j] for S in S_matrices])

    # Create figure
    if show_phase:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)
    else:
        fig, ax1 = plt.subplots(figsize=(figsize[0], figsize[1]/2))
        ax2 = None

    # Magnitude plot
    colors = plt.cm.tab10(np.linspace(0, 1, len(params)))
    for idx, ((i, j), color) in enumerate(zip(params, colors)):
        label = f"S{i+1}{j+1}"
        mag_db = 20 * np.log10(np.abs(S_data[label]) + 1e-30)
        ax1.plot(freqs_scaled, mag_db, label=label, color=color, linewidth=1.5)

    ax1.set_ylabel("Magnitude (dB)")
    ax1.legend(loc='best', ncol=min(4, len(params)))
    ax1.grid(True, alpha=0.3)
    ax1.set_title(title)

    # Phase plot
    if show_phase:
        for idx, ((i, j), color) in enumerate(zip(params, colors)):
            label = f"S{i+1}{j+1}"
            phase_deg = np.angle(S_data[label], deg=True)
            ax2.plot(freqs_scaled, phase_deg, label=label, color=color, linewidth=1.5)

        ax2.set_ylabel("Phase (degrees)")
        ax2.set_xlabel(f"Frequency ({freq_unit})")
        ax2.set_ylim(-180, 180)
        ax2.grid(True, alpha=0.3)
    else:
        ax1.set_xlabel(f"Frequency ({freq_unit})")

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')

    if show:
        plt.show()

    return fig


def plot_sparams_smith(
    S_matrices: List[np.ndarray],
    params: Optional[List[Tuple[int, int]]] = None,
    freqs: Optional[np.ndarray] = None,
    title: str = "Smith Chart",
    figsize: Tuple[float, float] = (8, 8),
    save: Optional[str] = None,
    show: bool = True,
) -> "plt.Figure":
    """
    Plot S-parameters on a Smith chart.

    Args:
        S_matrices: List of S-parameter matrices (complex, N x N)
        params: Which S-parameters to plot as (i, j) tuples.
                Default: S11 only
        freqs: Optional frequencies for colormap
        title: Plot title
        figsize: Figure size in inches
        save: Path to save figure
        show: Whether to display the figure

    Returns:
        matplotlib Figure object
    """
    _check_matplotlib()

    if params is None:
        params = [(0, 0)]  # S11 by default

    fig, ax = plt.subplots(figsize=figsize)

    # Draw Smith chart grid
    _draw_smith_chart(ax)

    # Plot S-parameters
    colors = plt.cm.tab10(np.linspace(0, 1, len(params)))
    for (i, j), color in zip(params, colors):
        label = f"S{i+1}{j+1}"
        S_values = np.array([S[i, j] for S in S_matrices])

        if freqs is not None:
            # Color by frequency
            scatter = ax.scatter(
                S_values.real, S_values.imag,
                c=freqs/1e9, cmap='viridis', s=20, label=label
            )
            plt.colorbar(scatter, ax=ax, label='Frequency (GHz)')
        else:
            ax.plot(S_values.real, S_values.imag, 'o-',
                   color=color, markersize=4, linewidth=1, label=label)

    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)
    ax.set_aspect('equal')
    ax.set_xlabel("Real")
    ax.set_ylabel("Imaginary")
    ax.set_title(title)
    ax.legend(loc='upper right')

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')

    if show:
        plt.show()

    return fig


def _draw_smith_chart(ax):
    """Draw Smith chart grid lines."""
    # Unit circle
    circle = Circle((0, 0), 1, fill=False, color='black', linewidth=1)
    ax.add_patch(circle)

    # Constant resistance circles
    r_values = [0, 0.2, 0.5, 1, 2, 5]
    for r in r_values:
        if r == 0:
            # Imaginary axis portion
            theta = np.linspace(-np.pi/2, np.pi/2, 100)
            x = np.cos(theta)
            y = np.sin(theta)
            ax.plot(x, y, 'gray', linewidth=0.5, alpha=0.5)
        else:
            center = r / (r + 1)
            radius = 1 / (r + 1)
            theta = np.linspace(0, 2*np.pi, 100)
            x = center + radius * np.cos(theta)
            y = radius * np.sin(theta)
            # Only plot portion inside unit circle
            mask = x**2 + y**2 <= 1.01
            x[~mask] = np.nan
            ax.plot(x, y, 'gray', linewidth=0.5, alpha=0.5)

    # Constant reactance arcs
    x_values = [0.2, 0.5, 1, 2, 5, -0.2, -0.5, -1, -2, -5]
    for x_val in x_values:
        if x_val == 0:
            ax.axhline(y=0, color='gray', linewidth=0.5, alpha=0.5)
        else:
            center_y = 1 / x_val
            radius = abs(1 / x_val)
            theta = np.linspace(0, 2*np.pi, 200)
            x = 1 + radius * np.cos(theta)
            y = center_y + radius * np.sin(theta)
            # Only plot portion inside unit circle
            mask = x**2 + y**2 <= 1.01
            x[~mask] = np.nan
            ax.plot(x, y, 'gray', linewidth=0.5, alpha=0.5)


def plot_sparams_polar(
    S_matrices: List[np.ndarray],
    params: Optional[List[Tuple[int, int]]] = None,
    freqs: Optional[np.ndarray] = None,
    title: str = "S-Parameters (Polar)",
    figsize: Tuple[float, float] = (8, 8),
    save: Optional[str] = None,
    show: bool = True,
) -> "plt.Figure":
    """
    Plot S-parameters in polar form (magnitude vs angle).

    Args:
        S_matrices: List of S-parameter matrices
        params: Which S-parameters to plot (default: S11)
        freqs: Optional frequencies for labeling
        title: Plot title
        figsize: Figure size
        save: Path to save figure
        show: Whether to display

    Returns:
        matplotlib Figure object
    """
    _check_matplotlib()

    if params is None:
        params = [(0, 0)]

    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': 'polar'})

    colors = plt.cm.tab10(np.linspace(0, 1, len(params)))
    for (i, j), color in zip(params, colors):
        label = f"S{i+1}{j+1}"
        S_values = np.array([S[i, j] for S in S_matrices])
        mag = np.abs(S_values)
        phase = np.angle(S_values)

        ax.plot(phase, mag, 'o-', color=color, markersize=4, linewidth=1, label=label)

    ax.set_title(title, pad=15)
    ax.legend(loc='upper right')
    ax.set_rlim(0, 1.1)

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')

    if show:
        plt.show()

    return fig


def plot_return_loss(
    freqs: Union[List[float], np.ndarray],
    S_matrices: List[np.ndarray],
    ports: Optional[List[int]] = None,
    title: str = "Return Loss",
    figsize: Tuple[float, float] = (10, 5),
    save: Optional[str] = None,
    show: bool = True,
    freq_unit: str = "GHz",
) -> "plt.Figure":
    """
    Plot return loss (negative of |Sii| in dB) vs frequency.

    Args:
        freqs: Frequencies in Hz
        S_matrices: S-parameter matrices
        ports: Which ports to show (1-indexed), default all
        title: Plot title
        figsize: Figure size
        save: Path to save figure
        show: Whether to display
        freq_unit: Frequency unit for display

    Returns:
        matplotlib Figure object
    """
    _check_matplotlib()

    freqs = np.asarray(freqs)
    n_ports = S_matrices[0].shape[0]

    if ports is None:
        ports = list(range(1, n_ports + 1))

    freq_scales = {"Hz": 1, "kHz": 1e3, "MHz": 1e6, "GHz": 1e9}
    freq_scale = freq_scales.get(freq_unit, 1e9)
    freqs_scaled = freqs / freq_scale

    fig, ax = plt.subplots(figsize=figsize)

    colors = plt.cm.tab10(np.linspace(0, 1, len(ports)))
    for port, color in zip(ports, colors):
        idx = port - 1  # Convert to 0-indexed
        s_ii = np.array([S[idx, idx] for S in S_matrices])
        rl_db = -20 * np.log10(np.abs(s_ii) + 1e-30)
        ax.plot(freqs_scaled, rl_db, label=f"Port {port}", color=color, linewidth=1.5)

    ax.set_xlabel(f"Frequency ({freq_unit})")
    ax.set_ylabel("Return Loss (dB)")
    ax.set_title(title)
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    # Common reference lines
    ax.axhline(y=10, color='gray', linestyle='--', alpha=0.5, label='10 dB')
    ax.axhline(y=20, color='gray', linestyle=':', alpha=0.5, label='20 dB')

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')

    if show:
        plt.show()

    return fig


def plot_insertion_loss(
    freqs: Union[List[float], np.ndarray],
    S_matrices: List[np.ndarray],
    port_pairs: Optional[List[Tuple[int, int]]] = None,
    title: str = "Insertion Loss",
    figsize: Tuple[float, float] = (10, 5),
    save: Optional[str] = None,
    show: bool = True,
    freq_unit: str = "GHz",
) -> "plt.Figure":
    """
    Plot insertion loss (|Sij| in dB) vs frequency for i != j.

    Args:
        freqs: Frequencies in Hz
        S_matrices: S-parameter matrices
        port_pairs: List of (input_port, output_port) tuples (1-indexed)
                   Default: S21 for 2-port network
        title: Plot title
        figsize: Figure size
        save: Path to save figure
        show: Whether to display
        freq_unit: Frequency unit for display

    Returns:
        matplotlib Figure object
    """
    _check_matplotlib()

    freqs = np.asarray(freqs)
    n_ports = S_matrices[0].shape[0]

    if port_pairs is None:
        if n_ports >= 2:
            port_pairs = [(1, 2)]  # S21 by default
        else:
            raise ValueError("Need at least 2 ports for insertion loss")

    freq_scales = {"Hz": 1, "kHz": 1e3, "MHz": 1e6, "GHz": 1e9}
    freq_scale = freq_scales.get(freq_unit, 1e9)
    freqs_scaled = freqs / freq_scale

    fig, ax = plt.subplots(figsize=figsize)

    colors = plt.cm.tab10(np.linspace(0, 1, len(port_pairs)))
    for (in_port, out_port), color in zip(port_pairs, colors):
        i_idx = out_port - 1
        j_idx = in_port - 1
        s_ij = np.array([S[i_idx, j_idx] for S in S_matrices])
        il_db = 20 * np.log10(np.abs(s_ij) + 1e-30)
        ax.plot(freqs_scaled, il_db,
               label=f"S{out_port}{in_port}", color=color, linewidth=1.5)

    ax.set_xlabel(f"Frequency ({freq_unit})")
    ax.set_ylabel("Insertion Loss (dB)")
    ax.set_title(title)
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    # Reference line at 0 dB (ideal transmission)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')

    if show:
        plt.show()

    return fig

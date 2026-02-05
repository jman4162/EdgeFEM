"""
Radiation pattern plotting functions.

Provides 2D polar/rectangular plots and 3D spherical visualization
for far-field antenna patterns.
"""

import numpy as np
from typing import Optional, List, Union, Tuple

try:
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    import plotly.graph_objects as go
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False


def _check_matplotlib():
    if not HAS_MATPLOTLIB:
        raise ImportError(
            "matplotlib is required for plotting. "
            "Install with: pip install matplotlib"
        )


def _check_plotly():
    if not HAS_PLOTLY:
        raise ImportError(
            "plotly is required for interactive 3D plots. "
            "Install with: pip install plotly"
        )


def plot_pattern_2d(
    pattern: Union[List, np.ndarray],
    title: str = "Radiation Pattern",
    ylabel: str = "Normalized Pattern (dB)",
    xlabel: str = "Angle (degrees)",
    min_db: float = -40.0,
    figsize: Tuple[float, float] = (8, 5),
    save: Optional[str] = None,
    show: bool = True,
    ax: Optional["plt.Axes"] = None,
    **kwargs
) -> Optional["plt.Figure"]:
    """
    Plot 2D radiation pattern in rectangular coordinates.

    Args:
        pattern: List of NTFPoint2D or array with columns [theta_deg, magnitude]
                 If NTFPoint2D list, extracts theta_deg and total E-field magnitude
        title: Plot title
        ylabel: Y-axis label
        xlabel: X-axis label
        min_db: Minimum dB value for clipping (default -40)
        figsize: Figure size in inches (width, height)
        save: Path to save figure (None = don't save)
        show: Whether to display the figure
        ax: Existing axes to plot on (None = create new figure)
        **kwargs: Additional arguments passed to plt.plot()

    Returns:
        Figure object if ax was None, otherwise None
    """
    _check_matplotlib()

    # Convert pattern to arrays
    if hasattr(pattern[0], 'theta_deg'):
        # List of NTFPoint2D
        theta = np.array([p.theta_deg for p in pattern])
        mag = np.array([
            np.sqrt(np.abs(p.e_theta)**2 + np.abs(p.e_phi)**2)
            for p in pattern
        ])
    else:
        # Assume numpy array
        pattern = np.asarray(pattern)
        theta = pattern[:, 0]
        mag = pattern[:, 1]

    # Normalize and convert to dB
    max_mag = np.max(mag)
    if max_mag > 0:
        mag_db = 20 * np.log10(mag / max_mag + 1e-30)
    else:
        mag_db = np.full_like(mag, min_db)

    # Clip to minimum dB
    mag_db = np.clip(mag_db, min_db, 0)

    # Create figure if needed
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
        created_fig = True
    else:
        fig = ax.figure

    # Plot
    ax.plot(theta, mag_db, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(min_db, 0)

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')

    if show and created_fig:
        plt.show()

    return fig if created_fig else None


def plot_pattern_polar(
    pattern: Union[List, np.ndarray],
    title: str = "Radiation Pattern",
    min_db: float = -40.0,
    figsize: Tuple[float, float] = (8, 8),
    save: Optional[str] = None,
    show: bool = True,
    ax: Optional["plt.Axes"] = None,
    **kwargs
) -> Optional["plt.Figure"]:
    """
    Plot 2D radiation pattern in polar coordinates.

    Args:
        pattern: List of NTFPoint2D or array with columns [theta_deg, magnitude]
        title: Plot title
        min_db: Minimum dB value for normalization (default -40)
        figsize: Figure size in inches (width, height)
        save: Path to save figure (None = don't save)
        show: Whether to display the figure
        ax: Existing polar axes to plot on (None = create new)
        **kwargs: Additional arguments passed to plt.plot()

    Returns:
        Figure object if ax was None, otherwise None
    """
    _check_matplotlib()

    # Convert pattern to arrays
    if hasattr(pattern[0], 'theta_deg'):
        theta = np.array([p.theta_deg for p in pattern])
        mag = np.array([
            np.sqrt(np.abs(p.e_theta)**2 + np.abs(p.e_phi)**2)
            for p in pattern
        ])
    else:
        pattern = np.asarray(pattern)
        theta = pattern[:, 0]
        mag = pattern[:, 1]

    # Convert theta to radians
    theta_rad = np.deg2rad(theta)

    # Normalize and convert to dB
    max_mag = np.max(mag)
    if max_mag > 0:
        mag_db = 20 * np.log10(mag / max_mag + 1e-30)
    else:
        mag_db = np.full_like(mag, min_db)

    # Shift so minimum is 0 (for polar plot radius)
    mag_normalized = mag_db - min_db
    mag_normalized = np.clip(mag_normalized, 0, -min_db)

    # Create figure if needed
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': 'polar'})
        created_fig = True
    else:
        fig = ax.figure

    # Plot
    ax.plot(theta_rad, mag_normalized, **kwargs)
    ax.set_title(title, pad=15)

    # Set radial ticks in dB
    r_ticks = np.linspace(0, -min_db, 5)
    ax.set_rticks(r_ticks)
    ax.set_yticklabels([f'{int(t + min_db)}' for t in r_ticks])
    ax.set_rlabel_position(45)

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')

    if show and created_fig:
        plt.show()

    return fig if created_fig else None


def plot_pattern_3d(
    pattern,
    title: str = "3D Radiation Pattern",
    colorscale: str = "Viridis",
    min_db: float = -40.0,
    interactive: bool = True,
    save: Optional[str] = None,
    show: bool = True,
    figsize: Tuple[float, float] = (10, 8),
) -> Optional[Union["plt.Figure", "go.Figure"]]:
    """
    Plot 3D radiation pattern on a spherical surface.

    For interactive=True, uses Plotly for a rotatable 3D view.
    For interactive=False, uses matplotlib for a static view.

    Args:
        pattern: FFPattern3D object from stratton_chu_3d()
        title: Plot title
        colorscale: Color scale name (Viridis, Jet, Hot, etc.)
        min_db: Minimum dB value for color scale (default -40)
        interactive: Use Plotly for interactive 3D (True) or matplotlib (False)
        save: Path to save figure (html for plotly, png for matplotlib)
        show: Whether to display the figure
        figsize: Figure size in inches (for matplotlib only)

    Returns:
        Plotly Figure or matplotlib Figure object
    """
    # Get pattern data
    theta = pattern.theta_grid
    phi = pattern.phi_grid
    pattern_db = pattern.pattern_dB()

    # Clip to minimum
    pattern_db = np.clip(pattern_db, min_db, 0)

    # Convert pattern to radius (normalized 0-1)
    radius = (pattern_db - min_db) / (-min_db)

    # Convert spherical to Cartesian coordinates
    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)

    if interactive:
        _check_plotly()

        # Create Plotly surface
        fig = go.Figure(data=[go.Surface(
            x=x, y=y, z=z,
            surfacecolor=pattern_db,
            colorscale=colorscale,
            cmin=min_db,
            cmax=0,
            colorbar=dict(title='dB', ticksuffix=' dB'),
        )])

        fig.update_layout(
            title=title,
            scene=dict(
                xaxis_title='X',
                yaxis_title='Y',
                zaxis_title='Z',
                aspectmode='data',
            ),
            width=int(figsize[0] * 80),
            height=int(figsize[1] * 80),
        )

        if save:
            if save.endswith('.html'):
                fig.write_html(save)
            else:
                fig.write_image(save)

        if show:
            fig.show()

        return fig

    else:
        _check_matplotlib()
        from mpl_toolkits.mplot3d import Axes3D

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')

        # Plot surface with colormap
        norm = mcolors.Normalize(vmin=min_db, vmax=0)
        cmap = plt.get_cmap(colorscale.lower() if colorscale != 'Viridis' else 'viridis')
        colors = cmap(norm(pattern_db))

        surf = ax.plot_surface(x, y, z, facecolors=colors, shade=False)

        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, shrink=0.5, aspect=10)
        cbar.set_label('Pattern (dB)')

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(title)

        # Equal aspect ratio
        max_range = np.max([x.max() - x.min(), y.max() - y.min(), z.max() - z.min()])
        mid_x = (x.max() + x.min()) / 2
        mid_y = (y.max() + y.min()) / 2
        mid_z = (z.max() + z.min()) / 2
        ax.set_xlim(mid_x - max_range/2, mid_x + max_range/2)
        ax.set_ylim(mid_y - max_range/2, mid_y + max_range/2)
        ax.set_zlim(mid_z - max_range/2, mid_z + max_range/2)

        plt.tight_layout()

        if save:
            fig.savefig(save, dpi=150, bbox_inches='tight')

        if show:
            plt.show()

        return fig


def plot_pattern_3d_cuts(
    pattern,
    title: str = "Pattern Cuts",
    min_db: float = -40.0,
    figsize: Tuple[float, float] = (12, 5),
    save: Optional[str] = None,
    show: bool = True,
) -> "plt.Figure":
    """
    Plot E-plane and H-plane cuts from a 3D pattern.

    Args:
        pattern: FFPattern3D object from stratton_chu_3d()
        title: Overall plot title
        min_db: Minimum dB value for y-axis (default -40)
        figsize: Figure size in inches (width, height)
        save: Path to save figure (None = don't save)
        show: Whether to display the figure

    Returns:
        matplotlib Figure object
    """
    _check_matplotlib()

    pattern_db = pattern.pattern_dB()
    theta = pattern.theta_grid[:, 0]  # First column (theta varies along rows)
    theta_deg = np.rad2deg(theta)

    # Find phi=0 and phi=90 columns
    phi_values = pattern.phi_grid[0, :]  # First row (phi varies along cols)
    phi_deg = np.rad2deg(phi_values)

    # Find indices closest to phi=0 and phi=90
    idx_phi0 = np.argmin(np.abs(phi_deg))
    idx_phi90 = np.argmin(np.abs(phi_deg - 90))

    e_plane = pattern_db[:, idx_phi0]
    h_plane = pattern_db[:, idx_phi90]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # E-plane (phi=0)
    ax1.plot(theta_deg, e_plane, 'b-', linewidth=2)
    ax1.set_xlabel('Theta (degrees)')
    ax1.set_ylabel('Pattern (dB)')
    ax1.set_title(f'E-plane (phi={phi_deg[idx_phi0]:.1f})')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(min_db, 0)
    ax1.set_xlim(0, 180)

    # H-plane (phi=90)
    ax2.plot(theta_deg, h_plane, 'r-', linewidth=2)
    ax2.set_xlabel('Theta (degrees)')
    ax2.set_ylabel('Pattern (dB)')
    ax2.set_title(f'H-plane (phi={phi_deg[idx_phi90]:.1f})')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(min_db, 0)
    ax2.set_xlim(0, 180)

    fig.suptitle(title)
    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')

    if show:
        plt.show()

    return fig

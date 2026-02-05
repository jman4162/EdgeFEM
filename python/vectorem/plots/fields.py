"""
Field visualization plotting functions.

Provides 2D slice plots for electric and magnetic fields
computed from FEM solutions.
"""

import numpy as np
from typing import Optional, Tuple, Union

try:
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def _check_matplotlib():
    if not HAS_MATPLOTLIB:
        raise ImportError(
            "matplotlib is required for plotting. "
            "Install with: pip install matplotlib"
        )


def plot_field_slice(
    x: np.ndarray,
    y: np.ndarray,
    field: np.ndarray,
    title: str = "Field Distribution",
    xlabel: str = "X (mm)",
    ylabel: str = "Y (mm)",
    cmap: str = "RdBu_r",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    show_colorbar: bool = True,
    colorbar_label: str = "E-field (V/m)",
    figsize: Tuple[float, float] = (10, 8),
    save: Optional[str] = None,
    show: bool = True,
    ax: Optional["plt.Axes"] = None,
) -> Optional["plt.Figure"]:
    """
    Plot a 2D field slice as a filled contour plot.

    Args:
        x: X coordinates (1D array or 2D meshgrid)
        y: Y coordinates (1D array or 2D meshgrid)
        field: Field values (2D array, same shape as meshgrid)
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        cmap: Colormap name
        vmin, vmax: Color scale limits (auto if None)
        show_colorbar: Whether to show colorbar
        colorbar_label: Label for colorbar
        figsize: Figure size
        save: Path to save figure
        show: Whether to display
        ax: Existing axes to use (None = create new)

    Returns:
        Figure object if ax was None
    """
    _check_matplotlib()

    # Create meshgrid if needed
    if x.ndim == 1:
        X, Y = np.meshgrid(x, y)
    else:
        X, Y = x, y

    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
        created_fig = True
    else:
        fig = ax.figure

    # Plot filled contours
    cf = ax.contourf(X, Y, field, levels=50, cmap=cmap, vmin=vmin, vmax=vmax)

    if show_colorbar:
        cbar = plt.colorbar(cf, ax=ax)
        cbar.set_label(colorbar_label)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_aspect('equal')

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')

    if show and created_fig:
        plt.show()

    return fig if created_fig else None


def plot_field_magnitude(
    positions: np.ndarray,
    field_values: np.ndarray,
    title: str = "Field Magnitude",
    xlabel: str = "X (mm)",
    ylabel: str = "Y (mm)",
    zlabel: str = "Z (mm)",
    cmap: str = "hot",
    figsize: Tuple[float, float] = (10, 8),
    save: Optional[str] = None,
    show: bool = True,
    projection: str = "xy",
    log_scale: bool = False,
) -> "plt.Figure":
    """
    Plot field magnitude at scattered points.

    Useful for visualizing field values at FEM element centroids
    without requiring structured grid interpolation.

    Args:
        positions: Point positions, shape (N, 3)
        field_values: Field magnitudes at each point, shape (N,)
        title: Plot title
        xlabel, ylabel, zlabel: Axis labels
        cmap: Colormap name
        figsize: Figure size
        save: Path to save figure
        show: Whether to display
        projection: Which 2D projection ("xy", "xz", "yz", or "3d")
        log_scale: Use log scale for field values

    Returns:
        matplotlib Figure object
    """
    _check_matplotlib()

    positions = np.asarray(positions)
    field_values = np.asarray(field_values)

    # Convert to display units (mm)
    pos_mm = positions * 1000

    if log_scale:
        field_plot = np.log10(field_values + 1e-30)
        cbar_label = "log10(|E|)"
    else:
        field_plot = field_values
        cbar_label = "|E| (V/m)"

    if projection == "3d":
        from mpl_toolkits.mplot3d import Axes3D

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')

        sc = ax.scatter(pos_mm[:, 0], pos_mm[:, 1], pos_mm[:, 2],
                       c=field_plot, cmap=cmap, s=10, alpha=0.7)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)

    else:
        fig, ax = plt.subplots(figsize=figsize)

        if projection == "xy":
            x_data = pos_mm[:, 0]
            y_data = pos_mm[:, 1]
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
        elif projection == "xz":
            x_data = pos_mm[:, 0]
            y_data = pos_mm[:, 2]
            ax.set_xlabel(xlabel)
            ax.set_ylabel(zlabel)
        else:  # yz
            x_data = pos_mm[:, 1]
            y_data = pos_mm[:, 2]
            ax.set_xlabel(ylabel)
            ax.set_ylabel(zlabel)

        sc = ax.scatter(x_data, y_data, c=field_plot, cmap=cmap, s=10, alpha=0.7)
        ax.set_aspect('equal')

    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label(cbar_label)

    ax.set_title(title)

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')

    if show:
        plt.show()

    return fig


def plot_field_vector(
    positions: np.ndarray,
    field_vectors: np.ndarray,
    title: str = "Field Vectors",
    xlabel: str = "X (mm)",
    ylabel: str = "Y (mm)",
    figsize: Tuple[float, float] = (10, 8),
    save: Optional[str] = None,
    show: bool = True,
    projection: str = "xy",
    scale: float = 1.0,
    color_by_magnitude: bool = True,
) -> "plt.Figure":
    """
    Plot field vectors as quiver plot.

    Args:
        positions: Point positions, shape (N, 3)
        field_vectors: Field vectors at each point, shape (N, 3) complex
        title: Plot title
        xlabel, ylabel: Axis labels
        figsize: Figure size
        save: Path to save figure
        show: Whether to display
        projection: Which 2D projection ("xy", "xz", "yz")
        scale: Arrow scaling factor
        color_by_magnitude: Color arrows by field magnitude

    Returns:
        matplotlib Figure object
    """
    _check_matplotlib()

    positions = np.asarray(positions)
    field_vectors = np.asarray(field_vectors)

    # Use real part of complex field
    if np.iscomplexobj(field_vectors):
        field_real = field_vectors.real
    else:
        field_real = field_vectors

    # Convert positions to mm
    pos_mm = positions * 1000

    fig, ax = plt.subplots(figsize=figsize)

    if projection == "xy":
        x = pos_mm[:, 0]
        y = pos_mm[:, 1]
        u = field_real[:, 0]
        v = field_real[:, 1]
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    elif projection == "xz":
        x = pos_mm[:, 0]
        y = pos_mm[:, 2]
        u = field_real[:, 0]
        v = field_real[:, 2]
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Z (mm)")
    else:  # yz
        x = pos_mm[:, 1]
        y = pos_mm[:, 2]
        u = field_real[:, 1]
        v = field_real[:, 2]
        ax.set_xlabel(ylabel)
        ax.set_ylabel("Z (mm)")

    if color_by_magnitude:
        mag = np.sqrt(u**2 + v**2)
        quiver = ax.quiver(x, y, u, v, mag, cmap='hot', scale=scale,
                          scale_units='inches', alpha=0.8)
        cbar = plt.colorbar(quiver, ax=ax)
        cbar.set_label("|E| (V/m)")
    else:
        ax.quiver(x, y, u, v, scale=scale, scale_units='inches', alpha=0.8)

    ax.set_title(title)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')

    if show:
        plt.show()

    return fig

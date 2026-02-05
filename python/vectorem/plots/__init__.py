"""
VectorEM Plotting Module

Matplotlib-based plotting utilities for common EM visualization tasks.
All functions are designed to produce publication-quality figures.

Example usage:
    import vectorem as em
    import vectorem.plots as vp

    # Plot 2D radiation pattern
    pattern = em.stratton_chu_2d(...)
    vp.plot_pattern_2d(pattern, title="E-plane", save="pattern.png")

    # Plot S-parameters vs frequency
    vp.plot_sparams_vs_freq(freqs, S_list, save="sweep.png")

    # Plot coupling matrix heatmap
    vp.plot_coupling_matrix(S, save="coupling.png")

    # Plot 3D radiation pattern (requires plotly for interactive)
    pattern_3d = em.stratton_chu_3d(...)
    vp.plot_pattern_3d(pattern_3d, interactive=True)

Dependencies:
    - matplotlib (required)
    - numpy (required)
    - plotly (optional, for 3D interactive plots)
"""

from .patterns import (
    plot_pattern_2d,
    plot_pattern_polar,
    plot_pattern_3d,
    plot_pattern_3d_cuts,
)

from .sparams import (
    plot_sparams_vs_freq,
    plot_sparams_smith,
    plot_sparams_polar,
    plot_return_loss,
    plot_insertion_loss,
)

from .coupling import (
    plot_coupling_matrix,
    plot_coupling_vs_distance,
    plot_active_impedance_scan,
)

from .fields import (
    plot_field_slice,
    plot_field_magnitude,
)

__all__ = [
    # Pattern plots
    "plot_pattern_2d",
    "plot_pattern_polar",
    "plot_pattern_3d",
    "plot_pattern_3d_cuts",
    # S-parameter plots
    "plot_sparams_vs_freq",
    "plot_sparams_smith",
    "plot_sparams_polar",
    "plot_return_loss",
    "plot_insertion_loss",
    # Coupling plots
    "plot_coupling_matrix",
    "plot_coupling_vs_distance",
    "plot_active_impedance_scan",
    # Field plots
    "plot_field_slice",
    "plot_field_magnitude",
]

"""
High-Level Design Classes for VectorEM

This module provides convenient wrappers around the low-level VectorEM API
for common electromagnetic structures.

Available designs:
    - RectWaveguideDesign: Rectangular waveguide S-parameter analysis
    - PatchAntennaDesign: Microstrip patch antenna simulation
    - UnitCellDesign: Periodic unit cell for metasurface/FSS analysis

Example - Rectangular Waveguide:
    from vectorem.designs import RectWaveguideDesign

    # Create WR-90 waveguide (X-band)
    wg = RectWaveguideDesign(a=22.86e-3, b=10.16e-3, length=50e-3)

    # Generate mesh and compute S-parameters
    wg.generate_mesh(density=10)  # 10 elements per wavelength
    S = wg.sparams_at_freq(10e9)

    # Frequency sweep
    freqs, S_list = wg.frequency_sweep(8e9, 12e9, n_points=11)

    # Export to Touchstone
    wg.export_touchstone("waveguide.s2p")

Example - Patch Antenna:
    from vectorem.designs import PatchAntennaDesign

    # Create 2.4 GHz patch on FR-4
    patch = PatchAntennaDesign(
        patch_length=29e-3,
        patch_width=38e-3,
        substrate_height=1.6e-3,
        substrate_eps_r=4.4,
    )
    patch.set_probe_feed()  # Probe feed at optimal location
    patch.generate_mesh(density=15)

    # Compute input impedance and pattern
    Z_in = patch.input_impedance(2.4e9)
    pattern = patch.radiation_pattern(2.4e9)
    D = patch.directivity(2.4e9)

Example - Unit Cell:
    from vectorem.designs import UnitCellDesign

    # Create 5mm x 5mm unit cell
    cell = UnitCellDesign(
        period_x=5e-3,
        period_y=5e-3,
        substrate_height=1e-3,
        substrate_eps_r=3.5,
    )
    cell.add_patch(width=4e-3, length=4e-3)
    cell.generate_mesh(density=20)

    # Reflection coefficient at normal incidence
    R, T = cell.reflection_transmission(10e9)

    # Surface impedance
    Zs = cell.surface_impedance(10e9)
"""

from .waveguide import RectWaveguideDesign
from .patch_antenna import PatchAntennaDesign
from .unit_cell import UnitCellDesign

__all__ = ["RectWaveguideDesign", "PatchAntennaDesign", "UnitCellDesign"]

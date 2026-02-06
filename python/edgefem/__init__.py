"""
EdgeFEM - 3D FEM Electromagnetic Simulator

A research-grade finite-element electromagnetics solver for metasurface,
metamaterial, and phased array unit cell RF modeling (100 kHz–110 GHz).

Example usage:
    import edgefem as em

    # Load mesh and run simulation
    mesh = em.load_gmsh("waveguide.msh")
    params = em.MaxwellParams()
    params.omega = 2 * 3.14159 * 10e9  # 10 GHz

    bc = em.build_edge_pec(mesh, 1)  # PEC on surface tag 1

    S = em.calculate_sparams_eigenmode(mesh, params, bc, ports)

    # Export to Touchstone
    em.write_touchstone_nport("output.s2p", [10e9], [S])

For high-level workflows, see edgefem.designs submodule:
    from edgefem.designs import RectWaveguideDesign

    wg = RectWaveguideDesign(a=22.86e-3, b=10.16e-3, length=50e-3)
    S = wg.sparams_at_freq(10e9)
"""

# Import everything from pyedgefem for convenience
from pyedgefem import *

# Import high-level designs submodule
from . import designs

# Import plotting module (optional - requires matplotlib)
try:
    from . import plots
except ImportError:
    plots = None  # matplotlib not installed

__version__ = "1.0.0"
__all__ = [
    # Re-exported from pyedgefem (documented for IDE support)
    "Mesh", "Node", "Element", "Edge", "ElemType",
    "load_gmsh",
    "BC", "build_scalar_pec", "build_edge_pec",
    "PhysicalTagInfo", "list_physical_tags", "has_surface_tag", "has_volume_tag",
    "MaxwellParams", "MaxwellAssembly", "PortABCType",
    "assemble_maxwell", "assemble_maxwell_periodic",
    "calculate_sparams", "calculate_sparams_periodic", "calculate_sparams_eigenmode",
    "normalize_port_weights",
    "WavePort", "PortMode", "ModePolarization", "RectWaveguidePort",
    "build_wave_port", "solve_port_eigens", "solve_te10_mode",
    "extract_surface_mesh", "populate_te10_field",
    "compute_te_eigenvector", "build_wave_port_from_eigenvector",
    "SolveOptions", "SolveResult", "solve_linear",
    "CouplingMatrix", "compute_coupling", "mutual_coupling_db",
    "active_reflection", "is_passive", "check_passivity", "is_lossless", "is_reciprocal",
    "PeriodicBC", "PeriodicPair", "build_periodic_pairs", "validate_periodic_bc",
    "set_floquet_phase", "floquet_phase_from_angle",
    "ActiveImpedanceResult", "uniform_scan_excitation", "compute_active_impedance",
    "active_impedance_scan", "active_vswr", "find_scan_blindness",
    "CouplingCoefficient", "CouplingStats",
    "extract_all_couplings", "extract_nearest_neighbor_coupling",
    "coupling_vs_separation", "compute_coupling_stats",
    "FFPoint2D", "EmbeddedPattern", "ArrayPatternConfig",
    "compute_embedded_patterns", "compute_single_embedded_pattern",
    "to_db", "normalize_pattern",
    "ArrayPackage", "export_array_package_json",
    "export_patterns_csv", "export_coupling_csv", "export_scan_csv",
    "NTFPoint2D", "stratton_chu_2d", "write_pattern_csv",
    "FFPattern3D", "stratton_chu_3d", "compute_directivity", "compute_max_gain",
    "compute_hpbw", "write_pattern_3d_csv", "write_pattern_3d_vtk",
    "VTKExportOptions", "export_fields_vtk", "export_fields_vtk_legacy",
    "TouchstoneFormat", "TouchstoneOptions",
    "write_touchstone", "write_touchstone_nport", "touchstone_extension",
    "MeshQualityReport", "validate_mesh",
    # Dispersive materials
    "materials",
    # High-level API
    "designs",
    "plots",
]

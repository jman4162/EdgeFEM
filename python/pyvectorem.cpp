#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "vectorem/bc.hpp"
#include "vectorem/coupling.hpp"
#include "vectorem/fem.hpp"
#include "vectorem/io/touchstone.hpp"
#include "vectorem/mesh.hpp"
#include "vectorem/mesh_quality.hpp"
#include "vectorem/periodic.hpp"
#include "vectorem/ports/port_eigensolve.hpp"
#include "vectorem/ports/wave_port.hpp"
#include "vectorem/solver.hpp"
#include "vectorem/maxwell.hpp"
#include "vectorem/array/active_impedance.hpp"
#include "vectorem/array/coupling_extract.hpp"

namespace py = pybind11;
using namespace vectorem;

PYBIND11_MAKE_OPAQUE(SpMatC);
PYBIND11_MAKE_OPAQUE(VecC);

PYBIND11_MODULE(pyvectorem, m) {
  m.doc() = "Python bindings for the VectorEM 3D FEM solver.";

  // Node structure
  py::class_<Node>(m, "Node", "Mesh node with position.")
      .def(py::init<>())
      .def_readonly("id", &Node::id, "Node ID")
      .def_readonly("xyz", &Node::xyz, "Node position (x, y, z)");

  // Element structure
  py::class_<Element>(m, "Element", "Mesh element (triangle or tetrahedron).")
      .def(py::init<>())
      .def_readonly("id", &Element::id, "Element ID")
      .def_readonly("type", &Element::type, "Element type")
      .def_readonly("conn", &Element::conn, "Node connectivity")
      .def_readonly("phys", &Element::phys, "Physical group tag");

  // Element type enum
  py::enum_<ElemType>(m, "ElemType")
      .value("Tri3", ElemType::Tri3)
      .value("Tet4", ElemType::Tet4);

  py::class_<Mesh>(m, "Mesh", "Represents a 3D mesh.")
      .def(py::init<>())
      .def_readonly("nodes", &Mesh::nodes, "List of mesh nodes.")
      .def_readonly("tets", &Mesh::tets, "List of tetrahedra.")
      .def_readonly("tris", &Mesh::tris, "List of triangles.")
      .def_readonly("edges", &Mesh::edges, "List of edges.")
      .def("num_nodes", [](const Mesh &m) { return m.nodes.size(); })
      .def("num_tets", [](const Mesh &m) { return m.tets.size(); })
      .def("num_tris", [](const Mesh &m) { return m.tris.size(); })
      .def("num_edges", [](const Mesh &m) { return m.edges.size(); });

  // Edge structure
  py::class_<Edge>(m, "Edge", "Mesh edge connecting two nodes.")
      .def(py::init<>())
      .def_readonly("n0", &Edge::n0, "First node ID")
      .def_readonly("n1", &Edge::n1, "Second node ID");

  m.def("load_gmsh", &load_gmsh_v2, "Loads a mesh from a Gmsh v2 .msh file.");

  py::class_<BC>(m, "BC", "Boundary condition definitions.")
      .def(py::init<>());

  m.def("build_scalar_pec", &build_scalar_pec, "Build scalar PEC BC",
        py::arg("mesh"), py::arg("pec_tag"));
  m.def("build_edge_pec", &build_edge_pec, "Build edge-based PEC BC",
        py::arg("mesh"), py::arg("pec_tag"));

  py::class_<ScalarHelmholtzParams>(m, "ScalarHelmholtzParams", "Parameters for the scalar Helmholtz equation.")
      .def(py::init<>())
      .def_readwrite("k0", &ScalarHelmholtzParams::k0, "Wavenumber")
      .def_readwrite("eps_r", &ScalarHelmholtzParams::eps_r, "Relative permittivity");

  py::class_<SpMatC>(m, "SpMatC", "Sparse matrix with complex coefficients.");
  py::class_<VecC>(m, "VecC", "Vector with complex coefficients.");

  py::class_<ScalarAssembly>(m, "ScalarAssembly", "Assembled scalar system (A, b).")
      .def_property_readonly(
          "A", [](ScalarAssembly &s) -> SpMatC & { return s.A; },
          py::return_value_policy::reference_internal, "System matrix A.")
      .def_property_readonly(
          "b", [](ScalarAssembly &s) -> VecC & { return s.b; },
          py::return_value_policy::reference_internal, "RHS vector b.");

  m.def("assemble_scalar_helmholtz", &assemble_scalar_helmholtz,
        "Assembles the scalar Helmholtz system.", py::arg("mesh"), py::arg("params"),
        py::arg("bc"));

  py::class_<MaxwellParams>(m, "MaxwellParams", "Parameters for Maxwell's equations.")
      .def(py::init<>())
      .def_readwrite("omega", &MaxwellParams::omega, "Angular frequency (rad/s).")
      .def_readwrite("eps_r", &MaxwellParams::eps_r, "Complex relative permittivity.")
      .def_readwrite("mu_r", &MaxwellParams::mu_r, "Complex relative permeability.")
      .def_readwrite("pml_sigma", &MaxwellParams::pml_sigma, "PML conductivity.")
      .def_readwrite("pml_regions", &MaxwellParams::pml_regions, "Set of physical region tags for PML.")
      .def_readwrite("use_abc", &MaxwellParams::use_abc, "Flag to enable absorbing boundary condition.");

  py::class_<MaxwellAssembly>(m, "MaxwellAssembly", "Assembled Maxwell system (A, b).")
      .def_property_readonly(
          "A", [](MaxwellAssembly &s) -> SpMatC & { return s.A; },
          py::return_value_policy::reference_internal, "System matrix A.")
      .def_property_readonly(
          "b", [](MaxwellAssembly &s) -> VecC & { return s.b; },
          py::return_value_policy::reference_internal, "RHS vector b.");

  py::enum_<ModePolarization>(m, "ModePolarization")
      .value("TE", ModePolarization::TE)
      .value("TM", ModePolarization::TM);

  py::class_<PortMode>(m, "PortMode", "Modal solution for a port cross-section.")
      .def(py::init<>())
      .def_readwrite("pol", &PortMode::pol)
      .def_readwrite("fc", &PortMode::fc)
      .def_readwrite("kc", &PortMode::kc)
      .def_readwrite("omega", &PortMode::omega)
      .def_readwrite("eps", &PortMode::eps)
      .def_readwrite("mu", &PortMode::mu)
      .def_readwrite("beta", &PortMode::beta)
      .def_readwrite("Z0", &PortMode::Z0)
      .def_readwrite("field", &PortMode::field);

  py::class_<PortSurfaceMesh>(m, "PortSurfaceMesh", "Surface mesh extracted from a volume port.")
      .def(py::init<>())
      .def_property_readonly(
          "mesh",
          [](PortSurfaceMesh &s) -> Mesh & { return s.mesh; },
          py::return_value_policy::reference_internal)
      .def_readonly("volume_tri_indices", &PortSurfaceMesh::volume_tri_indices);

  py::class_<WavePort>(m, "WavePort", "Wave port definition with modal data.")
      .def(py::init<>())
      .def_readwrite("surface_tag", &WavePort::surface_tag)
      .def_readwrite("mode", &WavePort::mode)
      .def_readwrite("edges", &WavePort::edges)
      .def_readwrite("weights", &WavePort::weights);

  m.def("extract_surface_mesh", &extract_surface_mesh,
        "Extract a surface mesh for a port.", py::arg("mesh"),
        py::arg("surface_tag"));

  m.def("build_wave_port", &build_wave_port, "Project a modal field onto port edges.",
        py::arg("volume_mesh"), py::arg("surface"), py::arg("mode"));

  m.def("assemble_maxwell", &assemble_maxwell, "Assembles the 3D Maxwell system.",
        py::arg("mesh"), py::arg("params"), py::arg("bc"),
        py::arg("ports") = std::vector<WavePort>(),
        py::arg("active_port_idx") = -1);

  m.def("calculate_sparams", &calculate_sparams,
        "Calculates S-parameters for a set of wave ports.",
        py::arg("mesh"), py::arg("params"), py::arg("bc"), py::arg("ports"));

  py::class_<SolveOptions>(m, "SolveOptions", "Options for the linear solver.")
      .def(py::init<>())
      .def_readwrite("use_bicgstab", &SolveOptions::use_bicgstab, "Use BiCGSTAB solver.")
      .def_readwrite("tolerance", &SolveOptions::tolerance, "Convergence tolerance.")
      .def_readwrite("max_iterations", &SolveOptions::max_iterations, "Maximum iterations.");

  py::class_<SolveResult>(m, "SolveResult", "Results from a linear solve.")
      .def_readonly("method", &SolveResult::method, "Solver method used.")
      .def_readonly("iters", &SolveResult::iters, "Number of iterations.")
      .def_readonly("residual", &SolveResult::residual, "Final residual error.")
      .def_readonly("converged", &SolveResult::converged, "Whether the solver converged.")
      .def_readonly("error_message", &SolveResult::error_message, "Error message if not converged.")
      .def_property_readonly(
          "x", [](SolveResult &s) -> VecC & { return s.x; },
          py::return_value_policy::reference_internal, "Solution vector x.");

  m.def("solve_linear", &solve_linear, "Solves a complex linear system Ax=b.",
        py::arg("A"), py::arg("b"), py::arg("options") = SolveOptions());

  py::class_<SParams2>(m, "SParams2", "2-port S-parameter data.")
      .def(py::init<>())
      .def_readwrite("s11", &SParams2::s11)
      .def_readwrite("s21", &SParams2::s21)
      .def_readwrite("s12", &SParams2::s12)
      .def_readwrite("s22", &SParams2::s22);

  m.def("write_touchstone", &write_touchstone, "Writes S-parameters to a Touchstone file.",
        py::arg("path"), py::arg("freq"), py::arg("data"));

  // Touchstone format enum
  py::enum_<TouchstoneFormat>(m, "TouchstoneFormat")
      .value("RI", TouchstoneFormat::RI, "Real/Imaginary format")
      .value("MA", TouchstoneFormat::MA, "Magnitude/Angle format")
      .value("DB", TouchstoneFormat::DB, "dB/Angle format");

  // Touchstone options
  py::class_<TouchstoneOptions>(m, "TouchstoneOptions", "Options for Touchstone export.")
      .def(py::init<>())
      .def_readwrite("format", &TouchstoneOptions::format, "Output format (RI, MA, or DB).")
      .def_readwrite("z0", &TouchstoneOptions::z0, "Reference impedance in ohms.");

  m.def("write_touchstone_nport", &write_touchstone_nport,
        "Writes N-port S-parameters to a Touchstone file.",
        py::arg("path"), py::arg("freq"), py::arg("S_matrices"),
        py::arg("opts") = TouchstoneOptions());

  m.def("touchstone_extension", &touchstone_extension,
        "Get appropriate Touchstone file extension for given port count.",
        py::arg("num_ports"));

  // Mesh quality validation
  py::class_<MeshQualityReport>(m, "MeshQualityReport", "Mesh quality validation report.")
      .def(py::init<>())
      .def_readonly("num_degenerate_tets", &MeshQualityReport::num_degenerate_tets)
      .def_readonly("num_inverted_tets", &MeshQualityReport::num_inverted_tets)
      .def_readonly("min_volume", &MeshQualityReport::min_volume)
      .def_readonly("max_aspect_ratio", &MeshQualityReport::max_aspect_ratio)
      .def_readonly("bad_tet_indices", &MeshQualityReport::bad_tet_indices)
      .def("is_valid", &MeshQualityReport::is_valid, "Returns true if mesh has no bad elements.");

  m.def("validate_mesh", &validate_mesh, "Validate mesh quality.",
        py::arg("mesh"));

  // Coupling matrix
  py::class_<CouplingMatrix>(m, "CouplingMatrix", "Coupling matrix representation.")
      .def(py::init<>())
      .def_readonly("Z", &CouplingMatrix::Z, "Impedance matrix.")
      .def_readonly("Y", &CouplingMatrix::Y, "Admittance matrix.")
      .def_readonly("S", &CouplingMatrix::S, "Scattering matrix.")
      .def_readonly("frequency", &CouplingMatrix::frequency, "Operating frequency (Hz).")
      .def_readonly("z0", &CouplingMatrix::z0, "Reference impedance (ohms).")
      .def_readonly("num_ports", &CouplingMatrix::num_ports, "Number of ports.");

  m.def("compute_coupling", &compute_coupling,
        "Compute impedance/admittance matrices from S-parameters.",
        py::arg("S"), py::arg("freq"), py::arg("z0") = 50.0);

  m.def("mutual_coupling_db", &mutual_coupling_db,
        "Convert S-parameters to mutual coupling in dB.",
        py::arg("cm"));

  m.def("active_reflection", &active_reflection,
        "Compute active reflection coefficients for given excitation.",
        py::arg("S"), py::arg("excitation"));

  m.def("is_passive", &is_passive,
        "Check if network is passive.",
        py::arg("S"), py::arg("tolerance") = 1e-6);

  m.def("is_lossless", &is_lossless,
        "Check if network is lossless.",
        py::arg("S"), py::arg("tolerance") = 1e-6);

  m.def("is_reciprocal", &is_reciprocal,
        "Check if network is reciprocal.",
        py::arg("S"), py::arg("tolerance") = 1e-6);

  // Periodic boundary conditions
  py::class_<PeriodicPair>(m, "PeriodicPair", "Edge pair for periodic BC.")
      .def(py::init<>())
      .def_readwrite("master_edge", &PeriodicPair::master_edge)
      .def_readwrite("slave_edge", &PeriodicPair::slave_edge)
      .def_readwrite("master_orient", &PeriodicPair::master_orient)
      .def_readwrite("slave_orient", &PeriodicPair::slave_orient)
      .def_readwrite("translation", &PeriodicPair::translation);

  py::class_<PeriodicBC>(m, "PeriodicBC", "Periodic boundary condition specification.")
      .def(py::init<>())
      .def_readwrite("pairs", &PeriodicBC::pairs)
      .def_readwrite("period_vector", &PeriodicBC::period_vector)
      .def_readwrite("phase_shift", &PeriodicBC::phase_shift);

  m.def("build_periodic_pairs", &build_periodic_pairs,
        "Build periodic edge pairs from mesh surface tags.",
        py::arg("mesh"), py::arg("master_tag"), py::arg("slave_tag"),
        py::arg("period_vector"), py::arg("tolerance") = 1e-9);

  m.def("validate_periodic_bc", &validate_periodic_bc,
        "Validate periodic BC pairs.",
        py::arg("mesh"), py::arg("pbc"));

  m.def("set_floquet_phase", &set_floquet_phase,
        "Set Floquet phase shift for transverse wave vector.",
        py::arg("pbc"), py::arg("k_transverse"));

  m.def("floquet_phase_from_angle", &floquet_phase_from_angle,
        "Compute Floquet phase from incidence angles.",
        py::arg("period_vector"), py::arg("theta"), py::arg("phi"), py::arg("k0"));

  m.def("count_surface_edges", &count_surface_edges,
        "Count edges on surface with given tag.",
        py::arg("mesh"), py::arg("surface_tag"));

  m.def("assemble_maxwell_periodic", &assemble_maxwell_periodic,
        "Assemble Maxwell system with periodic BC.",
        py::arg("mesh"), py::arg("params"), py::arg("bc"), py::arg("pbc"),
        py::arg("ports") = std::vector<WavePort>(),
        py::arg("active_port_idx") = -1);

  m.def("calculate_sparams_periodic", &calculate_sparams_periodic,
        "Calculate S-parameters with periodic BC.",
        py::arg("mesh"), py::arg("params"), py::arg("bc"), py::arg("pbc"),
        py::arg("ports"));

  // Active impedance
  py::class_<ActiveImpedanceResult>(m, "ActiveImpedanceResult",
                                     "Active impedance calculation result.")
      .def(py::init<>())
      .def_readwrite("theta", &ActiveImpedanceResult::theta)
      .def_readwrite("phi", &ActiveImpedanceResult::phi)
      .def_readwrite("Z_active", &ActiveImpedanceResult::Z_active)
      .def_readwrite("Gamma_active", &ActiveImpedanceResult::Gamma_active)
      .def_readwrite("scan_loss_db", &ActiveImpedanceResult::scan_loss_db);

  m.def("uniform_scan_excitation", &uniform_scan_excitation,
        "Generate uniform phase excitation for scanning.",
        py::arg("num_elements"), py::arg("positions"),
        py::arg("theta"), py::arg("phi"), py::arg("k0"));

  m.def("compute_active_impedance", &compute_active_impedance,
        "Compute active impedance for given excitation.",
        py::arg("coupling"), py::arg("excitation"));

  m.def("active_impedance_scan", &active_impedance_scan,
        "Compute active impedance over scan angle sweep.",
        py::arg("coupling"), py::arg("positions"),
        py::arg("theta_rad"), py::arg("phi_rad"), py::arg("k0"));

  m.def("active_vswr", &active_vswr,
        "Compute VSWR from active reflection coefficient.",
        py::arg("Gamma_active"));

  // Coupling extraction
  py::class_<CouplingCoefficient>(m, "CouplingCoefficient",
                                   "Coupling coefficient between ports.")
      .def(py::init<>())
      .def_readonly("port_i", &CouplingCoefficient::port_i)
      .def_readonly("port_j", &CouplingCoefficient::port_j)
      .def_readonly("S_ij", &CouplingCoefficient::S_ij)
      .def_readonly("magnitude_db", &CouplingCoefficient::magnitude_db)
      .def_readonly("phase_deg", &CouplingCoefficient::phase_deg)
      .def_readonly("distance", &CouplingCoefficient::distance);

  m.def("extract_all_couplings", &extract_all_couplings,
        "Extract all coupling coefficients from S-matrix.",
        py::arg("S"), py::arg("positions") = std::vector<Eigen::Vector3d>());

  py::class_<CouplingStats>(m, "CouplingStats", "Coupling statistics.")
      .def(py::init<>())
      .def_readonly("max_coupling_db", &CouplingStats::max_coupling_db)
      .def_readonly("min_coupling_db", &CouplingStats::min_coupling_db)
      .def_readonly("mean_coupling_db", &CouplingStats::mean_coupling_db)
      .def_readonly("std_coupling_db", &CouplingStats::std_coupling_db)
      .def_readonly("max_coupling_i", &CouplingStats::max_coupling_i)
      .def_readonly("max_coupling_j", &CouplingStats::max_coupling_j);

  m.def("compute_coupling_stats", &compute_coupling_stats,
        "Compute coupling statistics on S-matrix.",
        py::arg("S"));
}

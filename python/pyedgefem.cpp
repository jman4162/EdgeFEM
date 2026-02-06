#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <iostream>

#include "edgefem/array/active_impedance.hpp"
#include "edgefem/array/coupling_extract.hpp"
#include "edgefem/array/embedded_pattern.hpp"
#include "edgefem/bc.hpp"
#include "edgefem/coupling.hpp"
#include "edgefem/fem.hpp"
#include "edgefem/io/array_export.hpp"
#include "edgefem/io/touchstone.hpp"
#include "edgefem/io/vtk_export.hpp"
#include "edgefem/materials/dispersive.hpp"
#include "edgefem/maxwell.hpp"
#include "edgefem/mesh.hpp"
#include "edgefem/mesh_quality.hpp"
#include "edgefem/periodic.hpp"
#include "edgefem/ports/port_eigensolve.hpp"
#include "edgefem/ports/wave_port.hpp"
#include "edgefem/post/ntf.hpp"
#include "edgefem/solver.hpp"

namespace py = pybind11;
using namespace edgefem;

PYBIND11_MAKE_OPAQUE(SpMatC);
PYBIND11_MAKE_OPAQUE(VecC);

PYBIND11_MODULE(pyedgefem, m) {
  m.doc() = "Python bindings for the EdgeFEM 3D FEM solver.";

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
      .def_readonly("phys", &Element::phys, "Physical group tag")
      .def_readonly("edges", &Element::edges, "Edge indices")
      .def_readonly("edge_orient", &Element::edge_orient, "Edge orientations");

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
      .def_readonly("nodeIndex", &Mesh::nodeIndex, "Map from node ID to index.")
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
      .def(py::init<>())
      .def_property_readonly("dirichlet_nodes",
                             [](const BC &b) {
                               return std::vector<int>(
                                   b.dirichlet_nodes.begin(),
                                   b.dirichlet_nodes.end());
                             })
      .def_property_readonly("dirichlet_edges", [](const BC &b) {
        return std::vector<int>(b.dirichlet_edges.begin(),
                                b.dirichlet_edges.end());
      });

  m.def("build_scalar_pec", &build_scalar_pec,
        "Build scalar PEC BC from mesh surface tag.\n"
        "Prints a warning if no nodes are found (check tag value).",
        py::arg("mesh"), py::arg("pec_tag"));
  m.def("build_edge_pec", &build_edge_pec,
        "Build edge-based PEC BC from mesh surface tag.\n"
        "Prints a warning if no edges are found (check tag value).",
        py::arg("mesh"), py::arg("pec_tag"));

  // Physical tag helpers
  py::class_<PhysicalTagInfo>(m, "PhysicalTagInfo",
                              "Information about physical tags in a mesh.")
      .def(py::init<>())
      .def_readonly("volume_tags", &PhysicalTagInfo::volume_tags,
                    "Tags from tetrahedra (3D regions)")
      .def_readonly("surface_tags", &PhysicalTagInfo::surface_tags,
                    "Tags from triangles (2D surfaces)");

  m.def("list_physical_tags", &list_physical_tags,
        "List all unique physical tags in the mesh.\n"
        "Returns PhysicalTagInfo with volume_tags and surface_tags.",
        py::arg("mesh"));

  m.def("has_surface_tag", &has_surface_tag,
        "Check if a surface tag exists in the mesh.", py::arg("mesh"),
        py::arg("tag"));

  m.def("has_volume_tag", &has_volume_tag,
        "Check if a volume tag exists in the mesh.", py::arg("mesh"),
        py::arg("tag"));

  py::class_<ScalarHelmholtzParams>(
      m, "ScalarHelmholtzParams",
      "Parameters for the scalar Helmholtz equation.")
      .def(py::init<>())
      .def_readwrite("k0", &ScalarHelmholtzParams::k0, "Wavenumber")
      .def_readwrite("eps_r", &ScalarHelmholtzParams::eps_r,
                     "Relative permittivity");

  py::class_<SpMatC>(m, "SpMatC", "Sparse matrix with complex coefficients.")
      .def_property_readonly(
          "shape",
          [](const SpMatC &m) { return py::make_tuple(m.rows(), m.cols()); })
      .def("coeff",
           [](const SpMatC &m, int row, int col) { return m.coeff(row, col); })
      .def("to_dense",
           [](const SpMatC &m) {
             Eigen::MatrixXcd dense = m;
             py::array_t<std::complex<double>> result({m.rows(), m.cols()});
             auto r = result.mutable_unchecked<2>();
             for (Eigen::Index i = 0; i < m.rows(); ++i) {
               for (Eigen::Index j = 0; j < m.cols(); ++j) {
                 r(i, j) = dense(i, j);
               }
             }
             return result;
           })
      .def("nnz", [](const SpMatC &m) { return m.nonZeros(); });

  py::class_<VecC>(m, "VecC", "Vector with complex coefficients.")
      .def("__len__", [](const VecC &v) { return v.size(); })
      .def("__getitem__",
           [](const VecC &v, int i) {
             if (i < 0)
               i += v.size();
             if (i < 0 || i >= v.size())
               throw py::index_error();
             return v(i);
           })
      .def("__getitem__",
           [](const VecC &v, const std::vector<int> &indices) {
             py::array_t<std::complex<double>> result(indices.size());
             auto r = result.mutable_unchecked<1>();
             for (size_t k = 0; k < indices.size(); ++k) {
               int i = indices[k];
               if (i < 0 || i >= v.size())
                 throw py::index_error();
               r(k) = v(i);
             }
             return result;
           })
      .def("to_numpy", [](const VecC &v) {
        py::array_t<std::complex<double>> result(v.size());
        auto r = result.mutable_unchecked<1>();
        for (Eigen::Index i = 0; i < v.size(); ++i) {
          r(i) = v(i);
        }
        return result;
      });

  py::class_<ScalarAssembly>(m, "ScalarAssembly",
                             "Assembled scalar system (A, b).")
      .def_property_readonly(
          "A", [](ScalarAssembly &s) -> SpMatC & { return s.A; },
          py::return_value_policy::reference_internal, "System matrix A.")
      .def_property_readonly(
          "b", [](ScalarAssembly &s) -> VecC & { return s.b; },
          py::return_value_policy::reference_internal, "RHS vector b.");

  m.def("assemble_scalar_helmholtz", &assemble_scalar_helmholtz,
        "Assembles the scalar Helmholtz system.", py::arg("mesh"),
        py::arg("params"), py::arg("bc"));

  py::enum_<PortABCType>(
      m, "PortABCType", "ABC formulation variants for port boundary treatment.")
      .value("None", PortABCType::None, "No port ABC")
      .value("Beta", PortABCType::Beta, "j*beta (propagation constant)")
      .value("BetaNorm", PortABCType::BetaNorm, "j*beta/k0 (dimensionless)")
      .value("ImpedanceMatch", PortABCType::ImpedanceMatch,
             "j*beta*sqrt(Z0/eta0)")
      .value("ModalAdmittance", PortABCType::ModalAdmittance,
             "j*omega*eps0/Z0");

  py::class_<MaxwellParams>(m, "MaxwellParams",
                            "Parameters for Maxwell's equations.")
      .def(py::init<>())
      .def_readwrite("omega", &MaxwellParams::omega,
                     "Angular frequency (rad/s).")
      .def_readwrite("eps_r", &MaxwellParams::eps_r,
                     "Complex relative permittivity (global default).")
      .def_readwrite("mu_r", &MaxwellParams::mu_r,
                     "Complex relative permeability (global default).")
      .def_readwrite(
          "eps_r_regions", &MaxwellParams::eps_r_regions,
          "Region-based permittivity: maps physical tag to complex eps_r.\n"
          "Elements with unassigned tags use global eps_r.")
      .def_readwrite(
          "mu_r_regions", &MaxwellParams::mu_r_regions,
          "Region-based permeability: maps physical tag to complex mu_r.\n"
          "Elements with unassigned tags use global mu_r.")
      .def_property(
          "eps_models", [](const MaxwellParams &p) { return p.eps_models; },
          [](MaxwellParams &p,
             const std::unordered_map<
                 int, std::shared_ptr<materials::DispersiveMaterial>> &m) {
            p.eps_models = m;
          },
          "Dispersive permittivity models: maps physical tag to "
          "DispersiveMaterial.\n"
          "Takes precedence over eps_r_regions. Evaluated at current omega.")
      .def(
          "set_eps_model",
          [](MaxwellParams &p, int tag,
             std::shared_ptr<materials::DispersiveMaterial> model) {
            p.eps_models[tag] = model;
          },
          py::arg("phys_tag"), py::arg("model"),
          "Set a dispersive permittivity model for a physical region.")
      .def_property(
          "mu_models", [](const MaxwellParams &p) { return p.mu_models; },
          [](MaxwellParams &p,
             const std::unordered_map<
                 int, std::shared_ptr<materials::DispersiveMaterial>> &m) {
            p.mu_models = m;
          },
          "Dispersive permeability models: maps physical tag to "
          "DispersiveMaterial.\n"
          "Takes precedence over mu_r_regions. Evaluated at current omega.")
      .def(
          "set_mu_model",
          [](MaxwellParams &p, int tag,
             std::shared_ptr<materials::DispersiveMaterial> model) {
            p.mu_models[tag] = model;
          },
          py::arg("phys_tag"), py::arg("model"),
          "Set a dispersive permeability model for a physical region.")
      .def("get_eps_r",
           py::overload_cast<int>(&MaxwellParams::get_eps_r, py::const_),
           "Get effective permittivity for a physical region tag "
           "(non-dispersive).",
           py::arg("phys_tag"))
      .def(
          "get_eps_r",
          py::overload_cast<int, double>(&MaxwellParams::get_eps_r, py::const_),
          "Get effective permittivity for a physical region tag at given "
          "omega.\n"
          "Evaluates dispersive models if present.",
          py::arg("phys_tag"), py::arg("omega"))
      .def("get_mu_r",
           py::overload_cast<int>(&MaxwellParams::get_mu_r, py::const_),
           "Get effective permeability for a physical region tag "
           "(non-dispersive).",
           py::arg("phys_tag"))
      .def("get_mu_r",
           py::overload_cast<int, double>(&MaxwellParams::get_mu_r, py::const_),
           "Get effective permeability for a physical region tag at given "
           "omega.\n"
           "Evaluates dispersive models if present.",
           py::arg("phys_tag"), py::arg("omega"))
      .def_readwrite("pml_sigma", &MaxwellParams::pml_sigma,
                     "PML conductivity.")
      .def_readwrite("pml_regions", &MaxwellParams::pml_regions,
                     "Set of physical region tags for PML.")
      .def_readwrite("use_abc", &MaxwellParams::use_abc,
                     "Flag to enable absorbing boundary condition.")
      .def_readwrite("use_port_abc", &MaxwellParams::use_port_abc,
                     "Flag to enable port-specific ABC.")
      .def_readwrite("port_abc_type", &MaxwellParams::port_abc_type,
                     "Port ABC formulation type.")
      .def_readwrite("port_weight_scale", &MaxwellParams::port_weight_scale,
                     "Port weight scaling factor.")
      .def_readwrite("port_abc_scale", &MaxwellParams::port_abc_scale,
                     "ABC coefficient scale (0.5 optimal).")
      .def_readwrite("use_eigenmode_excitation",
                     &MaxwellParams::use_eigenmode_excitation,
                     "Use eigenmode excitation.");

  py::class_<MaxwellAssembly>(m, "MaxwellAssembly",
                              "Assembled Maxwell system (A, b).")
      .def_property_readonly(
          "A", [](MaxwellAssembly &s) -> SpMatC & { return s.A; },
          py::return_value_policy::reference_internal, "System matrix A.")
      .def_property_readonly(
          "b", [](MaxwellAssembly &s) -> VecC & { return s.b; },
          py::return_value_policy::reference_internal, "RHS vector b.");

  py::enum_<ModePolarization>(m, "ModePolarization")
      .value("TE", ModePolarization::TE)
      .value("TM", ModePolarization::TM);

  py::class_<PortMode>(m, "PortMode",
                       "Modal solution for a port cross-section.")
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

  py::class_<RectWaveguidePort>(m, "RectWaveguidePort",
                                "Rectangular waveguide port dimensions.")
      .def(py::init<>())
      .def(py::init([](double a, double b) {
             RectWaveguidePort p;
             p.a = a;
             p.b = b;
             return p;
           }),
           py::arg("a"), py::arg("b"))
      .def_readwrite("a", &RectWaveguidePort::a, "Width in meters")
      .def_readwrite("b", &RectWaveguidePort::b, "Height in meters");

  m.def("solve_port_eigens", &solve_port_eigens,
        "Solve 2D eigenmode problem on port cross-section.", py::arg("mesh"),
        py::arg("num_modes"), py::arg("omega"), py::arg("eps_r"),
        py::arg("mu_r"), py::arg("pol"));

  m.def("solve_te10_mode", &solve_te10_mode,
        "Compute analytical TE10 mode for rectangular waveguide.",
        py::arg("port"), py::arg("freq"));

  m.def("straight_waveguide_sparams", &straight_waveguide_sparams,
        "Compute analytical S-parameters for straight waveguide section.",
        py::arg("port"), py::arg("length"), py::arg("freq"));

  py::class_<PortSurfaceMesh>(m, "PortSurfaceMesh",
                              "Surface mesh extracted from a volume port.")
      .def(py::init<>())
      .def_property_readonly(
          "mesh", [](PortSurfaceMesh &s) -> Mesh & { return s.mesh; },
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

  m.def("build_wave_port", &build_wave_port,
        "Project a modal field onto port edges.", py::arg("volume_mesh"),
        py::arg("surface"), py::arg("mode"));

  m.def("populate_te10_field", &populate_te10_field,
        "Populate TE10 mode field on a surface mesh using analytical formula.",
        py::arg("surface"), py::arg("port"), py::arg("mode"));

  m.def("assemble_maxwell", &assemble_maxwell,
        "Assembles the 3D Maxwell system.", py::arg("mesh"), py::arg("params"),
        py::arg("bc"), py::arg("ports") = std::vector<WavePort>(),
        py::arg("active_port_idx") = -1);

  m.def("calculate_sparams", &calculate_sparams,
        "Calculates S-parameters for a set of wave ports.", py::arg("mesh"),
        py::arg("params"), py::arg("bc"), py::arg("ports"));

  m.def("normalize_port_weights", &normalize_port_weights,
        "Normalize port weights so w^H A^-1 w = Z0 for proper S-parameter "
        "extraction.",
        py::arg("mesh"), py::arg("params"), py::arg("bc"), py::arg("ports"));

  m.def("calculate_sparams_eigenmode", &calculate_sparams_eigenmode,
        "Calculate S-parameters using direct eigenmode excitation (more "
        "accurate).",
        py::arg("mesh"), py::arg("params"), py::arg("bc"), py::arg("ports"));

  m.def(
      "compute_te_eigenvector", &compute_te_eigenvector,
      "Compute 3D FEM eigenvector for a TE mode with given cutoff wavenumber.",
      py::arg("mesh"), py::arg("pec_edges"), py::arg("target_kc_sq"));

  m.def("build_wave_port_from_eigenvector", &build_wave_port_from_eigenvector,
        "Build wave port using 3D FEM eigenvector as weights.", py::arg("mesh"),
        py::arg("surface"), py::arg("eigenvector"), py::arg("mode"),
        py::arg("pec_edges"));

  py::class_<SolveOptions>(m, "SolveOptions", "Options for the linear solver.")
      .def(py::init<>())
      .def_readwrite("use_bicgstab", &SolveOptions::use_bicgstab,
                     "Use BiCGSTAB solver.")
      .def_readwrite("tolerance", &SolveOptions::tolerance,
                     "Convergence tolerance.")
      .def_readwrite("max_iterations", &SolveOptions::max_iterations,
                     "Maximum iterations.")
      .def_readwrite("verbose", &SolveOptions::verbose,
                     "Print progress to stderr.");

  py::class_<SolveResult>(m, "SolveResult", "Results from a linear solve.")
      .def_readonly("method", &SolveResult::method, "Solver method used.")
      .def_readonly("iters", &SolveResult::iters, "Number of iterations.")
      .def_readonly("residual", &SolveResult::residual, "Final residual error.")
      .def_readonly("converged", &SolveResult::converged,
                    "Whether the solver converged.")
      .def_readonly("error_message", &SolveResult::error_message,
                    "Error message if not converged.")
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

  m.def("write_touchstone", &write_touchstone,
        "Writes S-parameters to a Touchstone file.", py::arg("path"),
        py::arg("freq"), py::arg("data"));

  // Touchstone format enum
  py::enum_<TouchstoneFormat>(m, "TouchstoneFormat")
      .value("RI", TouchstoneFormat::RI, "Real/Imaginary format")
      .value("MA", TouchstoneFormat::MA, "Magnitude/Angle format")
      .value("DB", TouchstoneFormat::DB, "dB/Angle format");

  // Touchstone options
  py::class_<TouchstoneOptions>(m, "TouchstoneOptions",
                                "Options for Touchstone export.")
      .def(py::init<>())
      .def_readwrite("format", &TouchstoneOptions::format,
                     "Output format (RI, MA, or DB).")
      .def_readwrite("z0", &TouchstoneOptions::z0,
                     "Reference impedance in ohms.");

  m.def("write_touchstone_nport", &write_touchstone_nport,
        "Writes N-port S-parameters to a Touchstone file.", py::arg("path"),
        py::arg("freq"), py::arg("S_matrices"),
        py::arg("opts") = TouchstoneOptions());

  m.def("touchstone_extension", &touchstone_extension,
        "Get appropriate Touchstone file extension for given port count.",
        py::arg("num_ports"));

  // Mesh quality validation
  py::class_<MeshQualityReport>(m, "MeshQualityReport",
                                "Mesh quality validation report.")
      .def(py::init<>())
      .def_readonly("num_degenerate_tets",
                    &MeshQualityReport::num_degenerate_tets)
      .def_readonly("num_inverted_tets", &MeshQualityReport::num_inverted_tets)
      .def_readonly("min_volume", &MeshQualityReport::min_volume)
      .def_readonly("max_aspect_ratio", &MeshQualityReport::max_aspect_ratio)
      .def_readonly("bad_tet_indices", &MeshQualityReport::bad_tet_indices)
      .def("is_valid", &MeshQualityReport::is_valid,
           "Returns true if mesh has no bad elements.");

  m.def("validate_mesh", &validate_mesh, "Validate mesh quality.",
        py::arg("mesh"));

  // Coupling matrix
  py::class_<CouplingMatrix>(m, "CouplingMatrix",
                             "Coupling matrix representation.")
      .def(py::init<>())
      .def_readonly("Z", &CouplingMatrix::Z, "Impedance matrix.")
      .def_readonly("Y", &CouplingMatrix::Y, "Admittance matrix.")
      .def_readonly("S", &CouplingMatrix::S, "Scattering matrix.")
      .def_readonly("frequency", &CouplingMatrix::frequency,
                    "Operating frequency (Hz).")
      .def_readonly("z0", &CouplingMatrix::z0, "Reference impedance (ohms).")
      .def_readonly("num_ports", &CouplingMatrix::num_ports,
                    "Number of ports.");

  m.def("compute_coupling", &compute_coupling,
        "Compute impedance/admittance matrices from S-parameters.",
        py::arg("S"), py::arg("freq"), py::arg("z0") = 50.0);

  m.def("mutual_coupling_db", &mutual_coupling_db,
        "Convert S-parameters to mutual coupling in dB.", py::arg("cm"));

  m.def("active_reflection", &active_reflection,
        "Compute active reflection coefficients for given excitation.",
        py::arg("S"), py::arg("excitation"));

  m.def("is_passive", &is_passive, "Check if network is passive.", py::arg("S"),
        py::arg("tolerance") = 1e-6);

  // check_passivity: convenience function that prints warning if not passive
  m.def(
      "check_passivity",
      [](const Eigen::MatrixXcd &S, double tolerance) {
        bool passive = is_passive(S, tolerance);
        if (!passive) {
          // Compute max singular value to show how much it violates passivity
          Eigen::JacobiSVD<Eigen::MatrixXcd> svd(S);
          double max_sv = svd.singularValues()(0);
          std::cerr << "WARNING: S-parameters violate passivity! "
                    << "Max singular value = " << max_sv
                    << " (should be <= " << (1.0 + tolerance) << ")\n";
          // Also check individual elements for 2-port case
          if (S.rows() == 2 && S.cols() == 2) {
            double s11_mag = std::abs(S(0, 0));
            double s21_mag = std::abs(S(1, 0));
            double s12_mag = std::abs(S(0, 1));
            double s22_mag = std::abs(S(1, 1));
            double power_sum1 = s11_mag * s11_mag + s21_mag * s21_mag;
            double power_sum2 = s12_mag * s12_mag + s22_mag * s22_mag;
            std::cerr << "  |S11|=" << s11_mag << ", |S21|=" << s21_mag
                      << ", |S11|^2+|S21|^2=" << power_sum1 << "\n";
            std::cerr << "  |S12|=" << s12_mag << ", |S22|=" << s22_mag
                      << ", |S12|^2+|S22|^2=" << power_sum2 << "\n";
          }
        }
        return passive;
      },
      "Check passivity and print warning if violated.", py::arg("S"),
      py::arg("tolerance") = 1e-6);

  m.def("is_lossless", &is_lossless, "Check if network is lossless.",
        py::arg("S"), py::arg("tolerance") = 1e-6);

  m.def("is_reciprocal", &is_reciprocal, "Check if network is reciprocal.",
        py::arg("S"), py::arg("tolerance") = 1e-6);

  // Periodic boundary conditions
  py::class_<PeriodicPair>(m, "PeriodicPair", "Edge pair for periodic BC.")
      .def(py::init<>())
      .def_readwrite("master_edge", &PeriodicPair::master_edge)
      .def_readwrite("slave_edge", &PeriodicPair::slave_edge)
      .def_readwrite("master_orient", &PeriodicPair::master_orient)
      .def_readwrite("slave_orient", &PeriodicPair::slave_orient)
      .def_readwrite("translation", &PeriodicPair::translation);

  py::class_<PeriodicBC>(m, "PeriodicBC",
                         "Periodic boundary condition specification.")
      .def(py::init<>())
      .def_readwrite("pairs", &PeriodicBC::pairs)
      .def_readwrite("period_vector", &PeriodicBC::period_vector)
      .def_readwrite("phase_shift", &PeriodicBC::phase_shift);

  m.def("build_periodic_pairs", &build_periodic_pairs,
        "Build periodic edge pairs from mesh surface tags.", py::arg("mesh"),
        py::arg("master_tag"), py::arg("slave_tag"), py::arg("period_vector"),
        py::arg("tolerance") = 1e-9);

  m.def("validate_periodic_bc", &validate_periodic_bc,
        "Validate periodic BC pairs.", py::arg("mesh"), py::arg("pbc"));

  m.def("set_floquet_phase", &set_floquet_phase,
        "Set Floquet phase shift for transverse wave vector.", py::arg("pbc"),
        py::arg("k_transverse"));

  m.def("floquet_phase_from_angle", &floquet_phase_from_angle,
        "Compute Floquet phase from incidence angles.",
        py::arg("period_vector"), py::arg("theta"), py::arg("phi"),
        py::arg("k0"));

  m.def("count_surface_edges", &count_surface_edges,
        "Count edges on surface with given tag.", py::arg("mesh"),
        py::arg("surface_tag"));

  m.def("assemble_maxwell_periodic", &assemble_maxwell_periodic,
        "Assemble Maxwell system with periodic BC.", py::arg("mesh"),
        py::arg("params"), py::arg("bc"), py::arg("pbc"),
        py::arg("ports") = std::vector<WavePort>(),
        py::arg("active_port_idx") = -1);

  m.def("calculate_sparams_periodic", &calculate_sparams_periodic,
        "Calculate S-parameters with periodic BC.", py::arg("mesh"),
        py::arg("params"), py::arg("bc"), py::arg("pbc"), py::arg("ports"));

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
        py::arg("num_elements"), py::arg("positions"), py::arg("theta"),
        py::arg("phi"), py::arg("k0"));

  m.def("compute_active_impedance", &compute_active_impedance,
        "Compute active impedance for given excitation.", py::arg("coupling"),
        py::arg("excitation"));

  m.def("active_impedance_scan", &active_impedance_scan,
        "Compute active impedance over scan angle sweep.", py::arg("coupling"),
        py::arg("positions"), py::arg("theta_rad"), py::arg("phi_rad"),
        py::arg("k0"));

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
        "Extract all coupling coefficients from S-matrix.", py::arg("S"),
        py::arg("positions") = std::vector<Eigen::Vector3d>());

  py::class_<CouplingStats>(m, "CouplingStats", "Coupling statistics.")
      .def(py::init<>())
      .def_readonly("max_coupling_db", &CouplingStats::max_coupling_db)
      .def_readonly("min_coupling_db", &CouplingStats::min_coupling_db)
      .def_readonly("mean_coupling_db", &CouplingStats::mean_coupling_db)
      .def_readonly("std_coupling_db", &CouplingStats::std_coupling_db)
      .def_readonly("max_coupling_i", &CouplingStats::max_coupling_i)
      .def_readonly("max_coupling_j", &CouplingStats::max_coupling_j);

  m.def("compute_coupling_stats", &compute_coupling_stats,
        "Compute coupling statistics on S-matrix.", py::arg("S"));

  // ============================================================
  // Additional Array Analysis Functions
  // ============================================================

  // find_scan_blindness
  m.def(
      "find_scan_blindness", &find_scan_blindness,
      "Find scan blindness angles where active VSWR exceeds threshold.\n"
      "Returns a map of element index to vector of blindness angles (radians).",
      py::arg("scan_results"), py::arg("vswr_threshold") = 10.0);

  // Coupling vs separation analysis
  m.def("coupling_vs_separation", &coupling_vs_separation,
        "Compute coupling magnitude vs. element separation distance.\n"
        "Returns map from distance bin (mm) to average coupling (dB).",
        py::arg("S"), py::arg("positions"));

  // Extract nearest neighbor coupling
  m.def("extract_nearest_neighbor_coupling", &extract_nearest_neighbor_coupling,
        "Extract mutual coupling between nearest neighbors only.", py::arg("S"),
        py::arg("positions"), py::arg("max_neighbors") = 4);

  // ============================================================
  // Embedded Pattern Computation
  // ============================================================

  // FFPoint2D structure (from embedded_pattern.hpp)
  py::class_<FFPoint2D>(m, "FFPoint2D", "Far-field point in 2D pattern cut.")
      .def(py::init<>())
      .def_readwrite("theta", &FFPoint2D::theta, "Angle in radians")
      .def_readwrite("magnitude", &FFPoint2D::magnitude,
                     "Field magnitude (linear)")
      .def_readwrite("phase", &FFPoint2D::phase, "Phase in radians");

  // EmbeddedPattern structure
  py::class_<EmbeddedPattern>(m, "EmbeddedPattern",
                              "Embedded element pattern for one port.")
      .def(py::init<>())
      .def_readwrite("port_index", &EmbeddedPattern::port_index,
                     "Index of the active port")
      .def_readwrite("pattern_phi0", &EmbeddedPattern::pattern_phi0,
                     "E-plane cut (phi=0)")
      .def_readwrite("pattern_phi90", &EmbeddedPattern::pattern_phi90,
                     "H-plane cut (phi=90)")
      .def_readwrite("frequency", &EmbeddedPattern::frequency,
                     "Operating frequency (Hz)")
      .def_readwrite("element_position", &EmbeddedPattern::element_position,
                     "Element position");

  // ArrayPatternConfig structure
  py::class_<ArrayPatternConfig>(
      m, "ArrayPatternConfig",
      "Configuration for embedded pattern computation.")
      .def(py::init<>())
      .def_readwrite("theta_rad", &ArrayPatternConfig::theta_rad,
                     "Theta angles for pattern cuts")
      .def_readwrite("phi_e_plane", &ArrayPatternConfig::phi_e_plane,
                     "Phi angle for E-plane (usually 0)")
      .def_readwrite("phi_h_plane", &ArrayPatternConfig::phi_h_plane,
                     "Phi angle for H-plane (usually pi/2)")
      .def_readwrite("k0", &ArrayPatternConfig::k0,
                     "Free-space wavenumber (rad/m)")
      .def_readwrite("huygens_radius", &ArrayPatternConfig::huygens_radius,
                     "Radius for Huygens surface (0=auto)");

  // Compute embedded patterns (all ports)
  m.def("compute_embedded_patterns", &compute_embedded_patterns,
        "Compute embedded element patterns for all ports.\n"
        "Includes mutual coupling effects from other elements.",
        py::arg("mesh"), py::arg("params"), py::arg("bc"), py::arg("ports"),
        py::arg("config"));

  // Compute single embedded pattern
  m.def("compute_single_embedded_pattern", &compute_single_embedded_pattern,
        "Compute embedded pattern for a single port.\n"
        "Other ports are properly loaded.",
        py::arg("mesh"), py::arg("params"), py::arg("bc"), py::arg("ports"),
        py::arg("active_port"), py::arg("config"));

  // Utility functions
  m.def("to_db", &to_db, "Convert linear magnitude to dB.", py::arg("linear"));

  m.def("normalize_pattern", &normalize_pattern,
        "Normalize pattern to peak value of 1.0.", py::arg("pattern"));

  // ============================================================
  // Array Export Functions
  // ============================================================

  // ArrayPackage structure
  py::class_<ArrayPackage>(m, "ArrayPackage",
                           "Complete array characterization package.")
      .def(py::init<>())
      .def_readwrite("model_name", &ArrayPackage::model_name, "Model name")
      .def_readwrite("element_positions", &ArrayPackage::element_positions,
                     "Element positions")
      .def_readwrite("frequencies", &ArrayPackage::frequencies,
                     "Frequencies (Hz)")
      .def_readwrite("S_matrices", &ArrayPackage::S_matrices,
                     "S-matrices per frequency")
      .def_readwrite("embedded_patterns", &ArrayPackage::embedded_patterns,
                     "Embedded patterns [freq][port]")
      .def_readwrite("scan_data", &ArrayPackage::scan_data, "Scan data")
      .def_readwrite("reference_impedance", &ArrayPackage::reference_impedance,
                     "Reference impedance (ohms)");

  // Export functions
  m.def("export_array_package_json", &export_array_package_json,
        "Export complete array package to JSON format.\n"
        "Compatible with Phased-Array-Antenna-Model.",
        py::arg("path"), py::arg("pkg"));

  m.def("export_patterns_csv", &export_patterns_csv,
        "Export embedded element patterns to CSV.\n"
        "Columns: port_index, theta_deg, phi0_mag, phi0_phase, phi90_mag, "
        "phi90_phase",
        py::arg("path"), py::arg("patterns"));

  m.def("export_coupling_csv", &export_coupling_csv,
        "Export coupling matrix to CSV.\n"
        "Columns: port_i, port_j, S_real, S_imag, magnitude_db, phase_deg, "
        "distance",
        py::arg("path"), py::arg("S"),
        py::arg("positions") = std::vector<Eigen::Vector3d>());

  m.def("export_scan_csv", &export_scan_csv,
        "Export scan angle sweep results to CSV.\n"
        "Columns: theta_deg, phi_deg, element_idx, Gamma_real, Gamma_imag, "
        "Z_real, Z_imag, VSWR",
        py::arg("path"), py::arg("scan_results"));

  // ============================================================
  // Near-to-Far Field (Stratton-Chu)
  // ============================================================

  // NTFPoint2D from ntf.hpp (renamed from FFPoint2D to avoid conflict)
  py::class_<NTFPoint2D>(
      m, "NTFPoint2D",
      "Far-field point from near-to-far transformation.\n"
      "Contains E_theta and E_phi complex far-field components.")
      .def(py::init<>())
      .def_readwrite("theta_deg", &NTFPoint2D::theta_deg,
                     "Polar angle in degrees")
      .def_readwrite("e_theta", &NTFPoint2D::e_theta,
                     "E_theta far-field component")
      .def_readwrite("e_phi", &NTFPoint2D::e_phi, "E_phi far-field component");

  m.def("stratton_chu_2d", &stratton_chu_2d,
        "Compute far-field pattern using Stratton-Chu integral.\n"
        "The surface is represented by sample points with normals, E/H fields, "
        "and areas.",
        py::arg("r"), py::arg("n"), py::arg("E"), py::arg("H"), py::arg("area"),
        py::arg("theta_rad"), py::arg("phi_rad"), py::arg("k0"));

  m.def("write_pattern_csv", &write_pattern_csv,
        "Write a 2D NTF pattern to CSV for visualization.", py::arg("path"),
        py::arg("phi_deg"), py::arg("pattern"));

  // ============================================================
  // 3D Far-Field Patterns
  // ============================================================

  py::class_<FFPattern3D>(
      m, "FFPattern3D",
      "3D far-field pattern over spherical (theta, phi) grid.\n"
      "Matrices are Ntheta x Nphi where theta varies along rows, phi along "
      "columns.")
      .def(py::init<>())
      .def_readonly("theta_grid", &FFPattern3D::theta_grid,
                    "Theta angles in radians (Ntheta x Nphi)")
      .def_readonly("phi_grid", &FFPattern3D::phi_grid,
                    "Phi angles in radians (Ntheta x Nphi)")
      .def_readonly("E_theta", &FFPattern3D::E_theta,
                    "E_theta far-field component (complex)")
      .def_readonly("E_phi", &FFPattern3D::E_phi,
                    "E_phi far-field component (complex)")
      .def("total_magnitude", &FFPattern3D::total_magnitude,
           "Get total electric field magnitude at each point")
      .def("power_pattern", &FFPattern3D::power_pattern,
           "Get total power pattern (|E_theta|^2 + |E_phi|^2)")
      .def("pattern_dB", &FFPattern3D::pattern_dB,
           "Get pattern in dB relative to maximum (normalized)");

  m.def("stratton_chu_3d", &stratton_chu_3d,
        "Compute full 3D far-field pattern using Stratton-Chu integral.\n"
        "Returns pattern sampled over spherical (theta, phi) grid.\n\n"
        "Args:\n"
        "    r: Surface sample point positions\n"
        "    n: Outward surface normals at each point\n"
        "    E: Tangential electric field at each point\n"
        "    H: Tangential magnetic field at each point\n"
        "    area: Surface patch area at each point\n"
        "    theta_rad: Vector of theta angles (0 to pi)\n"
        "    phi_rad: Vector of phi angles (0 to 2*pi)\n"
        "    k0: Free-space wavenumber (2*pi/lambda)",
        py::arg("r"), py::arg("n"), py::arg("E"), py::arg("H"), py::arg("area"),
        py::arg("theta_rad"), py::arg("phi_rad"), py::arg("k0"));

  m.def("compute_directivity", &compute_directivity,
        "Compute directivity from 3D far-field pattern.\n"
        "D = 4*pi * U_max / P_rad where U is radiation intensity.\n\n"
        "Returns: Directivity (dimensionless, linear scale)",
        py::arg("pattern"));

  m.def("compute_max_gain", &compute_max_gain,
        "Compute maximum gain from 3D pattern with radiation efficiency.\n"
        "G_max = efficiency * D_max\n\n"
        "Args:\n"
        "    pattern: 3D far-field pattern\n"
        "    efficiency: Radiation efficiency (0 to 1, default 1.0)\n\n"
        "Returns: Maximum gain (dimensionless, linear scale)",
        py::arg("pattern"), py::arg("efficiency") = 1.0);

  m.def("compute_hpbw", &compute_hpbw,
        "Compute half-power beamwidth in E-plane and H-plane.\n\n"
        "Returns: Tuple of (E-plane HPBW, H-plane HPBW) in degrees",
        py::arg("pattern"));

  m.def("write_pattern_3d_csv", &write_pattern_3d_csv,
        "Export 3D pattern to CSV file.\n"
        "Columns: theta_deg, phi_deg, E_theta_mag, E_theta_phase, E_phi_mag, "
        "E_phi_phase, total_dB",
        py::arg("path"), py::arg("pattern"));

  m.def("write_pattern_3d_vtk", &write_pattern_3d_vtk,
        "Export 3D pattern to VTK file for ParaView visualization.\n"
        "Creates a spherical surface colored by pattern magnitude.",
        py::arg("path"), py::arg("pattern"), py::arg("radius") = 1.0);

  // ============================================================
  // VTK Field Export for ParaView Visualization
  // ============================================================

  py::class_<VTKExportOptions>(m, "VTKExportOptions",
                               "Options for VTK field export.")
      .def(py::init<>())
      .def_readwrite("export_magnitude", &VTKExportOptions::export_magnitude,
                     "Export E-field magnitude (default: true)")
      .def_readwrite("export_phase", &VTKExportOptions::export_phase,
                     "Export E-field phase in radians (default: true)")
      .def_readwrite("export_real", &VTKExportOptions::export_real,
                     "Export E-field real part (default: false)")
      .def_readwrite("export_imag", &VTKExportOptions::export_imag,
                     "Export E-field imaginary part (default: false)")
      .def_readwrite("export_vector", &VTKExportOptions::export_vector,
                     "Export E-field as 3D vector (default: true)");

  m.def("export_fields_vtk", &export_fields_vtk,
        "Export E-field solution to VTK XML Unstructured Grid (.vtu) format.\n"
        "The field is computed at element centroids from edge solution.\n"
        "Open the output file in ParaView to visualize fields.",
        py::arg("mesh"), py::arg("solution"), py::arg("filename"),
        py::arg("options") = VTKExportOptions());

  m.def("export_fields_vtk_legacy", &export_fields_vtk_legacy,
        "Export E-field solution to legacy VTK format (.vtk).\n"
        "Compatible with older visualization tools.",
        py::arg("mesh"), py::arg("solution"), py::arg("filename"),
        py::arg("options") = VTKExportOptions());

  // ============================================================
  // Dispersive Materials Submodule
  // ============================================================

  py::module_ materials = m.def_submodule(
      "materials",
      "Frequency-dependent (dispersive) material models.\n\n"
      "Supported models:\n"
      "  - DebyeMaterial: Polar dielectrics with single relaxation time\n"
      "  - LorentzMaterial: Resonant materials with oscillator poles\n"
      "  - DrudeMaterial: Free-electron metals and plasmas\n"
      "  - DrudeLorentzMaterial: Combined model for realistic metals");

  using namespace edgefem::materials;

  // Base class (for type hints and polymorphism)
  py::class_<DispersiveMaterial, std::shared_ptr<DispersiveMaterial>>(
      materials, "DispersiveMaterial",
      "Base class for dispersive material models.\n"
      "Use subclasses DebyeMaterial, LorentzMaterial, or DrudeMaterial.")
      .def("eval_eps", &DispersiveMaterial::eval_eps,
           "Evaluate complex permittivity at angular frequency omega (rad/s).",
           py::arg("omega"))
      .def("eval_mu", &DispersiveMaterial::eval_mu,
           "Evaluate complex permeability at angular frequency omega (rad/s).\n"
           "Returns 1.0 for non-magnetic materials (default).",
           py::arg("omega"));

  // Debye relaxation model
  py::class_<DebyeMaterial, DispersiveMaterial, std::shared_ptr<DebyeMaterial>>(
      materials, "DebyeMaterial",
      "Debye relaxation model for polar dielectrics.\n\n"
      "Formula: eps(omega) = eps_inf + (eps_s - eps_inf) / (1 + "
      "j*omega*tau)\n\n"
      "Example (water at 20C):\n"
      "    water = DebyeMaterial(eps_static=80.1, eps_inf=4.9, tau=9.3e-12)")
      .def(py::init<double, double, double>(), py::arg("eps_static"),
           py::arg("eps_inf"), py::arg("tau"),
           "Create Debye material.\n\n"
           "Args:\n"
           "    eps_static: Static (DC) relative permittivity\n"
           "    eps_inf: High-frequency relative permittivity\n"
           "    tau: Relaxation time in seconds (must be > 0)")
      .def("eval_eps", &DebyeMaterial::eval_eps, py::arg("omega"))
      .def_property_readonly("eps_static", &DebyeMaterial::eps_static,
                             "Static (DC) permittivity")
      .def_property_readonly("eps_inf", &DebyeMaterial::eps_inf,
                             "High-frequency permittivity")
      .def_property_readonly("tau", &DebyeMaterial::tau,
                             "Relaxation time constant (seconds)");

  // Lorentz oscillator model
  py::class_<LorentzMaterial, DispersiveMaterial,
             std::shared_ptr<LorentzMaterial>>(
      materials, "LorentzMaterial",
      "Multi-pole Lorentz oscillator model for resonant materials.\n\n"
      "Formula: eps(omega) = eps_inf + Sum_k [ delta_eps_k * omega0_k^2 /\n"
      "                                       (omega0_k^2 - omega^2 + "
      "j*gamma_k*omega) ]\n\n"
      "Example:\n"
      "    mat = LorentzMaterial(eps_inf=1.0)\n"
      "    mat.add_pole(delta_eps=5.0, omega0=2*pi*10e9, gamma=2*pi*0.5e9)")
      .def(py::init<>(), "Create Lorentz material with eps_inf=1.0")
      .def(py::init<double>(), py::arg("eps_inf"),
           "Create Lorentz material with specified eps_inf")
      .def("add_pole", &LorentzMaterial::add_pole, py::arg("delta_eps"),
           py::arg("omega0"), py::arg("gamma"),
           "Add a Lorentz pole.\n\n"
           "Args:\n"
           "    delta_eps: Oscillator strength (permittivity contribution)\n"
           "    omega0: Resonance angular frequency (rad/s, must be > 0)\n"
           "    gamma: Damping factor (rad/s, must be >= 0)")
      .def("eval_eps", &LorentzMaterial::eval_eps, py::arg("omega"))
      .def_property_readonly("eps_inf", &LorentzMaterial::eps_inf,
                             "High-frequency permittivity")
      .def_property_readonly("num_poles", &LorentzMaterial::num_poles,
                             "Number of Lorentz poles");

  // Drude free-electron model
  py::class_<DrudeMaterial, DispersiveMaterial, std::shared_ptr<DrudeMaterial>>(
      materials, "DrudeMaterial",
      "Drude free-electron model for metals and plasmas.\n\n"
      "Formula: eps(omega) = 1 - omega_p^2 / (omega^2 + j*gamma*omega)\n\n"
      "Example (gold-like):\n"
      "    gold = DrudeMaterial(omega_p=1.37e16, gamma=4.05e13)")
      .def(py::init<double, double>(), py::arg("omega_p"), py::arg("gamma"),
           "Create Drude material.\n\n"
           "Args:\n"
           "    omega_p: Plasma frequency in rad/s (must be > 0)\n"
           "    gamma: Collision/damping frequency in rad/s (must be >= 0)")
      .def("eval_eps", &DrudeMaterial::eval_eps, py::arg("omega"))
      .def_property_readonly("omega_p", &DrudeMaterial::omega_p,
                             "Plasma frequency (rad/s)")
      .def_property_readonly("gamma", &DrudeMaterial::gamma,
                             "Collision/damping frequency (rad/s)");

  // Drude-Lorentz combined model
  py::class_<DrudeLorentzMaterial, DispersiveMaterial,
             std::shared_ptr<DrudeLorentzMaterial>>(
      materials, "DrudeLorentzMaterial",
      "Combined Drude-Lorentz model for realistic metals.\n\n"
      "Combines free-electron (Drude) and bound-electron (Lorentz) "
      "contributions.\n\n"
      "Example:\n"
      "    gold = DrudeLorentzMaterial(eps_inf=1.0, omega_p=1.37e16, "
      "gamma_d=4.05e13)\n"
      "    gold.add_lorentz_pole(delta_eps=2.0, omega0=4e15, gamma=1e14)")
      .def(py::init<double, double, double>(), py::arg("eps_inf"),
           py::arg("omega_p"), py::arg("gamma_d"),
           "Create Drude-Lorentz material.\n\n"
           "Args:\n"
           "    eps_inf: High-frequency permittivity (background)\n"
           "    omega_p: Plasma frequency in rad/s\n"
           "    gamma_d: Drude damping in rad/s")
      .def("add_lorentz_pole", &DrudeLorentzMaterial::add_lorentz_pole,
           py::arg("delta_eps"), py::arg("omega0"), py::arg("gamma"),
           "Add a Lorentz pole for interband transitions.")
      .def("eval_eps", &DrudeLorentzMaterial::eval_eps, py::arg("omega"))
      .def_property_readonly("eps_inf", &DrudeLorentzMaterial::eps_inf)
      .def_property_readonly("omega_p", &DrudeLorentzMaterial::omega_p)
      .def_property_readonly("gamma_d", &DrudeLorentzMaterial::gamma_d);
}

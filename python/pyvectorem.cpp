#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "vectorem/bc.hpp"
#include "vectorem/fem.hpp"
#include "vectorem/io/touchstone.hpp"
#include "vectorem/mesh.hpp"
#include "vectorem/ports/port_eigensolve.hpp"
#include "vectorem/ports/wave_port.hpp"
#include "vectorem/solver.hpp"
#include "vectorem/maxwell.hpp"

namespace py = pybind11;
using namespace vectorem;

PYBIND11_MAKE_OPAQUE(SpMatC);
PYBIND11_MAKE_OPAQUE(VecC);

PYBIND11_MODULE(pyvectorem, m) {
  m.doc() = "Python bindings for the VectorEM 3D FEM solver.";

  py::class_<Mesh>(m, "Mesh", "Represents a 3D mesh.")
      .def(py::init<>())
      .def_readonly("nodes", &Mesh::nodes, "List of mesh nodes.")
      .def_readonly("tets", &Mesh::tets)
      .def_readonly("tris", &Mesh::tris);

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
      .def_readwrite("use_bicgstab", &SolveOptions::use_bicgstab, "Use BiCGSTAB solver.");

  py::class_<SolveResult>(m, "SolveResult", "Results from a linear solve.")
      .def_readonly("method", &SolveResult::method, "Solver method used.")
      .def_readonly("iters", &SolveResult::iters, "Number of iterations.")
      .def_readonly("residual", &SolveResult::residual, "Final residual error.")
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
}

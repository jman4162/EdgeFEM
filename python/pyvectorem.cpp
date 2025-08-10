#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "vectorem/bc.hpp"
#include "vectorem/fem.hpp"
#include "vectorem/io/touchstone.hpp"
#include "vectorem/mesh.hpp"
#include "vectorem/ports/port_eigensolve.hpp"
#include "vectorem/solver.hpp"
#include "vectorem/maxwell.hpp"

namespace py = pybind11;
using namespace vectorem;

PYBIND11_MAKE_OPAQUE(SpMatC);
PYBIND11_MAKE_OPAQUE(VecC);

PYBIND11_MODULE(pyvectorem, m) {
  m.doc() = "pybind11 bindings for VectorEM";

  py::class_<Mesh>(m, "Mesh")
      .def(py::init<>())
      .def_readonly("nodes", &Mesh::nodes)
      .def_readonly("tets", &Mesh::tets)
      .def_readonly("tris", &Mesh::tris);

  m.def("load_gmsh", &load_gmsh_v2, "Load Gmsh v2 mesh");

  py::class_<BC>(m, "BC").def(py::init<>());
  m.def("build_scalar_pec", &build_scalar_pec, "Build PEC BC", py::arg("mesh"),
        py::arg("pec_tag"));
  m.def("build_edge_pec", &build_edge_pec, "Build edge PEC BC", py::arg("mesh"),
        py::arg("pec_tag"));

  py::class_<ScalarHelmholtzParams>(m, "ScalarHelmholtzParams")
      .def(py::init<>())
      .def_readwrite("k0", &ScalarHelmholtzParams::k0)
      .def_readwrite("eps_r", &ScalarHelmholtzParams::eps_r);

  py::class_<SpMatC>(m, "SpMatC");
  py::class_<VecC>(m, "VecC");

  py::class_<ScalarAssembly>(m, "ScalarAssembly")
      .def_property_readonly(
          "A", [](ScalarAssembly &s) -> SpMatC & { return s.A; },
          py::return_value_policy::reference_internal)
      .def_property_readonly(
          "b", [](ScalarAssembly &s) -> VecC & { return s.b; },
          py::return_value_policy::reference_internal);

  m.def("assemble_scalar_helmholtz", &assemble_scalar_helmholtz,
        "Assemble scalar Helmholtz", py::arg("mesh"), py::arg("params"),
        py::arg("bc"));

  py::class_<MaxwellParams>(m, "MaxwellParams")
      .def(py::init<>())
      .def_readwrite("omega", &MaxwellParams::omega)
      .def_readwrite("eps_r", &MaxwellParams::eps_r)
      .def_readwrite("mu_r", &MaxwellParams::mu_r)
      .def_readwrite("pml_sigma", &MaxwellParams::pml_sigma)
      .def_readwrite("pml_regions", &MaxwellParams::pml_regions)
      .def_readwrite("use_abc", &MaxwellParams::use_abc);

  py::class_<MaxwellAssembly>(m, "MaxwellAssembly")
      .def_property_readonly(
          "A", [](MaxwellAssembly &s) -> SpMatC & { return s.A; },
          py::return_value_policy::reference_internal)
      .def_property_readonly(
          "b", [](MaxwellAssembly &s) -> VecC & { return s.b; },
          py::return_value_policy::reference_internal);

  m.def("assemble_maxwell", &assemble_maxwell, "Assemble Maxwell", py::arg("mesh"),
        py::arg("params"), py::arg("bc"));

  py::class_<SolveOptions>(m, "SolveOptions")
      .def(py::init<>())
      .def_readwrite("use_bicgstab", &SolveOptions::use_bicgstab);

  py::class_<SolveResult>(m, "SolveResult")
      .def_readonly("method", &SolveResult::method)
      .def_readonly("iters", &SolveResult::iters)
      .def_readonly("residual", &SolveResult::residual)
      .def_property_readonly(
          "x", [](SolveResult &s) -> VecC & { return s.x; },
          py::return_value_policy::reference_internal);

  m.def("solve_linear", &solve_linear, "Solve linear system", py::arg("A"),
        py::arg("b"), py::arg("options") = SolveOptions());

  py::class_<SParams2>(m, "SParams2")
      .def(py::init<>())
      .def_readwrite("s11", &SParams2::s11)
      .def_readwrite("s21", &SParams2::s21)
      .def_readwrite("s12", &SParams2::s12)
      .def_readwrite("s22", &SParams2::s22);

  m.def("write_touchstone", &write_touchstone, "Write Touchstone file",
        py::arg("path"), py::arg("freq"), py::arg("data"));
}

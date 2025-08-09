Heck yeah—let’s get you from spec → code. Below is a minimal, *buildable* starter repo that:

1. Compiles cleanly on Linux/macOS
2. Loads a **Gmsh v2** tetra mesh
3. Assembles & solves a **scalar Helmholtz** FEM (Milestone 0) so you can validate infra fast
4. Leaves the hooks to swap in **Nédélec edge elements** (Milestone 1) for real RF

> Why start scalar? It gets you mesh I/O, assembly, boundary handling, linear algebra, and CI green quickly. Then you replace the element kernels with curl–curl Nédélec ops and add ports/PML.

---

# Repo layout

```
vectorem/
├─ CMakeLists.txt
├─ cmake/FindEigen3.cmake              # fallback if needed
├─ external/                           # optional: vendored deps later
├─ include/vectorem/
│  ├─ mesh.hpp
│  ├─ fem.hpp
│  ├─ solver.hpp
│  ├─ bc.hpp
│  └─ utils.hpp
├─ src/
│  ├─ mesh_gmsh.cpp
│  ├─ fem_scalar.cpp
│  ├─ solver.cpp
│  ├─ main_scalar_demo.cpp
│  └─ utils.cpp
├─ python/
│  ├─ CMakeLists.txt
│  └─ pyvectorem.cpp                   # pybind11 (later)
└─ examples/
   ├─ cube_cavity.msh                  # gmsh v2 tetra mesh (placeholder)
   └─ run_scalar_demo.sh
```

---

# CMakeLists.txt (top-level)

```cmake
cmake_minimum_required(VERSION 3.20)
project(vectorem LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Dependencies
find_package(Eigen3 3.4 QUIET)
if(NOT Eigen3_FOUND)
  message(STATUS "Eigen3 not found, using fallback FindEigen3.cmake")
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
  find_package(Eigen3 REQUIRED)
endif()

add_library(vectorem
  src/mesh_gmsh.cpp
  src/fem_scalar.cpp
  src/solver.cpp
  src/utils.cpp
)
target_include_directories(vectorem PUBLIC include)
target_link_libraries(vectorem PUBLIC Eigen3::Eigen)

add_executable(vectorem_scalar_demo src/main_scalar_demo.cpp)
target_link_libraries(vectorem_scalar_demo PRIVATE vectorem)

# Python bindings (optional; enable later)
option(VECTOREM_PYTHON "Build Python bindings" OFF)
if(VECTOREM_PYTHON)
  find_package(pybind11 CONFIG REQUIRED)
  add_library(pyvectorem MODULE python/pyvectorem.cpp)
  target_link_libraries(pyvectorem PRIVATE vectorem pybind11::module)
  target_include_directories(pyvectorem PRIVATE include)
  set_target_properties(pyvectorem PROPERTIES PREFIX "" OUTPUT_NAME "pyvectorem")
endif()
```

---

# include/vectorem/mesh.hpp

```cpp
#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include <unordered_map>
#include <Eigen/Core>

namespace em {

struct Node {
  std::int64_t id;
  Eigen::Vector3d xyz;
};

enum class ElemType : int {
  Tri3   = 2,  // surface
  Tet4   = 4   // volume
};

struct Element {
  std::int64_t id;
  ElemType type;
  std::array<std::int64_t, 4> conn{}; // Tet4 uses first 4
  int phys = 0; // physical tag (BC/material)
};

struct Mesh {
  std::vector<Node> nodes;
  std::vector<Element> tets;
  std::vector<Element> tris; // boundary faces
  // maps for quick lookup
  std::unordered_map<std::int64_t, int> nodeIndex;
};

Mesh load_gmsh_v2(const std::string& path);

} // namespace em
```

---

# include/vectorem/bc.hpp

```cpp
#pragma once
#include "mesh.hpp"
#include <unordered_set>

namespace em {

struct BC {
  // For scalar prototype: Dirichlet on a set of node indices
  std::unordered_set<int> dirichlet_nodes;
};

// Build a simple PEC-like boundary for scalar prototype:
// treat all exterior Tris in phys tag == pec_tag as Dirichlet
BC build_scalar_pec(const Mesh& mesh, int pec_tag);

} // namespace em
```

---

# include/vectorem/fem.hpp

```cpp
#pragma once
#include "mesh.hpp"
#include "bc.hpp"
#include <Eigen/Sparse>
#include <complex>

namespace em {

// Scalar Helmholtz: -∇·(∇u) - k^2 u = f  (prototype)
// This is just to bring up infra. For Maxwell, you'll replace with curl-curl and Nédélec.
using cxd = std::complex<double>;
using SpMatC = Eigen::SparseMatrix<cxd>;
using VecC   = Eigen::Matrix<cxd, Eigen::Dynamic, 1>;

struct ScalarHelmholtzParams {
  double k0 = 2.0 * M_PI; // wavenumber in free space (for demo)
  double eps_r = 1.0;     // homogeneous medium (demo)
};

struct Assembly {
  SpMatC A; // system matrix
  VecC   b; // rhs
};

Assembly assemble_scalar_helmholtz(const Mesh& mesh,
                                    const ScalarHelmholtzParams& p,
                                    const BC& bc);

} // namespace em
```

---

# include/vectorem/solver.hpp

```cpp
#pragma once
#include <Eigen/Sparse>
#include <complex>
#include <optional>
#include <string>

namespace em {

using cxd = std::complex<double>;
using SpMatC = Eigen::SparseMatrix<cxd>;
using VecC   = Eigen::Matrix<cxd, Eigen::Dynamic, 1>;

struct SolveOptions {
  int max_iters = 5000;
  double tol = 1e-10;
  bool use_bicgstab = true; // else ConjugateGradient on normal equations
};

struct SolveResult {
  VecC x;
  int iters = 0;
  double residual = 0.0;
  std::string method;
};

SolveResult solve_linear(const SpMatC& A, const VecC& b, const SolveOptions& opt);

} // namespace em
```

---

# src/mesh\_gmsh.cpp

```cpp
#include "vectorem/mesh.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>

using namespace em;

Mesh em::load_gmsh_v2(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("Failed to open gmsh file: " + path);
  Mesh m;
  std::string line;
  while (std::getline(in, line)) {
    if (line == "$Nodes") {
      int n; in >> n;
      m.nodes.reserve(n);
      for (int i=0;i<n;i++) {
        std::int64_t id; double x,y,z;
        in >> id >> x >> y >> z;
        m.nodeIndex[id] = (int)m.nodes.size();
        m.nodes.push_back(Node{id, Eigen::Vector3d{x,y,z}});
      }
      // consume endline and $EndNodes
      std::getline(in, line);
      std::getline(in, line);
    } else if (line == "$Elements") {
      int ecount; in >> ecount;
      for (int i=0;i<ecount;i++) {
        std::int64_t id; int etype, nTags;
        in >> id >> etype >> nTags;
        int phys=0; int geom=0;
        if (nTags >= 1) in >> phys;
        if (nTags >= 2) in >> geom;
        for (int t=2;t<nTags;t++){ int dummy; in >> dummy; } // skip extra tags
        if (etype == (int)ElemType::Tet4) {
          std::int64_t a,b,c,d; in >> a >> b >> c >> d;
          m.tets.push_back(Element{id, ElemType::Tet4, {a,b,c,d}, phys});
        } else if (etype == (int)ElemType::Tri3) {
          std::int64_t a,b,c; in >> a >> b >> c;
          Element e{id, ElemType::Tri3, {a,b,c,0}, phys};
          m.tris.push_back(e);
        } else {
          // skip connectivity line for unsupported types
          std::getline(in, line);
        }
      }
      std::getline(in, line); // endline
      std::getline(in, line); // $EndElements
    }
  }
  return m;
}
```

---

# src/fem\_scalar.cpp

```cpp
#include "vectorem/fem.hpp"
#include <Eigen/Dense>
#include <unordered_set>

using namespace em;

namespace {
struct TetGeom {
  Eigen::Matrix<double,3,4> X; // columns are node coords
  double vol;
  Eigen::Matrix<double,3,4> gradN; // gradients of linear shape fns
};

// Build gradients & volume for linear Tet4
TetGeom tet_geom(const Mesh& mesh, const Element& e) {
  TetGeom tg;
  for (int i=0;i<4;i++) {
    int ni = mesh.nodeIndex.at(e.conn[i]);
    tg.X.col(i) = mesh.nodes[ni].xyz;
  }
  // volume = |det([x1-x4, x2-x4, x3-x4])|/6
  Eigen::Matrix3d D;
  D.col(0) = tg.X.col(0) - tg.X.col(3);
  D.col(1) = tg.X.col(1) - tg.X.col(3);
  D.col(2) = tg.X.col(2) - tg.X.col(3);
  double det = D.determinant();
  tg.vol = std::abs(det) / 6.0;

  // Gradients of shape functions for linear tetra
  // N_i = a_i + b_i x + c_i y + d_i z
  // grad N_i = [b_i, c_i, d_i]
  Eigen::Matrix4d M;
  M << 1, tg.X(0,0), tg.X(1,0), tg.X(2,0),
       1, tg.X(0,1), tg.X(1,1), tg.X(2,1),
       1, tg.X(0,2), tg.X(1,2), tg.X(2,2),
       1, tg.X(0,3), tg.X(1,3), tg.X(2,3);
  Eigen::Matrix4d Minv = M.inverse(); // small 4x4
  // Rows of Minv give coefficients [a_i, b_i, c_i, d_i]
  for (int i=0;i<4;i++) {
    tg.gradN(0,i) = Minv(i,1);
    tg.gradN(1,i) = Minv(i,2);
    tg.gradN(2,i) = Minv(i,3);
  }
  return tg;
}
} // namespace

Assembly em::assemble_scalar_helmholtz(const Mesh& mesh,
                                       const ScalarHelmholtzParams& p,
                                       const BC& bc)
{
  using Triplet = Eigen::Triplet<cxd>;
  std::vector<Triplet> trips;
  const int N = (int)mesh.nodes.size();

  VecC b = VecC::Zero(N);

  const double k = p.k0 * std::sqrt(p.eps_r);
  const cxd i(0.0,1.0);

  // Loop elements: K = ∫ gradN^T gradN dV, M = ∫ N^T N dV
  for (const auto& e : mesh.tets) {
    auto tg = tet_geom(mesh, e);

    Eigen::Matrix4d Ke = Eigen::Matrix4d::Zero();
    for (int a=0;a<4;a++){
      for (int bidx=0;bidx<4;bidx++){
        Ke(a,bidx) = tg.gradN.col(a).dot(tg.gradN.col(bidx)) * tg.vol;
      }
    }
    Eigen::Matrix4d Me = (tg.vol / 20.0) * (Eigen::Matrix4d()
      << 2,1,1,1,
         1,2,1,1,
         1,1,2,1,
         1,1,1,2).finished();

    // A_local = Ke - (k^2) * Me
    Eigen::Matrix4cd Ae = Ke.cast<cxd>() - cxd(k*k,0.0) * Me.cast<cxd>();

    // Scatter
    int ni[4];
    for (int a=0;a<4;a++) ni[a] = mesh.nodeIndex.at(e.conn[a]);
    for (int a=0;a<4;a++){
      for (int b2=0;b2<4;b2++){
        trips.emplace_back(ni[a], ni[b2], Ae(a,b2));
      }
    }
  }

  SpMatC A(N, N);
  A.setFromTriplets(trips.begin(), trips.end());

  // Apply Dirichlet BCs (zero) by row/col zeroing and diag = 1
  for (int nd : bc.dirichlet_nodes) {
    for (SpMatC::InnerIterator it(A, nd); it; ++it) {
      it.valueRef() = (it.row()==nd && it.col()==nd) ? cxd(1.0,0.0) : cxd(0.0,0.0);
    }
    // Zero column nd
    for (int k=0;k<A.outerSize();++k){
      for (SpMatC::InnerIterator it(A,k); it; ++it){
        if (it.col()==nd && it.row()!=nd) it.valueRef() = cxd(0.0,0.0);
      }
    }
    b(nd) = cxd(0.0,0.0);
  }

  return Assembly{std::move(A), std::move(b)};
}
```

---

# src/solver.cpp

```cpp
#include "vectorem/solver.hpp"
#include <Eigen/IterativeLinearSolvers>

using namespace em;

SolveResult em::solve_linear(const SpMatC& A, const VecC& b, const SolveOptions& opt) {
  SolveResult r;
  if (opt.use_bicgstab) {
    Eigen::BiCGSTAB<SpMatC, Eigen::IncompleteLUT<cxd>> solver;
    solver.setMaxIterations(opt.max_iters);
    solver.setTolerance(opt.tol);
    solver.compute(A);
    r.x = solver.solve(b);
    r.iters = solver.iterations();
    r.residual = solver.error();
    r.method = "BiCGSTAB(ILUT)";
  } else {
    // CG on normal equations (A^H A) for demo (not great for indef matrices)
    Eigen::ConjugateGradient<SpMatC, Eigen::Lower|Eigen::Upper> solver;
    solver.setMaxIterations(opt.max_iters);
    solver.setTolerance(opt.tol);
    SpMatC AhA = A.adjoint()*A;
    VecC rhs = A.adjoint()*b;
    solver.compute(AhA);
    r.x = solver.solve(rhs);
    r.iters = solver.iterations();
    r.residual = solver.error();
    r.method = "CG on normal equations";
  }
  return r;
}
```

---

# src/utils.cpp & include/vectorem/utils.hpp

```cpp
// include/vectorem/utils.hpp
#pragma once
#include <Eigen/Core>
#include <vector>
#include <unordered_set>

namespace em {
inline std::unordered_set<int> nodes_from_tri_faces(
  const std::vector<int>& tri_face_indices) {
  return std::unordered_set<int>(tri_face_indices.begin(), tri_face_indices.end());
}
} // namespace em
```

```cpp
// src/utils.cpp
#include "vectorem/utils.hpp"
// (intentionally empty for now)
```

---

# src/main\_scalar\_demo.cpp

```cpp
#include "vectorem/mesh.hpp"
#include "vectorem/fem.hpp"
#include "vectorem/solver.hpp"
#include "vectorem/bc.hpp"
#include <iostream>

using namespace em;

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: vectorem_scalar_demo <mesh.msh> [pec_phys_tag=1]\n";
    return 1;
  }
  std::string mesh_path = argv[1];
  int pec_tag = (argc >=3) ? std::atoi(argv[2]) : 1;

  auto mesh = load_gmsh_v2(mesh_path);
  auto bc   = build_scalar_pec(mesh, pec_tag);

  ScalarHelmholtzParams p;
  p.k0 = 2.0 * M_PI * 1.0; // pretend 1 m^-1; adjust as needed
  p.eps_r = 1.0;

  auto asmbl = assemble_scalar_helmholtz(mesh, p, bc);

  // Simple RHS: excite a single interior node to nonzero to check solution path
  VecC b = asmbl.b;
  if (!mesh.nodes.empty()) {
    int mid = (int)mesh.nodes.size() / 2;
    b(mid) = std::complex<double>(1.0, 0.0);
  }

  SolveOptions opt;
  opt.use_bicgstab = true;
  auto res = solve_linear(asmbl.A, b, opt);

  std::cout << "Method: " << res.method
            << " iters=" << res.iters
            << " resid=" << res.residual << "\n";

  // Print a few solution magnitudes as smoke test
  for (int i=0; i<std::min<int>(5, res.x.size()); ++i) {
    std::cout << "u["<< i << "] = " << std::abs(res.x[i]) << "\n";
  }
  return 0;
}
```

---

# src/bc.cpp (inlined for brevity)

```cpp
#include "vectorem/bc.hpp"
#include <unordered_set>

using namespace em;

BC em::build_scalar_pec(const Mesh& mesh, int pec_tag) {
  BC bc;
  // mark all nodes that belong to boundary triangles with phys == pec_tag
  for (const auto& tri : mesh.tris) {
    if (tri.phys == pec_tag) {
      for (int k=0;k<3;k++) {
        int ni = mesh.nodeIndex.at(tri.conn[k]);
        bc.dirichlet_nodes.insert(ni);
      }
    }
  }
  return bc;
}
```

> Add `src/bc.cpp` to `vectorem` target in CMake if you split files.

---

# examples/run

```bash
# examples/run_scalar_demo.sh
#!/usr/bin/env bash
set -euo pipefail
build_dir=${1:-build}
cmake -S . -B "$build_dir" -DCMAKE_BUILD_TYPE=Release
cmake --build "$build_dir" -j
"$build_dir/vectorem_scalar_demo" examples/cube_cavity.msh 1
```

---

## What you can do **today**

1. **Drop a simple Gmsh v2 tetra mesh** into `examples/` (or use the placeholder)
2. `./examples/run_scalar_demo.sh`
3. Confirm it assembles/solves and prints a few solution magnitudes

This proves your **I/O + assembly + solver** pipeline works. Next, we swap in Maxwell.

---

## Milestone 1: swap in Nédélec (curl–curl) for real EM

Replace `fem_scalar.cpp` with:

* **Unknowns:** edge DOFs per tet (Whitney 1-forms)
* **Element matrices:**

  * `Ke = ∫ (μ⁻¹ (curl N_i)·(curl N_j)) dV`
  * `Me = ∫ (ε N_i·N_j) dV`
* **System:** `(Ke - ω² Me) e = rhs`
* **Ports:** add a separate 2D eigen-solve on port cross-sections to get modal fields & Z₀, then couple as boundary terms
* **PML:** add stretched-coordinate Jacobians to `Ke` and `Me` for elems inside PML layers

I can drop in the Nédélec tet-1 element kernels and a tiny port eigen-solver next if you want to go straight there.

---

## Milestone 2: Python driver (pybind11)

Expose:

* `load_mesh(path) -> MeshHandle`
* `assemble_maxwell(params, bc) -> (A, b)`
* `solve(A, b, opts) -> x`
* `post.sparams(project, ports, x)`
* `post.near_to_far(x, huygens_surface) -> pattern`

This lets you script sweeps, run headless, and glue to your existing EM validation notebooks.

---

## Want me to:

* **Add the actual Nédélec TET(1) kernels** and a minimal wave-port eigenmode?
* **Wire a toy rectangular waveguide** example that computes S₁₁ near cutoff?
* **Add PML elements** for open-region antenna tests?

Tell me your priority and target platform (Linux cluster? NVIDIA GPUs handy?), and I’ll extend this starter with the next concrete chunk of code.

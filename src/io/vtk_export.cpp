#include "edgefem/io/vtk_export.hpp"
#include "edgefem/edge_basis.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace edgefem {

namespace {

/// Compute E-field at a point within a tetrahedron from edge DOFs.
/// Uses Whitney (Nédélec) edge basis interpolation.
Eigen::Vector3cd compute_field_at_point(const Mesh &mesh, const Element &tet,
                                        const VecC &solution,
                                        const Eigen::Vector3d &point) {
  // Get node positions for this tet
  std::array<Eigen::Vector3d, 4> X;
  for (int i = 0; i < 4; ++i) {
    int idx = mesh.nodeIndex.at(tet.conn[i]);
    X[i] = mesh.nodes[idx].xyz;
  }

  // Compute barycentric coordinates
  Eigen::Vector4d lambda = compute_barycentric(X, point);

  // Sum contributions from all 6 edges
  Eigen::Vector3cd E = Eigen::Vector3cd::Zero();

  // Edge connectivity: (n0,n1) for local nodes
  constexpr int edge_nodes[6][2] = {{0, 1}, {0, 2}, {0, 3},
                                    {1, 2}, {1, 3}, {2, 3}};

  for (int e = 0; e < 6; ++e) {
    int n0 = edge_nodes[e][0];
    int n1 = edge_nodes[e][1];

    // Whitney basis: W_e = lambda_n0 * grad(lambda_n1) - lambda_n1 * grad(lambda_n0)
    Eigen::Vector3d grad_lambda0 = compute_grad_lambda(X, n0);
    Eigen::Vector3d grad_lambda1 = compute_grad_lambda(X, n1);

    Eigen::Vector3d W_e = lambda(n0) * grad_lambda1 - lambda(n1) * grad_lambda0;

    // Get solution coefficient for this edge (with orientation)
    int global_edge = tet.edges[e];
    int orient = tet.edge_orient[e];
    std::complex<double> coeff = solution(global_edge) * static_cast<double>(orient);

    E += coeff * W_e;
  }

  return E;
}

/// Compute centroid of tetrahedron
Eigen::Vector3d compute_centroid(const Mesh &mesh, const Element &tet) {
  Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
  for (int i = 0; i < 4; ++i) {
    int idx = mesh.nodeIndex.at(tet.conn[i]);
    centroid += mesh.nodes[idx].xyz;
  }
  return centroid / 4.0;
}

} // namespace

void export_fields_vtk(const Mesh &mesh, const VecC &solution,
                       const std::string &filename,
                       const VTKExportOptions &options) {
  std::ofstream f(filename);
  if (!f.is_open()) {
    throw std::runtime_error("Cannot open file for VTK export: " + filename);
  }

  const size_t num_tets = mesh.tets.size();
  const size_t num_nodes = mesh.nodes.size();

  // Compute fields at element centroids
  std::vector<Eigen::Vector3cd> E_fields(num_tets);
  std::vector<Eigen::Vector3d> centroids(num_tets);

  for (size_t i = 0; i < num_tets; ++i) {
    centroids[i] = compute_centroid(mesh, mesh.tets[i]);
    E_fields[i] = compute_field_at_point(mesh, mesh.tets[i], solution, centroids[i]);
  }

  // Write VTK XML header
  f << "<?xml version=\"1.0\"?>\n";
  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  f << "  <UnstructuredGrid>\n";
  f << "    <Piece NumberOfPoints=\"" << num_nodes << "\" NumberOfCells=\"" << num_tets << "\">\n";

  // Point data (at cell centroids would require CellData, but we use original nodes)
  f << "      <Points>\n";
  f << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const auto &node : mesh.nodes) {
    f << "          " << std::scientific << std::setprecision(8)
      << node.xyz[0] << " " << node.xyz[1] << " " << node.xyz[2] << "\n";
  }
  f << "        </DataArray>\n";
  f << "      </Points>\n";

  // Cell connectivity
  f << "      <Cells>\n";
  f << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for (const auto &tet : mesh.tets) {
    // Map from global node IDs to local indices
    f << "          ";
    for (int i = 0; i < 4; ++i) {
      auto it = mesh.nodeIndex.find(tet.conn[i]);
      if (it != mesh.nodeIndex.end()) {
        f << it->second << " ";
      }
    }
    f << "\n";
  }
  f << "        </DataArray>\n";

  f << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  f << "          ";
  for (size_t i = 1; i <= num_tets; ++i) {
    f << (i * 4) << " ";
    if (i % 20 == 0) f << "\n          ";
  }
  f << "\n        </DataArray>\n";

  f << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  f << "          ";
  for (size_t i = 0; i < num_tets; ++i) {
    f << "10 ";  // VTK_TETRA = 10
    if ((i + 1) % 30 == 0) f << "\n          ";
  }
  f << "\n        </DataArray>\n";
  f << "      </Cells>\n";

  // Cell data (fields at centroids)
  f << "      <CellData>\n";

  if (options.export_magnitude) {
    f << "        <DataArray type=\"Float64\" Name=\"E_magnitude\" format=\"ascii\">\n";
    for (size_t i = 0; i < num_tets; ++i) {
      double mag = E_fields[i].norm();
      f << "          " << std::scientific << mag << "\n";
    }
    f << "        </DataArray>\n";
  }

  if (options.export_phase) {
    f << "        <DataArray type=\"Float64\" Name=\"E_phase_rad\" format=\"ascii\">\n";
    for (size_t i = 0; i < num_tets; ++i) {
      // Use phase of dominant component
      double max_mag = 0.0;
      double phase = 0.0;
      for (int j = 0; j < 3; ++j) {
        double m = std::abs(E_fields[i](j));
        if (m > max_mag) {
          max_mag = m;
          phase = std::arg(E_fields[i](j));
        }
      }
      f << "          " << std::scientific << phase << "\n";
    }
    f << "        </DataArray>\n";
  }

  if (options.export_real) {
    f << "        <DataArray type=\"Float64\" Name=\"E_real\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (size_t i = 0; i < num_tets; ++i) {
      f << "          " << std::scientific
        << E_fields[i](0).real() << " "
        << E_fields[i](1).real() << " "
        << E_fields[i](2).real() << "\n";
    }
    f << "        </DataArray>\n";
  }

  if (options.export_imag) {
    f << "        <DataArray type=\"Float64\" Name=\"E_imag\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (size_t i = 0; i < num_tets; ++i) {
      f << "          " << std::scientific
        << E_fields[i](0).imag() << " "
        << E_fields[i](1).imag() << " "
        << E_fields[i](2).imag() << "\n";
    }
    f << "        </DataArray>\n";
  }

  if (options.export_vector) {
    // Export real part as 3D vector for glyph visualization
    f << "        <DataArray type=\"Float64\" Name=\"E_vector\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (size_t i = 0; i < num_tets; ++i) {
      f << "          " << std::scientific
        << E_fields[i](0).real() << " "
        << E_fields[i](1).real() << " "
        << E_fields[i](2).real() << "\n";
    }
    f << "        </DataArray>\n";
  }

  // Also export physical group tag for region identification
  f << "        <DataArray type=\"Int32\" Name=\"PhysicalTag\" format=\"ascii\">\n";
  f << "          ";
  for (size_t i = 0; i < num_tets; ++i) {
    f << mesh.tets[i].phys << " ";
    if ((i + 1) % 20 == 0) f << "\n          ";
  }
  f << "\n        </DataArray>\n";

  f << "      </CellData>\n";

  f << "    </Piece>\n";
  f << "  </UnstructuredGrid>\n";
  f << "</VTKFile>\n";
}

void export_fields_vtk_legacy(const Mesh &mesh, const VecC &solution,
                              const std::string &filename,
                              const VTKExportOptions &options) {
  std::ofstream f(filename);
  if (!f.is_open()) {
    throw std::runtime_error("Cannot open file for VTK export: " + filename);
  }

  const size_t num_tets = mesh.tets.size();
  const size_t num_nodes = mesh.nodes.size();

  // Compute fields at element centroids
  std::vector<Eigen::Vector3cd> E_fields(num_tets);
  for (size_t i = 0; i < num_tets; ++i) {
    Eigen::Vector3d centroid = compute_centroid(mesh, mesh.tets[i]);
    E_fields[i] = compute_field_at_point(mesh, mesh.tets[i], solution, centroid);
  }

  // Write legacy VTK header
  f << "# vtk DataFile Version 3.0\n";
  f << "EdgeFEM Field Export\n";
  f << "ASCII\n";
  f << "DATASET UNSTRUCTURED_GRID\n";

  // Points
  f << "POINTS " << num_nodes << " double\n";
  for (const auto &node : mesh.nodes) {
    f << std::scientific << std::setprecision(8)
      << node.xyz[0] << " " << node.xyz[1] << " " << node.xyz[2] << "\n";
  }

  // Cells
  f << "CELLS " << num_tets << " " << (num_tets * 5) << "\n";
  for (const auto &tet : mesh.tets) {
    f << "4";
    for (int i = 0; i < 4; ++i) {
      auto it = mesh.nodeIndex.find(tet.conn[i]);
      if (it != mesh.nodeIndex.end()) {
        f << " " << it->second;
      }
    }
    f << "\n";
  }

  f << "CELL_TYPES " << num_tets << "\n";
  for (size_t i = 0; i < num_tets; ++i) {
    f << "10\n";  // VTK_TETRA = 10
  }

  // Cell data
  f << "CELL_DATA " << num_tets << "\n";

  if (options.export_magnitude) {
    f << "SCALARS E_magnitude double 1\n";
    f << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < num_tets; ++i) {
      f << std::scientific << E_fields[i].norm() << "\n";
    }
  }

  if (options.export_vector) {
    f << "VECTORS E_vector double\n";
    for (size_t i = 0; i < num_tets; ++i) {
      f << std::scientific
        << E_fields[i](0).real() << " "
        << E_fields[i](1).real() << " "
        << E_fields[i](2).real() << "\n";
    }
  }
}

} // namespace edgefem

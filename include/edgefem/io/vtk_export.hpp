#pragma once

#include <string>

#include "edgefem/fem.hpp"
#include "edgefem/mesh.hpp"

namespace edgefem {

/// Options for VTK field export.
struct VTKExportOptions {
  /// Export E-field magnitude (default: true)
  bool export_magnitude = true;

  /// Export E-field phase in radians (default: true)
  bool export_phase = true;

  /// Export E-field real part (default: false)
  bool export_real = false;

  /// Export E-field imaginary part (default: false)
  bool export_imag = false;

  /// Export E-field as 3D vector (real parts) (default: true)
  bool export_vector = true;
};

/// Export E-field solution to VTK XML Unstructured Grid (.vtu) format.
/// The field is computed at element centroids from edge solution.
///
/// @param mesh The FEM mesh
/// @param solution The edge-element solution vector (complex)
/// @param filename Output filename (should end in .vtu)
/// @param options Export options for selecting which fields to include
void export_fields_vtk(const Mesh &mesh, const VecC &solution,
                       const std::string &filename,
                       const VTKExportOptions &options = VTKExportOptions());

/// Export E-field solution to legacy VTK format (.vtk).
/// Compatible with older visualization tools.
///
/// @param mesh The FEM mesh
/// @param solution The edge-element solution vector (complex)
/// @param filename Output filename (should end in .vtk)
/// @param options Export options
void export_fields_vtk_legacy(const Mesh &mesh, const VecC &solution,
                              const std::string &filename,
                              const VTKExportOptions &options = VTKExportOptions());

} // namespace edgefem

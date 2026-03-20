// Coaxial Transmission Line (50 Ohm)
// Air-filled: Z0 = (1/(2π)) * sqrt(μ/ε) * ln(b/a) = 50 Ω
// For Z0 = 50Ω in air: ln(b/a) = 2π * 50 / 377 ≈ 0.833
// => b/a ≈ 2.3
//
// Design: a = 1mm (inner radius), b = 2.3mm (outer radius), L = 30mm

lc_inner = 0.0005;  // Finer mesh near inner conductor
lc_outer = 0.001;   // Coarser mesh elsewhere

// Dimensions in meters
a = 0.001;      // Inner conductor radius (1 mm)
b = 0.0023;     // Outer conductor radius (2.3 mm)
L = 0.030;      // Length (30 mm)

// Create inner conductor boundary circle at z=0
Point(1) = {0, 0, 0, lc_inner};       // Center
Point(2) = {a, 0, 0, lc_inner};       // +x
Point(3) = {0, a, 0, lc_inner};       // +y
Point(4) = {-a, 0, 0, lc_inner};      // -x
Point(5) = {0, -a, 0, lc_inner};      // -y

// Create outer conductor boundary circle at z=0
Point(6) = {b, 0, 0, lc_outer};       // +x
Point(7) = {0, b, 0, lc_outer};       // +y
Point(8) = {-b, 0, 0, lc_outer};      // -x
Point(9) = {0, -b, 0, lc_outer};      // -y

// Inner circle arcs
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Outer circle arcs
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

// Create curve loops
Curve Loop(1) = {1, 2, 3, 4};   // Inner circle
Curve Loop(2) = {5, 6, 7, 8};   // Outer circle

// Create annular cross-section (dielectric region)
Plane Surface(1) = {2, 1};      // Outer minus inner = annulus

// Extrude to create coaxial line
Extrude {0, 0, L} { Surface{1}; }

// Physical groups
// Note: After extrusion, surface IDs change. These need to be verified with Gmsh GUI
// Typical assignment:
// - Original surface (1) becomes Port 1
// - Extruded end surface becomes Port 2
// - Lateral surfaces are walls

// Inner conductor PEC (lateral surface of inner cylinder)
Physical Surface("InnerPEC", 1) = {25, 29, 33, 37};

// Outer conductor PEC (lateral surface of outer cylinder)
Physical Surface("OuterPEC", 2) = {26, 30, 34, 38};

// Port 1 (z=0, annular region)
Physical Surface("Port1", 3) = {1};

// Port 2 (z=L, annular region)
Physical Surface("Port2", 4) = {50};

// Dielectric volume (air)
Physical Volume("Air", 100) = {1};

// Mesh settings
Mesh.Algorithm = 6;    // Frontal-Delaunay for 2D
Mesh.Algorithm3D = 4;  // Frontal Delaunay for 3D
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

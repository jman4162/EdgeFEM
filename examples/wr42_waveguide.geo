// Rectangular Waveguide (WR-42) - Ka-band
// a = 10.668 mm, b = 4.318 mm, L = 30 mm
// TE10 cutoff: fc = c/(2a) ≈ 14.05 GHz
// Operating band: 18-26.5 GHz

lc = 0.0015;  // Mesh size ~1.5mm (lambda/10 at 20 GHz)

// Waveguide dimensions in meters
a = 0.010668;  // 10.668 mm
b = 0.004318;  // 4.318 mm
L = 0.030;     // 30 mm

// Cross-section at z=0
Point(1) = {0, 0, 0, lc};
Point(2) = {a, 0, 0, lc};
Point(3) = {a, b, 0, lc};
Point(4) = {0, b, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Extrude to create waveguide
Extrude {0, 0, L} { Surface{1}; }

// Physical groups
Physical Surface("PEC", 1) = {13, 17, 21, 25};  // Walls
Physical Surface("Port1", 2) = {1};              // Input (z=0)
Physical Surface("Port2", 3) = {26};             // Output (z=L)
Physical Volume("Air", 100) = {1};

// Mesh settings
Mesh.Algorithm3D = 4;  // Frontal Delaunay
Mesh.Optimize = 1;

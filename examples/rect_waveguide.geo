// Rectangular Waveguide (WR-90)
// a = 22.86 mm, b = 10.16 mm, L = 50 mm
// TE10 cutoff: fc = c/(2a) ≈ 6.56 GHz

lc = 0.003;  // Mesh size ~3mm

// Waveguide dimensions in meters
a = 0.02286;
b = 0.01016;
L = 0.05;

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

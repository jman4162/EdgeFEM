// Gmsh geometry for a WR-90 waveguide cross-section
// Dimensions are in meters.

a = 0.02286; // Width
b = 0.01016; // Height
lc = 0.002;  // Mesh characteristic length

Point(1) = {0, 0, 0, lc};
Point(2) = {a, 0, 0, lc};
Point(3) = {a, b, 0, lc};
Point(4) = {0, b, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Define physical groups for boundaries and the surface
Physical Curve("pec", 1) = {1, 2, 3, 4};
Physical Surface("domain", 2) = {1};

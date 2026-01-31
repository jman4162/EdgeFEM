// unit_cell.geo
// Periodic unit cell for phased array element simulation
// Square unit cell with patch element on dielectric substrate

// Unit cell parameters
cell_size = 0.015;      // 15mm cell size (~lambda/2 at 10 GHz)
substrate_h = 0.00157;  // 1.57mm substrate (typical FR4 thickness)
patch_w = 0.0085;       // 8.5mm patch width
patch_l = 0.0095;       // 9.5mm patch length
air_h = 0.010;          // 10mm air region above patch

// Mesh size
lc = cell_size / 8;     // ~8 elements per wavelength

// Ground plane (z=0)
Point(1) = {0, 0, 0, lc};
Point(2) = {cell_size, 0, 0, lc};
Point(3) = {cell_size, cell_size, 0, lc};
Point(4) = {0, cell_size, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Substrate top (z = substrate_h)
Point(5) = {0, 0, substrate_h, lc};
Point(6) = {cell_size, 0, substrate_h, lc};
Point(7) = {cell_size, cell_size, substrate_h, lc};
Point(8) = {0, cell_size, substrate_h, lc};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// Patch corners
offset_x = (cell_size - patch_w) / 2;
offset_y = (cell_size - patch_l) / 2;

Point(9) = {offset_x, offset_y, substrate_h, lc/2};
Point(10) = {offset_x + patch_w, offset_y, substrate_h, lc/2};
Point(11) = {offset_x + patch_w, offset_y + patch_l, substrate_h, lc/2};
Point(12) = {offset_x, offset_y + patch_l, substrate_h, lc/2};

Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 9};

// Air top (z = substrate_h + air_h)
Point(13) = {0, 0, substrate_h + air_h, lc};
Point(14) = {cell_size, 0, substrate_h + air_h, lc};
Point(15) = {cell_size, cell_size, substrate_h + air_h, lc};
Point(16) = {0, cell_size, substrate_h + air_h, lc};

Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 13};

Curve Loop(2) = {13, 14, 15, 16};
Plane Surface(2) = {2};

// Vertical lines for substrate
Line(17) = {1, 5};
Line(18) = {2, 6};
Line(19) = {3, 7};
Line(20) = {4, 8};

// Vertical lines for air
Line(21) = {5, 13};
Line(22) = {6, 14};
Line(23) = {7, 15};
Line(24) = {8, 16};

// Substrate side surfaces
Curve Loop(3) = {1, 18, -5, -17};
Plane Surface(3) = {3};
Curve Loop(4) = {2, 19, -6, -18};
Plane Surface(4) = {4};
Curve Loop(5) = {3, 20, -7, -19};
Plane Surface(5) = {5};
Curve Loop(6) = {4, 17, -8, -20};
Plane Surface(6) = {6};

// Substrate top surface (with patch hole)
Curve Loop(7) = {5, 6, 7, 8};
Curve Loop(8) = {9, 10, 11, 12};
Plane Surface(7) = {7, 8};

// Patch surface
Curve Loop(9) = {9, 10, 11, 12};
Plane Surface(8) = {9};

// Air side surfaces
Curve Loop(10) = {5, 22, -13, -21};
Plane Surface(9) = {10};
Curve Loop(11) = {6, 23, -14, -22};
Plane Surface(10) = {11};
Curve Loop(12) = {7, 24, -15, -23};
Plane Surface(11) = {12};
Curve Loop(13) = {8, 21, -16, -24};
Plane Surface(12) = {13};

// Substrate volume
Surface Loop(1) = {1, 3, 4, 5, 6, 7, 8};
Volume(1) = {1};

// Air volume
Surface Loop(2) = {7, 8, 9, 10, 11, 12, 2};
Volume(2) = {2};

// Physical groups
Physical Surface("Ground", 1) = {1};        // Ground plane (PEC)
Physical Surface("Patch", 2) = {8};         // Patch (PEC)
Physical Surface("FloquetTop", 3) = {2};    // Floquet port (top)

// Periodic pairs: x-direction
Physical Surface("MasterX", 4) = {6, 12};   // x=0 faces (substrate + air)
Physical Surface("SlaveX", 5) = {4, 10};    // x=cell_size faces

// Periodic pairs: y-direction
Physical Surface("MasterY", 6) = {3, 9};    // y=0 faces (substrate + air)
Physical Surface("SlaveY", 7) = {5, 11};    // y=cell_size faces

Physical Volume("Substrate", 100) = {1};    // Dielectric (eps_r = 4.4)
Physical Volume("Air", 101) = {2};          // Air (eps_r = 1.0)

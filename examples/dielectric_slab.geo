// Dielectric Slab Unit Cell
// For Fresnel coefficient validation with periodic boundary conditions
//
// Geometry: Unit cell with dielectric slab at center
// - Slab thickness: d = 3mm
// - Slab permittivity: εr = 4.0 (glass-like)
// - Cell size: 10mm x 10mm x 30mm (z-propagation)
// - Slab centered at z = 15mm
//
// Fresnel coefficients at normal incidence:
//   r = (n1 - n2)/(n1 + n2) = (1 - 2)/(1 + 2) = -1/3
//   t = 2*n1/(n1 + n2) = 2/3
//
// With multiple reflections (Fabry-Pérot), the total transmission depends on
// the slab electrical thickness.

lc = 0.002;  // 2mm mesh size (~lambda/15 at 10 GHz in air)

// Cell dimensions
Lx = 0.010;   // 10mm
Ly = 0.010;   // 10mm
Lz = 0.030;   // 30mm

// Slab parameters
d_slab = 0.003;      // 3mm thickness
z_slab_start = 0.0135;  // Slab starts at 13.5mm
z_slab_end = 0.0165;    // Slab ends at 16.5mm

// Create unit cell geometry
// Bottom face (z = 0)
Point(1) = {0, 0, 0, lc};
Point(2) = {Lx, 0, 0, lc};
Point(3) = {Lx, Ly, 0, lc};
Point(4) = {0, Ly, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Extrude to slab bottom (air region 1)
out1[] = Extrude {0, 0, z_slab_start} { Surface{1}; };

// Extrude slab region
out2[] = Extrude {0, 0, d_slab} { Surface{out1[0]}; };

// Extrude to top (air region 2)
out3[] = Extrude {0, 0, Lz - z_slab_end} { Surface{out2[0]}; };

// Physical groups
// Port 1 (bottom, z=0)
Physical Surface("Port1", 1) = {1};

// Port 2 (top, z=Lz)
Physical Surface("Port2", 2) = {out3[0]};

// Periodic faces (x-direction)
Physical Surface("PeriodicX_minus", 3) = {out1[2], out2[2], out3[2]};
Physical Surface("PeriodicX_plus", 4) = {out1[4], out2[4], out3[4]};

// Periodic faces (y-direction)
Physical Surface("PeriodicY_minus", 5) = {out1[5], out2[5], out3[5]};
Physical Surface("PeriodicY_plus", 6) = {out1[3], out2[3], out3[3]};

// Volumes
Physical Volume("Air1", 100) = {out1[1]};       // Air below slab
Physical Volume("Dielectric", 101) = {out2[1]}; // Dielectric slab
Physical Volume("Air2", 102) = {out3[1]};       // Air above slab

// Mesh settings
Mesh.Algorithm = 6;
Mesh.Algorithm3D = 4;
Mesh.Optimize = 1;

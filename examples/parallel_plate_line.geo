// Air-filled parallel-plate line section for the lumped-port tests
// (tests/test_lumped_port_physics.cpp).
//
// Box 0 <= x <= L, 0 <= y <= w, 0 <= z <= h. The plates are z = 0 and
// z = h. With PEC plates and PMC (natural boundary) side walls y = 0 and
// y = w, the line carries an exact TEM mode, E = z_hat E0 exp(-j k0 x),
// with characteristic impedance Zc = eta0 * h / w. Choosing
// w = eta0 * h / 50 gives Zc = 50 ohm.
//
// Physical tags:
//   1  plates (z = 0 and z = h)
//   2  end face x = 0 (port 1)
//   3  end face x = L (port 2)
//   4  side walls (y = 0 and y = w)
//   10 air volume
//
// Regenerate: gmsh -3 parallel_plate_line.geo -format msh2 -o parallel_plate_line.msh

SetFactory("OpenCASCADE");

h = 1.0e-3;
w = 376.730313668 * h / 50.0;  // 7.5346 mm, Zc = 50 ohm
L = 20.0e-3;
lc = 0.75e-3;
e = 1.0e-7;

Box(1) = {0, 0, 0, L, w, h};

plates() = Surface In BoundingBox{-e, -e, -e, L + e, w + e, e};
plates() += Surface In BoundingBox{-e, -e, h - e, L + e, w + e, h + e};
port1() = Surface In BoundingBox{-e, -e, -e, e, w + e, h + e};
port2() = Surface In BoundingBox{L - e, -e, -e, L + e, w + e, h + e};
sides() = Surface In BoundingBox{-e, -e, -e, L + e, e, h + e};
sides() += Surface In BoundingBox{-e, w - e, -e, L + e, w + e, h + e};

Physical Surface("plates", 1) = {plates()};
Physical Surface("port1", 2) = {port1()};
Physical Surface("port2", 3) = {port2()};
Physical Surface("sides", 4) = {sides()};
Physical Volume("air", 10) = {1};

Mesh.CharacteristicLengthMax = lc;
Mesh.CharacteristicLengthMin = lc / 2;

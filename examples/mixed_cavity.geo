// Gmsh geometry for a cube cavity with mixed PEC/PMC walls.
SetFactory("OpenCASCADE");
lc = 0.5;
Box(1) = {0, 0, 0, 1, 1, 1};

// Wall at x=0 is PEC (tag 1)
// Wall at x=1 is PMC (tag 2)
// Other walls are also PEC (tag 1)
Physical Surface("pec", 1) = {2, 3, 4, 5, 6};
Physical Surface("pmc", 2) = {1};
Physical Volume("vacuum", 1) = {1};

Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

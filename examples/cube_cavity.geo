SetFactory("OpenCASCADE");
lc = 0.5;
Box(1) = {0, 0, 0, 1, 1, 1};
Physical Surface(1) = {1,2,3,4,5,6};
Physical Volume(1) = {1};
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;


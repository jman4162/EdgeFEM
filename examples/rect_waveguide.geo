SetFactory("OpenCASCADE");
// WR-90 rectangular waveguide dimensions
a = 0.02286; // width (m)
b = 0.01016; // height (m)
L = 0.05;    // length (m)

Rectangle(1) = {0, 0, 0, a, b};
Extrude {0, 0, L} {
  Surface{1}; Layers{1}; Recombine;
}

Physical Surface("walls") = {1,3,4,5};
Physical Surface("ports") = {2,6};
Physical Volume("wg") = {1};

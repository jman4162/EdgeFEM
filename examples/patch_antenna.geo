SetFactory("OpenCASCADE");
// Units: meters
h=0.0016; W=0.03723426118288438; L=0.028809290261854397; inset=0.009999999999999998; feedW=0.00308;
board_x=0.08; board_y=0.07; air=0.0306;
Rectangle(1) = {-board_x/2, -board_y/2, 0, board_x, board_y};
Rectangle(2) = {-W/2, -L/2, h, W, L};
Rectangle(3) = {-feedW/2, -board_y/2, h, feedW, board_y/2 - L/2 + inset};
Rectangle(4) = {-feedW/2, -L/2, h, feedW, inset};
tmp1[] = BooleanDifference{ Surface{2}; Delete; }{ Surface{4}; Delete; };
tmp2[] = BooleanUnion{ Surface{tmp1[0]}; Delete; }{ Surface{3}; Delete; };
Rectangle(5) = {-(board_x/2+air), -(board_y/2+air), -air, board_x+2*air, board_y+2*air};
Extrude {0,0,h} { Surface{1}; }
Extrude {0,0,h+air} { Surface{5}; }

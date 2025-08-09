SetFactory("OpenCASCADE");
// Units: meters
h=0.0016; W=0.03723426118288438; L=0.028809290261854397; inset=0.004; feedW=0.00308;
board_x=0.08; board_y=0.07; air=0.0306;
board = Rectangle(1, -board_x/2, -board_y/2, 0, board_x, board_y);
patch = Rectangle(2, -W/2, -L/2, h, W, L);
feed = Rectangle(3, -feedW/2, -board_y/2, h, feedW, board_y/2 - L/2 + inset);
slot = Rectangle(4, -feedW/2, -L/2, h, feedW, inset);
BooleanDifference{ Surface{patch}; Delete; }{ Surface{slot}; Delete; }
BooleanUnion{ Surface{patch}; Delete; }{ Surface{feed}; Delete; }
airbox = Rectangle(5, -(board_x/2+air), -(board_y/2+air), -air, board_x+2*air, board_y+2*air);
Extrude {0,0,h} { Surface{board}; }
Extrude {0,0,0} { Surface{patch}; }
Extrude {0,0,0} { Surface{feed}; }
Extrude {0,0,h+air} { Surface{airbox}; }

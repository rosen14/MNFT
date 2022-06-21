N = 20 + 1;
L = 1;
dx = L/N;

Point(1) = {0, 0, 0};
Point(2) = {L, 0, 0};
Point(3) = {L, L, 0};
Point(4) = {0, L, 0};
Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};

Line Loop(1) = {4, 3, 2, 1};

Plane Surface(1) = {1};

Transfinite Line{4,3,2,1} = N Using Progression 1.0;

Transfinite Surface "*";
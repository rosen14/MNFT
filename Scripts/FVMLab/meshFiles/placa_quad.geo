dx = 0.05;
L = 1;
N = L/dx;

Point(1) = {-L, -L, 0};
Point(2) = {L, -L, 0};
Point(3) = {L, L, 0};
Point(4) = {-L, L, 0};
Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};

Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Transfinite Line{1,2,3,4} = N Using Progression 1.0;

Transfinite Surface "*";
Recombine Surface "*";

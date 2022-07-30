N = 3 + 1;
L = 1;
dx = L/N;

Point(1) = {0, 0, 0};
Point(2) = {L, 0, 0};
Point(3) = {L, 0.1*dx, 0};
Point(4) = {0, 0.1*dx, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Transfinite Line{1,3} = N;

Transfinite Surface "*";
Recombine Surface "*";
// Target element size of outer boundary
lc2 = 1.25;

// Bottom left coords
xL = 0;
yL = -5;

// Top right coords
xR = 10;
yR = 5;

// Outline points
p1 = newp; Point(p1) = {xL, yL, 0, lc2};
p2 = newp; Point(p2) = {xL, yR, 0, lc2};
p3 = newp; Point(p3) = {xR, yR, 0, lc2};
p4 = newp; Point(p4) = {xR, yL, 0, lc2};

// Outline lines
l1 = newl; Line(l1) = {p1, p2};
l2 = newl; Line(l2) = {p2, p3};
l3 = newl; Line(l3) = {p3, p4};
l4 = newl; Line(l4) = {p4, p1};

Curve Loop(1) = {l1, l2, l3, l4};

// Surface of mesh
Plane Surface(1) = {1};

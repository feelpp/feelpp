lr = 0.1;
Point(1) = {0, 0, 0, lr};
Point(2) = {1, 0, 0, lr};
Line(1) = {1, 2};

Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
Line{1};Layers{16};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
Line{2};Layers{16};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
Line{5};Layers{16};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
Line{8};Layers{16};
}

Extrude {0, 0, 5} {
Surface{4};Layers{20};Recombine;
}
Extrude {0, 0, 5} {
Surface{7};Layers{20};Recombine;
}
Extrude {0, 0, 5} {
Surface{10};Layers{20};Recombine;
}
Extrude {0, 0, 5} {
Surface{13};Layers{20};Recombine;
}
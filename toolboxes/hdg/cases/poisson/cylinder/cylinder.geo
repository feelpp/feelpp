h = 0.1;
r = 0.5;   // radius
L = 5;     // lenght
nbCut = 5; // number of part

// radius
Point(1) = {0, 0, 0, h};
line[] = Extrude{r,0,0} { Point{1};};

// for each part of the cylinder
For k In {1:nbCut}
    // if first part, extrude radius for a surface
    If( k == 1 )
        surf0[] = Extrude{0,0,L/nbCut} { Line{line[1]};};
        surf[0] = surf0[1];
    EndIf
    // else translate first surface
    If( k != 1 )
        surf[] = Translate {0, 0, (k-1)*L/nbCut} { Duplicata{ Surface{surf0[1]}; } } ;
    EndIf

    // rotate the surface by Pi/2
    out1[] = Extrude {{0,0,1},{0,0,0},Pi/2} {
                    Surface{surf[0]};
                    };
    out2[] = Extrude {{0,0,1},{0,0,0},Pi/2} {
                    Surface{out1[0]};
                    };
    out3[] = Extrude {{0,0,1},{0,0,0},Pi/2} {
                    Surface{out2[0]};
                    };
    out4[] = Extrude {{0,0,1},{0,0,0},Pi/2} {
                    Surface{out3[0]};
                    };

    // keep the indices of surfaces
    If( k == 1 )
        inflow[] = {out1[2], out2[2], out3[2], out4[2]};
    EndIf
    If( k == nbCut )
        outflow[] = {out1[4], out2[4], out3[4], out4[4]};
    EndIf
    If( k < nbCut )
        cut[] += {out1[4], out2[4], out3[4], out4[4]};
    EndIf
    wall[] += {out1[3],out2[3],out3[3],out4[3]};
    cutH1[] += {out1[0],out3[0]};
    cutH2[] += {out2[0],out4[0]};
    vol[] += {out1[1],out2[1],out3[1],out4[1]};
EndFor

Physical Surface("inlet") = inflow[];
Physical Surface("outlet") = outflow[];
Physical Surface("wall") = wall[];
// each surface between the parts
For k In {0:nbCut-2}
    Physical Surface(Sprintf("cut%g", k+1)) = {cut[k*4+0],cut[k*4+1],cut[k*4+2],cut[k*4+3]};
EndFor
// surface x=0 and y=0
Physical Surface(Sprintf("cut%g", nbCut)) = cutH1[];
Physical Surface(Sprintf("cut%g", nbCut+1)) = cutH2[];

Physical Volume("omega") = vol[];

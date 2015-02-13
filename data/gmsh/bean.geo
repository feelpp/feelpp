/*************************************************************************[ LIC ]
  This file is part of the Feel library

  Author(s): Guillaume Dolle <gdolle@unistra.fr>
       Date: 2014-08-01

  Copyright (C) 2014 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
********************************************************************************/
/*************************************************************************[ MAN ]
  This file describes a geometry for a bean shape. It is parametrized by
  two radius `r1`, `r2`, and the bean length `d`. You can choose the topological
  dimension (2D or 3D) filling the integer parameter `dim` (=2 or =3).
********************************************************************************/


////////////////////////////////////////////////////////////////////////////////
// Mesh size
h=0.5;
// Topological dimension
dim=2;
// Bean length
d=1;
// Bean left radius
r1=0.4;
// Bean right radius
r2=0.5;

// IMPORTANT NOTE: Sometimes 3D mesh generation crash for some shape!
//                 (For example d=1,r1=0.5,r2=0.6)

////////////////////////////////////////////////////////////////////////////////


//------------------------------------------------------------------------------
// Constraints.
//------------------------------------------------------------------------------

If( dim > 3 || dim < 2 )
    Error("dim should be equal to 2 or 3!");
    Abort;
EndIf

// Check GMSH version (> 2.8.4)
If( GMSH_MAJOR_VERSION < 2 )
    Error("This geometry requires GMSH version > 2.8.5");
EndIf

If( GMSH_MINOR_VERSION < 8 )
    Error("This geometry requires GMSH version > 2.8.5");
EndIf

If( GMSH_PATCH_VERSION < 5 )
    Error("This geometry requires GMSH version > 2.8.5");
EndIf

//------------------------------------------------------------------------------
// DRAW2DPOINTS creates points for the 2D plane (x,y).
//------------------------------------------------------------------------------
Function Draw2DPoints
    Printf("-- Creating 2D Points");
    p=newp; Point(p) = { 0, 0, 0, h };
    p=newp; Point(p) = { -d/2, 0, 0, h };
    p=newp; Point(p) = { d/2, 0, 0, h };

    // Left circle
    ctrl=r1/2;
    x=-d/2; y=r1;
    p=newp; pleftlist2d[]+=p; Point(p) = { x+ctrl, y, 0, h };
    p=newp; pleftlist2d[]+=p; Po2d[]+=p; Point(p) = { x, y, 0, h };
    p=newp; pleftlist2d[]+=p; Po2d[]+=p; Point(p) = { x-ctrl, y, 0, h };

    x=-d/2-r1; y=0;
    p=newp; pleftlist2d[]+=p; Point(p) = { x, y+ctrl, 0, h };
    p=newp; pleftlist2d[]+=p; Point(p) = { x, y, 0, h };
    p=newp; pleftlist2d[]+=p; Point(p) = { x, y-ctrl, 0, h};

    x=-d/2; y=-r1;
    p=newp; pleftlist2d[]+=p; Point(p) = { x-ctrl, y, 0, h };
    p=newp; pleftlist2d[]+=p; Point(p) = { x, y, 0, h };
    p=newp; pleftlist2d[]+=p; Point(p) = { x+ctrl, y, 0, h };

    // Bottom mid point
    p=newp; point2dlist[]+=p;
    Point(p) = { 0, 0.2*d*Sin(-Pi/2), 0, h };

    // Right circle
    ctrl=r2/2;
    x=d/2; y=-r2;
    p=newp; prightlist2d[]+=p; Point(p) = { x-ctrl, y, 0, h };
    p=newp; prightlist2d[]+=p; Point(p) = { x, y, 0, h };
    p=newp; prightlist2d[]+=p; Point(p) = { x+ctrl, y, 0, h };

    x=d/2+r2; y=0;
    p=newp; prightlist2d[]+=p; Point(p) = { x, y-ctrl, 0, h };
    p=newp; prightlist2d[]+=p; Point(p) = { x, y, 0, h };
    p=newp; prightlist2d[]+=p; Point(p) = { x, y+ctrl, 0, h };

    x=d/2; y=r2;
    p=newp; prightlist2d[]+=p; Point(p) = { x+ctrl, y, 0, h };
    p=newp; prightlist2d[]+=p; Point(p) = { x, y, 0, h };
    p=newp; prightlist2d[]+=p; Point(p) = { x-ctrl, y, 0, h };

    // Top mid point
    R=Sqrt((d/2)*(d/2)+((r1+r2)/2)*((r1+r2)/2));
    p=newp; point2dlist[]+=p;
    Point(p) = { 0, R*Sin(Pi/6), 0, h };
Return

//------------------------------------------------------------------------------
// DRAW2DLINES creates all lines for the 2D model.
//------------------------------------------------------------------------------
Function Draw2DLines
    Printf("-- Creating 2D Lines");
    l2d[]={};
    l=newl; l2d[]+=l; Bezier(l)={5,6,7,8};
    l=newl; l2d[]+=l; Bezier(l)={8,9,10,11};
    l=newl; l2d[]+=l; BSpline(l)={11,12,13,14,15};
    l=newl; l2d[]+=l; Bezier(l)={15,16,17,18};
    l=newl; l2d[]+=l; Bezier(l)={18,19,20,21};
    l=newl; l2d[]+=l; BSpline(l)={21,22,23,4,5};
Return

//------------------------------------------------------------------------------
// DRAW2DSURFACES creates the surface for the 2D model.
//------------------------------------------------------------------------------
Function Draw2DSurfaces
    Printf("-- Creating 2D Surfaces");
    l=newll; Line Loop(l)=l2d[];
    s=news; Plane Surface(s) = {l};
Return

//------------------------------------------------------------------------------
// ADD2DPHYSICALS gives a name for the whole domain and its boundary (2D case).
//------------------------------------------------------------------------------
Function Add2DPhysicals
    Printf("-- Add 2D Physicals");
    Physical Line( "Wall" ) = l2d[];
    Physical Surface( "Omega" ) = { s };
    Recombine Surface {s};
Return

//------------------------------------------------------------------------------
// DRAW3DPOINTS Create points for the third coordinates (z).
//------------------------------------------------------------------------------
Function Draw3DPoints
    Printf("-- Creating 3D Points");
    // CENTER
    p=newp; Point(p) = { 0, 0, -(r1+r2)/2, h };
    p=newp; Point(p) = { 0, 0, (r1+r2)/2, h };
    // LEFT BALL.
    ctrl=r1/2;
    // Middle
    p=newp; Point(p) = { -d/2, 0, r1, h };
    p=newp; Point(p) = { -d/2, ctrl, r1, h };
    p=newp; Point(p) = { -d/2, -ctrl, r1, h };
    p=newp; Point(p) = { -d/2-ctrl, 0, r1, h };
    p=newp; Point(p) = { -d/2+ctrl, 0, r1, h };

    p=newp; Point(p) = { -d/2, 0, -r1, h };
    p=newp; Point(p) = { -d/2, ctrl, -r1, h };
    p=newp; Point(p) = { -d/2, -ctrl, -r1, h };
    p=newp; Point(p) = { -d/2-ctrl, 0, -r1, h };
    p=newp; Point(p) = { -d/2+ctrl, 0, -r1, h };

    // Top
    x=-d/2; y=r1;
    p=newp; pleftlist3d[]+=p; Point(p) = { x, y, ctrl, h };
    p=newp; pleftlist3d[]+=p; Point(p) = { x, y, -ctrl, h };
    // Left
    ctrl=r1/2;
    x=-d/2-r1; y=0;
    p=newp; pleftlist3d[]+=p; Point(p) = { x, y, ctrl, h };
    p=newp; pleftlist3d[]+=p; Point(p) = { x, y, -ctrl, h };
    // Bottom
    x=-d/2; y=-r1;
    p=newp; pleftlist3d[]+=p; Point(p) = { x, y, ctrl, h };
    p=newp; pleftlist3d[]+=p; Point(p) = { x, y, -ctrl, h };

    // RIGHT BALL
    ctrl=r2/2;
    p=newp; Point(p) = { d/2, 0, r2, h };
    p=newp; Point(p) = { d/2, ctrl, r2, h };
    p=newp; Point(p) = { d/2, -ctrl, r2, h };
    p=newp; Point(p) = { d/2+ctrl, 0, r2, h };
    p=newp; Point(p) = { d/2-ctrl, 0, r2, h };

    p=newp; Point(p) = { d/2, 0, -r2, h };
    p=newp; Point(p) = { d/2, ctrl, -r2, h };
    p=newp; Point(p) = { d/2, -ctrl, -r2, h };
    p=newp; Point(p) = { d/2+ctrl, 0, -r2, h };
    p=newp; Point(p) = { d/2-ctrl, 0, -r2, h };
    // Top
    x=d/2; y=r2;
    p=newp; pleftlist3d[]+=p; Point(p) = { x, y, ctrl, h };
    p=newp; pleftlist3d[]+=p; Point(p) = { x, y, -ctrl, h };
    // Right
    x=d/2+r2; y=0;
    p=newp; pleftlist3d[]+=p; Point(p) = { x, y, ctrl, h };
    p=newp; pleftlist3d[]+=p; Point(p) = { x, y, -ctrl, h };
    // Bottom
    x=d/2; y=-r2;
    p=newp; pleftlist3d[]+=p; Point(p) = { x, y, ctrl, h };
    p=newp; pleftlist3d[]+=p; Point(p) = { x, y, -ctrl, h };
Return


//------------------------------------------------------------------------------
// DRAW3DLINES creates all lines for the 3D model.
//------------------------------------------------------------------------------
Function Draw3DLines
    Printf("-- Creating 3D Lines");
    // Horizontal.
    l=newl; Bezier(l)={26,29,38,8};
    l=newl; Bezier(l)={8,39,34,31};
    l=newl; BSpline(l)={26,30,25,46,42};
    l=newl; Bezier(l)={42,45,54,18};
    l=newl; Bezier(l)={18,55,50,47};
    l=newl; BSpline(l)={47,51,24,35,31};
    // Right ring.
    l=newl; Bezier(l)={21,53,48,47};
    l=newl; Bezier(l)={47,49,57,15};
    l=newl; Bezier(l)={15,56,44,42};
    l=newl; Bezier(l)={42,43,52,21};
    // Left ring.
    l=newl; Bezier(l)={5,36,27,26};
    l=newl; Bezier(l)={26,28,40,11};
    l=newl; Bezier(l)={11,41,33,31};
    l=newl; Bezier(l)={31,32,37,5};
Return

//------------------------------------------------------------------------------
// DRAW3DSURFACES creates all surfaces for the 3D model.
//------------------------------------------------------------------------------
Function Draw3DSurfaces
    Printf("-- Creating 3D Surfaces");
    ll3d[]={};
    l=newll; ll3d[]+=l; Line Loop(l) = {17, 7, -1};
    l=newll; ll3d[]+=l; Line Loop(l) = {1, 8, 20};
    l=newll; ll3d[]+=l; Line Loop(l) = {2, 19, -8};
    l=newll; ll3d[]+=l; Line Loop(l) = {2, -18, 7};
    l=newll; ll3d[]+=l; Line Loop(l) = {9, 16, 6, 17};
    l=newll; ll3d[]+=l; Line Loop(l) = {18, 3, 15, -9};
    l=newll; ll3d[]+=l; Line Loop(l) = {19, -12, 14, -3};
    l=newll; ll3d[]+=l; Line Loop(l) = {20, -6, 13, 12};
    l=newll; ll3d[]+=l; Line Loop(l) = {14, 4, 11};
    l=newll; ll3d[]+=l; Line Loop(l) = {13, -11, 5};
    l=newll; ll3d[]+=l; Line Loop(l) = {5, -16, 10};
    l=newll; ll3d[]+=l; Line Loop(l) = {10, -4, 15};

    N=#ll3d[];
    surf3d[]={};
    For k In {0:N-1}
       s=news; surf3d[]+=s; Ruled Surface(s) = {ll3d[k]};
    EndFor
Return

//------------------------------------------------------------------------------
// DRAW3DVOLUMES creates the volume for the 3D model.
//------------------------------------------------------------------------------
Function Draw3DVolumes
    Printf("-- Creating 3D Volumes");
    s=news; Surface Loop(s) = surf3d[];
    v=newv; Volume(v)={s};
Return

//------------------------------------------------------------------------------
// ADD3DPHYSICALS gives a name for the domain and its boundary (3D case).
//------------------------------------------------------------------------------
Function Add3DPhysicals
    Printf("-- Add 2D Physicals");
    Physical Surface( "Wall" ) = surf3d[];
    Physical Volume( "Omega" ) = { v };
Return

//------------------------------------------------------------------------------
// MAIN.
//------------------------------------------------------------------------------
Call Draw2DPoints;
Call Draw2DLines;
If( dim==2 )
    Call Draw2DSurfaces;
    Call Add2DPhysicals;
EndIf
If( dim==3 )
    Call Draw3DPoints;
    Call Draw3DLines;
    Call Draw3DSurfaces;
    Call Draw3DVolumes;
    Call Add3DPhysicals;
EndIf

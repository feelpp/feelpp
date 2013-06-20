# -*- mode: octave -*-
# load function with feel++ dof numbering
uf;
F1=F;
# load function with gmsh dof numbering
ug;
F2=F;
# load correspondance between feel++ and gmsh
# first column Feel++ and second column gmsh
M=load("feelpp2msh.m");
# renumber F2(gmsh) so that it matches feel++ numbering
F3=F2(M(:,2));

# load A.m, the matrix is named from filename 'A'
A;
# load b.m, the vector is named from filename 'b'
b;

u=A\b;

# check that matlab solve and feel++ solve do the same thing
assert(max(abs(u-F1)),0,1e-13);
assert(max(abs(F3-F1)),0,1e-13);

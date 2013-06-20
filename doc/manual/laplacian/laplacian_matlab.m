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

assert(max(abs(F3-F1)),0,1e-13);

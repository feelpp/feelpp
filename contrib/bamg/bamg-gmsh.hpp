Mesh2 *bamg2msh( bamg::Triangles* tTh,bool renumbering);
bamg::Triangles * msh2bamg(const Mesh2 & Th,double cutoffradian,long * reqedgeslab,int nreqedgeslab);
bamg::Triangles * msh2bamg(const Mesh2 & Th,double cutoffradian, 
                           int  nbdfv, int * ndfv,int  nbdfe, int * ndfe,
			   long * reqedgeslab,int nreqedgeslab);
Mesh2 *Bamg(Mesh2 *Thh, double * args,double *mm11,double *mm12,double *mm22);

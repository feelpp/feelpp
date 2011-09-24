// -*- Mode : c++ -*-
//
// SUMMARY  :
// USAGE    :
// ORG      :
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
//

/*

 This file is part of Freefem++

 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.


 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
#include <iostream>
using namespace std;
#include <stdio.h>
#include <string.h>


#include "Mesh2.h"
#include "RNM.hpp"
#include "Mesh2d.hpp"
#include "bamg-gmsh.hpp"

#include <set>
void ExecError(const char *s) { cout << s << endl;exit(1);}


Mesh2 *bamg2msh( bamg::Triangles* tTh,bool renumbering)
{
  using  bamg::Triangles;
  using  bamg::Triangle;
  using  bamg::Vertex;
  using  bamg::TriangleAdjacent;
  using  bamg::EdgesVertexTriangle;
  using  bamg::Int4;
  using  bamg::VerticesOfTriangularEdge;

  bamg::Triangles & th (*tTh);
  tTh->ReNumberingTheTriangleBySubDomain(!renumbering);//  just compress
  //tTh->NbRef++;
  Int4  i,j,k=0;
  int nv  =  tTh->nbv;
  int nt  =   tTh->nbt - tTh->NbOutT;
  int neb =   tTh->nbe;

  int nbcrakev = 0;
  tTh->ReMakeTriangleContainingTheVertex();
  Triangle2 * t =  new Triangle2[nt]  ;
  Seg * b_e = new Seg[neb];

  Vertex2 vbase;
  Vertex2 *vb(&vbase);
  if (verbosity>5)
    cout << "  -- Before cracking mesh:  Nb Triangles = " << nt << " Nb of Vertices " << nv << endl;
  for (int iv=0;iv<th.nbv;iv++) // vertex
    {
      // cout << iv << " : " ;
      const Vertex & v(th[iv]);
      int kk=0; // nb cracked
      int kc=0;
      int kkk =0; // nb triangle  with same number
      Triangle * tbegin = v.t;
      Vertex2 * vv = vb+iv;
      int i  = v.vint;
      assert(tbegin && (i >= 0 ) && (i <3));
      // turn around the vertex v
      TriangleAdjacent ta(tbegin,EdgesVertexTriangle[i][0]);// previous edge
      int k=0;
      do {
        int kv = VerticesOfTriangularEdge[ta][1];
        k++;
        Triangle * tt (ta);
        assert( &v == & (*  tt)[kv] );
        if ( ta.Cracked() )
          {   // cout << " || "    ;
            if ( kk == 0) tbegin=ta,kkk=0;  //  begin by a cracked edge  => restart
            if (  kkk ) { kc =1;vv = vb +  nv++;  kkk = 0; } // new vertex if use
            kk++;
            // number of cracked edge view
          }
        if ( tt->link ) { // if good triangles store the value
          int it = th.Number(tt);
          assert(it < nt);
          //int iiv=vv-vb;
          t[it](kv) = vv;
          /*
          cout << it << " " << kv << " "<< iiv  << endl;
          if (&th(it)[kv] != &th[iiv])
             cout << it << " " << kv << " "<< iiv << " != " << th.Number(th(it)[kv]) << endl ;
          */
          kkk++;
        } else if (kk) { // crack + boundary
          if (  kkk ) { kc =1;vv = vb +  nv++;  kkk = 0; } // new vertex if use
        }

        ta = Next(ta).Adj();
      } while ( (tbegin != ta));
      assert(k);
      if (kc)  nbcrakev++;
    }
  Vertex2 * v = new Vertex2[nv];
  //  set the vertices --
  for (i=0;i<nt;i++)
    {
      for (j=0;j<3;j++)
        {
          assert( t[i](j) );
          int k = t[i](j) - vb;
          t[i](j) = v+ k;
          assert(k>=0 && k < nv);
          Vertex & thv(th(i)[j]);
          v[k].x    =  thv.r.x;
          v[k].y    =  thv.r.y;
          v[k].lab  =  thv.ref();
        }
    }
  // warning in cracked edges
  // construction of the edges --

  if (nbcrakev && verbosity>2)
    cout << "  -- Nb of craked vertices = " << nbcrakev << " Nb of created vertices " << nv - th.nbv << endl;


  for (i=0;i<tTh->nbe;i++)
    {
        int ii[]={(int)tTh->Number(tTh->edges[i][0]),(int)tTh->Number(tTh->edges[i][1])};
      assert(ii[0]>=0 && ii[0] <nv);
      assert(ii[1]>=0 && ii[1] <nv);
      b_e[i].init(v,ii,tTh->edges[i].ref);
    }
  Int4 *reft = new Int4[tTh->nbt];
  //Int4 nbref =
  tTh->ConsRefTriangle(reft);
  for( i=0,k=0;i<tTh->nbt;i++)
    if(tTh->triangles[i].link)
      {

        R2 A(t[k][0]),B(t[k][1]),C(t[k][2]);
        t[k].area = (( B-A)^(C-A))*0.5 ;
        t[k].lab = tTh->subdomains[reft[i]].ref;  // a faire
        assert(k == i);
        k++;
      }
  delete [] reft;
  assert ( nt == k);
  tTh->ReMakeTriangleContainingTheVertex();

  if (verbosity)
    cout << "  --  mesh:  Nb of Triangles = "  << setw(6) <<  nt << ", Nb of Vertices " << nv << endl;

  {
    Mesh2 *m = new Mesh2(nv,nt,neb,v,t,b_e);
    //    if (renumbering) m->renum();
    //< m->MakeQuadTree();
    return m;
  }
}



bamg::Triangles * msh2bamg(const Mesh2 & Th,double cutoffradian,long * reqedgeslab,int nreqedgeslab)

{
    using namespace bamg;
    Triangles *Tn=new Triangles(Th.nv);
    Tn->nbv = Th.nv;
    Tn->nbt = Th.nt;
    Tn->nbe = Th.nbe;
    Tn->name= new char[strlen("msh2bamg")+1];
    strcpy(Tn->name,"msh2bamg");
    //  Tn->triangles = new Triangle [Tn->nbtx];
    assert(Tn->triangles);
    //  Tn->vertices = new Vertex [Tn->nbvx];
    //  Tn->ordre = new (Vertex* [Tn->nbvx]);
    Tn->edges = new Edge [Th.nbe];

    Int4 i;
    Metric Mid(1.);
    for (i = 0; i < Th.nv; i++)
      {
	Tn->vertices[i].r.x = Th(i).x;
	Tn->vertices[i].r.y = Th(i).y;
	Tn->vertices[i].m=Mid;
	Tn->vertices[i].ReferenceNumber = Th(i).lab;
      }

    //  Int4 i1 [nbt],i2 [nbt],i3 [nbt];
    for (i = 0; i < Th.nt; i++)
      {
	int i1 = Th(Th[i][0]);
	int i2 = Th(Th[i][1]);
	int i3 = Th(Th[i][2]);
	Tn->triangles[i]= Triangle( Tn,i1 ,i2 ,i3 );
	Tn->triangles[i].color = Th[i].lab;
      }
    //  Real8 cutoffradian = -1;
    // add code   un change boundary part ...  frev 2009 JYU FH
    set<int> labreq;
    if(nreqedgeslab && verbosity) cout << " label of required edges " ;
    for (int i=0; i <nreqedgeslab;++i)
      {
	if(verbosity)
	    cout << " " << reqedgeslab[i];
	labreq.insert(reqedgeslab[i]);
      }
    bamg::GeometricalEdge paszero;  // add JYU    fevr 2009   for  required edge ....
    if(nreqedgeslab && verbosity) cout << endl;
    int k=0;
    for (i = 0; i < Th.nbe; i++)
      {
	Tn->edges[i].v[0] = Tn->vertices + Th(Th.be(i)[0]);
	Tn->edges[i].v[1] = Tn->vertices + Th(Th.be(i)[1]);
	Tn->edges[i].ref = Th.be(i).lab;
	Tn->edges[i].on = 0;
	if( labreq.find( Tn->edges[i].ref) != labreq.end())
	  {
	    k++;
	    Tn->edges[i].on = &paszero;
	  }

      }
    if(verbosity)cout << "  number of required edges : "<< k << endl;


    Tn->ConsGeometry(cutoffradian);
    Tn->Gh.AfterRead();
    Tn->SetIntCoor();
    Tn->FillHoleInMesh();
    return Tn;
}


bamg::Triangles * msh2bamg(const Mesh2 & Th,double cutoffradian,
                           int  nbdfv, int * ndfv,int  nbdfe, int * ndfe,
			   long * reqedgeslab,int nreqedgeslab)
{
    using namespace bamg;
    Triangles *Tn=new Triangles(Th.nv);
    KN<int> equiedges(Th.nbe);
    for(int i=0;i<Th.nbe;i++)
	equiedges[i]=2*i;
    if(nbdfe !=0 )
      {
	KN<int>  kk(Th.nbe),kn(Th.nbe);
	kk=0;
	for(int i=0;i<Th.nbe;i++)
	  {
	    int df=ndfe[i];
	    kk[df]++;
	    if(kk[df]==1) kn[df]=i;
	    else {
		int k=kn[df],sens=0;
		int di0=ndfv[Th(Th.be(i)[0])];
		int di1=ndfv[Th(Th.be(i)[1])];
		int dk0=ndfv[Th(Th.be(k)[0])];
		int dk1=ndfv[Th(Th.be(k)[1])];
		if ((di0==dk0) &&(di1==dk1) ) sens=0;
		else if ((di1==dk0) &&(di0==dk1) ) sens=1;
		else  {
		    cout << "Error in periodic mesh " << di0 << " " << di1 << " <=> " << dk0 << " " << dk1 << endl;
		    ExecError("bug periodic mesh in ??? ");
		}
		equiedges[i]=2*k+sens;

	    }
	  }

      }; // a faire pour les maillages periodique

    Tn->nbv = Th.nv;
    Tn->nbt = Th.nt;
    Tn->nbe = Th.nbe;
    Tn->name= new char[strlen("msh2bamg")+1];
    strcpy(Tn->name,"msh2bamg");
    //  Tn->triangles = new Triangle [Tn->nbtx];
    assert(Tn->triangles);
    //  Tn->vertices = new Vertex [Tn->nbvx];
    //  Tn->ordre = new (Vertex* [Tn->nbvx]);
    Tn->edges = new Edge [Th.nbe];

    Int4 i;
    Metric Mid(1.);
    for (i = 0; i < Th.nv; i++)
      {
	Tn->vertices[i].r.x = Th(i).x;
	Tn->vertices[i].r.y = Th(i).y;
	Tn->vertices[i].ReferenceNumber = Th(i).lab;
	Tn->vertices[i].m=Mid;
      }

    //  Int4 i1 [nbt],i2 [nbt],i3 [nbt];
    for (i = 0; i < Th.nt; i++)
      {
	int i1 = Th(Th[i][0]);
	int i2 = Th(Th[i][1]);
	int i3 = Th(Th[i][2]);
	Tn->triangles[i]= Triangle( Tn,i1 ,i2 ,i3 );
	Tn->triangles[i].color = Th[i].lab;
      }

    // add code   un change boundary part ...  frev 2009 JYU FH
    set<int> labreq;
    if(nreqedgeslab && verbosity) cout << " label of required edges " ;
    for (int i=0; i <nreqedgeslab;++i)
      {
	if(verbosity)
	    cout << " " << reqedgeslab[i];
	labreq.insert(reqedgeslab[i]);
      }
    bamg::GeometricalEdge paszero;  // add JYU    fevr 2009   for  required edge ....
    if(nreqedgeslab && verbosity) cout << endl;
    int k=0;

    for (i = 0; i < Th.nbe; i++)
      {
	Tn->edges[i].v[0] = Tn->vertices + Th(Th.be(i)[0]);
	Tn->edges[i].v[1] = Tn->vertices + Th(Th.be(i)[1]);
	Tn->edges[i].ref = Th.be(i).lab;
	Tn->edges[i].on = 0;
	if( labreq.find( Tn->edges[i].ref) != labreq.end())
	  {
	    k++;
	    Tn->edges[i].on = &paszero;
	  }
      }
    //  Real8 cutoffradian = -1;
    Tn->ConsGeometry(cutoffradian,equiedges);
    Tn->Gh.AfterRead();
    Tn->SetIntCoor();
    Tn->FillHoleInMesh();
    return Tn;
}


template<class T> T arg(int i,double *args,const T & d)
{
    return args[i]<= -1.e100  ? d: (T) args[i];
}
Mesh2 *Bamg(Mesh2 *Thh, double * args,double *mm11,double *mm12,double *mm22, bool initialMesh)
{
  using namespace bamg;
  using RNM::Min;
  using RNM::Max;
  using RNM::Abs;


  Real8 err         = arg(2,args,0.01);    // coef in the metric
  Real8  errg       = Min(arg(3,args,0.01),err);
  long nbsx         = Max(100L,arg(4,args,900000L));
  long nbsmooth     =  arg(5,args,3L);
  long nbjacobi     =  arg(6,args,0L) ;              // if increased will be more smooth
  const Real8 raison = arg(7,args,1.8); // 1.8
  const Real8 omega =  arg(8,args,1.0) ;
  bool iso          =   arg(9,args,false);
  bool AbsError     =   arg(10,args,true);
  Real8 CutOff      =   arg(11,args, 1.0e-6);
  verbosity         =   arg(12,args, (long) verbosity);
  bool inq          =   arg(13,args,false);
  bool SplitEdgeWith2Boundary =  arg(14,args,true);
  double maxsubdiv             = Max(Min( arg(15,args,10.0),10.0),0.1);
  double anisomax              =   Max((double) arg(16,args,1.0e6),1.0);
  bool rescaling               = arg(17,args,true) ;
  bool KeepBackVertices        =   arg(18,args,true) ;
  int givenmetric              =  arg(19,args,false) ;
  double powerM                = arg(20,args,1.0) ;
  double cutoffradian          = arg(21,args,-1.0)* bamg::Pi/180. ;
  bool split                    = arg(22,args,false) ;
  bool nomeshgeneration         = arg(23,args,false) ;

  //  printf("raison = %g %g\n",raison,args[7]);

  //   the 24th param is metrix  and is store at compilation time
  //  const E_Array * expmetrix = dynamic_cast<const E_Array *>(nargs[24]);
  //   the 25th param is periodic and it store at compilation time
  // in nbcperiodic,periodic  variable
  //  KN<long> reqedges0;
  //  list of label of required edges , for no adapattion on this part of the boundary.
  //    KN<long> reqedges ( nargs[26] ? GetAny< KN_<long> >( (*nargs[26])(stack) ): (KN_<long>)reqedges0);

  KN<long> reqedges(Thh->nbe);
  for (int i=0;i<Thh->nbe;i++)reqedges[i]=Thh->be(i).lab;

  if(reqedges.N() && verbosity)
      cout << " reqedges labels "  << reqedges << endl;
  //  KN<double> *mm11=0, *mm12=0,* mm22=0;

  // using MeshPoint;
  //using Mesh;
  ffassert(Thh);
    Triangles * oTh =0;
    int nbcperiodic=0;
  /*   A change ...
  if (nbcperiodic) {
    KN<int> ndfv(Thh->nv);
    KN<int> ndfe(Thh->nbe);
    int nbdfv=0,nbdfe=0;
    BuildPeriodic(nbcperiodic,periodic,*Thh,stack,nbdfv,ndfv,nbdfe,ndfe);
    oTh = msh2bamg(*Thh,cutoffradian,nbdfv,ndfv,nbdfe,ndfe,reqedges,reqedges.N());
  }
  else */
   oTh = msh2bamg(*Thh,cutoffradian,reqedges,reqedges.N());
  Triangles &Th(*oTh);

  //  printf("COUCOUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU\n");
  //  Th.Write("toto.mesh",bamg::Triangles::AutoMesh);


  bool mtx=true;
  KN_<double> m11(mm11,Thh->nv);
  KN_<double> m12(mm12,Thh->nv);
  KN_<double> m22(mm22,Thh->nv);

  Real8 hmax = 0.3*Th.MaximalHmax(); // final largest edge
  Real8 hmin = Th.MinimalHmin();        // final smallest edge


  Real8 coef =1;                // a priori don't touch
  // gestion des arguments
  hmin              = Max(hmin, arg(0,args,hmin));
  hmax              = Min(hmax,arg(1,args,hmax));


  if (iso)  anisomax=1;
  if (verbosity>2)
    {
      cout << endl  << endl;
      cout << " \t\t ## adapt :   nbsx = " << nbsx << ", err = " << err ;
      cout << ", hmin = " << hmin << ", hmax = " << hmax <<endl;
      cout << " \t\t    ratio  = " << raison << ", nbsmooth = " << nbsmooth ;
      cout << ", omega = " << omega <<  ", coef = " << coef << ", iso = " << iso << endl;
      cout << " \t\t    AbsError =" << AbsError << ", CutOff = " << CutOff << ", nbjacobi = " << nbjacobi <<endl;
      cout << " \t\t    maxsubdiv = " << maxsubdiv << " splitpbedge = " << SplitEdgeWith2Boundary  <<endl;
      cout << " \t\t    anisomax = " << anisomax << ", rescaling = " << rescaling << ", power = " << powerM
           << ", KeepBackvertices = " << KeepBackVertices << " IsMetric = " << givenmetric
      << endl << endl ;
    }

 //
 Th.ReMakeTriangleContainingTheVertex();
    /*
  //MeshPoint* mp(MeshPointStack(stack));

  Int4 i,iv;
  int ksol =0;
  for (i=0;i<nbsol;i++)
    ksol += typesol[i]+1; //  marche en 2d

  double * lessol = new double [Th.nbv*ksol];
  double *ss = lessol;
  // be careful because renum --
  // the triangle was no renum
  for ( iv=0;iv<Th.nbv;iv++)
    Th[iv].color=1; // color
 for (Int4  it = 0; it < Thh->nt; it++)
    for (Int4  jt = 0; jt < 3; jt++)
      {
        bamg::Vertex & v= Th(it)[jt];
	//   const Vertex & vf = (*Thh)[it][jt];
        if (&v && v.color)
          {
            v.color =0; // uncolor
            mp->setP(Thh ,it,jt);

            ss = lessol + ksol* Th.Number(v);
            for (int j =0; j < ksol; j++)
              *ss++= GetAny<double>( (*sol[j])(stack) );

          }
      }
  mp->unset();
     */
  // computation of the metric ---
  // better thing -> create keyword in the language
  //    a faire F Hecht .
  Metric Mhmax(hmax);
  for (int iv=0;iv<Th.nbv;iv++)
    Th[iv].m = Mhmax;

   if (mtx)
     for (int iv=0;iv<Th.nbv;iv++) {
       //      if ( Max(m11[iv],m12[iv],m22[iv]) > hmax)
       Th[iv].m.IntersectWith(MetricAnIso(m11[iv],m12[iv],m22[iv]));
     }
  /*
  if ( givenmetric)
    if (ksol == 1)
      {
        for (Int4  iv = 0,k=0; iv < Th.nbv ; iv++)
          Th[iv].m.IntersectWith(Metric(lessol[k++]));
      }
    else if (ksol == 3)
      {
        for (Int4  iv = 0,k=0; iv < Th.nbv ; iv++, k += 3)
          {
	    Metric MM(lessol[k],lessol[k+1],lessol[k+2]);
	    MatVVP2x2 vp(MM);
	    vp.Abs();
	    Th[iv].m.IntersectWith(vp);
          }
      }
    else
      lgerror("Adapt mesh: ksol is wrong, IsMetric  and ksol != 1 or 3");
  else
  Th.IntersectConsMetric(lessol,nbsol,typesol,hmin,hmax,sqrt(err)*coef,anisomax,AbsError?0.0:CutOff,nbjacobi,rescaling,powerM,0);

    delete [] lessol;
   */
  //Th.IntersectGeomMetric(errg,iso);

  Th.SmoothMetric(raison);
  Th.MaxSubDivision(maxsubdiv);
  Th.BoundAnisotropy(anisomax);
  // end of metric's computation
   if (mtx)
    for (int iv=0;iv<Th.nbv;iv++)
      {
       m11[iv] = Th[iv].m.a11  ;
       m22[iv] = Th[iv].m.a22 ;
       m12[iv] = Th[iv].m.a21;
      }

  Triangles* nTh = 0;

   //   Th.Write("toto.msh",bamg::Triangles::AutoMesh);

  if (initialMesh){
    nTh= new Triangles(nbsx,Th.Gh);
  }
  else {
    nTh= new Triangles(nbsx,Th,KeepBackVertices); // Adaption is here
  }

  //  nTh->Write("tata.mesh",bamg::Triangles::AutoMesh);


  if (split)
    nTh->SplitElement(1); // modif FH mai 2009 (thank J-M Mirebeau) : Th ->nTh

  if(SplitEdgeWith2Boundary)
    nTh->SplitInternalEdgeWithBorderVertices();
  if(verbosity>3)
    nTh->ShowHistogram();
  if (nbsmooth)
    nTh->SmoothingVertex(nbsmooth,omega);
  if(verbosity>2 && nbsmooth)
    nTh->ShowHistogram();
  if(verbosity>0)
      nTh->ShowRegulaty()  ;

  inq=0;
  Metric M(hmax);
  for (int iv=0;iv < Th.nbv;iv++)
    Th[iv].m = M;

  Mesh2 * g=  bamg2msh(nTh,true);

  delete nTh;
  delete oTh;
 //  Add2StackOfPtr2FreeRC(stack,g);// 07/2008  FH

  return g;

}



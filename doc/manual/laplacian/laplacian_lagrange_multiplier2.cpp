// -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

#include <feel/feel.hpp>
#include <feel/feelalg/vectorblock.hpp>

namespace Feel
{

boost::shared_ptr<Mesh<Simplex<2> > >
createMeshLaplacianLM()
{
    typedef Mesh<Simplex<2> > mesh_type;
    double meshSize = option("gmsh.hsize").as<double>();
    GeoTool::Node x1(-2,-1);
    GeoTool::Node x2(2,1);
    GeoTool::Rectangle R( meshSize ,"OMEGA",x1,x2);
    //R.setMarker(_type="line",_name="Paroi",_marker3=true);
    R.setMarker(_type="line",_name="gamma",_markerAll=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);
    GeoTool::Node x3(0,0); //center
    GeoTool::Node x4(1);//majorRadiusParam
    GeoTool::Node x5(0.7);//minorRadiusParam
    GeoTool::Node x6(0.1);//penautRadiusParam
    GeoTool::Peanut P( meshSize ,"OMEGA2",x3,x4,x5,x6);
    P.setMarker(_type="line",_name="peanut",_markerAll=true);
    P.setMarker(_type="surface",_name="Omega2",_markerAll=true);
    auto mesh = (R-P).createMesh(_mesh=new mesh_type,_name="mymesh.msh");

    return mesh;
}


}


void runLaplacianLMV1()
{
    using namespace Feel;

    auto mesh = createMeshLaplacianLM();
    auto submesh = createSubmesh(mesh,boundaryfaces(mesh));

    auto Vh1 = Pch<2>(mesh);
    auto Vh2 = Pch<2>(submesh);
    auto u1 = Vh1->elementPtr();
    auto u2 = Vh2->elementPtr();

    BlocksBaseGraphCSR myblockGraph(2,2);
    myblockGraph(0,0) = stencil( _test=Vh1,_trial=Vh1, _diag_is_nonzero=false, _close=false)->graph();
    myblockGraph(1,0) = stencil( _test=Vh2,_trial=Vh1, _diag_is_nonzero=false, _close=false)->graph();
    myblockGraph(0,1) = stencil( _test=Vh1,_trial=Vh2, _diag_is_nonzero=false, _close=false)->graph();
    myblockGraph(1,1) = stencil( _test=Vh2,_trial=Vh2, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    auto A = backend()->newBlockMatrix(_block=myblockGraph);

    BlocksBaseVector<double> myblockVec(2);
    myblockVec(0,0) = backend()->newVector( Vh1 );
    myblockVec(1,0) = backend()->newVector( Vh2 );
    auto F = backend()->newBlockVector(_block=myblockVec, _copy_values=false);

    BlocksBaseVector<double> myblockVecSol(2);
    myblockVecSol(0,0) = u1;
    myblockVecSol(1,0) = u2;
    auto U = backend()->newBlockVector(_block=myblockVecSol, _copy_values=false);

    //auto g = vf::cos(vf::Px()*8*M_PI)*vf::sin(vf::Px()*8*M_PI);
    auto g = vf::cst(0.);

    form2( _trial=Vh1, _test=Vh1,_matrix=A)
        += integrate(_range=elements(mesh),
                     _expr=gradt(u1)*trans(grad(u1)) );

    form2( _trial=Vh2, _test=Vh1 ,_matrix=A,
           _rowstart=0, _colstart=Vh1->nLocalDofWithGhost() )
        += integrate( _range=elements(submesh),
                      _expr=idt(u2)*id(u1) );

    form2( _trial=Vh1, _test=Vh2 ,_matrix=A,
           _rowstart=Vh1->nLocalDofWithGhost(), _colstart=0 )
        += integrate( _range=elements(submesh),
                      _expr=idt(u1)*id(u2) );

    form1( _test=Vh1, _vector=F )
        = integrate(_range=elements(mesh),
                    _expr=id(u1));

    form1( _test=Vh2, _vector=F,
           _rowstart=Vh1->nLocalDofWithGhost() )
        += integrate(_range=elements(submesh),
                    _expr=g*id(u2));

    backend(_rebuild=true)->solve( _matrix=A, _rhs=F, _solution=U );

    myblockVecSol.localize(U);

    auto e = exporter( _mesh=mesh,_name="exportV1" );
    e->add( "u1V1", *u1 );
    e->save();
}

void runLaplacianLMV2()
{
    using namespace Feel;

    auto mesh = createMeshLaplacianLM();
    auto submesh = createSubmesh(mesh,boundaryfaces(mesh));

    typedef Mesh<Simplex<2> > mesh_type;
    typedef FunctionSpace<meshes<mesh_type,mesh_type::trace_mesh_type>, bases<Lagrange<2, Scalar>,Lagrange<2,Scalar> > > space_type;
    auto Vh = space_type::New(fusion::make_vector(mesh,submesh));

    auto U = Vh->element();
    auto u1 = U.element<0>();
    auto u2 = U.element<1>();

    auto g = vf::cst(0.);

    auto a = form2( _trial=Vh, _test=Vh);
    a = integrate(_range=elements(mesh),
                  _expr=gradt(u1)*trans(grad(u1)) );
    a += integrate( _range=elements(submesh),
                    _expr=idt(u2)*id(u1) );
    a += integrate( _range=elements(submesh),
                    _expr=idt(u1)*id(u2) );
    a += integrate(_range=elements(mesh),
                    _expr=id(u1));
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=id(u1));
    l += integrate(_range=elements(submesh),
                    _expr=g*id(u2));

    a.solve(_rhs=l,_solution=U);

    auto e = exporter( _mesh=mesh, _name="exportV2" );
    e->add( "u1V2", u1 );
    e->save();
}


int main(int argc, char**argv )
{
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="laplacian_lm2",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    runLaplacianLMV1();
    //runLaplacianLMV2();

    return 0;
}

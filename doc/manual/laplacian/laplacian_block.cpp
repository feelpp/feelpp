// -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

#include <feel/feel.hpp>

int main(int argc, char**argv )
{
    //# marker1 #
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="laplacian_block",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    GeoTool::Rectangle R( option(_name="mesh2d.hsize").as<double>(),"Rectangle",
                          GeoTool::Node(-1,-1),
                          GeoTool::Node( 1, 1) );
    R.setMarker(_type="line",_name="BoundaryA",_marker1=true);
    R.setMarker(_type="line",_name="BoundaryB",_marker2=true,_marker3=true,_marker4=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);
    auto mesh = R.createMesh(_mesh=new Mesh<Simplex<2> >,
                             _name="meshrect",
                             _hmax=option(_name="mesh2d.hsize").as<double>() );

    auto P0d = Pdh<0>(mesh);
    double r=0.2;
    auto proj = vf::project(_space=P0d,_range=elements(mesh),
                            _expr=vf::chi( (Px()*Px()+Py()*Py()) < r*r ) );
    mesh->updateMarker3( proj );
    auto submesh = createSubmesh( mesh, marked3elements(mesh,1) );
    //saveGMSHMesh(_mesh=submesh,_filename="mysubmesh.msh");

    auto Vh1 = Pch<1>( mesh );
    auto Vh2 = Pch<1>( submesh );
    auto u1 = Vh1->elementPtr();
    auto u2 = Vh2->elementPtr();

    BlocksBaseGraphCSR myblockGraph(2,2);
    myblockGraph(0,0) = stencil(_test=Vh1,_trial=Vh1, _diag_is_nonzero=false)->graph();
    myblockGraph(0,1) = stencil(_test=Vh1,_trial=Vh2, _diag_is_nonzero=false)->graph();
    myblockGraph(1,0) = stencil(_test=Vh2,_trial=Vh1, _diag_is_nonzero=false)->graph();
    myblockGraph(1,1) = stencil(_test=Vh2,_trial=Vh2, /*_pattern=(size_type)Pattern::ZERO,*/ _diag_is_nonzero=false)->graph();
    auto A = backend()->newBlockMatrix(_block=myblockGraph);

    BlocksBaseVector<double> myblockVec(2);
    myblockVec(0,0) = backend()->newVector( Vh1 );
    myblockVec(1,0) = backend()->newVector( Vh2 );
    auto F = backend()->newBlockVector(_block=myblockVec, _copy_values=false);

    BlocksBaseVector<double> myblockVecSol(2);
    myblockVecSol(0,0) = u1;
    myblockVecSol(1,0) = u2;
    auto U = backend()->newBlockVector(_block=myblockVecSol, _copy_values=false);

    form2( _trial=Vh1, _test=Vh1 ,_matrix=A )
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

    form2( _trial=Vh1, _test=Vh1 ,_matrix=A )
        +=on(_range=boundaryfaces(mesh), _rhs=F, _element=*u1,
        _expr=constant(0.) );

    backend()->solve( _matrix=A, _rhs=F, _solution=U );

    myblockVecSol.localize(U);

    auto e = exporter( _mesh=mesh );
    e->add( "u1", *u1 );
    e->save();

    return 0;
}

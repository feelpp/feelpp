// -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

#include <feel/feel.hpp>
#include <feel/feelalg/vectorblock.hpp>

namespace Feel
{

template <uint16_type OrderGeo>
boost::shared_ptr<Mesh<Simplex<2,OrderGeo> > >
createMeshStokesDirichletLM( mpl::int_<2> /**/ )
{
    typedef Mesh<Simplex<2,OrderGeo> > mesh_type;
    double meshSize = option("gmsh.hsize").as<double>();

    GeoTool::Node x1(0,-0.5);
    GeoTool::Special_1b R( meshSize,"unrectangle",x1);
    R.setMarker(_type="line",_name="inlet",_marker2=true);
    R.setMarker(_type="line",_name="wall",_marker1=true);
    R.setMarker(_type="line",_name="outlet",_marker3=true);
    R.setMarker(_type="surface",_name="OmegaFluide",_markerAll=true);

    auto mesh = R.createMesh(_mesh=new mesh_type,_name="mymesh2d.msh");

    return mesh;
}


template <uint16_type OrderGeo>
boost::shared_ptr<Mesh<Simplex<3,OrderGeo> > >
createMeshStokesDirichletLM( mpl::int_<3> /**/ )
{
    typedef Mesh<Simplex<3,OrderGeo> > mesh_type;
    double meshSize = option("gmsh.hsize").as<double>();

    GeoTool::Node Centre(0,0,0.205);
    GeoTool::Node Rayon( 0.205);
    GeoTool::Node Dir(1,0,0);
    GeoTool::Node Lg(2.5,0,0);
    GeoTool::Cylindre C( meshSize,"Cyl",Centre,Dir,Rayon,Lg);
    C.setMarker(_type="surface",_name="inlet",_marker1=true);
    C.setMarker(_type="surface",_name="outlet",_marker2=true);
    C.setMarker(_type="surface",_name="wall",_marker3=true);
    C.setMarker(_type="volume",_name="OmegaFluid",_markerAll=true);

    GeoTool::Special3D_1 S( meshSize/3.,"shapeStruct");
    S.setMarker(_type="surface",_name="wall",_markerAll=true);
    //S.setMarker(_type="surface",_name="fsiwall",_marker1=true);
    //S.setMarker(_type="surface",_name="wallcylinder",_marker2=true);
    S.setMarker(_type="volume",_name="OmegaStruct",_markerAll=true);

    auto mesh = (C-S).createMesh(_mesh=new mesh_type,
                                 _name="mymesh3d.msh",
                                 _hmax=meshSize );

    return mesh;
}

decltype( (Py()-1)*(Py()+1)*N() )
inletVelocityExpr( mpl::int_<2> /**/ )
{
    return (Py()-1)*(Py()+1)*N();
}

decltype( 1.5*2*(4./0.1681)*( pow(Py()-0 ,2) + pow(Pz()-0.205 ,2) - cst(0.205*0.205) )*N() )
inletVelocityExpr( mpl::int_<3> /**/ )
{
    return 1.5*2*(4./0.1681)*( pow(Py()-0 ,2) + pow(Pz()-0.205 ,2) - cst(0.205*0.205) )*N();
}


template <uint16_type Dim,uint16_type OrderGeo>
void runStokesDirichletLM()
{
    using namespace Feel;

    std::string configstr = (boost::format("%1%dGeo%2%")%Dim %OrderGeo).str();

    auto mesh = createMeshStokesDirichletLM<OrderGeo>( mpl::int_<Dim>() );
    std::list<std::string> listMarker{"inlet","wall"};
    auto submesh = createSubmesh(mesh,markedfaces(mesh,listMarker));

    auto Vh1 = THch<OrderGeo>(mesh);
    auto Vh2 = Pchv<2>(submesh);

    if (Environment::worldComm().isMasterRank())
    {
        std::cout << "mesh->numGlobalElements() "<<mesh->numGlobalElements() << std::endl;
        std::cout << "submesh->numGlobalElements() "<<submesh->numGlobalElements() << std::endl;
        std::cout << "Vh1->nDof() "<<Vh1->nDof() << std::endl;
        std::cout << "Vh2->nDof() "<<Vh2->nDof() << std::endl;
    }

    auto U = Vh1->elementPtr();
    auto u = U->template element<0>();
    auto p = U->template element<1>();
    auto lambda = Vh2->elementPtr();
    BlocksBaseGraphCSR myblockGraph(2,2);

    myblockGraph(0,0) = stencil( _test=Vh1,_trial=Vh1, _diag_is_nonzero=false, _close=false)->graph();
    myblockGraph(1,0) = stencil( _test=Vh2,_trial=Vh1, _diag_is_nonzero=false, _close=false)->graph();
    myblockGraph(0,1) = stencil( _test=Vh1,_trial=Vh2, _diag_is_nonzero=false, _close=false)->graph();
    //myblockGraph(1,1) = stencil( _test=Vh2,_trial=Vh2, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    auto A = backend()->newBlockMatrix(_block=myblockGraph);

    BlocksBaseVector<double> myblockVec(2);
    myblockVec(0,0) = backend()->newVector( Vh1 );
    myblockVec(1,0) = backend()->newVector( Vh2 );
    auto F = backend()->newBlockVector(_block=myblockVec, _copy_values=false);

    BlocksBaseVector<double> myblockVecSol(2);
    myblockVecSol(0,0) = U;
    myblockVecSol(1,0) = lambda;
    auto UVec = backend()->newBlockVector(_block=myblockVecSol, _copy_values=false);

    auto deft = gradt( u );
    auto def = grad( u );

    auto u_in = inletVelocityExpr( mpl::int_<Dim>() );

    double mu=1;
    auto stokes = form2( _test=Vh1, _trial=Vh1, _matrix=A );
    stokes += integrate( _range=elements( mesh ),
                         _expr=mu*inner( deft,def ) );
    stokes += integrate( _range=elements( mesh ),
                         _expr=-div(u)*idt(p) + divt(u)*id(p) );

    form2( _trial=Vh2, _test=Vh1 ,_matrix=A,
           _rowstart=0, _colstart=Vh1->nLocalDofWithGhost() )
        += integrate( //_range=elements(submesh),
                     _range=markedfaces(mesh,listMarker),
                     _expr=inner(idt(lambda),id(u)) );

    form2( _trial=Vh1, _test=Vh2 ,_matrix=A,
           _rowstart=Vh1->nLocalDofWithGhost(), _colstart=0 )
        += integrate( _range=elements(submesh),
                      //_range=markedfaces(mesh,listMarker),
                      _expr=inner(idt(u),id(lambda)) );

    form1( _test=Vh2, _vector=F,
           _rowstart=Vh1->nLocalDofWithGhost() )
        += integrate(//_range=markedelements(submesh,"inlet"),
                     _range=markedfaces(mesh,"inlet"),
                     _expr=inner(u_in,id(lambda) ) );

    backend(_rebuild=true)->solve( _matrix=A, _rhs=F, _solution=UVec );

    double normL1_A = A->l1Norm();
    double normL2_F = F->l2Norm();
    double normL2_U = UVec->l2Norm();

    if (Environment::worldComm().isMasterRank())
        std::cout << " normL1_A " << std::setprecision( 9 ) << normL1_A
                  << " normL2_F " << std::setprecision( 9 ) << normL2_F
                  << " normL2_U " << std::setprecision( 9 ) << normL2_U
                  << std::endl;

    myblockVecSol.localize(UVec);

    if (OrderGeo==1)
    {
        auto e = exporter( _mesh=mesh,_name="export"+configstr );
        e->add( "u"+configstr, u );
        e->add( "p"+configstr, p );
        e->save();
    }
    else
    {
        auto meshVisu = lagrangeP1(_space=Vh1->template functionSpace<0>())->mesh();
        auto XhVisuVel = Pchv<1>(meshVisu);
        auto XhVisuPressure = Pch<1>(meshVisu);
        auto opIVisuVel = opInterpolation(_domainSpace=Vh1->template functionSpace<0>(),
                                          _imageSpace=XhVisuVel,
                                          _type=InterpolationNonConforme(false,true,false) );
        auto uVisu = opIVisuVel->operator()(u);
        auto opIVisuPressure = opInterpolation(_domainSpace=Vh1->template functionSpace<1>(),
                                               _imageSpace=XhVisuPressure,
                                               _type=InterpolationNonConforme(false,true,false) );
        auto pVisu = opIVisuPressure->operator()(p);
        auto e = exporter( _mesh=meshVisu,_name="export"+configstr );
        e->add( "uVisu"+configstr, uVisu );
        e->add( "pVisu"+configstr, pVisu );
        e->save();

    }
}

} // namespace Feel

int main(int argc, char**argv )
{
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="stokes_dirichletlm",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));


    //runStokesDirichletLM<2,1>();
    runStokesDirichletLM<3,1>();

    return 0;
}

// -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

#include <feel/feel.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feelvf/print.hpp>
#include <feel/feeldiscr/pch.hpp>

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

    auto mesh = C.createMesh(_mesh=new mesh_type,
                             _name="mymesh3d",
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
    std::list<std::string> presslm{"inlet","outlet"};
    auto submesh = createSubmesh(mesh,markedfaces(mesh,presslm));

    auto Vh1 = THch<OrderGeo>(mesh);
    auto Vh21 = Pch<2,double, PointSetEquiSpaced,Mesh<Simplex<2,1,3>>,0>(submesh);
    auto Vh22 = Pch<2,double, PointSetEquiSpaced,Mesh<Simplex<2,1,3>>,1>(submesh);

    if (Environment::worldComm().isMasterRank())
    {
        std::cout << "mesh->numGlobalElements() "<<mesh->numGlobalElements() << std::endl;
        std::cout << "submesh->numGlobalElements() "<<submesh->numGlobalElements() << std::endl;
        std::cout << "Vh1->nDof() "<<Vh1->nDof() << std::endl;
        std::cout << "Vh21->nDof() "<<Vh21->nDof() << std::endl;
        std::cout << "Vh22->nDof() "<<Vh22->nDof() << std::endl;
    }

    auto U = Vh1->elementPtr();
    auto u = U->template element<0>();
    auto p = U->template element<1>();
    auto lambda1 = Vh21->elementPtr();
    auto lambda2 = Vh22->elementPtr();
    BlocksBaseGraphCSR myblockGraph(3,3);

    myblockGraph(0,0) = stencil( _test=Vh1,_trial=Vh1, _diag_is_nonzero=false, _close=false)->graph();
    myblockGraph(1,0) = stencil( _test=Vh21,_trial=Vh1, _diag_is_nonzero=false, _close=false)->graph();
    myblockGraph(0,1) = stencil( _test=Vh1,_trial=Vh21, _diag_is_nonzero=false, _close=false)->graph();
    myblockGraph(1,1) = stencil( _test=Vh21,_trial=Vh21, _diag_is_nonzero=true, _close=false)->graph();
    myblockGraph(2,0) = stencil( _test=Vh22, _trial=Vh1, _diag_is_nonzero=false, _close=false)->graph();
    myblockGraph(0,2) = stencil( _test=Vh1,_trial=Vh22, _diag_is_nonzero=false, _close=false)->graph();
    myblockGraph(2,2) = stencil( _test=Vh22,_trial=Vh22, _diag_is_nonzero=true, _close=false)->graph();

    auto A = backend()->newBlockMatrix(_block=myblockGraph);

    auto a_01 = form2( _trial=Vh21, _test=Vh1 ,_matrix=A,
                       _rowstart=0, _colstart=Vh1->nLocalDofWithGhost() );
    auto a_10 = form2( _test=Vh21, _trial=Vh1 ,_matrix=A,
                       _rowstart=Vh1->nLocalDofWithGhost(), _colstart=0 );
    auto a_11=form2( _test=Vh21, _trial=Vh21 ,_matrix=A,
                     _rowstart=Vh1->nLocalDofWithGhost(),
                     _colstart=Vh1->nLocalDofWithGhost() );
    auto a_02 = form2( _trial=Vh22, _test=Vh1 ,_matrix=A,
                       _rowstart=0, _colstart=Vh1->nLocalDofWithGhost()+Vh21->nLocalDofWithGhost() );
    auto a_20 = form2( _test=Vh22, _trial=Vh1 ,_matrix=A,
                       _rowstart=Vh1->nLocalDofWithGhost()+Vh21->nLocalDofWithGhost(), _colstart=0 );
    auto a_22=form2( _test=Vh22, _trial=Vh22 ,_matrix=A,
                     _rowstart=Vh1->nLocalDofWithGhost()+Vh21->nLocalDofWithGhost(),
                     _colstart=Vh1->nLocalDofWithGhost()+Vh21->nLocalDofWithGhost() );


    BlocksBaseVector<double> myblockVec(3);
    myblockVec(0,0) = backend()->newVector( Vh1 );
    myblockVec(1,0) = backend()->newVector( Vh21 );
    myblockVec(2,0) = backend()->newVector( Vh22 );
    auto F = backend()->newBlockVector(_block=myblockVec, _copy_values=false);

    BlocksBaseVector<double> myblockVecSol(3);
    myblockVecSol(0,0) = U;
    myblockVecSol(1,0) = lambda1;
    myblockVecSol(2,0) = lambda2;
    auto UVec = backend()->newBlockVector(_block=myblockVecSol, _copy_values=false);

    auto deft = sym(gradt( u ));
    auto def = sym(grad( u ));

    auto u_in = inletVelocityExpr( mpl::int_<Dim>() );

    double mu=1;
    auto stokes = form2( _test=Vh1, _trial=Vh1, _matrix=A );
    stokes += integrate( _range=elements( mesh ),
                         _expr=2*mu*inner( deft,def ) );
    stokes += integrate( _range=elements( mesh ),
                         _expr=-div(u)*idt(p) - divt(u)*id(p) );
    // total stress tensor (trial)
    auto SigmaNt = -idt( p )*N()+2*mu*deft*N();
    // total stress tensor (test)
    auto SigmaN = -id( p )*N()+2*mu*def*N();
    stokes +=integrate( _range=markedfaces( mesh, "wall" ), _expr=- trans(SigmaNt)*id(u) - trans(SigmaN)*idt(u) + 100*trans(idt(u))*id(u)/hFace());

    std::string inlet( "inlet" ), outlet( "outlet" );
    auto alpha = 1./sqrt(1-Nz()*Nz());
    auto C = alpha*mat<3,2>( cst(0.), Ny(), cst(0.), -Nx(), cst(1.), cst(0.) );
    auto lagt=vec(idt(lambda1),idt(lambda2));
    auto lag=vec(id(lambda1),id(lambda2));
    auto Clag1 = alpha*vec( cst(0.), cst(0.), id(lambda1) );
    auto Clag2 = alpha*vec( id(lambda2)*Ny(), -id(lambda2)*Nx(), cst(0.) );
    auto Clag1t = alpha*vec( cst(0.), cst(0.), idt(lambda1) );
    auto Clag2t = alpha*vec( idt(lambda2)*Ny(), -idt(lambda2)*Nx(), cst(0.) );

    //stokes +=integrate( markedfaces( mesh,inlet ), -trans(id(v))*(C*lagt));
    std::cout << "Vh21/Vh22 trial\n";
    for( auto bdy : { inlet, outlet } )
    {
        a_01 +=integrate( markedfaces( mesh, bdy ),
                          -trans(cross(id(u),N()))*(Clag1t) );


        a_02 +=integrate( markedfaces( mesh, bdy ),
                          -trans(cross(id(u),N()))*(Clag2t) );
    }


    auto Clag1TT = vec( cst(0.), cst(0.), id(lambda1) );
    auto Clag2inlet = vec( cst(0.), id(lambda2), cst(0.) );
    auto Clag2outlet = vec( cst(0.), -id(lambda2), cst(0.) );


    std::cout << "Vh21 test " << inlet << "\n";



    a_10   +=integrate( markedelements( submesh, inlet  ),
                     -trans(cross(idt(u),vec(cst(-1.0),cst(0.),cst(0.))))*Clag1TT);

    std::cout << "Vh21 test " << outlet << "\n";

    a_10 +=integrate( markedelements( submesh, outlet ),
                      -trans(cross(idt(u),vec(cst(1.0),cst(0.),cst(0.))))*Clag1TT);

    std::cout << "Vh22 test " << inlet << "\n";

    a_20 +=integrate( markedelements( submesh, inlet ),
                      -trans(cross(idt(u),vec(cst(-1.0),cst(0.),cst(0.))))*(Clag2inlet));

    std::cout << "Vh22 test " << outlet << "\n";
    form2( _test=Vh22, _trial=Vh1 ,_matrix=A,
               _rowstart=Vh1->nLocalDofWithGhost()+Vh21->nLocalDofWithGhost(), _colstart=0 )
        +=integrate( markedelements( submesh, outlet ),
                     -trans(cross(idt(u),vec(cst(1.0),cst(0.),cst(0.))))*(Clag2outlet));


    std::cout << "diag Vh1\n";
    form2( _test=Vh21, _trial=Vh21 ,_matrix=A,
           _rowstart=Vh1->nLocalDofWithGhost(), _colstart=Vh1->nLocalDofWithGhost() )
        +=integrate( elements( submesh ), doption("eps-lag")*idt(lambda1)*id(lambda1) );

    std::cout << "diag Vh2\n";
    form2( _test=Vh22, _trial=Vh22 ,_matrix=A,
           _rowstart=Vh1->nLocalDofWithGhost()+Vh21->nLocalDofWithGhost(),
           _colstart=Vh1->nLocalDofWithGhost()+Vh21->nLocalDofWithGhost() )
        +=integrate( elements( submesh ), doption("eps-lag")*idt(lambda2)*id(lambda2));


    form1( _test=Vh1, _vector=F,_rowstart=0 )
        += integrate(
            _range=markedfaces(mesh,"inlet"),
            _expr= -10*trans(N())*id(u) );
    form1( _test=Vh1, _vector=F,_rowstart=0 )
        += integrate(
            _range=markedfaces(mesh,"outlet"),
            _expr= -trans(N())*id(u) );

    a_11        += on( _range=boundaryfaces(submesh), _rhs=F, _element=*lambda1, _expr=cst(0.));

    a_22
        += on( _range=boundaryfaces(submesh), _rhs=F, _element=*lambda2, _expr=cst(0.));

    A->printMatlab("A.m");

    //stokes += on ( _range=markedfaces(mesh,"wall"), _rhs=F, _element=u, _expr=cst(0.));
    F->printMatlab("F.m");
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
                     _about=about(_name="stokes_pressurelm",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));


    //runStokesDirichletLM<2,1>();
    runStokesDirichletLM<3,1>();

    return 0;
}

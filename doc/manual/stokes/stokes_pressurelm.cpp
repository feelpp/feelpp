// -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

#include <feel/feel.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feelvf/print.hpp>
#include <feel/feeldiscr/pch.hpp>

namespace Feel
{
template <uint16_type Dim,uint16_type OrderGeo>
void runStokesDirichletLM()
{
    using namespace Feel;

    std::string configstr = (boost::format("%1%dGeo%2%")%Dim %OrderGeo).str();

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<3>> );
    std::list<std::string> listMarker{"inlet","wall"};
    std::list<std::string> presslm{"inlet","outlet"};
    auto submesh = createSubmesh(_mesh=mesh,_range=markedfaces(mesh,presslm));

    auto Vh1 = THch<OrderGeo>(mesh);
    auto Vh21 = Pch<2,double, PointSetEquiSpaced,Mesh<Simplex<2,1,3>>,0>(submesh);
    auto Vh22 = Pch<2,double, PointSetEquiSpaced,Mesh<Simplex<2,1,3>>,1>(submesh);

    cout << "mesh->numGlobalElements() "<<mesh->numGlobalElements() << std::endl;
    cout << "submesh->numGlobalElements() "<<submesh->numGlobalElements() << std::endl;
    cout << "Vh1->nDof() "<<Vh1->nDof() << std::endl;
    cout << "Vh21->nDof() "<<Vh21->nDof() << std::endl;
    cout << "Vh22->nDof() "<<Vh22->nDof() << std::endl;

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
                       _rowstart=0, _colstart=2 );
    auto a_10 = form2( _test=Vh21, _trial=Vh1 ,_matrix=A,
                       _rowstart=2, _colstart=0 );
    auto a_11=form2( _test=Vh21, _trial=Vh21 ,_matrix=A,
                     _rowstart=2,
                     _colstart=2 );
    auto a_02 = form2( _trial=Vh22, _test=Vh1 ,_matrix=A,
                       _rowstart=0, _colstart=3 );
    auto a_20 = form2( _test=Vh22, _trial=Vh1 ,_matrix=A,
                       _rowstart=3, _colstart=0 );
    auto a_22=form2( _test=Vh22, _trial=Vh22 ,_matrix=A,
                     _rowstart=3,
                     _colstart=3 );


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

    cout << ".. assembler" << std::endl;
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
    cout << ".. block(0,1) block(1,0) done";
    for( auto bdy : { inlet, outlet } )
    {
        a_01 +=integrate( _range=markedfaces( mesh, bdy ),
                          _expr=-trans(cross(id(u),N()))(0,2)*idt(lambda1)*alpha);
                          //-trans(cross(id(u),N()))*(Clag1t) );
        a_10 +=integrate( _range=markedfaces( mesh, bdy ),
                          _expr=-trans(cross(idt(u),N()))(0,2)*id(lambda1)*alpha);
        a_10 += on( _range=boundaryfaces(submesh), _rhs=F, _element=*lambda1, _expr=cst(0.));
            //-trans(cross(idt(u),N()))*(Clag1) );

    }
    cout << ".. block(0,2) block(2,0) done";
    for( auto bdy : { inlet, outlet } )
    {
        a_02 +=integrate( _range=markedfaces( mesh, bdy ),
                          _expr=
                          -trans(cross(id(u),N()))(0,0)*alpha*idt(lambda2)*Ny()
                          +trans(cross(id(u),N()))(0,1)*alpha*idt(lambda2)*Nx());
        
                          //-trans(cross(id(u),N()))*(Clag2t) );
        a_20 +=integrate( _range=markedfaces( mesh, bdy ),
                          _expr=
                          -trans(cross(idt(u),N()))(0,0)*alpha*id(lambda2)*Ny()
                          +trans(cross(idt(u),N()))(0,1)*alpha*id(lambda2)*Nx());
        a_20 += on( _range=boundaryfaces(submesh), _rhs=F, _element=*lambda2, _expr=cst(0.));
            //-trans(cross(idt(u),N()))*(Clag2) );
    }
    cout << ".. forcing terms on inlet" << std::endl;
    form1( _test=Vh1, _vector=F,_rowstart=0 )
        += integrate(
            _range=markedfaces(mesh,"inlet"),
            _expr= -10*trans(N())*id(u) );
    cout << ".. forcing terms on outlet" << std::endl;
    form1( _test=Vh1, _vector=F,_rowstart=0 )
        += integrate(
            _range=markedfaces(mesh,"outlet"),
            _expr= -trans(N())*id(u) );
    cout << ".. boundary conditions for lagrange multiplier 1" << std::endl;
    a_11        += on( _range=boundaryfaces(submesh), _rhs=F, _element=*lambda1, _expr=cst(0.));

    cout << ".. boundary conditions for lagrange multiplier 2" << std::endl;
    a_22
        += on( _range=boundaryfaces(submesh), _rhs=F, _element=*lambda2, _expr=cst(0.));

    cout << ".. solver" << std::endl;
    backend(_rebuild=true)->solve( _matrix=A, _rhs=F, _solution=UVec );
    myblockVecSol.localize(UVec);

    auto e1 = exporter( _mesh=mesh,_name="export"+configstr );
    if ( (soption( "functions.u" ) != "0") && (soption( "functions.p" ) != "0") )
    {
        auto Ue = Vh1->element();
        auto ue = Ue.template element<0>();
        auto pe = Ue.template element<1>();
        ue.on( _range=elements(mesh), _expr=expr<3,1>(soption( "functions.u" )));
        pe.on( _range=elements(mesh), _expr=expr(soption( "functions.p" )));
        e1->add("ue"+configstr, ue );
        e1->add("pe"+configstr, pe );
        cout << "||u-ue||_L2 = " << normL2( _range=elements(mesh), _expr=idv(u)-idv(ue)) << std::endl;
        cout << "||p-pe||_L2 = " << normL2( _range=elements(mesh), _expr=idv(p)-idv(pe)) << std::endl;
    }
    if (OrderGeo==1)
    {
        cout << ".. exporter" << std::endl;
        e1->add( "u"+configstr, u );
        e1->add( "p"+configstr, p );
        e1->save();

        auto es = exporter( _mesh=submesh,_name="subexport"+configstr );
        es->add( "lambda1"+configstr, *lambda1 );
        es->add( "lambda2"+configstr, *lambda2 );
        es->save();
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

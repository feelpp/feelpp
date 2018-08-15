/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feel.hpp>

namespace Feel
{
template <typename SpaceType>
NullSpace<double> qsNullSpace( SpaceType const& space, mpl::int_<2> /**/ )
{
    auto mode1 = space->element( oneX() );
    auto mode2 = space->element( oneY() );
    auto mode3 = space->element( vec(Py(),-Px()) );
    NullSpace<double> userNullSpace( { mode1,mode2,mode3 } );
    return userNullSpace;
}
template <typename SpaceType>
NullSpace<double> qsNullSpace( SpaceType const& space, mpl::int_<3> /**/ )
{
    auto mode1 = space->element( oneX() );
    auto mode2 = space->element( oneY() );
    auto mode3 = space->element( oneZ() );
    auto mode4 = space->element( vec(Py(),-Px(),cst(0.)) );
    auto mode5 = space->element( vec(-Pz(),cst(0.),Px()) );
    auto mode6 = space->element( vec(cst(0.),Pz(),-Py()) );
    NullSpace<double> userNullSpace( { mode1,mode2,mode3,mode4,mode5,mode6 } );
    return userNullSpace;
}

}

int main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description laplacianoptions( "Elasticity options" );
	laplacianoptions.add_options()
        ( "E", po::value<double>()->default_value( 1.0e6 ), "Young modulus" )
        ( "nu", po::value<double>()->default_value( 0.3 ), "Poisson ratio" )
        ( "no-solve", po::value<bool>()->default_value( false ), "No solve" )
        ( "weakdir", po::value<bool>()->default_value( false ), "use weak dirichlet" )
        ( "gamma", po::value<double>()->default_value( 100 ), "penalisation term" )
        ( "nullspace", po::value<bool>()->default_value( false ), "add null space" )
		;

	Environment env( _argc=argc, _argv=argv,
                   _desc=laplacianoptions,
                   _about=about(_name="qs_elasticity",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));

    tic();
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,1>>);
    toc("loadMesh");

    tic();
    auto Vh = Pchv<1>( mesh );
    toc("Vh");

    auto u = Vh->element("u");
    auto v = Vh->element("v");
    auto nu = doption(_name="nu");
    auto E = doption(_name="E");
    auto lambda = E*nu/( (1+nu)*(1-2*nu) );
    auto mu = E/(2*(1+nu));
    auto deft = sym(gradt(u));
    auto def = sym(grad(u));
    auto Id = eye<FEELPP_DIM,FEELPP_DIM>();
    auto sigmat = lambda*trace(deft)*Id + 2*mu*deft;
    auto sigma = lambda*trace(def)*Id + 2*mu*def;
    auto f = expr<FEELPP_DIM,1>( soption(_name="functions.f"), "f" );
    auto g = expr<FEELPP_DIM,1>( soption(_name="functions.g"), "g" );

    tic();
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=inner(f,id(v)));
    toc("l");

    tic();
    auto a = form2( _trial=Vh, _test=Vh);
    a = integrate(_range=elements(mesh),
                  _expr=inner( sigmat, grad(v) ) );

    if ( boption(_name="weakdir") )
    {
        double penaldir = doption(_name="gamma");
        a += integrate(_range=markedfaces(mesh,"Dirichlet"),
                       _expr=-inner(sigmat*N(),id(u)) + inner(-sigma*N()+std::max(2*mu,lambda)*penaldir*id(u)/hFace(),idt(u)) );

    }
    else
    {
        a+=on(_range=markedfaces(mesh,"Dirichlet"), _rhs=l, _element=u, _expr=g );
    }
    toc("a");

    //! solve the linear system, find u s.t. a(u,v)=l(v) for all v
    if ( !boption( "no-solve" ) )
    {
        tic();
        std::shared_ptr<NullSpace<double> > myNullSpace( new NullSpace<double>(backend(),qsNullSpace(Vh,mpl::int_<FEELPP_DIM>())) );
        backend()->attachNearNullSpace( myNullSpace );
        if ( boption(_name="nullspace") )
            backend()->attachNearNullSpace( myNullSpace );

        a.solve(_rhs=l,_solution=u);
        toc("a.solve");
    }

    tic();
    auto e = exporter( _mesh=mesh );
    e->addRegions();
    e->add( "u", u );
    e->save();
    toc("Exporter");
    return 0;
}

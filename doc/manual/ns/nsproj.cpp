// -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
#include <feel/feel.hpp>

using namespace Feel;


inline
po::options_description
makeOptions()
{
    po::options_description NSProjoptions( "NSproj options" );
    NSProjoptions.add_options()
    ( "dt", po::value<double>()->default_value( 0.01 ), "time step" )
    ( "mu", po::value<double>()->default_value( 1 ), "viscosity" )
    ( "dp", po::value<double>()->default_value( 1 ), "pressure difference" )
    ( "Niter", po::value<int>()->default_value( 1 ), "time iterations number" )
    ( "geo.file", po::value<std::string>()->default_value( "tube.geo" ), "geo file" )
    ;
    return NSProjoptions.add( feel_options().add(backend_options("vitesse").add(backend_options("pression"))) );
}

int main(int argc, char**argv )
{
	using namespace Feel;

    typedef Mesh<Simplex<2> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

	Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="NSProj",
                                  _author="Mourad Ismail",
                                  _email="mourad.ismail@ujf-grenoble.fr"));
    
    backend_ptrtype backend_V( backend_type::build( env.vm() , "vitesse") );
    backend_ptrtype backend_P( backend_type::build( env.vm() , "pression") );

    auto mu = env.vm(_name="mu").as<double>() ;
    auto dt = env.vm(_name="dt").as<double>() ;
    auto Niter = env.vm(_name="Niter").as<int>() ;
    auto dp = env.vm(_name="dp").as<double>() ;

    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                _desc=geo(
                                          _filename=env.vm(_name="geo.file").as<std::string>(),
                                          _h=env.vm(_name="mesh2d.hsize").as<double>() 
                                          )
                                );
    auto Vh = Pchv<2>( mesh );

    auto UTn  = Vh->element( "(u1,u2)" );
    auto UTn1 = Vh->element( "(u1,u2)" );
    auto Un1  = Vh->element( "(u1,u2)" );
    auto V = Vh->element( "(v1,v2)" );

    auto Dvit = backend_V->newMatrix( _test=Vh, _trial=Vh  );
    auto  aVit = form2( _trial=Vh, _test=Vh, _matrix=Dvit );

    auto Ph = Pch<1>( mesh );
    auto Dpre = backend_P->newMatrix( _test=Ph, _trial=Ph  );    
    auto aPre = form2( _trial=Ph, _test=Ph, _matrix=Dpre );

    auto pn1 = Ph->element( "p" );
    auto pn  = Ph->element( "p" );
    auto pnm1  = Ph->element( "p" );
    auto q   = Ph->element( "q" );

    //    auto poiseuille = vec( 4*0.3*Py()*(0.41-Py())/(0.41*0.41),cst(0.) );
    auto poiseuille = vec( 4*Py()*(0.41-Py()),cst(0.) );

    auto fn1 = vec( cst(0.),cst(0.) );
    auto Fvit = backend_V->newVector( Vh );
    auto lVit = form1( _test=Vh, _vector=Fvit  );


    auto Fpr = backend_P->newVector( Ph );
    auto lPre = form1( _test=Ph, _vector=Fpr );

    export_ptrtype exporter( export_type::New() );

    for(int i=0; i<Niter; i++)
        {
            /*
            auto t = i*dt;
            auto uex = vec( sin(Px())*sin(Py()+t),cos(Px())*cos(Py()+t) );
            auto pex = cos(Px())*sin(Py()+t);
            auto fex = vec( sin(Px())*cos(Py()+t) + sin(Px())*sin(Py()+t),
                            -cos(Px())*sin(Py()+t) + 3*cos(Px())*cos(Py()+t)
                            );
            auto fn1 = fex;

            UTn = vf::project(Vh,elements(mesh),uex );
            pn = vf::project(Ph,elements(mesh),pex );
            pnm1 = pn;
            */
            lVit = integrate(_range=elements(mesh),
                             _expr=
                             dt*inner(id(V),fn1) 
                             +inner(id(V),idv(UTn)) 
                             -2*dt* inner(trans(gradv(pn)),id(V)) 
                             + dt* inner(trans(gradv(pnm1)),id(V)) 
                             );

            aVit = integrate(_range=elements(mesh),
                             _expr=
                             inner(idt(UTn1),id(V)) 
                             + dt*mu*inner(gradt(UTn1),grad(V))
                               + dt*inner((gradt(UTn1)*idv(UTn)),id(V)) // ( a verifier UgradU )
                             );

            /*
            aVit+=on(_range=boundaryfaces(mesh), _rhs=lVit, _element=UTn1,
                     _expr=uex );
            */
            aVit+=on(_range=markedfaces(mesh,"wall"), _rhs=lVit, _element=UTn1,
                   _expr=vec(cst(0.),cst(0.)) );
            aVit+=on(_range=markedfaces(mesh,"wall"), _rhs=lVit, _element=UTn1,
                   _expr=vec(cst(0.),cst(0.)) );
            aVit+=on(_range=markedfaces(mesh,"cylinder"), _rhs=lVit, _element=UTn1,
                   _expr=vec(cst(0.),cst(0.)) );

            aVit+=on(_range=markedfaces(mesh,"inlet"), _rhs=lVit, _element=UTn1,
                     _expr=poiseuille );
            //            aVit+=on(_range=markedfaces(mesh,"outlet"), _rhs=lVit, _element=UTn1,
            //                   _expr=vec(cst(0.),cst(0.)) );

            backend_V->solve( _matrix=Dvit, _solution=UTn1, _rhs=Fvit );
            

            std::cout << " End velocity resolution " << std::endl;
            std::cout << " Begin pressure resolution " << std::endl;
            // Pressure

            lPre = integrate(_range=elements(mesh),
                             _expr=
                             -(1./dt)*id(q)*divv(UTn1)
                             + inner(gradv(pn),grad(q))
                             );
            
            aPre = integrate(_range=elements(mesh),
                             _expr=inner(gradt(pn1),grad(q))
                             );
            //            aPre+=on(_range=markedfaces(mesh,"inlet"), _rhs=lPre, _element=pn1,
            //                     _expr=cst(dp) );
            aPre+=on(_range=markedfaces(mesh,"outlet"), _rhs=lPre, _element=pn1,
                     _expr=cst(0.) );
            /*
            aPre+=on(_range=boundaryfaces(mesh), _rhs=lPre, _element=pn1,
                     _expr=pex );
            */
            
            backend_P->solve( _matrix=Dpre, _solution=pn1, _rhs=Fpr );
            
            auto gradPn1Proj = vf::project(Vh,elements(mesh), trans(gradv(pn1)));
            auto gradPnProj = vf::project(Vh,elements(mesh), trans(gradv(pn)));
            Un1 = UTn1 - dt*gradPn1Proj + dt*gradPnProj;            

            auto divUn = vf::project(Ph,elements(mesh), divv(Un1) );
            auto divUTn = vf::project(Ph,elements(mesh), divv(UTn1) );
            /*
            auto uexproj = vf::project(Vh,elements(mesh),uex );
            auto pexproj = vf::project(Ph,elements(mesh),pex );

            auto divUex = vf::project(Ph,elements(mesh), divv(uexproj) );

            auto udiff = uexproj - Un1;
            auto pdiff = pexproj - pn;
            */
            exporter->step( i )->setMesh( mesh );
            exporter->step( i )->add( "p", pn1 );
            exporter->step( i )->add( "u", Un1 );
            //            exporter->step( i )->add( "pdiff", pdiff );
            //            exporter->step( i )->add( "udiff", udiff );
            exporter->step( i )->add( "uT", UTn1 );
            exporter->step( i )->add( "divu", divUn );
            exporter->step( i )->add( "divuT", divUTn );
            //            exporter->step( i )->add( "divuex", divUex );
            exporter->save();
            
            UTn = UTn1;
            pnm1 = pn;
            pn = pn1;
        }

}




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
        ;
    return NSProjoptions.add( feel_options().add(backend_options("velocity").add(backend_options("pressure"))) );
}

int main(int argc, char**argv )
{

    typedef Mesh<Simplex<2> > mesh_type;

	Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="nsproj",
                                  _author="Mourad Ismail",
                                  _email="mourad.ismail@ujf-grenoble.fr"));

    auto mu = option(_name="mu").as<double>() ;
    auto dt = option(_name="dt").as<double>() ;
    auto Niter = option(_name="Niter").as<int>() ;
    auto dp = option(_name="dp").as<double>() ;

    auto mesh = loadMesh( _mesh=new mesh_type );

    auto Vh = Pchv<2>( mesh );

    auto UTn  = Vh->element( "(u1,u2)" );
    auto UTn1 = Vh->element( "(u1,u2)" );
    auto Un1  = Vh->element( "(u1,u2)" );
    auto V = Vh->element( "(v1,v2)" );

    auto aVit = form2( _trial=Vh, _test=Vh );

    auto Ph = Pch<1>( mesh );
    auto aPre = form2( _trial=Ph, _test=Ph );

    auto pn1 = Ph->element( "p" );
    auto pn  = Ph->element( "p" );
    auto pnm1  = Ph->element( "p" );
    auto q   = Ph->element( "q" );

    auto poiseuille = vec( 4*0.3*Py()*(0.41-Py())/(0.41*0.41),cst(0.) );

    auto fn1 = vec( cst(0.),cst(0.) );

    auto lVit = form1( _test=Vh  );
    auto lPre = form1( _test=Ph );

    auto exp = exporter( _mesh=mesh );

    node_type aa(2);
    aa[0]=0.15;
    aa[1]=0.2;
    node_type bb(2);
    bb[0]=0.25;
    bb[1]=0.2;

    for(int i=0; i<Niter; i++)
    {
        LOG(INFO) << "Iteration " << i << "/" << Niter << " Time = " << i*dt << "s";

        LOG(INFO) << "Velocity...";
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
                         + dt*inner((gradt(UTn1)*idv(UTn)),id(V))
            );

        aVit+=on(_range=markedfaces(mesh,"wall"), _rhs=lVit, _element=UTn1,
                 _expr=vec(cst(0.),cst(0.)) );
        aVit+=on(_range=markedfaces(mesh,"cylinder"), _rhs=lVit, _element=UTn1,
                 _expr=vec(cst(0.),cst(0.)) );

        aVit+=on(_range=markedfaces(mesh,"inlet"), _rhs=lVit, _element=UTn1,
                 _expr=poiseuille );
        //            aVit+=on(_range=markedfaces(mesh,"outlet"), _rhs=lVit, _element=UTn1,
        //                   _expr=vec(cst(0.),cst(0.)) );
        aVit.solve( _solution=UTn1, _rhs=lVit, _name="velocity" );
        LOG(INFO) << "Velocity problem done.";
        // Pressure
        LOG(INFO) << "Pressure...";
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

        aPre.solve( _solution=pn1, _rhs=lPre, _name="pressure" );
        LOG(INFO) << "Pressure problem done.";
        LOG(INFO) << "Velocity correction...";
        auto gradPn1Proj = vf::project(_space=Vh,_range=elements(mesh), _expr=trans(gradv(pn1)));
        auto gradPnProj = vf::project(_space=Vh,_range=elements(mesh), _expr=trans(gradv(pn)));
        //Un1 = UTn1 - dt*gradPn1Proj + dt*gradPnProj;
        Un1 = vf::project(_space=Vh,_range=elements(mesh), _expr=idv(UTn1) + dt*trans(gradv(pn)-gradv(pn1) ) );

        auto divUn = vf::project(_space=Ph,_range=elements(mesh), _expr=divv(Un1) );
        LOG(INFO) << "Velocity correction done.";

        double time = i*dt;

        if ( exp->doExport() )
        {
            exp->step( time )->setMesh( mesh );
            exp->step( time )->addScalar( "viscosity", mu, true);
            exp->step( time )->add( "p", pn1 );
            exp->step( time )->add( "u", Un1 );
            exp->step( time )->add( "uT", UTn1 );
            exp->step( time )->add( "divu", divUn );
            exp->save();
        }


        UTn = UTn1;
        pnm1 = pn;
        pn = pn1;


        LOG(INFO) << "----> pn(aa)-pn(bb) =  " << pn1(aa)(0,0,0)-pn1(bb)(0,0,0);
        LOG(INFO) << "----> Un(aa) =  " << Un1(aa)(0,0,0) << " , " << Un1(aa)(1,0,0);
        LOG(INFO) << "----> Un(bb) =  " << Un1(bb)(0,0,0) << " , " << Un1(bb)(1,0,0);
    }

}

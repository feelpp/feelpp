#ifndef PENALISATION_IMPL
#define PENALISATION_IMPL

#include "penalisation.hpp"

using namespace Feel;

template <int Dim>
Penalisation<Dim>::Penalisation():
    Application(),
    H1( 2. ),
    H2( 2. ),
    L1( 4. ),
    E( 2. ),
/*TODO : in the mesh, call the boundary H1, H2 ...
  and calculate here the lenghts and do not enter it in the code !
*/
    radius(       doption("radius"      ) ),
    epsilon(      doption("epsilon"     ) ),
    epsilonpress( doption("epsilonpress") ),
    nu(           doption("nu"          ) ),
    Q(            doption("Q"           ) ),
    Qtop(         doption("RQ"          )*Q ),
    Qbot( Q-Qtop ),
    dt(     doption("DT"    ) ),
    Ylimit( doption("Ylimit") )
{
    LOG(INFO)<<"Dimension : "<<Dim<<"\n";

    if ( this->vm().count( "help" ) )
    {
        LOG(INFO) << this->optionsDescription() << "\n";
        exit( 0 );
    }

    std::string OutFolder=soption("OutFolder");

    if ( !OutFolder.empty() )
        this->changeRepository( boost::format( OutFolder ) );

    mesh = loadMesh( _mesh = new mesh_type );
    mesh_visu = mesh;

    mesh_visu = mesh;

    exporter =  Exporter<mesh_type>::New();

    Xh=space_type::New( mesh );
    Yh=space_carac_type::New( mesh );
    P0h=space_p0_type::New( mesh );
    LOG(INFO) << "mesh.elements:" << mesh->numElements() << "\n";
    LOG(INFO) << "P0h.nLocalDof:" << P0h->nLocalDof() << "\n";
    LOG(INFO) << "Xi.nLocalDof:" << Yh->nLocalDof() << "\n";

    t=0;
    iter=0;
    Tfinal=doption("Tfinal");
    // if we are in test mode then do only one (or a small multiple of) time step
    if ( Environment::vm().count( "test" ) )
        Tfinal =  Environment::vm( _name="test" ).template as<int>()*dt;
    xp=doption("x0");
    yp=doption("y0");
    zp=doption("z0");

    M_backend= backend(_name="stokes_backend" );

    U=Xh->elementPtr();
    LOG(INFO)<<"U.size = "<<U->size()<<"\n";
    LOG(INFO)<<"U.constenair.size = "<<( U->container() ).size()<<"\n";

    LOG(INFO)<<"out of constructor"<<"\n";
}//Penalisation


template <int Dim>
void Penalisation<Dim>::initStokes()
{
    boost::timer local_chrono;
    auto V = Xh->element();

    element_veloc_type u = U->template element<0>();
    element_veloc_type v = V.template element<0>();

    element_pressure_type p= U->template element<1>();
    element_pressure_type q= V.template element<1>();

    element_lag_type lambda= U->template element<2>();

    auto deft = sym( gradt( u ) );
    auto def  = sym( grad( v ) );

    LOG(INFO)<<"check mesh\n";
    LOG(INFO)<<"H1 = "<<integrate( markedfaces( mesh, "Inflow" ), vf::cst( 1. ) ).evaluate()( 0,0 )<<"\n";
    LOG(INFO)<<"Wall = "<<integrate( markedfaces( mesh, "Wall" ), vf::cst( 1. ) ).evaluate()( 0,0 )<<"\n";
    LOG(INFO)<<"H2 = "<<integrate( markedfaces( mesh, "Outflowtop" ), vf::cst( 1. ) ).evaluate()( 0,0 )<<"\n";

    D=M_backend->newMatrix( _test=Xh, _trial=Xh );
    auto a = form2( _test=Xh, _trial=Xh, _matrix=D );

    local_chrono.restart();

    //try not use F
    F = M_backend->newVector( Xh );

    form1( Xh, _vector=F, _init=true ) =
        integrate( elements( mesh ), trans( vf::one()-vf::one() ) * id( v ) ) ;

    F->close();
    // F->printMatlab("F.m");
    LOG(INFO)<<"F closed"<<"\n";
    LOG(INFO)<<"assemblage F  : "<<local_chrono.elapsed()<<" s"<<"\n";

    local_chrono.restart();

    C = M_backend->newMatrix( Xh, Xh );

    //    def(u):def(v) = trace ( def(u)* trans( def(v) ) )
    form2( Xh, Xh, _matrix=C, _init=true )= integrate( elements( mesh ),
                                    nu * vf::trace( deft * trans( def ) ) );
    form2( Xh, Xh, _matrix=C )+= integrate( elements( mesh ),
                                    - div( v ) * idt( p ) +id( q ) * divt( u ) );
    form2( Xh, Xh, _matrix=C )+= integrate( elements( mesh ),
                                    - epsilonpress * idt( p ) * id( q ) );
    form2( Xh, Xh, _matrix=C )+= integrate( elements( mesh ),idt( lambda )*id( q ) );
    form2( Xh, Xh, _matrix=C )+= integrate( elements( mesh ),id( lambda )*idt( q ) );
    C->close();

    LOG(INFO)<<"Fin de l'assemblage statique (C) temps : "<<
             local_chrono.elapsed()<<" s"<<"\n";

    local_chrono.restart();

} //initStokes


template <int Dim>
void Penalisation<Dim>::stokes()
{
    boost::timer local_chrono;
    auto u = U->template element<0>();

    auto V = Xh->element();
    auto v = V.template element<0>();

    auto deft = sym( gradt( u ) );
    auto def  = sym( grad( v ) );

    local_chrono.restart();
    D=M_backend->newMatrix( Xh, Xh );
    form2( Xh, Xh, _matrix=D )= integrate( marked3elements( mesh,1 ),
                                   idv( carac ) * nu * trace( def*trans( deft ) ) / epsilon );
    LOG(INFO)<<"fin assemblage : "<<local_chrono.elapsed()<<" s"<<"\n";

    local_chrono.restart();
    D->addMatrix( 1.0, C );
    LOG(INFO)<<"fin copie ajout D : "<<local_chrono.elapsed()<<" s"<<"\n";

    D->close();

    addCL();

    Feel::backend(_name="stokes_backend",_rebuild=true)->solve( _matrix=D,_solution=U,_rhs=F );

    LOG(INFO)<<"fin resolution : "<<local_chrono.elapsed()<<" s"<<"\n";

    double Q0, Q1, Q2;
    Q0 = integrate( markedfaces( mesh,mesh->markerName( "Inflow" ) ),trans( idv( u ) )*vf::N() ).evaluate()( 0,0 );
    Q1 = integrate( markedfaces( mesh,mesh->markerName( "Outflowtop" ) ),trans( idv( u ) )*vf::N() ).evaluate()( 0,0 );
    Q2 = integrate( markedfaces( mesh,mesh->markerName( "Outflowbottom" ) ),trans( idv( u ) )*vf::N() ).evaluate()( 0,0 );

    LOG(INFO)<<"Q0 = "<<Q0<<"\n";
    LOG(INFO)<<"Q1 = "<<Q1<<"\n";
    LOG(INFO)<<"Q2 = "<<Q2<<"\n";
    LOG(INFO)<<"div(u)="<<integrate( elements( mesh ),
                                     divv( u ) ).evaluate()( 0,0 )<<"\n";

}//stokes


template <int Dim>
void Penalisation<Dim>::addCL()
{
    LOG(INFO) << "adding boundary conditions...\n";
    auto u = U->template element<0>();

    if ( Dim==2 )
    {
        auto inflow = ( 6*Q /( pow( H1,3 ) ) ) * ( Py()+H1/2 )*( Py()-H1/2 );
        auto outflowtop = ( 6*Qtop / ( pow( H2,3 ) ) ) * ( Px()-L1 )*( Px()-( L1+H2 ) );
        auto outflowbottom = ( 6*Qbot /( pow( H2,3 ) ) ) * ( Px()-L1 )*( Px()-( L1+H2 ) );

        form2( Xh, Xh, _matrix=D )+=
            on( markedfaces( mesh, "Wall" ),
                u, F, cst( 0. ) * vf::N()  );
        form2( Xh, Xh, _matrix=D )+=
            on( markedfaces( mesh, "Inflow" ),
                u, F, inflow * vf::N() );
        form2( Xh, Xh, _matrix=D )+=
            on( markedfaces( mesh, "Outflowtop" ),
                u, F, - outflowtop * vf::N() );
        form2( Xh, Xh, _matrix=D )+=
            on( markedfaces( mesh, "Outflowbottom" ),
                u, F,  - outflowbottom * vf::N() );
    }

    else
    {

        auto inflow = ( 6*Q /( pow( H1,3 ) ) ) * ( Py()+H1/2 )*( Py()-H1/2 ) * Pz() * ( Pz()-E )  ;
        auto outflowtop = ( 6*Qtop / ( pow( H2,3 ) ) ) * ( Px()-L1 )*( Px()-( L1+H2 ) ) * Pz() * ( Pz()-E ) ;
        auto outflowbottom = ( 6*Qbot /( pow( H2,3 ) ) ) * ( Px()-L1 )*( Px()-( L1+H2 ) ) * Pz() * ( Pz()-E ) ;

        form2( Xh, Xh, _matrix=D )+=
            on( markedfaces( mesh, "Wall" ),
                u, F, cst( 0. ) * vf::N()  );
        form2( Xh, Xh, _matrix=D )+=
            on( markedfaces( mesh, "Inflow" ),
                u, F, inflow * vf::N() );
        form2( Xh, Xh, _matrix=D )+=
            on( markedfaces( mesh, "Outflowtop" ),
                u, F, - outflowtop * vf::N() );
        form2( Xh, Xh, _matrix=D )+=
            on( markedfaces( mesh, "Outflowbottom" ),
                u, F,  - outflowbottom * vf::N() );
    }
    LOG(INFO) << "adding boundary conditions done.\n";
}


template <int Dim>
void Penalisation<Dim>::updateChi()
{
    LOG(INFO) << "update chi...\n";
    auto carac_expr = vf::chi( radius*radius >
                               ( vf::Px()-xp )*( vf::Px()-xp ) +
                               ( vf::Py()-yp )*( vf::Py()-yp ) +
                               ( vf::Pz()-zp )*( vf::Pz()-zp ) );
    carac=vf::project( Yh, elements( mesh ),carac_expr );
    LOG(INFO) << "caracteristic function built...\n";
    p0.zero();
    p0 = integrate( _range=elements( mesh ),_expr=idv( carac ), _quad=_Q<6>() ).broken( P0h );
    LOG(INFO) << "p0 built...\n";
    google::FlushLogFiles(google::INFO);
#if 0
    for ( auto it=p0.begin(),en=p0.end(); it != en; ++it )
    {
        if ( math::abs( *it ) > 1e-10 ) *it = 1;

        else *it = 0;
    }
#else
    p0 = vf::project( _space=P0h, _range=elements(mesh), _expr=chi(abs(idv(p0)) > 1e-10) );
#endif
    LOG(INFO) << "update marker3 built...\n";
    google::FlushLogFiles(google::INFO);
    mesh->updateMarker3( p0 );
    LOG(INFO) << "updated marker3...\n";
    google::FlushLogFiles(google::INFO);
    LOG(INFO) << "number of marked 3 elements(1): " << nelements( marked3elements( mesh, 1 ) ) << "\n";
    google::FlushLogFiles(google::INFO);
    LOG(INFO) << "number of marked 3 elements(0): " << nelements( marked3elements( mesh, 0 ) ) << "\n";
    google::FlushLogFiles(google::INFO);
    LOG(INFO) << "updateChi: p0.size=  " << p0.size() << "\n";
    google::FlushLogFiles(google::INFO);
    LOG(INFO) << "updateChi: carac.size=  " << carac.size() << "\n";
    google::FlushLogFiles(google::INFO);

    double aire = integrate( marked3elements( mesh, 1 ), idv( carac ) ).evaluate()( 0,0 );
    LOG(INFO)<<"updateChi: aire = "<<aire << "  (exact: " << M_PI*radius*radius << ")" <<"\n";

    google::FlushLogFiles(google::INFO);
#if 0
    p0.printMatlab( "p0" );
    carac.printMatlab( "ca" );
#endif
}//updateChi


template <int Dim>
void Penalisation<Dim>::updatePosition()
{
    element_veloc_type u = U->template element<0>();

    double aire = integrate( marked3elements( mesh,1 ),
                             idv( carac ) ).evaluate()( 0,0 );
    LOG(INFO)<<"aire = "<<aire<<"\n";

    auto Velo = integrate( marked3elements( mesh,1 ),
                           idv( carac )*idv( u ) ).evaluate() / aire ;

    //TO DO : try to go to order 2
    Vx = Velo( 0,0 );
    Vy = Velo( 1,0 );

    if ( Dim == 3 )
        Vz = Velo( 2,0 );

    else
        Vz = 0;

    LOG(INFO)<<"Vx, Vy, Vz \n"<<Vx<<" "<<Vy<<" "<<Vz<<"\n";

    //Euler
    xp +=  dt * Vx;
    yp +=  dt * Vy;
    zp +=  dt * Vz;

    LOG(INFO)<<"x y z \n"<<xp<<" "<<yp<<" "<<zp<<"\n";

}//updatePosition


template <int Dim>
void Penalisation<Dim>::exportResults()
{
    exporter->step( t )->setMesh( mesh_visu );
    exporter->step( t )->add( "elements", p0 );
    exporter->step( t )->add( "chi", carac );
    exporter->step( t )->add( "u", U->template element<0>() );
    exporter->step( t )->add( "p", U->template element<1>() );
    exporter->save();
}//exportResults


template <int Dim>
void Penalisation<Dim>::run()
{
    carac=Yh->element();
    p0=P0h->element();
    std::ofstream TrajFile;
    TrajFile.open( "Trajectory" );
    TrajFile<<"# time  x  y  z  Vx  Vy  Vz"<<"\n";

    chrono.restart();
    this->initStokes();
    LOG(INFO)<<"Init_Stokes : "<<chrono.elapsed()<<" s"<<"\n";
    google::FlushLogFiles(google::INFO);
    updateChi();
    exportResults();

    for ( t=dt; t<Tfinal; t+=dt )
    {
        chrono.restart();
        // need to rebuild backend (it stores information about chi !)
        // backend= backend_type::build(this->vm(), "stokes_backend");
        // LOG(INFO)<<"reinit backend : "<<chrono.elapsed()<<" s"<<"\n";

        iter++;
        LOG(INFO)<<"===========================================================================================\n";
        LOG(INFO)<<"stokes, iter "<<iter<<"\n";
        LOG(INFO)<<"t= "<<t<<"\n";
        chrono.restart();
        stokes();
        LOG(INFO)<<"stokes : "<<chrono.elapsed()<<" s"<<"\n";
        updatePosition();
        updateChi();

        exportResults();
        TrajFile<<std::left<<std::setw( 15 )<<t
                <<std::left<<std::setw( 15 )<<xp
                <<std::left<<std::setw( 15 )<<yp
                <<std::left<<std::setw( 15 )<<zp
                <<std::left<<std::setw( 15 )<<Vx
                <<std::left<<std::setw( 15 )<<Vy
                <<std::left<<std::setw( 15 )<<Vz
                <<std::endl;

        if (  ( Ylimit!=0 ) && ( std::abs( yp )>Ylimit )  )
            t=Tfinal+1;
    }//for time

    TrajFile.close();
}//run



#endif // PENALISATION_IMPL

#ifndef PENALISATION
#define PENALISATION

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
    radius( this->vm()["radius"].template as<double>() ),
    epsilon( this->vm()["epsilon"].template as<double>() ),
    epsilonpress( this->vm()["epsilonpress"].template as<double>() ),
    nu( this->vm()["nu"].template as<double>() ),
    Q( this->vm()["Q"].template as<double>() ),
    Qtop( this->vm()["RQ"].template as<double>()*Q ),
    Qbot( Q-Qtop ),
    dt( this->vm()["DT"].template as<double>() ),
    Ylimit( this->vm()["Ylimit"].template as<double>() )
{
    std::cout<<"Dimension : "<<Dim<<std::endl;

    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        exit( 0 );
    }

    std::string OutFolder=this->vm()["OutFolder"].template as<std::string>();

    if ( !OutFolder.empty() )
        this->changeRepository( boost::format( OutFolder ) );

    std::string File_Mesh= this->vm()["ImportMeshFromFile"].template as<std::string>();
    std::cout<<"chargement maillage"<<std::endl;
#if 0
    mesh = loadGMSHMesh( _mesh=new mesh_type,
                         _filename=File_Mesh,
                         _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );
#else
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=geo( _filename=File_Mesh,_dim=2,_h=this->vm()["hsize"].template as<double>() ),
                           _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                           _partitions=this->comm().size() );

#if 0
    mesh_visu = createGMSHMesh( _mesh=new mesh_type,
                                _desc=geo( _filename=File_Mesh,_dim=2,_h=this->vm()["hsize"].template as<double>()/2 ),
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );
#else
    mesh_visu = mesh;
#endif
#endif
    exporter.reset( Exporter<mesh_type>::New( this->vm(),this->about().appName() ) );

    Xh=space_type::New( mesh );
    Yh=space_carac_type::New( mesh );
    P0h=space_p0_type::New( mesh );
    std::cout << "mesh.elements:" << mesh->numElements() << "\n";
    std::cout << "P0h.nLocalDof:" << P0h->nLocalDof() << "\n";
    std::cout << "Xi.nLocalDof:" << Yh->nLocalDof() << "\n";

    t=0;
    iter=0;
    Tfinal=this->vm()["Tfinal"].template as<double>();
    xp=this->vm()["x0"].template as<double>();
    yp=this->vm()["y0"].template as<double>();
    zp=this->vm()["z0"].template as<double>();

    backend= backend_type::build( this->vm(), "stokes_backend" );

    U=Xh->element();
    std::cout<<"U.size = "<<U.size()<<std::endl;
    std::cout<<"U.constenair.size = "<<( U.container() ).size()<<std::endl;

    std::cout<<"out of constructor"<<std::endl;
}//Penalisation


template <int Dim>
void Penalisation<Dim>::initStokes()
{
    boost::timer local_chrono;
    auto V = Xh->element();

    element_veloc_type u = U.template element<0>();
    element_veloc_type v = V.template element<0>();

    element_pressure_type p= U.template element<1>();
    element_pressure_type q= V.template element<1>();

    element_lag_type lambda= U.template element<2>();

    auto deft = sym( gradt( u ) );
    auto def  = sym( grad( v ) );

    std::cout<<"check mesh\n";
    std::cout<<"H1 = "<<integrate( markedfaces( mesh, "Inflow" ), vf::cst( 1. ) ).evaluate()( 0,0 )<<std::endl;
    std::cout<<"Wall = "<<integrate( markedfaces( mesh, "Wall" ), vf::cst( 1. ) ).evaluate()( 0,0 )<<std::endl;
    std::cout<<"H2 = "<<integrate( markedfaces( mesh, "Outflowtop" ), vf::cst( 1. ) ).evaluate()( 0,0 )<<std::endl;

    D=backend->newMatrix( Xh, Xh );
    form2( Xh, Xh, D, _init=true );

    local_chrono.restart();

    //try not use F
    F = backend->newVector( Xh );

    form1( Xh, F, _init=true ) =
        integrate( elements( mesh ), trans( vf::one()-vf::one() ) * id( v ) ) ;

    F->close();
    // F->printMatlab("F.m");
    std::cout<<"F closed"<<std::endl;
    std::cout<<"assemblage F  : "<<local_chrono.elapsed()<<" s"<<std::endl;

    local_chrono.restart();

    C = backend->newMatrix( Xh, Xh );

    //    def(u):def(v) = trace ( def(u)* trans( def(v) ) )
    form2( Xh, Xh, C, _init=true )= integrate( elements( mesh ),
                                    nu * vf::trace( deft * trans( def ) ) );
    form2( Xh, Xh, C )+= integrate( elements( mesh ),
                                    - div( v ) * idt( p ) +id( q ) * divt( u ) );
    form2( Xh, Xh, C )+= integrate( elements( mesh ),
                                    - epsilonpress * idt( p ) * id( q ) );
    form2( Xh, Xh, C )+= integrate( elements( mesh ),idt( lambda )*id( q ) );
    form2( Xh, Xh, C )+= integrate( elements( mesh ),id( lambda )*idt( q ) );
    C->close();

    std::cout<<"Fin de l'assemblage statique (C) temps : "<<
             local_chrono.elapsed()<<" s"<<std::endl;

    local_chrono.restart();

} //initStokes


template <int Dim>
void Penalisation<Dim>::stokes()
{
    boost::timer local_chrono;
    auto u = U.template element<0>();

    auto V = Xh->element();
    auto v = V.template element<0>();

    auto deft = sym( gradt( u ) );
    auto def  = sym( grad( v ) );

    local_chrono.restart();
    D=backend->newMatrix( Xh, Xh );
    form2( Xh, Xh, D )= integrate( marked3elements( mesh,1 ),
                                   idv( carac ) * nu * trace( def*trans( deft ) ) / epsilon );
    std::cout<<"fin assemblage : "<<local_chrono.elapsed()<<" s"<<std::endl;

    local_chrono.restart();
    D->addMatrix( 1.0, C );
    std::cout<<"fin copie ajout D : "<<local_chrono.elapsed()<<" s"<<std::endl;

    D->close();

    addCL();

    backend->solve( _matrix=D,_solution=U,_rhs=F );

    std::cout<<"fin resolution : "<<local_chrono.elapsed()<<" s"<<std::endl;

    double Q0, Q1, Q2;
    Q0 = integrate( markedfaces( mesh,mesh->markerName( "Inflow" ) ),trans( idv( u ) )*vf::N() ).evaluate()( 0,0 );
    Q1 = integrate( markedfaces( mesh,mesh->markerName( "Outflowtop" ) ),trans( idv( u ) )*vf::N() ).evaluate()( 0,0 );
    Q2 = integrate( markedfaces( mesh,mesh->markerName( "Outflowbottom" ) ),trans( idv( u ) )*vf::N() ).evaluate()( 0,0 );

    std::cout<<"Q0 = "<<Q0<<std::endl;
    std::cout<<"Q1 = "<<Q1<<std::endl;
    std::cout<<"Q2 = "<<Q2<<std::endl;
    std::cout<<"div(u)="<<integrate( elements( mesh ),
                                     divv( u ) ).evaluate()( 0,0 )<<std::endl;

}//stokes


template <int Dim>
void Penalisation<Dim>::addCL()
{
    auto u = U.template element<0>();

    if ( Dim==2 )
    {
        auto inflow = ( 6*Q /( pow( H1,3 ) ) ) * ( Py()+H1/2 )*( Py()-H1/2 );
        auto outflowtop = ( 6*Qtop / ( pow( H2,3 ) ) ) * ( Px()-L1 )*( Px()-( L1+H2 ) );
        auto outflowbottom = ( 6*Qbot /( pow( H2,3 ) ) ) * ( Px()-L1 )*( Px()-( L1+H2 ) );

        form2( Xh, Xh, D )+=
            on( markedfaces( mesh, "Wall" ),
                u, F, cst( 0. ) * vf::N()  );
        form2( Xh, Xh, D )+=
            on( markedfaces( mesh, "Inflow" ),
                u, F, inflow * vf::N() );
        form2( Xh, Xh, D )+=
            on( markedfaces( mesh, "Outflowtop" ),
                u, F, - outflowtop * vf::N() );
        form2( Xh, Xh, D )+=
            on( markedfaces( mesh, "Outflowbottom" ),
                u, F,  - outflowbottom * vf::N() );
    }

    else
    {

        auto inflow = ( 6*Q /( pow( H1,3 ) ) ) * ( Py()+H1/2 )*( Py()-H1/2 ) * Pz() * ( Pz()-E )  ;
        auto outflowtop = ( 6*Qtop / ( pow( H2,3 ) ) ) * ( Px()-L1 )*( Px()-( L1+H2 ) ) * Pz() * ( Pz()-E ) ;
        auto outflowbottom = ( 6*Qbot /( pow( H2,3 ) ) ) * ( Px()-L1 )*( Px()-( L1+H2 ) ) * Pz() * ( Pz()-E ) ;

        form2( Xh, Xh, D )+=
            on( markedfaces( mesh, "Wall" ),
                u, F, cst( 0. ) * vf::N()  );
        form2( Xh, Xh, D )+=
            on( markedfaces( mesh, "Inflow" ),
                u, F, inflow * vf::N() );
        form2( Xh, Xh, D )+=
            on( markedfaces( mesh, "Outflowtop" ),
                u, F, - outflowtop * vf::N() );
        form2( Xh, Xh, D )+=
            on( markedfaces( mesh, "Outflowbottom" ),
                u, F,  - outflowbottom * vf::N() );
    }

}


template <int Dim>
void Penalisation<Dim>::updateChi()
{
    auto carac_expr = vf::chi( radius*radius >
                               ( vf::Px()-xp )*( vf::Px()-xp ) +
                               ( vf::Py()-yp )*( vf::Py()-yp ) +
                               ( vf::Pz()-zp )*( vf::Pz()-zp ) );
    carac=vf::project( Yh, elements( mesh ),carac_expr );
    p0.zero();
    p0 = integrate( _range=elements( mesh ),_expr=idv( carac ), _quad=_Q<6>() ).broken( P0h );

    for ( auto it=p0.begin(),en=p0.end(); it != en; ++it )
    {
        if ( math::abs( *it ) > 1e-10 ) *it = 1;

        else *it = 0;
    }

    mesh->updateMarker3( p0 );
    double aire = integrate( marked3elements( mesh,1 ),
                             idv( carac ) ).evaluate()( 0,0 );
    std::cout<<"updateChi: aire = "<<aire << "  (exact: " << M_PI*radius*radius << ")" <<std::endl;

    std::cout << "number of marked 3 elements: " << std::distance( marked3elements( mesh, 1 ).get<1>(),
              marked3elements( mesh, 1 ).get<2>() ) << "\n";
    std::cout << "updateChi: p0.size=  " << p0.size() << "\n";
    std::cout << "updateChi: carac.size=  " << carac.size() << "\n";
    std::cout << "updateChi: mesh.size=  " << mesh->numElements() << "\n";
    p0.printMatlab( "p0" );
    carac.printMatlab( "ca" );
}//updateChi


template <int Dim>
void Penalisation<Dim>::updatePosition()
{
    element_veloc_type u = U.template element<0>();

    double aire = integrate( marked3elements( mesh,1 ),
                             idv( carac ) ).evaluate()( 0,0 );
    std::cout<<"aire = "<<aire<<std::endl;

    auto Velo = integrate( marked3elements( mesh,1 ),
                           idv( carac )*idv( u ) ).evaluate() / aire ;

    //TO DO : try to go to order 2
    Vx = Velo( 0,0 );
    Vy = Velo( 1,0 );

    if ( Dim == 3 )
        Vz = Velo( 2,0 );

    else
        Vz = 0;

    std::cout<<"Vx, Vy, Vz \n"<<Vx<<" "<<Vy<<" "<<Vz<<std::endl;

    //Euler
    xp +=  dt * Vx;
    yp +=  dt * Vy;
    zp +=  dt * Vz;

    std::cout<<"x y z \n"<<xp<<" "<<yp<<" "<<zp<<std::endl;

}//updatePosition


template <int Dim>
void Penalisation<Dim>::exportResults()
{
    exporter->step( t )->setMesh( mesh_visu );
    exporter->step( t )->add( "elements", p0 );
    exporter->step( t )->add( "chi", carac );
    exporter->step( t )->add( "u", U.template element<0>() );
    exporter->step( t )->add( "p", U.template element<1>() );
    exporter->save();
}//exportResults


template <int Dim>
void Penalisation<Dim>::run()
{
    carac=Yh->element();
    p0=P0h->element();
    std::ofstream TrajFile;
    TrajFile.open( "Trajectory" );
    TrajFile<<"# time  x  y  z  Vx  Vy  Vz"<<std::endl;

    chrono.restart();
    this->initStokes();
    std::cout<<"Init_Stokes : "<<chrono.elapsed()<<" s"<<std::endl;

    updateChi();
    exportResults();

    for ( t=dt; t<Tfinal; t+=dt )
    {
        chrono.restart();
        // need to rebuild backend (it stores information about chi !)
        // backend= backend_type::build(this->vm(), "stokes_backend");
        // std::cout<<"reinit backend : "<<chrono.elapsed()<<" s"<<std::endl;

        iter++;
        std::cout<<"===========================================================================================\n";
        std::cout<<"stokes, iter "<<iter<<std::endl;
        std::cout<<"t= "<<t<<std::endl;
        chrono.restart();
        stokes();
        std::cout<<"stokes : "<<chrono.elapsed()<<" s"<<std::endl;
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



#endif

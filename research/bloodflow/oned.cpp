/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-11-08

  Copyright (C) 2007 Universite Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file oned.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-11-08
 */
#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/geotool.hpp>

#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelpoly/polynomialset.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feelmesh/elements.hpp>

#include <feel/feelvf/vf.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description onedoptions( "OneD options" );
    onedoptions.add_options()
    ( "dt", Feel::po::value<double>()->default_value( 1e-5 ), "time step" )
    ( "Tf", Feel::po::value<double>()->default_value( 1 ), "final time" )
    ( "Kr", Feel::po::value<double>()->default_value( 1 ), "friction parameter" )
    ( "rho", Feel::po::value<double>()->default_value( 1 ), "blood density" )
    ( "E", Feel::po::value<double>()->default_value( 1e6 ), "artery yound modulus" )
    ( "h0", Feel::po::value<double>()->default_value( 0.05 ), "artery thickness" )
    ( "A0", Feel::po::value<double>()->default_value( 1.0 ), "artery at reference state (pressure=0)" )
    ( "alpha", Feel::po::value<double>()->default_value( 1 ), "momentum-flux correction (coriolis coefficient)" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "mesh size" )
    ( "export", "export results(ensight, data file(1D)" )
    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return onedoptions.add( Feel::feel_options() ) ;
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "oned" ,
                           "oned" ,
                           "0.1",
                           "1D euler code for blood flow simulation",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2007 University Joseph Fourier Grenoble 1" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Tiziano Passerini", "developer", "tiziano.passerini@polimit.it", "" );
    return about;

}


namespace Feel
{
/**
 * Diffussion Advection Reaction Solver
 *
 * solve \f$-\epsilon \Delta u -\beta\cdot\nabla u + \mu u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma_{in}\f$
 */
class OneD
    :
    //public Application
public Application
{
    //typedef Application super;
    typedef Application super;
public:

    // -- TYPEDEFS --
    static const uint16_type Dim = 1;
    static const uint16_type Order = 1;
    static const uint16_type imOrder = 2*Order;

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<Backend<value_type> > backend_ptrtype;

    /*matrix*/
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Simplex<1, 1, 1> entity_type;
    typedef Mesh<entity_type > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptr_type;

    typedef Lagrange< 1, Scalar, Continuous> basis_type;

    typedef FunctionSpace<mesh_type, bases<basis_type>, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef space_type::element_type element_type;

    /*quadrature*/
    //typedef IM_PK<Dim, imOrder, value_type> im_type;
    //typedef IM<Dim, imOrder, value_type, Simplex> im_type;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef Exporter<mesh_type>::timeset_type timeset_type;

    OneD( int argc, char** argv, AboutData const& ad )
        :
        super( argc, argv, ad ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm(), "Export-oned" ) ),
        //timeSet( new timeset_type( "oned" ) ),
        timers(),
        stats()
    {
        LOG(INFO) << "[OneD] hsize = " << meshSize << "\n";
        //LOG(INFO) << "[OneD] export = " << this->vm().count("export") << "\n";

        //timeSet->setTimeIncrement( 1.0 );
        //exporter->addTimeSet( timeSet );
    }

    OneD( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm(), "Export-oned" ) ),
        //timeSet( new timeset_type( "oned" ) ),
        timers(),
        stats()
    {
        LOG(INFO) << "[OneD] hsize = " << meshSize << "\n";
        LOG(INFO) << "[OneD] export = " << this->vm().count( "export" ) << "\n";

        //timeSet->setTimeIncrement( 1.0 );
        //exporter->addTimeSet( timeSet );
    }

    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptr_type createMesh( double meshSize );

    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * solve system
     */
    void solve( sparse_matrix_ptrtype const& D, element_type& u, vector_ptrtype const& F );

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    template<typename f1_type, typename f2_type>
    void exportResults( double t,
                        f1_type& u,
                        f2_type& v );

private:

    backend_ptrtype M_backend;

    double meshSize;
    double bcCoeff;

    boost::shared_ptr<export_type> exporter;
    //export_type::timeset_ptrtype timeSet;

    std::map<std::string,std::pair<boost::timer,double> > timers;
    std::map<std::string,double> stats;
}; // OneD

OneD::mesh_ptr_type
OneD::createMesh( double meshSize )
{
    timers["mesh"].first.restart();


    GeoTool::Node x1( 0 );
    GeoTool::Node x2( 1 );
    GeoTool::Line L( meshSize,"hola",x1,x2 );
    L.setMarker( _type="point",_name="Left",_marker1=true );
    L.setMarker( _type="point",_name="Right",_marker2=true );
    L.setMarker( _type="line",_name="Omega",_markerAll=true );


    auto mesh = L.createMesh<mesh_type>( "hola" );
#if 0
    GmshHypercubeDomain<entity_type::nDim,entity_type::nOrder,entity_type::nRealDim,Simplex> td;
    td.setCharacteristicLength( meshSize );
    td.setX( std::make_pair( 0.0, 1.0 ) );
    std::string fname = td.generate( entity_type::name().c_str() );

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );
#endif
    timers["mesh"].second = timers["mesh"].first.elapsed();
    LOG(INFO) << "[timer] createMesh(): " << timers["mesh"].second << "\n";

    return mesh;
} // OneD::createMesh


void
OneD::run()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    //    int maxIter = 10.0/meshSize;
    using namespace Feel::vf;

    this->changeRepository( boost::format( "%1%/%2%/P%3%/h_%4%/" )
                            % this->about().appName()
                            % entity_type::name()
                            % Order
                            % this->vm()["hsize"].as<double>()
                          );
    /*
     * logs will be in <feel repo>/<app name>/<entity>/P<p>/h_<h>
     */
    this->setLogs();

    /*
     * First we create the mesh
     */
    mesh_ptr_type mesh = createMesh( meshSize );
    stats["nelt"] = mesh->elements().size();
    //const int Left = 1;
    //const int Right = 2;

    LOG(INFO) << "----------------1----------------\n";

    /*
     * The function space and some associate elements are then defined
     */
    timers["init"].first.restart();
    space_ptrtype Xh = space_type::New( mesh );
    //Xh->dof()->showMe();
    element_type An( Xh, "An" );
    element_type Qn( Xh, "Qn" );
    element_type Anp1( Xh, "Anp1" );
    element_type Qnp1( Xh, "Qnp1" );
    timers["init"].second = timers["init"].first.elapsed();
    stats["ndof"] = Xh->nDof();

    LOG(INFO) << "----------------2----------------\n";

    /*
     * a quadrature rule for numerical integration
     */
    //im_type im;

    value_type pi = 4.0*math::atan( 1.0 );

    value_type dt = this->vm()["dt"].as<value_type>();
    value_type Tf = this->vm()["Tf"].as<value_type>();

    value_type rho = this->vm()["rho"].as<value_type>();
    value_type E_coeff = this->vm()["E"].as<value_type>();
    value_type h0_coeff = this->vm()["h0"].as<value_type>();
    value_type beta_coeff = E_coeff*h0_coeff*math::sqrt( pi )/( 1-0.5*0.5 );
    value_type A0_coeff = this->vm()["A0"].as<value_type>();
    value_type alpha = this->vm()["alpha"].as<value_type>();
    value_type Kr = this->vm()["Kr"].as<value_type>();

    LOG(INFO) << "----------------3----------------\n";



    AUTO( E , constant( E_coeff ) );
    AUTO( h0 , constant( h0_coeff ) );
    AUTO( A0 , constant( A0_coeff ) );
    AUTO( beta , E*h0*math::sqrt( pi )/( 1-0.5*0.5 ) );


    element_type Beta( Xh, "Beta" );
    Beta = vf::project( Xh, elements( mesh ), beta );
    element_type Area0( Xh, "Area0" );
    Area0 = vf::project( Xh, elements( mesh ), A0 );

    element_type FA( Xh, "FA" );
    element_type FQ( Xh, "FQ" );
    vector_ptrtype rhsA( M_backend->newVector( Xh ) );
    vector_ptrtype rhsQ( M_backend->newVector( Xh ) );

    LOG(INFO) << "----------------4----------------\n";

    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );
    sparse_matrix_ptrtype MA( M_backend->newMatrix( Xh, Xh ) );
    sparse_matrix_ptrtype MQ( M_backend->newMatrix( Xh, Xh ) );

    form2( Xh, Xh, M, _init=true ) =
        integrate( elements( *mesh ), // im,
                   idt( Qn )*id( Qn )
                 );

    M->close();

    // initialization
    An = ublas::scalar_vector<double>( An.size(), A0_coeff );
    Qn = ublas::scalar_vector<double>( Qn.size(), 0 );

    LOG(INFO) << "----------------5----------------\n";

    int iteration = 1;

    // time loop
    for ( double t = dt; t < Tf; t+=dt, ++iteration )
    {
        LOG(INFO) << "============================================================\n";
        LOG(INFO) << "time = " << t << "\n";

        auto H = mat<2,2>( constant( 0 ), constant( 1 ),
                           -alpha*( ( idv( Qn )/idv( An ) )^( 2 ) )+beta*sqrt( idv( An ) )/( 2*rho*A0 ), 2*alpha*idv( Qn )/idv( An ) );
        auto F = vec( idv( Qn ), alpha*( idv( Qn )^( 2 ) )/idv( An ) + beta*pow( idv( An ),1.5 )/( 3*rho*A0 ) );
        auto Fx = trans( vec( constant( 1. ),constant( 0. ) ) )*F;
        auto Fy = trans( vec( constant( 0. ),constant( 1. ) ) )*F;

        //FA = vf::project( Xh, elements(*mesh), Fx );
        //FQ = vf::project( Xh, elements(*mesh), trans(vec(constant(0.),constant(1.)))*F);
        form1( Xh, rhsA )  = integrate( elements( mesh ), Fx*id( An ) );
        form1( Xh, rhsQ )  = integrate( elements( mesh ), Fy*id( Qn ) );
        this->solve( M, FA, rhsA );
        this->solve( M, FQ, rhsQ );

        auto B = vec( constant( 0 ),
                      Kr*idv( Qn )/idv( An )+
                      idv( An )/( A0*rho )*( 2/3*sqrt( idv( An ) )-sqrt( A0 ) )*gradv( Beta )-
                      beta*idv( An )/( rho*pow( A0,2 ) )*( 2/3*sqrt( idv( An ) )-sqrt( A0 ) )*gradv( Area0 ) );
        auto BU = mat<2,2>( constant( 0 ), constant( 0 ),
                            -Kr*idv( Qn )/pow( idv( An ),2 )+
                            1.0/( A0*rho )*( 2/3*sqrt( idv( An ) )-sqrt( A0 ) )*gradv( Beta )+
                            idv( An )/( A0*rho )*( 1/3*pow( idv( An ),-0.5 ) )*gradv( Beta )-

                            beta/( rho*pow( A0,2 ) )*( 2/3*sqrt( idv( An ) )-sqrt( A0 ) )*gradv( Area0 )-
                            beta*idv( An )/( rho*pow( A0,2 ) )*( 1/3*pow( idv( An ),-0.5 ) )*gradv( Area0 ), Kr/idv( An ) );
        auto Flw = F-constant( dt/2 )*H*B;
        auto Blw = B+constant( dt/2 )*BU*B;

        auto rhs = dt*Flw*grad( An )
                   - dt*dt/2*BU*vec( gradv( FA ),gradv( FQ ) )*id( An )
                   - dt*dt/2*H*vec( gradv( FA ),gradv( FQ ) )*grad( An )
                   - dt*Blw*id( Qn );



        timers["assembly"].first.restart();
        form1( Xh, rhsA )  = integrate( elements( *mesh ), idv( An )*id( An )+trans( vec( constant( 1. ),constant( 0. ) ) )*rhs );
        form1( Xh, rhsQ )  = integrate( elements( *mesh ), idv( Qn )*id( An )+trans( vec( constant( 0. ),constant( 1. ) ) )*rhs );

        timers["assembly"].second = timers["assembly"].first.elapsed();
        timers["assembly_F"].second = timers["assembly"].first.elapsed();


        timers["assembly"].first.restart();

        //size_type pattern = Pattern::COUPLED|Pattern::EXTENDED;
        size_type pattern = Pattern::COUPLED;
        form2( Xh, Xh, MA, _init=true, _pattern=pattern ) =
            integrate( elements( *mesh ),
                       idt( Qn )*id( Qn )
                     );

        MA->close();
        rhsA->close();
        form2( Xh, Xh, MQ, _init=true, _pattern=pattern ) =
            integrate( elements( *mesh ),
                       idt( Qn )*id( Qn )
                     );

        MQ->close();
        rhsQ->close();

        // 1- compute the eigenvalues and eigenvectors
        // define lambda and L
        ublas::vector<value_type> ubar( 2 );
        ubar[0] = Qn( 0.0 )/An( 0.0 );
        ubar[1] = Qn( 1.0 )/An( 1.0 );
        LOG(INFO) << "ubar(z=0)=" << ubar[0]<< "\n";
        LOG(INFO) << "ubar(z=L)=" << ubar[1]<< "\n";

        ublas::vector<value_type> c1( 2 );
        c1[0] = math::sqrt( beta_coeff/( 2*rho*A0_coeff ) )*math::pow( An( 0. ), .25 );
        c1[1] = math::sqrt( beta_coeff/( 2*rho*A0_coeff ) )*math::pow( An( 1. ), .25 );
        ublas::vector<value_type> calpha( 2 );
        calpha[0] = math::sqrt( c1[0]*c1[0]+ubar[0]*ubar[0]*alpha*( alpha-1 ) );
        calpha[1] = math::sqrt( c1[1]*c1[1]+ubar[1]*ubar[1]*alpha*( alpha-1 ) );
        ublas::matrix<double,ublas::row_major> L0( 2, 2 );
        L0( 0, 0 )= calpha[0]-alpha*ubar[0];
        L0( 0, 1 )= 1;
        L0( 1, 0 )= -calpha[0]-alpha*ubar[0];
        L0( 1, 1 )= 1;
        ublas::matrix<double,ublas::row_major> L1( 2, 2 );
        L1( 0, 0 )= calpha[1]-alpha*ubar[1];
        L1( 0, 1 )= 1;
        L1( 1, 0 )= -calpha[1]-alpha*ubar[1];
        L1( 1, 1 )= 1;
        LOG(INFO) << "L(z=0)=" << L0 << "]\n";
        LOG(INFO) << "L(z=L)=" << L1 << "]\n";

        ublas::vector<value_type> lambda0( 2 );
        lambda0[0] = alpha*ubar[0]+calpha[0];
        lambda0[1] = alpha*ubar[0]-calpha[0];
        ublas::vector<value_type> lambda1( 2 );
        lambda1[0] = alpha*ubar[1]+calpha[1];
        lambda1[1] = alpha*ubar[1]-calpha[1];

        LOG(INFO) << "Lambda(z=0)=[" << lambda0[0] << "," << lambda0[1] << "]\n";
        LOG(INFO) << "Lambda(z=L)=[" << lambda1[0] << "," << lambda1[1] << "]\n";


        // 2- compute the foot of the characteristic
        // solve for dz_i/dt = lambda_i
        ublas::vector<value_type> z0( 2 );
        z0[0] = 0.0 - dt*lambda0[0];
        z0[1] = 0.0 - dt*lambda0[1]; // this is inside (<0)
        ublas::vector<value_type> z1( 2 );
        z1[0] = 1.0 - dt*lambda1[0]; // this is inside(>0)
        z1[1] = 1.0 - dt*lambda1[1];
        LOG(INFO) << "foot(z=0)=[" << z0[1] << "\n";
        LOG(INFO) << "foot(z=K)=[" << z1[0] << "\n";


        // 3- solve the pseudo-characteristics variable ode
        // dZ_i(t,z_i)/dt + G_i(Z_1,Z_2) = 0 i=1,2
        // where G_i(Z1,Z2) = L B(U(zi))
        ublas::vector<value_type> Z0( 2 );
        Z0[0]=math::sin( 100*t ); // on Q
        Z0[1]= L0( 1,0 )*An( z0[1] )+L0( 1,1 )*Qn( z0[1] )- // l_2 Un^*
               dt*( L0( 1,0 )*0+ L0( 1,1 )*Kr*Qn( z0[1] )/An( z0[1] ) );

        ublas::vector<value_type> Z1( 2 );
        Z1[1]=L1( 1,0 )*A0_coeff; // absorbing condition (reference condition)
        Z1[0]= L1( 0,0 )*An( z1[0] )+L0( 0,1 )*Qn( z1[0] )- // l_1 Un^*
               dt*( L1( 0,0 )*0+ L1( 0,1 )*Kr*Qn( z1[0] )/An( z1[0] ) );

        LOG(INFO) << "Z(z=0)=" << Z0 << "\n";
        LOG(INFO) << "Z(z=L)=" << Z1 << "\n";

        // WARNING: we simplified the L*B(U^*) because Beta
        // and A0 are actually constant for the moment

        //aaAn(z0[1])/(A0_coeff*rho)*(2/3*sqrt(An(z0[1]))-sqrt(A0_coeff))*gradv(Beta)-
        //beta*idv(An)/(rho*pow(A0,2))*(2/3*sqrt(idv(An))-sqrt(A0))*gradv(Area0))

        // 4- define A and Q value to be imposed at z=0 and z=L by solving
        // a 2x2 system
        // solve L [A;Q] = [Z_1;Z_2]
        value_type Q_0=Z0[0];//(Z0[1]-L0(1,0)*Z0[0]/L0(0,0))/(L0(1,1)-L0(1,0)*L0(0,1)/L0(0,0));
        value_type A_0=Z0[1]/L0( 1,0 )-L0( 1,1 )*Q_0/L0( 1,0 ); //-L0(0,1)*Q_0/L0(0,0) + Z0[0]/L0(0,0);
        value_type Q_1=( Z1[1]-L1( 1,0 )*Z1[0]/L1( 0,0 ) )/( L1( 1,1 )-L1( 1,0 )*L1( 0,1 )/L1( 0,0 ) );
        value_type A_1=-L1( 0,1 )*Q_1/L1( 0,0 ) + Z1[0]/L1( 0,0 );

        LOG(INFO) << "A(z=0)=" << A_0 << " Q(z=0)=" << Q_0 << "\n";
        LOG(INFO) << "A(z=L)=" << A_1 << " Q(z=L)=" << Q_1 << "\n";

        form2( Xh, Xh, MA ) += on( markedfaces( mesh,"Left" ), Anp1, rhsA, constant( A_0 ) );
        form2( Xh, Xh, MA ) += on( markedfaces( mesh,"Right" ), Anp1, rhsA, constant( A_1 ) );
        form2( Xh, Xh, MQ ) += on( markedfaces( mesh,"Left" ), Qnp1, rhsQ, constant( Q_0 ) );
        form2( Xh, Xh, MQ ) += on( markedfaces( mesh,"Right" ), Qnp1, rhsQ, constant( Q_1 ) );
        timers["assembly"].second += timers["assembly"].first.elapsed();
        timers["assembly_D"].second += timers["assembly"].first.elapsed();

        if ( this->vm().count( "export-matlab" ) )
        {
            rhsA->printMatlab( "rhsA.m" );
            rhsQ->printMatlab( "rhsB.m" );
            MA->printMatlab( "MA.m" );
            MQ->printMatlab( "MQ.m" );

        }


        this->solve( MA, Anp1, rhsA );
        this->solve( MQ, Qnp1, rhsQ );

        {
            std::ostringstream fname_u;
            fname_u << "A-" << iteration << ".dat";
            std::ofstream ofs3( fname_u.str().c_str() );
            mesh_type::element_iterator it = Anp1.functionSpace()->mesh()->beginElementWithProcessId( Application::processId() );
            mesh_type::element_iterator en = Anp1.functionSpace()->mesh()->endElementWithProcessId( Application::processId() );
            Anp1.updateGlobalValues();

            for ( ; it!=en; ++it )
            {
                const int nLocalDof = 2;

                for ( size_type i = 0; i < nLocalDof; ++i )
                {
                    size_type dof0 = boost::get<0>( Anp1.functionSpace()->dof()->localToGlobal( it->id(), i ) );
                    ofs3 << std::setw( 5 ) << it->id() << " "
                         << std::setw( 5 ) << i << " "
                         << std::setw( 5 ) << dof0 << " "
                         << std::setw( 15 ) << Anp1.globalValue( dof0 ) << " ";
                    value_type a = it->point( 0 ).node()[0];
                    value_type b = it->point( 1 ).node()[0];

                    if ( i == 0 )
                        ofs3 << a;

                    else if ( i == 1 )
                        ofs3 <<  b;

                    else
                        ofs3 <<  a + ( i-1 )*( b-a )/( nLocalDof-1 );

                    ofs3 << "\n";

                }
            }

            ofs3.close();
        }
        {
            std::ostringstream fname_u;
            fname_u << "Q-" << iteration << ".dat";
            std::ofstream ofs3( fname_u.str().c_str() );
            mesh_type::element_iterator it = Qnp1.functionSpace()->mesh()->beginElementWithProcessId( Application::processId() );
            mesh_type::element_iterator en = Qnp1.functionSpace()->mesh()->endElementWithProcessId( Application::processId() );
            Qnp1.updateGlobalValues();

            for ( ; it!=en; ++it )
            {
                const int nLocalDof = 2;

                for ( size_type i = 0; i < nLocalDof; ++i )
                {
                    size_type dof0 = boost::get<0>( Qnp1.functionSpace()->dof()->localToGlobal( it->id(), i ) );
                    ofs3 << std::setw( 5 ) << it->id() << " "
                         << std::setw( 5 ) << i << " "
                         << std::setw( 5 ) << dof0 << " "
                         << std::setw( 15 ) << Qnp1.globalValue( dof0 ) << " ";
                    value_type a = it->point( 0 ).node()[0];
                    value_type b = it->point( 1 ).node()[0];

                    if ( i == 0 )
                        ofs3 << a;

                    else if ( i == 1 )
                        ofs3 <<  b;

                    else
                        ofs3 <<  a + ( i-1 )*( b-a )/( nLocalDof-1 );

                    ofs3 << "\n";

                }
            }

            ofs3.close();
        }
        this->exportResults( t, Anp1, Qnp1 );

        An = Anp1;
        Qn = Qnp1;

        LOG(INFO) << "[timer] run():     init: " << timers["init"].second << "\n";
        LOG(INFO) << "[timer] run(): assembly: " << timers["assembly"].second << "\n";
        LOG(INFO) << "[timer] run():     o D : " << timers["assembly_D"].second << "\n";
        LOG(INFO) << "[timer] run():     o F : " << timers["assembly_F"].second << "\n";
        LOG(INFO) << "[timer] run():     o M : " << timers["assembly_M"].second << "\n";
        LOG(INFO) << "[timer] run():     o L : " << timers["assembly_L"].second << "\n";
        LOG(INFO) << "[timer] run():     o i : " << timers["assembly_evaluate"].second << "\n";
        LOG(INFO) << "[timer] run():   solver: " << timers["solver"].second << "\n";
        LOG(INFO) << "[timer] run():   solver: " << timers["export"].second << "\n";
    }// time loop
} // OneD::run

void
OneD::solve( sparse_matrix_ptrtype const& D,
             element_type& u,
             vector_ptrtype const& F )
{
    timers["solver"].first.restart();

    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    M_backend->solve( D, D, U, F, false );
    u = *U;

    //LOG(INFO) << "u = " << u.container() << "\n";
    timers["solver"].second = timers["solver"].first.elapsed();
    LOG(INFO) << "[timer] solve(): " << timers["solver"].second << "\n";
} // OneD::solve

template<typename f1_type, typename f2_type>
void
OneD::exportResults( double t,
                     f1_type& U,
                     f2_type& V )
{
    timers["export"].first.restart();

    //timeset_type::step_ptrtype timeStep = timeSet->step( t );
    exporter->step( t )->setMesh( U.functionSpace()->mesh() );
    exporter->step( t )->add( "AA", U );
    exporter->step( t )->add( "QQ", V );
    exporter->save();

    timers["export"].second = timers["export"].first.elapsed();
    LOG(INFO) << "[timer] exportResults(): " << timers["export"].second << "\n";
} // OneD::export

const uint16_type OneD::Order;
const uint16_type OneD::Dim;
const uint16_type OneD::imOrder;
} // Feel




int
main( int argc, char** argv )
{
    using namespace Feel;

    /* define and run application */
    Feel::OneD oned( argc, argv, makeAbout(), makeOptions() );
    oned.run();
}






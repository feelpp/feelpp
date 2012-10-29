/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
             Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-09-11

  Copyright (C) 2006 EPFL
  Copyright (C) 2008 Université Joseph Fourier (Grenoble 1)

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
   \file splitting.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-09-11
 */

/* references
   - A. Ern and J.-L. Guermond: Theory and Practice of Finite Elements,
     Springer New York 2004 (splitting scheme)
   - E. Burman, M. Fernandez and P. Hansbo: Continuous interior penalty finite
     element method for Oseen's equations, SIAM J. Numer. Anal. Vol. 44, No. 3,
     pp. 1248-1274 (weak bc)
*/

#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporterensight.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelvf/vf.hpp>

#define VELOCITY_UPDATE 0

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description splittingoptions( "Splitting options" );
    splittingoptions.add_options()
    ( "dt", Feel::po::value<double>()->default_value( 0.1 ), "time step value" )
    ( "ft", Feel::po::value<double>()->default_value( 1 ), "Final time value" )
    ( "nu", Feel::po::value<double>()->default_value( 1 ), "viscosity value" )
    ( "ulid", Feel::po::value<double>()->default_value( 1 ), "lid velocity" )
    ( "peps", Feel::po::value<double>()->default_value( 0. ), "epsilon for pressure term" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.5 ), "first h value to start convergence" )
    ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
    ( "bctype", Feel::po::value<int>()->default_value( 0 ), "Dirichlet condition type(0=elimination,1=penalisation, 2=weak" )
    ( "divtol", Feel::po::value<double>()->default_value( 1.e10 ), "divergence tolerance" )

    ;
    Feel::po::options_description convdomoptions( "Convection dominated options" );
    convdomoptions.add_options()
    ( "supg", "Steamline Upwind Petrov Galerkin formulation" )
    ( "gals", "Galerkin Least Square formulation" )
    ( "dwg", "Douglas Wang formulation" )
    ( "ls", Feel::po::value<double>()->default_value( 0 ), "control symmetric coefficient, SUPG=0, GALS=1, DWG=-1" )
    ;
    return splittingoptions.add( convdomoptions ).add( Feel::feel_options() );
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "splitting" ,
                           "splitting" ,
                           "0.1",
                           "2D and 3D Cavity Problem in Splitting Approach",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2006 EPFL" );

    about.addAuthor( "Christoph Winkelmann", "developer", "christoph.winkelmann@epfl.ch", "" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}


namespace Feel
{
template<int Dim>
class Splitting : public Application
{
    typedef Application super;
public:

    // -- TYPEDEFS --
    static const uint16_type uOrder = 2;
    static const uint16_type pOrder = 1;

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<Backend<value_type> > backend_ptrtype;

    /* matrix */
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /* mesh */
    typedef Mesh<GeoEntity<Simplex<Dim, 1> > > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /* bases */
    typedef fusion::vector<fem::Lagrange<Dim, uOrder, Vectorial, Continuous, double> >
    basis_U_type;
    typedef fusion::vector<fem::Lagrange<Dim, pOrder, Scalar, Continuous, double> >
    basis_p_type;

    /* spaces */
    typedef FunctionSpace<mesh_type, basis_U_type, value_type> space_U_type;
    typedef FunctionSpace<mesh_type, basis_p_type, value_type> space_p_type;
    typedef boost::shared_ptr<space_U_type> space_U_ptrtype;
    typedef boost::shared_ptr<space_p_type> space_p_ptrtype;
    typedef typename space_U_type::element_type element_U_type;
    typedef typename space_p_type::element_type element_p_type;

    /* quadrature */
    typedef IM<Dim, 3*uOrder-1, value_type, Simplex> im_u_type;
    typedef IM<Dim, 2*pOrder,   value_type, Simplex> im_p_type;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef typename Exporter<mesh_type>::timeset_type timeset_type;

    Splitting( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        bcCoeff( this->vm()["bccoeff"].template as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm()["exporter"].template as<std::string>() )->setOptions( this->vm() ) ),
        timeSet( new timeset_type( "splitting" ) ),
        timers()
    {
        timeSet->setTimeIncrement( 1.0 );
        exporter->addTimeSet( timeSet );
        exporter->setPrefix( "splitting" );
    }



    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptrtype createMesh( double meshSize );


    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * solve non-symmetric system
     */
    void solve( sparse_matrix_ptrtype const& D, element_p_type& u, vector_ptrtype const& F, bool is_sym = false );
    void solve( sparse_matrix_ptrtype const& D, element_U_type& u, vector_ptrtype const& F, bool is_sym = false );

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( double time, element_U_type& U, element_p_type& p );

private:

    backend_ptrtype M_backend;

    double meshSize;
    double bcCoeff;
    int bctype;
    double nu;

    boost::shared_ptr<export_type> exporter;
    typename export_type::timeset_ptrtype timeSet;

    std::map<std::string,std::pair<boost::timer,double> > timers;
}; // Splitting

template<int Dim>
typename Splitting<Dim>::mesh_ptrtype
Splitting<Dim>::createMesh( double meshSize )
{
    timers["mesh"].first.restart();
    mesh_ptrtype mesh( new mesh_type );

    Gmsh __gmsh;
    __gmsh.setOrder( GMSH_ORDER_ONE );
    std::ostringstream ostr;
    std::ostringstream nameStr;
    std::string fname;

    switch ( Dim )
    {
    case 2:
        fname = __gmsh.generateSquare( "splitting2d", meshSize );
        break;

    case 3:
        fname = __gmsh.generateCube( "splitting3d", meshSize );
        break;

    default:
        std::ostringstream os;
        os << "invalid dimension: " << Dim;
        throw std::logic_error( os.str() );
    }

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );
    timers["mesh"].second = timers["mesh"].first.elapsed();
    LOG(INFO) << "[timer] createMesh(): " << timers["mesh"].second << "\n";
    return mesh;
} // Splitting::createMesh


template<int Dim>
void
Splitting<Dim>::run()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    this->changeRepository( boost::format( "%1%/h_%2%/Re_%3%" )
                            % this->about().appName()
                            % meshSize
                            % ( 1./nu ) );
    using namespace Feel::vf;

    /*
     * First we create the mesh
     */
    mesh_ptrtype mesh = createMesh( meshSize );

    /*
     * The function spaces and some associate elements are then defined
     */
    timers["init"].first.restart();
    space_U_ptrtype Uh = space_U_type::New( mesh );
    element_U_type U( Uh, "U" );
    element_U_type u( Uh, "u" );
    element_U_type v( Uh, "v" );
    element_U_type Un( Uh, "Un" );

    space_p_ptrtype Yh = space_p_type::New( mesh );
    //p = project( Yh, elements(mesh), constant(0.0) );
    element_p_type p( Yh, "p" );
    element_p_type q( Yh, "p" );
    element_p_type phi( Yh, "phi" );
    timers["init"].second = timers["init"].first.elapsed();

    /*
     * Construction of the right hand side
     *
     * \f$ f = \int_\Omega g * v \f$ where \f$ g \f$ is a vector
     * directed in the \f$ y \f$ direction.
     */
    vector_ptrtype F( M_backend->newVector( Uh ) );
    vector_ptrtype fp( M_backend->newVector( Yh ) );

    timers["assembly"].first.restart();
    timers["assembly"].second = timers["assembly"].first.elapsed();

    // --- Construction of regularized laplacian on pressure space Ap
    sparse_matrix_ptrtype Ap( M_backend->newMatrix( Uh, Uh ) );
    double peps = this->vm()["peps"].template as<double>();
    // Note: peps = 0 leads to a singular matrix, because there it contains
    //       a null space of dimension one (the constant pressures)
    //       However, this is not a problem for the Krylov solver, because
    //       the Krylov subspace built up will not contain the constant
    //       unless you start with it as an initial vector. The incomplete
    //       factorization preconditioner doesn't have a problem either, as
    //       it would encounter a zero pivot only at the last step, where the
    //       factorization would in fact be complete. So peps = 0 is safe.
    form2( Uh, Uh, Ap, _init=true ) = integrate( elements( mesh ), im_p_type(),
                                      gradt( p )*trans( grad( q ) )+ peps*idt( p )*id( p )
                                               );
    Ap->close();

    // --- Construction of velocity mass matrix Mu
#if VELOCITY_UPDATE
    sparse_matrix_ptrtype Mu( M_backend->newMatrix( Uh, Uh ) );
    form2( Uh, Uh, Mu, _init=true ) = integrate( elements( mesh ), im_u_type(),
                                      trans( idt( U ) )*id( U )
                                               );
    Mu->close();
#endif

    double dt = this->vm()["dt"].template as<double>();
    double time  = dt;
    double ulid = this->vm()["ulid"].template as<double>();

    // --- Construction of velocity diffusion-reaction matrix Lu
    sparse_matrix_ptrtype Lu( M_backend->newMatrix( Uh, Uh ) );
    form2( Uh, Uh, Lu, _init=true ) =
        integrate( elements( mesh ), im_u_type(),
                   nu*dt*trace( gradt( u )*trans( grad( u ) ) )+
                   trans( idt( u ) )*id( u ) );

    if ( bctype == 2 ) // weak bc
    {
        form2( Uh, Uh, Lu ) +=
            integrate( boundaryfaces( mesh ), im_u_type(),
                       nu*dt*( - trans( gradt( u )*N() )*id( u )
                               - trans( grad( u )*N() )*idt( u )
                               + bcCoeff/hFace() * trans( idt( u ) ) * id( u ) )
                     );
    }

    Lu->close();

    // --- Time loop
    for ( int iter = 0;
            time < this->vm()["ft"].template as<double>();
            ++iter, time += dt )
    {
        double divTol = this->vm()["divtol"].template as<double>();
        double divError = 2*divTol + 1;
        int subiter = 0;

        while ( divError > divTol )
        {
            ++subiter;
            timers["assembly"].first.restart();

            // --- update rhs for solve for intermediate ux and uy
            form1( Uh, F, _init=true ) = integrate( elements( mesh ), im_u_type(),
                                                    ( trans( idv( u ) )-dt*grad( p ) )*id( u )
                                                  );

            // --- Construction of convection operator on velocity
            //     space
            sparse_matrix_ptrtype Ac( M_backend->newMatrix( Uh, Uh ) );
            form2( Uh, Uh, Ac, _init=true ) =
                integrate( elements( mesh ), im_u_type(),
                           dt*trans( gradt( u )*idv( Un ) )*id( u ) );
            Ac->close();

            // --- Construction of convection-diffusion-reaction
            //     operators on velocity space Aux and Auy
            sparse_matrix_ptrtype Aux( M_backend->newMatrix( Uh, Uh ) );
#warning TODO
            //gmm::add(Lu.mat(), Ac.mat(),Aux.wmat());

            switch ( bctype )
            {
            case 0: // elimination
            case 1: // penalisation
                // use last arg 'false' to tell the form to
                // _not_ initialise its representation
                form2( Uh, Uh, Aux ) +=
                    on( markedfaces( mesh,10 ), u, F,
                        oneX() );
                break;

            case 2: // weak
                form1( Uh, F ) +=
                    integrate( markedfaces( mesh,20 ), im_u_type(),
                               dt*nu*trans( oneX() )* ( bcCoeff/hFace()*id( u )- grad( u )*N() ) );
                break;

            default: // wrong
                std::ostringstream os;
                os << "invalid bctype: " << bctype;
                throw std::logic_error( os.str() );
            }

            Aux->close();

            timers["assembly"].second += timers["assembly"].first.elapsed();

            /*
             * Solution phase
             */
            // --- solve for u
            this->solve( Aux, Un, F );


            // rhs for pressure increment
            form1( Yh, fp, _init=true ) = integrate( elements( mesh ), im_p_type(), -divv( Un )*id( p )/dt );

            // --- solve for pressure increment phi
            this->solve( Ap, phi, fp );

            divError = math::sqrt( integrate( elements( mesh ), im_p_type(),
                                              divv( Un )^2 ).evaluate()( 0,0 ) );
            std::cout << "[Splitting] ||div u||_2 = " << divError << std::endl;

            // --- update pressure
            p += phi;

#if VELOCITY_UPDATE
            // --- update velocity
            form1( Uh, F ) = integrate( elements( mesh ), im_u_type(),
                                        ( trans( idv( Un ) )-dt*gradv( phi ) )*id( U )
                                      );

            this->solve( Mu, Un, F );

            divError = std::sqrt( integrate( elements( mesh ), im_p_type(),
                                             divv( Un )^2 ).evaluate() );
            std::cout << "[Splitting] ||div u||_2 = " << divError << std::endl;
#endif

        } // inner loop

        U = Un;

        this->exportResults( time, U, p );

        LOG(INFO) << "[Splitting] t = " << time << ", " << subiter << " subiterations" << "\n";
    } // time loop

    LOG(INFO) << "[timer] run():     init: " << timers["init"].second << "\n";
    LOG(INFO) << "[timer] run(): assembly: " << timers["assembly"].second << "\n";

} // Splitting::run

template<int Dim>
void
Splitting<Dim>::solve( sparse_matrix_ptrtype const& D,
                       element_p_type & u,
                       vector_ptrtype const& F,
                       bool is_sym )
{
    timers["solver"].first.restart();

    //M_backend->set_symmetric( is_sym );

    vector_ptrtype U = M_backend->newVector( u.functionSpace()->map() );
    M_backend->solve( D, D, U, F );
    u = *U;

    timers["solver"].second = timers["solver"].first.elapsed();
    LOG(INFO) << "[timer] solveNonSym(): " << timers["solver"].second << "\n";
} // Splitting::solveNonSym

template<int Dim>
void
Splitting<Dim>::solve( sparse_matrix_ptrtype const& D,
                       element_U_type & u,
                       vector_ptrtype const& F,
                       bool is_sym )
{
    timers["solver"].first.restart();

    //M_backend->set_symmetric( is_sym );

    vector_ptrtype U = M_backend->newVector( u.functionSpace()->map() );
    M_backend->solve( D, D, U, F );
    u = *U;

    timers["solver"].second = timers["solver"].first.elapsed();
    LOG(INFO) << "[timer] solveNonSym(): " << timers["solver"].second << "\n";
} // Splitting::solveNonSym


template<int Dim>
void
Splitting<Dim>::exportResults( double time, element_U_type& U, element_p_type& p )
{
    timers["export"].first.restart();

    // -- EXPORT --
    if ( this->vm().count( "export" ) )
    {
        typename timeset_type::step_ptrtype timeStep = timeSet->step( time );
        timeStep->setMesh( U.functionSpace()->mesh() );
        timeStep->add( "u", U );
        timeStep->add( "p", p );
        exporter->save();
    } // export

    timers["export"].second = timers["export"].first.elapsed();
    LOG(INFO) << "[timer] exportResults(): " << timers["export"].second << "\n";
} // Splitting::export
} // Feel

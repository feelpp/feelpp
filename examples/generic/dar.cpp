/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-10-29

  Copyright (C) 2007-2008 University Joseph Fourier Grenoble 1

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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file dar.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-10-29
 */
#include <life/options.hpp>
#include <life/lifecore/application.hpp>
#include <life/lifecore/application.hpp>
#include <life/lifediscr/functionspace.hpp>
#include <life/lifepoly/im.hpp>

#include <life/lifefilters/importergmsh.hpp>
#include <life/lifefilters/exporter.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>
#include <life/lifepoly/polynomialset.hpp>

#include <life/lifealg/backend.hpp>

#include <life/lifemesh/elements.hpp>

#include <life/lifevf/vf.hpp>

inline
Life::po::options_description
makeOptions()
{
    Life::po::options_description daroptions("Dar options");
    daroptions.add_options()
        ("penal", Life::po::value<double>()->default_value( 0.5 ), "penalisation parameter")
        ("f", Life::po::value<double>()->default_value( 0 ), "forcing term")
        ("stab", Life::po::value<bool>()->default_value( true ), "true to enable stabilisation, false otherwise")
        ("bx", Life::po::value<double>()->default_value( 1.0 ), "convection X component")
        ("by", Life::po::value<double>()->default_value( 0.0 ), "convection Y component")
        ("bz", Life::po::value<double>()->default_value( 0.0 ), "convection Z component")
        ("mu", Life::po::value<double>()->default_value( 1.0 ), "reaction coefficient component")
        ("epsilon", Life::po::value<double>()->default_value( 1.0 ), "diffusion coefficient")
        ("hsize", Life::po::value<double>()->default_value( 0.5 ), "first h value to start convergence")
        ("bctype", Life::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet")
        ("on-sym", Life::po::value<bool>()->default_value( true ), "true: on() is symmetric, false otherwise")
        ("on-diag", Life::po::value<bool>()->default_value( true ), "true: on() keep diagonal value , false otherwise")

        ("bccoeff", Life::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions")
        ("export", "export results(ensight, data file(1D)")
        ("export-matlab", "export matrix and vectors in matlab" )
        ;
    return daroptions.add( Life::life_options() ) ;
}
inline
Life::AboutData
makeAbout()
{
    Life::AboutData about( "dar" ,
                           "dar" ,
                           "0.1",
                           "Dar equation on simplices or simplex products",
                           Life::AboutData::License_GPL,
                           "Copyright (c) 2007 University Joseph Fourier Grenoble 1");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
   return about;

}


namespace Life
{
/**
 * Diffussion Advection Reaction Solver
 *
 * solve \f$-\epsilon \Delta u -\beta\cdot\nabla u + \mu u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma_{in}\f$
 */
template<int Dim,
         int Order,
         typename Cont,
         template<uint16_type,uint16_type,uint16_type> class Entity>
class Dar : public Application

{
    typedef Application super;

public:

    // -- TYPEDEFS --
    static const uint16_type imOrder = 2*Order;

    typedef double value_type;


    typedef Application application_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Entity<Dim,1,Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptr_type;

    template<typename Conti = Cont>
    struct space
    {
        /*basis*/
        typedef fusion::vector<Lagrange<Order, Scalar> > basis_type;
        /*space*/
        typedef FunctionSpace<mesh_type, basis_type, Conti, value_type> type;
        typedef boost::shared_ptr<type> ptrtype;
        typedef typename type::element_type element_type;
        typedef typename element_type::template sub_element<0>::type element_0_type;
        typedef typename element_type::template sub_element<1>::type element_1_type;
    };
    /*quadrature*/
    //typedef IM_PK<Dim, imOrder, value_type> im_type;
    typedef IM<Dim, imOrder, value_type, Entity> im_type;

    /* export */
    typedef Exporter<mesh_type> export_type;

    Dar( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        bcCoeff( this->vm()["bccoeff"].template as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) ),
        timers(),
        stats()
    {
        Log() << "[Dar] hsize = " << meshSize << "\n";
        Log() << "[Dar] bccoeff = " << bcCoeff << "\n";
        Log() << "[Dar] export = " << this->vm().count("export") << "\n";


    }


    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptr_type createMesh( double meshSize );

    /**
     * alias for run()
     */
    void operator()()
    {
        run();
    }

    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * solve system
     */
    template<typename Mat, typename Vec1, typename Vec2>
    void solve( Mat const& D, Vec1& u, Vec2 const& F, bool is_sym );

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    template<typename f1_type, typename f2_type, typename f3_type>
    void exportResults( f1_type& u,
                        f2_type& v,
                        f3_type& e );

private:

    double meshSize;
    double bcCoeff;

    boost::shared_ptr<export_type> exporter;

    std::map<std::string,std::pair<boost::timer,double> > timers;
    std::map<std::string,double> stats;
}; // Dar

template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity>
typename Dar<Dim,Order,Cont,Entity>::mesh_ptr_type
Dar<Dim,Order,Cont,Entity>::createMesh( double meshSize )
{
    timers["mesh"].first.restart();
    mesh_ptr_type mesh( new mesh_type );

    GmshTensorizedDomain<entity_type::nDim,entity_type::nOrder,entity_type::nRealDim,Entity> td;
    td.setCharacteristicLength( meshSize );
    std::string fname = td.generate( entity_type::name().c_str() );

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );
    timers["mesh"].second = timers["mesh"].first.elapsed();
    Log() << "[timer] createMesh(): " << timers["mesh"].second << "\n";
    return mesh;
} // Dar::createMesh


template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Dar<Dim, Order, Cont, Entity>::run()
{
    if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }

    //    int maxIter = 10.0/meshSize;
    using namespace Life::vf;

    this->changeRepository( boost::format( "%1%/%2%/P%3%/h_%4%/" )
                            % this->about().appName()
                            % entity_type::name()
                            % Order
                            % this->vm()["hsize"].template as<double>()
                            );
    /*
     * logs will be in <life repo>/<app name>/<entity>/P<p>/h_<h>
     */
    this->setLogs();

    /*
     * First we create the mesh
     */
    mesh_ptr_type mesh = createMesh( meshSize );
    stats["nelt"] = mesh->elements().size();

    /*
     * The function space and some associate elements are then defined
     */
    timers["init"].first.restart();
    typename space<Cont>::ptrtype Xh = space<Cont>::type::New( mesh );
    //Xh->dof()->showMe();
    typename space<Cont>::element_type u( Xh, "u" );
    typename space<Cont>::element_type v( Xh, "v" );
    timers["init"].second = timers["init"].first.elapsed();
    stats["ndof"] = Xh->nDof();

    /*
     * a quadrature rule for numerical integration
     */
    im_type im;

    value_type penalisation = this->vm()["penal"].template as<value_type>();
    int bctype = this->vm()["bctype"].template as<int>();

    double beta_x = this->vm()["bx"].template as<value_type>();
    double beta_y = this->vm()["by"].template as<value_type>();
    double beta_z = this->vm()["bz"].template as<value_type>();
    if ( Dim == 1 )
        {
            beta_y = 0;
            beta_z = 0;
        }
    if ( Dim == 2 )
        {
            beta_z = 0;
        }
    value_type mu = this->vm()["mu"].template as<value_type>();
    value_type epsilon = this->vm()["epsilon"].template as<value_type>();
    bool stab = this->vm()["stab"].template as<bool>();
    value_type pi = 4.0*math::atan(1.0);

    Log() << "Data Summary:\n";
    Log() << "    beta = (" << beta_x << ", " << beta_y << ", " << beta_z << ");\n";
    Log() << "      mu = " << mu << "\n";
    Log() << " epsilon = " << epsilon << "\n";
    Log() << "    stab = " << stab << "\n";


    //AUTO( beta , vec(constant(beta_x),constant(beta_y),constant(beta_z)) );
    AUTO( beta , vec(constant(beta_x),constant(beta_y)) );
    AUTO( g ,    sin(pi*Px())*cos(pi*Py())*cos(pi*Pz()));
#if 0
    AUTO( grad_g, vec( +pi*cos(pi*Px())*cos(pi*Py())*cos(pi*Pz()),
                       -pi*sin(pi*Px())*sin(pi*Py())*cos(pi*Pz()),
                       -pi*sin(pi*Px())*cos(pi*Py())*sin(pi*Pz()) ) );
#else
    AUTO( grad_g, vec( +pi*cos(pi*Px())*cos(pi*Py())*cos(pi*Pz()),
                       -pi*sin(pi*Px())*sin(pi*Py())*cos(pi*Pz()) ) );
#endif
    AUTO( f , 2*(pi*pi)*epsilon*g+ trans(beta)*grad_g + mu*g );
    AUTO( delta, constant(1.0)/(1.0/h() + epsilon/(h()*h()) ) );
    AUTO( Aepsi, (-epsilon*trace(hess(v))+ grad(v)*beta + mu*id(v)) );
    AUTO( Aepsit, (-epsilon*trace(hesst(u))+ gradt(u)*beta + mu*idt(u)) );


    Xh->dof()->showMe();

    backend_ptrtype backend( backend_type::build( this->vm() ) );
    vector_ptrtype F( backend->newVector( Xh->map() ) );
    timers["assembly"].first.restart();
    form1( Xh, F, _init=true )  = integrate( elements(mesh), im, f*id(v) + stab*delta*f*Aepsi );
    Log() << "[dar] vector local assembly done\n";
    timers["assembly"].second = timers["assembly"].first.elapsed();
    timers["assembly_F"].second = timers["assembly"].first.elapsed();

    /*
     * Construction of the left hand side
     */
    sparse_matrix_ptrtype D( backend->newMatrix( Xh->map(), Xh->map() ) );
    timers["assembly"].first.restart();

    //size_type pattern = DOF_PATTERN_COUPLED|DOF_PATTERN_NEIGHBOR;
    size_type pattern = DOF_PATTERN_COUPLED;
    form2( Xh, Xh, D, _init=true, _pattern=pattern ) =
        integrate( elements(mesh), im,
                   epsilon*gradt(u)*trans(grad(v)) +(gradt(u)*beta)*id(v) + mu*idt(u)*id(v)
                   // stabilisation
                   + stab*delta*Aepsi*Aepsit
                   );
    Log() << "[dar] matrix local assembly done\n";
    D->close();
    F->close();
    Log() << "[dar] vector/matrix global assembly done\n";
    if ( this->vm().count( "export-matlab" ) )
        {
            F->printMatlab( "F.m" );
            D->printMatlab( "D.m" );
        }
    bool on_sym = this->vm()["on-sym"].template as<bool>();
    bool on_diag = this->vm()["on-diag"].template as<bool>();
    size_type on_op = ON_ELIMINATION;
    if ( on_sym )
        on_op |= ON_ELIMINATION_SYMMETRIC;
    if ( on_diag )
        on_op |= ON_ELIMINATION_KEEP_DIAGONAL;
    Log() << "On() operation : " << on_op << "\n";
    form2( Xh, Xh, D ) += on( boundaryfaces(mesh), u, F, g, on_op );

    Log() << "[dar] dirichlet condition applied\n";
    timers["assembly"].second += timers["assembly"].first.elapsed();
    timers["assembly_D"].second += timers["assembly"].first.elapsed();

    if ( this->vm().count( "export-matlab" ) )
        {
            F->printMatlab( "F_dir.m" );
            D->printMatlab( "D_dir.m" );
        }

    Log() << "[dar] starting solve for D\n";
    this->solve( D, u, F, ( bctype == 1 || !Cont::is_continuous ) );
    if ( this->vm().count( "export-matlab" ) )
        {

            u.printMatlab( "u.m" );
        }
    Log() << "[dar] solve for D done\n";

    typename space<Continuous>::ptrtype Xch = space<Continuous>::type::New( mesh );
    typename space<Continuous>::element_type uEx( Xch, "uEx" );
    sparse_matrix_ptrtype M( backend->newMatrix( Xch->map(), Xch->map() ) );
    vector_ptrtype L( backend->newVector( Xch->map() ) );
    timers["assembly"].first.restart();
    form2( Xch, Xch, M, _init=true, _pattern=pattern ) = integrate( elements( mesh ), im, idt(uEx)*id(uEx) );
    M->close();
    timers["assembly_M"].second += timers["assembly"].first.elapsed();

    form1( Xch, L, _init=true ) = integrate( elements( mesh ), im, g*id(uEx) );
    L->close();
    timers["assembly"].second += timers["assembly"].first.elapsed();
    timers["assembly_L"].second += timers["assembly"].first.elapsed();
    if ( this->vm().count( "export-matlab" ) )
        {
            M->printMatlab( "M.m" );
            L->printMatlab( "L.m" );

        }

    this->solve( M, uEx, L, true );
    if ( this->vm().count( "export-matlab" ) )
        {

            uEx.printMatlab( "uEx.m" );
        }
    timers["assembly"].first.restart();
    double error = integrate( elements(mesh), im, (idv(u)-g)*(idv(u)-g) ).evaluate()( 0, 0 );
    double global_error = 0;
    mpi::all_reduce( Application::comm(), error, global_error, std::plus<double>() );
    timers["assembly"].second += timers["assembly"].first.elapsed();
    timers["assembly_evaluate"].second += timers["assembly"].first.elapsed();

    Log() << "local  ||error||_0 =" << math::sqrt(error) << "\n";
    Log() << "global ||error||_0 = " << math::sqrt(global_error) << "\n";

    if ( Cont::is_continuous )
        this->exportResults( u, u, uEx );
    else
        {

            form1( Xch, L, _init=true ) = integrate( elements( mesh ), im, idv(u)*id(uEx) );
            typename space<Continuous>::element_type uc( Xch, "uc" );
            this->solve( M, uc, L, true );
            this->exportResults( u, uc, uEx );
        }

    Log() << "[timer] run():     init: " << timers["init"].second << "\n";
    Log() << "[timer] run(): assembly: " << timers["assembly"].second << "\n";
    Log() << "[timer] run():     o D : " << timers["assembly_D"].second << "\n";
    Log() << "[timer] run():     o F : " << timers["assembly_F"].second << "\n";
    Log() << "[timer] run():     o M : " << timers["assembly_M"].second << "\n";
    Log() << "[timer] run():     o L : " << timers["assembly_L"].second << "\n";
    Log() << "[timer] run():     o i : " << timers["assembly_evaluate"].second << "\n";
    Log() << "[timer] run():   solver: " << timers["solver"].second << "\n";
    Log() << "[timer] run():   solver: " << timers["export"].second << "\n";

} // Dar::run

template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity>
template<typename Mat, typename Vec1, typename Vec2>
void
Dar<Dim, Order, Cont, Entity>::solve( Mat const& D,
                                      Vec1& u,
                                      Vec2 const& F,
                                      bool is_sym )
{
    timers["solver"].first.restart();

    backend_ptrtype backend( backend_type::build( this->vm() ) );
    //backend.set_symmetric( is_sym );

    vector_ptrtype U( backend->newVector( u.map() ) );
    backend->solve( D, D, U, F, false );
    u = *U;

    timers["solver"].second = timers["solver"].first.elapsed();
    Log() << "[timer] solve(): " << timers["solver"].second << "\n";
} // Dar::solve

template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity>
template<typename f1_type, typename f2_type, typename f3_type>
void
Dar<Dim, Order, Cont, Entity>::exportResults( f1_type& U,
                                                    f2_type& V,
                                                    f3_type& E )
{
    timers["export"].first.restart();


    exporter->step(1)->setMesh( U.functionSpace()->mesh() );
    exporter->step(1)->add( "u", U );
    exporter->step(1)->add( "v", V );
    exporter->step(1)->add( "e", E );
    exporter->save();

    timers["export"].second = timers["export"].first.elapsed();
    Log() << "[timer] exportResults(): " << timers["export"].second << "\n";
} // Dar::export
} // Life




int
main( int argc, char** argv )
{
    using namespace Life;

    /* change parameters below */
    const int nDim = 2;
    const int nOrder = 1;
    typedef Continuous MyContinuity;
    //typedef Discontinuous MyContinuity;
    //typedef Life::Dar<nDim, nOrder, MyContinuity, SimplexProduct> dar_type;
    typedef Life::Dar<nDim, nOrder, MyContinuity, Simplex> dar_type;

    /* assertions handling */
    Life::Assert::setLog( "dar.assert");

    /* define and run application */
    dar_type dar( argc, argv, makeAbout(), makeOptions() );
    dar.run();
}





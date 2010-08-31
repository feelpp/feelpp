/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2006-11-23

  Copyright (C) 2006 University Joseph Fourier

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
   \file quad.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \author Benjamin Stamm <benjamin.stamm@epfl.ch>
  \date 2006-11-23
 */
#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelpoly/polynomialset.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feelmesh/elements.hpp>

#include <feel/feelvf/vf.hpp>


/**
 * Create the ring geometry
 */
std::pair<std::string,std::string> createRing( int Dim, double meshSize );


inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description advectionoptions("Advection options");
    advectionoptions.add_options()
        ("penal", Feel::po::value<double>()->default_value( 0.5 ), "penalisation parameter")
        ("f", Feel::po::value<double>()->default_value( 0 ), "forcing term")
//        ("g", Feel::po::value<double>()->default_value( 0 ), "boundary term")
        ("bx", Feel::po::value<double>()->default_value( 1.0 ), "convection X component")
        ("by", Feel::po::value<double>()->default_value( 0.0 ), "convection Y component")
        ("mu", Feel::po::value<double>()->default_value( 1.0 ), "reaction coefficient component")
        ("stiff", Feel::po::value<double>()->default_value( 1.0 ), "stiffness parameter of solution")
        ("ring", Feel::po::value<bool>()->default_value( 0 ), "0 = square computational domain, 1 = quarter of a ring as computational domain")
        ("hsize", Feel::po::value<double>()->default_value( 0.5 ), "first h value to start convergence")
        ("bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet")
        ("bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions")
        ("export", "export results(ensight, data file(1D)")

        ;
    return advectionoptions.add( Feel::feel_options() ) ;
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "advection" ,
                            "advection" ,
                            "0.2",
                            "nD(n=1,2,3)Advection equation on simplices or simplex products",
                            Feel::AboutData::License_GPL,
                            "Copyright (c) 2006 University Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    about.addAuthor("Benjamin Stamm", "developer", "benjamin.stamm@epfl.ch", "");
   return about;

}


namespace Feel
{
std::pair<std::string,std::string>
createRing( int Dim, double meshSize )
{
    std::ostringstream ostr;
    std::ostringstream nameStr;
//    std::string fname;//
    switch( Dim ) {
    case 2:
        ostr << "h=" << meshSize << ";\n"
             << "Point(1) = {0.1,0,0,h/2};\n"
             << "Point(2) = {1,0,0,h};\n"
             << "Point(3) = {0,1,0,h};\n"
             << "Point(4) = {0,0.1,0,h/2};\n"
             << "Point(5) = {0,0,0,h/2};\n"
             << "Line(1) = {1,2};\n"
             << "Circle(2) = {2,5,3};\n"
             << "Line(3) = {3,4};\n"
             << "Circle(4) = {4,5,1};\n"
             << "Line Loop(5) = {1,2,3,4};\n"
             << "Plane Surface(6) = {5};\n"
             << "Physical Line(10) = {1};\n"
             << "Physical Line(20) = {2};\n"
             << "Physical Line(30) = {3};\n"
             << "Physical Line(40) = {4};\n"
//             << "Physical Line(20) = {1,2,4};\n"
             << "Physical Surface(7) = {6};\n";
        nameStr << "ring." << meshSize;
//        fname = __gmsh.generateSquare( "advectiondg2d", meshSize );//
        break;
// To be added for 3D something like:
/*    case 3:
        ostr << "h=" << meshSize << ";\n"
             << "Point(1) = {-1,-1,-1,h};\n"
             << "Point(2) = {-1, 1,-1,h};\n"
             << "Point(3) = { 1, 1,-1,h};\n"
             << "Point(4) = { 1,-1,-1,h};\n"
             << "Line(1) = {1,4};\n"
             << "Line(2) = {4,3};\n"
             << "Line(3) = {3,2};\n"
             << "Line(4) = {2,1};\n"
             << "Line Loop(5) = {3,4,1,2};\n"
             << "Plane Surface(6) = {5};\n"
             << "Extrude Surface {6, {0,0,2}};\n"
             << "Physical Surface(10) = {15,23,6,28};\n"
             << "Physical Surface(20) = {19,27};\n"
             << "Surface Loop(31) = {28,15,-6,19,23,27};\n"
             << "Volume(1) = {31};\n"
             << "Physical Volume(2) = {1};\n";
        nameStr << "cube." << meshSize;
        break;*/
    default:
        std::ostringstream os;
        os << "invalid dimension: " << Dim;
        throw std::logic_error( os.str() );
    }
    return std::make_pair( nameStr.str(), ostr.str() );
}

/**
 * Advection Solver using discontinous approximation spaces
 *
 * solve \f$ -\beta\cdot\nabla u + \mu u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma_{in}\f$
 */
template<int Dim,
         int Order,
         typename Cont,
         template<uint16_type,uint16_type,uint16_type> class Entity>
class Advection
    :
        public Application
{
    typedef Application super;
public:

    // -- TYPEDEFS --
    static const uint16_type imOrder = 2*Order;

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Entity<Dim,1,Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

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
    typedef typename Exporter<mesh_type>::timeset_type timeset_type;

    Advection( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        bcCoeff( this->vm()["bccoeff"].template as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm()["exporter"].template as<std::string>() )->setOptions( this->vm() ) ),
        timeSet( new timeset_type( "advection" ) )
    {
        Log() << "[Advection] hsize = " << meshSize << "\n";
        Log() << "[Advection] bccoeff = " << bcCoeff << "\n";
        Log() << "[Advection] export = " << this->vm().count("export") << "\n";

        timeSet->setTimeIncrement( 1.0 );
        exporter->addTimeSet( timeSet );
        exporter->setPrefix( "advection" );
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
    typename export_type::timeset_ptrtype timeSet;
}; // Advection

template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity>
typename Advection<Dim,Order,Cont,Entity>::mesh_ptrtype
Advection<Dim,Order,Cont,Entity>::createMesh( double meshSize )
{
    mesh_ptrtype mesh( new mesh_type );

    std::ostringstream entity_str;
    if (entity_type::is_simplex_product )
        entity_str << "Hypercube";
    else
        entity_str << "Simplex";
    entity_str << "_"
               << entity_type::nDim
               << "_"
               << entity_type::nOrder;

    bool ring = this->vm()["ring"].template as<bool>();
    if (!ring)
	{
	    GmshHypercubeDomain<entity_type::nDim,entity_type::nOrder,entity_type::nRealDim,Entity> td;
	    td.setCharacteristicLength( meshSize );
	    ImporterGmsh<mesh_type> import( td.generate( entity_str.str() ) );
	    mesh->accept( import );
	}
    else
	{
	    Gmsh gmsh;
	    gmsh.setOrder( GMSH_ORDER_ONE );
	    std::string mesh_name, mesh_desc;
	    boost::tie( mesh_name, mesh_desc ) = createRing(Dim,meshSize);
	    std::string fname = gmsh.generate( mesh_name, mesh_desc );
	    ImporterGmsh<mesh_type> import( fname );
	    mesh->accept( import );
	}
    return mesh;
} // Advection::createMesh


template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Advection<Dim, Order, Cont, Entity>::run()
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
                            % this->vm()["hsize"].template as<double>()
                            );
    /*
     * logs will be in <feel repo>/<app name>/<entity>/P<p>/h_<h>
     */
    this->setLogs();

    /*
     * First we create the mesh
     */
    mesh_ptrtype mesh = createMesh( meshSize );

    /*
     * The function space and some associate elements are then defined
     */
    typename space<Cont>::ptrtype Xh = space<Cont>::type::New( mesh );
    //Xh->dof()->showMe();
    typename space<Cont>::element_type u( Xh, "u" );
    typename space<Cont>::element_type v( Xh, "v" );

    /*
     * a quadrature rule for numerical integration
     */
    im_type im;

    value_type penalisation = this->vm()["penal"].template as<value_type>();
    int bctype = this->vm()["bctype"].template as<int>();

    double beta_x = this->vm()["bx"].template as<value_type>();
    double beta_y = this->vm()["by"].template as<value_type>();
    value_type mu = this->vm()["mu"].template as<value_type>();
    value_type stiff = this->vm()["stiff"].template as<value_type>();
    bool ring = this->vm()["ring"].template as<bool>();

    //AUTO( r , sqrt(Px()*Px()+Py()*Py()) );
    //AUTO( r, norm2( P() ) );
    AUTO( r, sqrt(trans(P())*P()) );
    AUTO( beta , beta_x*(ring*Py()/r + !ring)*oneX()+beta_y*(-ring*Px()/r + !ring)*oneY() );
    //AUTO( beta , (ones<Dim,1>()) );
    AUTO( beta_N , (trans(N())*beta) );
    AUTO( beta_abs , abs(beta_N) );
    AUTO( beta_minus , constant(0.5)*(beta_abs-beta_N) );
    AUTO( beta_plus , constant(0.5)*(beta_abs+beta_N) );
    AUTO( g , ((!ring*exp(-mu*Px())*atan((Py()-	0.5)/stiff) + ring*exp(-mu*r*acos(Py()/r))*atan((r-0.5)/stiff)) ) );
    AUTO( f , ( constant(0.0) ) );

    backend_ptrtype backend( backend_type::build( this->vm() ) );
    vector_ptrtype F( backend->newVector( Xh ) );
    form1( Xh, F )  = integrate( elements(mesh), im, trans(f)*id(v) );
    if ( bctype == 1 || !Cont::is_continuous )
        form1( Xh, F, false ) += integrate( boundaryfaces(mesh), im, trans(beta_minus*g)*id(v) );


    /*
     * Construction of the left hand side
     */
    sparse_matrix_ptrtype D( backend->newMatrix( Xh, Xh ) );
    //size_type pattern = DOF_PATTERN_COUPLED|DOF_PATTERN_NEIGHBOR;
    size_type pattern = DOF_PATTERN_COUPLED;
    form2( Xh, Xh, D, _init=true, _pattern=pattern ) =
        integrate( elements(mesh), im,
                   // -(u,beta*grad(v))+(mu*u,v)-(u,div(beta)*v)
                   -trans(idt(u))*(grad(v)*beta) + mu*trans(idt(u))*id(v) //- idt(u)*id(v)*(dx(beta_x)+dy(beta_y))
                   );

    if ( !Cont::is_continuous )
        {
            form2( Xh, Xh, D ) +=integrate( internalfaces(mesh), im,
                                            // {beta u} . [v]
                                            //( trans(averaget(trans(beta)*idt(u))) * jump(trans(id(v))) )
                                            ( averaget(trans(beta)*idt(u)) * jump(id(v)) )
                                            // penal*[u] . [v]
                                            + penalisation*beta_abs*( trans(jumpt(trans(idt(u))))*jump(trans(id(v))) ) );
        }

    else // continuous case: stabilization by interior penalty
        form2( Xh, Xh, D ) +=integrate( internalfaces(mesh), im,
                                        // penal*[grad(u)] . [grad(v)]
                                        + penalisation*beta_abs*hFace()*hFace()*(trans(jumpt(gradt(u)))*jump(grad(v)) ) );

    if ( bctype == 1 || !Cont::is_continuous )
        {
            form2( Xh, Xh, D ) += integrate( boundaryfaces(mesh), im, beta_plus*trans(idt(u))*id(v) );

            D->close();
            F->close();
        }
    else if ( bctype == 0 )
        {
            D->close();
            F->close();
            form2( Xh, Xh, D, false ) += on( boundaryfaces(mesh), u, F, g );
        }

    F->printMatlab( "F" );
    D->printMatlab( "D" );

    this->solve( D, u, F, ( bctype == 1 || !Cont::is_continuous ) );

    typename space<Continuous>::ptrtype Xch = space<Continuous>::type::New( mesh );
    typename space<Continuous>::element_type uEx( Xch, "uEx" );
    sparse_matrix_ptrtype M( backend->newMatrix( Xch, Xch ) );
    form2( Xch, Xch, M ) = integrate( elements( mesh ), im, trans(idt(uEx))*id(uEx) );
    M->close();
    vector_ptrtype L( backend->newVector( Xch ) );
    form1( Xch, L ) = integrate( elements( mesh ), im, trans(g)*id(uEx) );
    this->solve( M, uEx, L, true );

    double error = integrate( elements(mesh), im, trans(idv(u)-g)*(idv(u)-g) ).evaluate()( 0, 0 );
    double global_error = 0;
    mpi::all_reduce( Application::comm(), error, global_error, std::plus<double>() );

    Log() << "local  ||error||_0 =" << math::sqrt(error) << "\n";
    Log() << "global ||error||_0 = " << math::sqrt(global_error) << "\n";

    if ( Cont::is_continuous )
        this->exportResults( u, u, uEx );
    else
        {

            form1( Xch, L ) = integrate( elements( mesh ), im, trans(idv(u))*id(uEx) );
            typename space<Continuous>::element_type uc( Xch, "uc" );
            this->solve( M, uc, L, true );
            this->exportResults( u, uc, uEx );
        }

} // Advection::run

template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity>
template<typename Mat, typename Vec1, typename Vec2>
void
Advection<Dim, Order, Cont, Entity>::solve( Mat const& D,
                                            Vec1& u,
                                            Vec2 const& F,
                                            bool is_sym )
{
    backend_ptrtype backend( backend_type::build( this->vm() ) );
    //backend.set_symmetric( is_sym );

    vector_ptrtype U( backend->newVector( u.functionSpace() ) );
    backend->solve( D, D, U, F, false );
    u = *U;

} // Advection::solve

template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity>
template<typename f1_type, typename f2_type, typename f3_type>
void
Advection<Dim, Order, Cont, Entity>::exportResults( f1_type& U,
                                                    f2_type& V,
                                                    f3_type& E )
{
    typename timeset_type::step_ptrtype timeStep = timeSet->step( 1.0 );
    timeStep->setMesh( U.functionSpace()->mesh() );
    timeStep->add( "u", U );
    timeStep->add( "v", V );
    timeStep->add( "e", E );
    exporter->save();

} // Advection::export
} // Feel




    int
    main( int argc, char** argv )
{
    using namespace Feel;

    /* change parameters below */
    const int nDim = 2;
    const int nOrder = 1;
    //typedef Continuous MyContinuity;
    typedef Discontinuous MyContinuity;
    typedef Feel::Advection<nDim, nOrder, MyContinuity, Hypercube> advection_type;

    /* assertions handling */
    Feel::Assert::setLog( "advection.assert");

    /* define and run application */
    advection_type advection( argc, argv, makeAbout(), makeOptions() );
    advection.run();
}





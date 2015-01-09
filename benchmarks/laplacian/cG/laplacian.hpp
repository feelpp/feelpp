/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4


  This file is part of the Feel library

  Author(s):
  Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  Perrimond Benoit <Benoit.Perrimond@bvra.e.ujf-grenoble.fr>
  Vincent Chabannes <vincent.chabannes@gmail.com>

  Copyright (C) 2008-2010 Université Joseph Fourier (Grenoble I)

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
   \file laplacian.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Perrimond Benoit <Benoit.Perrimond@bvra.e.ujf-grenoble.fr>
   \author Vincent Chabannes <vincent.chabannes@gmail.com>
   \date 2008-02-07
 */
#include <boost/timer.hpp>
#include <boost/shared_ptr.hpp>

#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelpoly/polynomialset.hpp>


#include <feel/feelvf/vf.hpp>
#include <fstream>
#include <sstream>

#include <feel/feelcore/applicationxml.hpp>
#include <feel/feelcore/xmlparser.hpp>

using namespace Feel;

inline
po::options_description
makeOptions()
{
    po::options_description laplacianoptions( "Laplacian options" );
    laplacianoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.5 ), "mesh size" )
    ( "shape", po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (hypercube, simplex, ellipsoid)" )
    ( "nu", po::value<double>()->default_value( 1 ), "coef diffusion" )
    ( "beta", po::value<double>()->default_value( 1 ), "coef reaction " )
    ( "gammabc", po::value<double>()->default_value( 80 ), "weak Dirichlet penalisation parameter " )
    ( "penal", Feel::po::value<double>()->default_value( 10 ), "jump penalisation parameter for dG" )
    ( "weak", Feel::po::value<bool>()->default_value( true ), "use weak dirichlet conditions" )
    ;
    return laplacianoptions.add( Feel::feel_options() );
}
inline
AboutData
makeAbout()
{
    AboutData about( "laplacian" ,
                     "laplacian" ,
                     "0.2",
                     "nD(n=1,2,3) Laplacian on simplices or simplex products",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2008-2010 Université Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Benoit Perrimond", "developer", "Benoit.Perrimond@bvra.e.ujf-grenoble.fr", "" );
    about.addAuthor( "Vincent Chabannes", "developer", "vincent.chabannes@gmail.com", "" );
    return about;

}


/**
 * Laplacian Solver using discontinous approximation spaces
 *
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 */
template<int Dim, int Order, int RDim = Dim, typename ContinuityType = Continuous, template<uint16_type,uint16_type,uint16_type> class Entity=Simplex>
class Laplacian
    :
public ApplicationXML
{
    typedef ApplicationXML super;
public:


    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Entity<Dim, 1,RDim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpace<mesh_type, bases<Lagrange<0, Scalar, Discontinuous> > > p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar, ContinuityType> > basis_type;
    typedef bases<Lagrange<Order+2, Scalar> > exact_basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    typedef FunctionSpace<mesh_type, exact_basis_type, value_type> exact_space_type;
    typedef boost::shared_ptr<exact_space_type> exact_space_ptrtype;
    typedef typename exact_space_type::element_type exact_element_type;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    Laplacian( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        backend( backend_type::build( this->vm() ) ),
        meshSize( doption("hsize") ),
        shape( soption("shape") ),

        M_use_weak_dirichlet( boption("weak") ),
        M_gammabc( doption("gammabc") ),

        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
    {
        if ( M_use_weak_dirichlet )
            LOG(INFO) << "use weak Dirichlet BC\n";

        if ( exporter->doExport() )
            LOG(INFO) << "export results to ensight format\n";

        Parameter h;

        if ( Dim == 1 )         //=== 1D ===
            if ( Order < 5 )
                h=Parameter( _name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.01:0.09:0.1" );

            else
                h=Parameter( _name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.08:0.09:0.4" );

        else if ( Dim == 2 )    //=== 2D ===
            if ( Order < 5 )
                h=Parameter( _name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.03:0.09:0.1" );

            else
                h=Parameter( _name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.06:0.09:0.1" );

        else

            //=== 3D ===
            switch ( Order )
            {
            case 1:
                h=Parameter( _name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.05:0.02:0.4" );
                break;

            case 2:
                h=Parameter( _name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.2:0.02:0.8" );
                break;

            case 3:
                h=Parameter( _name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.15:0.02:0.8" );
                break;

            case 4:
                h=Parameter( _name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.1:0.02:0.4" );
                break;

            case 5:
                h=Parameter( _name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.05:0.1:0.6" );
                break;
            }

        this->
        addParameter( Parameter( _name="cont",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( ContinuityType::is_continuous ).c_str() ) )
        .addParameter( Parameter( _name="dim",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( Dim  ).c_str() ) )
        .addParameter( Parameter( _name="order",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( Order  ).c_str() ) )
        .addParameter( Parameter( _name="beta",_type=CONT_ATTR,_latex="\\beta",_values="0.01:1:10" ) )
        .addParameter( Parameter( _name="nu",_type=CONT_ATTR,_latex="\\nu",_values="0.01:1:10" ) )
        .addParameter( h );

        std::vector<Parameter> depend;
        std::vector<std::string> funcs;
        depend.push_back( h );
        std::ostringstream oss;
        oss << "h**" << boost::lexical_cast<std::string>( Order + 1  ) ;
        funcs.push_back( oss.str() );
        oss.str( "" );
        std::vector<std::string> funcs2;
        oss << "h**" << boost::lexical_cast<std::string>( Order ) ;
        funcs2.push_back( oss.str() );

        this->
        addOutput( Output( _name="norm_L2",_latex="\\left\\| . \\right\\|_{L^2}",_dependencies=depend,_funcs=funcs ) )
        .addOutput( Output( _name="norm_H1",_latex="\\left\\| . \\right\\|_{H^1}",_dependencies=depend,_funcs=funcs2 ) );

    }

    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u, element_type& v );

private:

    backend_ptrtype backend;

    double meshSize;
    std::string shape;
    bool M_use_weak_dirichlet;
    double M_gammabc;

    export_ptrtype exporter;

}; // Laplacian

template<int Dim, int Order, int RDim, typename ContinuityType, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Laplacian<Dim, Order, RDim, ContinuityType, Entity>::run()
{
    boost::timer t1;

    this->addParameterValue( ContinuityType::is_continuous )
    .addParameterValue( Dim )
    .addParameterValue( Order )
    .addParameterValue( doption("beta") )
    .addParameterValue( doption("nu") )
    .addParameterValue( doption("hsize") );

    if ( this->preProcessing() == RUN_EXIT ) return;

    using namespace Feel::vf;

    /*
     * First we create the mesh
     */
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                                _usenames=true,
                                                _convex=( entity_type::is_hypercube )?"Hypercube":"Simplex",
                                                _shape=shape,
                                                _dim=Dim,
                                                _xmin=-1.,_ymin=-1.,_zmin=-1.,
                                                _h=meshSize ),
                                        _update=MESH_CHECK| MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_RENUMBER );
    LOG(INFO) << "mesh created in " << t1.elapsed() << "s\n";
    t1.restart();

    /*
     * The function space and some associate elements are then defined
     */
    auto Xh = space_type::New( mesh );
    auto u = Xh->element();
    auto v = Xh->element();
    LOG(INFO) << "[functionspace] Number of dof " << Xh->nLocalDof() << "\n";
    LOG(INFO) << "function space and elements created in " << t1.elapsed() << "s\n";
    t1.restart();

    exact_space_ptrtype Eh = exact_space_type::New( mesh );
    exact_element_type fproj( Eh, "f" );
    exact_element_type gproj( Eh, "g" );
    LOG(INFO) << "[functionspace] Number of dof " << Eh->nLocalDof() << "\n";
    LOG(INFO) << "function space and elements created in " << t1.elapsed() << "s\n";
    t1.restart();


    value_type nu = doption("nu");
    value_type beta = doption("beta");


    value_type pi = M_PI;
    auto g = sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() );
    auto gradg = trans( +pi*cos( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() )*unitX()+
                        -pi*sin( pi*Px() )*sin( pi*Py() )*cos( pi*Pz() )*unitY()+
                        -pi*sin( pi*Px() )*cos( pi*Py() )*sin( pi*Pz() )*unitZ() );
    auto f = ( pi*pi*Dim*nu+beta )*g;

    fproj = vf::project( Eh, elements( mesh ), f );
    gproj = vf::project( Eh, elements( mesh ), g );
    v = vf::project( Xh, elements( mesh ), g );

    // Construction of the right hand side

    auto F=  backend->newVector( Xh );


    form1( _test=Xh, _vector=F, _init=true ) =
        integrate( elements( mesh ), f*id( v ) )+
        integrate( markedfaces( mesh, "Neumann" ), nu*gradg*vf::N()*id( v ) );

    if ( M_use_weak_dirichlet || ( ContinuityType::is_continuous == false ) )
    {
        form1( Xh, F ) +=
            integrate( markedfaces( mesh,"Dirichlet" ),
                       g*( -nu*grad( v )*N()+M_gammabc*id( v )/hFace() ) );
    }

    F->close();
    LOG(INFO) << "F assembled in " << t1.elapsed() << "s\n";
    t1.restart();

    //Construction of the left hand side

    auto D = backend->newMatrix( Xh, Xh );


    size_type pattern = ( ContinuityType::is_continuous?Pattern::COUPLED:Pattern::COUPLED|Pattern::EXTENDED );
    form2( _trial=Xh, _test=Xh, _matrix=D, _init=true, _pattern=pattern );
    LOG(INFO) << "D initialized in " << t1.elapsed() << "s\n";
    t1.restart();

    form2( Xh, Xh, D ) +=
        integrate( elements( mesh ),
                   nu*( gradt( u )*trans( grad( v ) ) )
                   + beta*( idt( u )*id( v ) ) );
    LOG(INFO) << "D stiffness+mass assembled in " << t1.elapsed() << "s\n";
    t1.restart();

    if ( ContinuityType::is_continuous == false )
    {
        value_type penalisation = this->vm()["penal"].template as<value_type>();
        form2( Xh, Xh, D ) +=integrate( internalfaces( mesh ),
                                        // - {grad(u)} . [v]
                                        -averaget( gradt( u ) )*jump( id( v ) )
                                        // - [u] . {grad(v)}
                                        -average( grad( v ) )*jumpt( idt( u ) )
                                        // penal*[u] . [v]/h_face
                                        + penalisation* ( trans( jumpt( idt( u ) ) )*jump( id( v ) ) )/hFace()
                                      );
        LOG(INFO) << "D consistency and stabilization terms assembled in " << t1.elapsed() << "s\n";
        t1.restart();
    }

    if ( M_use_weak_dirichlet || ( ContinuityType::is_continuous == false ) )
    {

        form2( Xh, Xh, D ) += integrate( markedfaces( mesh,mesh->markerName( "Dirichlet" ) ),
                                         ( - nu*trans( id( v ) )*( gradt( u )*N() )
                                           - nu*trans( idt( u ) )*( grad( v )*N() )
                                           + M_gammabc*trans( idt( u ) )*id( v )/hFace() ) );
        LOG(INFO) << "D weak bc assembled in " << t1.elapsed() << "s\n";
        t1.restart();

    }

    D->close();

    LOG(INFO) << "D assembled in " << t1.elapsed() << "s\n";


    if ( ( M_use_weak_dirichlet == false )  && ContinuityType::is_continuous )
    {
        t1.restart();
        form2( Xh, Xh, D ) +=
            on( markedfaces( mesh, "Dirichlet" ), u, F, g );

        LOG(INFO) << "Strong Dirichlet assembled in " << t1.elapsed() << "s on faces " << mesh->markerName( "Dirichlet" ) << " \n";
    }

    t1.restart();

    backend_type::build( this->vm() )->solve( _matrix=D, _solution=u, _rhs=F );

    LOG(INFO) << "solve in " << t1.elapsed() << "s\n";
    t1.restart();

    double L2error2 =integrate( elements( mesh ),
                                ( idv( u )-g )*( idv( u )-g ) ).evaluate()( 0, 0 );
    double L2error =   math::sqrt( L2error2 );

    LOG(INFO) << "||error||_L2=" << L2error << "\n";
    LOG(INFO) << "L2 norm computed in " << t1.elapsed() << "s\n";
    t1.restart();


    double semiH1error2 =integrate( elements( mesh ),
                                    ( gradv( u )-gradg )*trans( gradv( u )-gradg ) ).evaluate()( 0, 0 ) ;

    LOG(INFO) << "semi H1 norm computed in " << t1.elapsed() << "s\n";
    t1.restart();

    double H1error =   math::sqrt( semiH1error2+L2error2 );


    LOG(INFO) << "||error||_H1=" << H1error << "\n";
    LOG(INFO) << "H1 norm computed in " << t1.elapsed() << "s\n";
    t1.restart();

    this->exportResults( u, v );


    this->addOutputValue( L2error ).addOutputValue( H1error );
    this->postProcessing();

} // Laplacian::run

template<int Dim, int Order, int RDim, typename ContinuityType, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Laplacian<Dim, Order, RDim, ContinuityType, Entity>::exportResults( element_type& U, element_type& v )
{
    if ( exporter->doExport() )
    {
        LOG(INFO) << "exportResults starts\n";

        exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );

        exporter->step( 0 )->add( "pid",
                                  regionProcess( boost::shared_ptr<p0_space_type>( new p0_space_type( U.functionSpace()->mesh() ) ) ) );
        exporter->step( 0 )->add( "u", U );
        exporter->step( 0 )->add( "exact", v );

        exporter->save();
    }
} // Laplacian::export


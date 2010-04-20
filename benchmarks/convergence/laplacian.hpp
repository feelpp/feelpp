/* -*- mode: c++ -*-


  This file is part of the Life library

  Author(s):
  Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
  Perrimond Benoit <Benoit.Perrimond@bvra.e.ujf-grenoble.fr>
  Vincent Chabannes <vincent.chabannes@gmail.com>

  Date: 2008-02-07

  Copyright (C) 2008 Universit� Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \author Perrimond Benoit <Benoit.Perrimond@bvra.e.ujf-grenoble.fr>
   \author Vincent Chabannes <vincent.chabannes@gmail.com>
   \date 2008-02-07
 */
#include <boost/timer.hpp>
#include <boost/shared_ptr.hpp>

#include <life/options.hpp>
#include <life/lifecore/life.hpp>

#include <life/lifealg/backend.hpp>

#include <life/lifediscr/functionspace.hpp>
#include <life/lifediscr/region.hpp>
#include <life/lifepoly/im.hpp>

#include <life/lifefilters/gmsh.hpp>
#include <life/lifefilters/exporter.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>
#include <life/lifepoly/polynomialset.hpp>


#include <life/lifevf/vf.hpp>
#include <fstream>
#include <sstream>

#include <life/lifecore/applicationxml.hpp>
#include <life/lifecore/xmlparser.hpp>

using namespace Life;

inline
po::options_description
makeOptions()
{
    po::options_description laplacianoptions("Laplacian options");
    laplacianoptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.5 ), "mesh size")
        ("nu", po::value<double>()->default_value( 1 ), "coef diffusion")
        ("beta", po::value<double>()->default_value( 1 ), "coef reaction " )
        ("gammabc", po::value<double>()->default_value( 20 ), "weak Dirichlet penalisation parameter " )

        ("weak", "use weak dirichlet conditions")
        ;
    return laplacianoptions.add( Life::life_options() );
}
inline
AboutData
makeAbout()
{
    AboutData about( "laplacian" ,
                     "laplacian" ,
                     "0.2",
                     "nD(n=1,2,3) Laplacian on simplices or simplex products",
                     Life::AboutData::License_GPL,
                     "Copyright (c) 2008 Universit� Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    about.addAuthor("Benoit Perrimond", "developer", "Benoit.Perrimond@bvra.e.ujf-grenoble.fr", "");
    about.addAuthor("Vincent Chabannes", "developer", "vincent.chabannes@gmail.com", "");
    return about;

}


/**
 * Laplacian Solver using discontinous approximation spaces
 *
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 */
template<int Dim, int Order, int RDim = Dim, template<uint16_type,uint16_type,uint16_type> class Entity=Simplex>
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

    typedef FunctionSpace<mesh_type, fusion::vector<Lagrange<0, Scalar> >, Discontinuous> p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;

    /*basis*/
    typedef fusion::vector<Lagrange<Order, Scalar> > basis_type;
    typedef fusion::vector<Lagrange<5, Scalar> > exact_basis_type;

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
        meshSize( this->vm()["hsize"].template as<double>() ),

        M_use_weak_dirichlet( this->vm().count( "weak" ) ),
        M_gammabc( this->vm()["gammabc"].template as<double>() ),

        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
    {
        if ( M_use_weak_dirichlet )
            Log() << "use weak Dirichlet BC\n";
        if ( exporter->doExport() )
            Log() << "export results to ensight format\n";
        Parameter h;
        if (Dim == 1)           //=== 1D ===
            if (Order < 5)
                h=Parameter(_name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.01:0.09:0.4" );
            else
                h=Parameter(_name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.08:0.09:0.4" );
        else if (Dim == 2)      //=== 2D ===
            if (Order < 5)
                h=Parameter(_name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.03:0.09:0.1" );
            else
                h=Parameter(_name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.06:0.09:0.1" );
        else
            //=== 3D ===
            switch( Order )
            {
            case 1:
                h=Parameter(_name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.05:0.02:0.08" );
                break;
            case 2:
                h=Parameter(_name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.1:0.02:0.2" );
                break;
            case 3:
                h=Parameter(_name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.15:0.02:0.25" );
                break;
            case 4:
            case 5:
                h=Parameter(_name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.2:0.1:0.25" );
                break;
            }
        this->
            addParameter( Parameter(_name="dim",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( Dim  ).c_str()) )
            .addParameter( Parameter(_name="order",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( Order  ).c_str() ) )
            .addParameter( Parameter(_name="beta",_type=CONT_ATTR,_latex="\\beta",_values="0.01:1:10") )
            .addParameter( Parameter(_name="nu",_type=CONT_ATTR,_latex="\\nu",_values="0.01:1:10") )
            .addParameter( h );

        std::vector<Parameter> depend;
        std::vector<std::string> funcs;
        depend.push_back(h);
        std::ostringstream oss;
        oss << "h**" << boost::lexical_cast<std::string>( Order + 1  ) ;
        funcs.push_back(oss.str());
        oss.str("");
        std::vector<std::string> funcs2;
        oss << "h**" << boost::lexical_cast<std::string>( Order ) ;
        funcs2.push_back(oss.str());

        this->
            addOutput( Output(_name="norm_L2",_latex="\\left\\| . \\right\\|_{L^2}",_dependencies=depend,_funcs=funcs) )
            .addOutput( Output(_name="norm_H1",_latex="\\left\\| . \\right\\|_{H^1}",_dependencies=depend,_funcs=funcs2) );

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
     * solve the system D u = F
     */
    void solve( sparse_matrix_ptrtype& D, element_type& u, vector_ptrtype& F );


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u, element_type& v );

private:

    backend_ptrtype backend;

    double meshSize;

    bool M_use_weak_dirichlet;
    double M_gammabc;

    export_ptrtype exporter;

}; // Laplacian

template<int Dim, int Order, int RDim, template<uint16_type,uint16_type,uint16_type> class Entity>
typename Laplacian<Dim,Order,RDim,Entity>::mesh_ptrtype
Laplacian<Dim,Order,RDim,Entity>::createMesh( double meshSize )
{
    mesh_ptrtype mesh( new mesh_type );

    GmshTensorizedDomain<entity_type::nDim,entity_type::nOrder,entity_type::nRealDim,Entity> td;
    td.setCharacteristicLength( meshSize );
    td.setX( std::make_pair( -1, 1 ) );
    if((Dim==1) && (RDim==2))  td.setY( std::make_pair( 1, 1 ) );
    if(Dim>=2)  td.setY( std::make_pair( -1, 1 ) );
    if((Dim==2) && (RDim==3))  td.setZ( std::make_pair( 1, 1 ) );
    if(Dim==3)  td.setZ( std::make_pair( -1, 1 ) );
    std::string fname = td.generate( entity_type::name().c_str() );

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );

    return mesh;
} // Laplacian::createMesh


template<int Dim, int Order, int RDim, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Laplacian<Dim, Order, RDim, Entity>::run()
{
    boost::timer t1;

    this->addParameterValue( Dim )
        .addParameterValue( Order )
        .addParameterValue( this->vm()["beta"].template as<double>() )
        .addParameterValue( this->vm()["nu"].template as<double>() )
        .addParameterValue( this->vm()["hsize"].template as<double>() );

    if (this->preProcessing() == RUN_EXIT) return;

    using namespace Life::vf;

    /*
     * First we create the mesh
     */
    mesh_ptrtype mesh = createMesh( meshSize );
    Log() << "mesh created in " << t1.elapsed() << "s\n"; t1.restart();

    /*
     * The function space and some associate elements are then defined
     */
    space_ptrtype Xh = space_type::New( mesh );
    element_type u( Xh, "u" );
    element_type v( Xh, "v" );
    Log() << "[functionspace] Number of dof " << Xh->nLocalDof() << "\n";
    Log() << "function space and elements created in " << t1.elapsed() << "s\n"; t1.restart();

    exact_space_ptrtype Eh = exact_space_type::New( mesh );
    exact_element_type fproj( Eh, "f" );
    exact_element_type gproj( Eh, "g" );
    Log() << "[functionspace] Number of dof " << Eh->nLocalDof() << "\n";
    Log() << "function space and elements created in " << t1.elapsed() << "s\n"; t1.restart();


    value_type nu = this->vm()["nu"].template as<double>();
    value_type beta = this->vm()["beta"].template as<double>();


    value_type pi = M_PI;
    AUTO( g, sin(pi*Px())*cos(pi*Py())*cos(pi*Pz()) );
    AUTO( f, (pi*pi*Dim*nu+beta)*g );
    AUTO( zf, 0*Px()+0*Py() );

    int tag1,tag2;
    if ( ( Dim == 1 ) || ( Dim == 2 ) )
        {
            tag1 = 1;
            tag2 = 3;
        }
    else if ( Dim == 3 )
        {
            tag1 = 15;
            tag2 = 23;
        }

    fproj = vf::project( Eh, elements(mesh), f );
    gproj = vf::project( Eh, elements(mesh), g );

    // Construction of the right hand side

    vector_ptrtype F( backend->newVector( Xh ) );


    form1( _test=Xh, _vector=F, _init=true ) =
        integrate( elements(mesh), idv(fproj)*id(v) );
    if ( M_use_weak_dirichlet )
        {
            form1( Xh, F ) +=
                integrate( markedfaces(mesh,tag1),
                           zf*(-nu*grad(v)*N()+M_gammabc*id(v)/hFace() ) );
            form1( Xh, F ) +=
                integrate( markedfaces(mesh,tag2),
                           zf*(-nu*grad(v)*N()+M_gammabc*id(v)/hFace() ) );

        }

    F->close();
    Log() << "F assembled in " << t1.elapsed() << "s\n"; t1.restart();

    //Construction of the left hand side

    sparse_matrix_ptrtype D( backend->newMatrix( Xh, Xh ) );


    form2( Xh, Xh, D, _init=true );
    Log() << "D initialized in " << t1.elapsed() << "s\n";t1.restart();

    form2( Xh, Xh, D ) +=
        integrate( elements(mesh),
                   nu*(gradt(u)*trans(grad(v)))
                   + beta*(idt(u)*id(v)) );
    Log() << "D stiffness+mass assembled in " << t1.elapsed() << "s\n";t1.restart();
    if ( M_use_weak_dirichlet )
        {

            form2( Xh, Xh, D ) += integrate( markedfaces(mesh,tag1),
                                             ( - nu*trans(id(v))*(gradt(u)*N())
                                               - nu*trans(idt(u))*(grad(v)*N())
                                               + M_gammabc*trans(idt(u))*id(v)/hFace()) );
            form2( Xh, Xh, D ) += integrate( markedfaces(mesh,tag2),
                                             ( - nu*trans(id(v))*(gradt(u)*N())
                                               - nu*trans(idt(u))*(grad(v)*N())
                                               + M_gammabc*trans(idt(u))*id(v)/hFace()) );
            Log() << "D weak bc assembled in " << t1.elapsed() << "s\n";t1.restart();

        }

    D->close();

    Log() << "D assembled in " << t1.elapsed() << "s\n";


    if ( ! M_use_weak_dirichlet )
        {
            t1.restart();
            form2( Xh, Xh, D ) +=
                on( markedfaces(mesh, tag1), u, F, g )+
                on( markedfaces(mesh, tag2), u, F, g );
            Log() << "Strong Dirichlet assembled in " << t1.elapsed() << "s on faces " << tag1 << " and " << tag2 << " \n";
        }

    t1.restart();

    this->solve( D, u, F );

    Log() << "solve in " << t1.elapsed() << "s\n";
    t1.restart();

    double L2error2 =integrate( elements(mesh),
                                (idv(u)-idv(gproj))*(idv(u)-idv(gproj))).evaluate()( 0, 0 );
    double L2error =   math::sqrt( L2error2 );

    Log() << "||error||_L2=" << L2error << "\n";
    Log() << "L2 norm computed in " << t1.elapsed() << "s\n";
    t1.restart();


    double semiH1error2 =integrate( elements(mesh),
                                    (gradv(u)-gradv(gproj))*trans(gradv(u)-gradv(gproj)) ).evaluate()( 0, 0 ) ;

    Log() << "semi H1 norm computed in " << t1.elapsed() << "s\n";
    t1.restart();

    double H1error =   math::sqrt( semiH1error2+L2error2 );


    Log() << "||error||_H1=" << H1error << "\n";
    Log() << "H1 norm computed in " << t1.elapsed() << "s\n";
    t1.restart();

    this->exportResults( u, v );


    this->addOutputValue( L2error ).addOutputValue( H1error );
    this->postProcessing();

} // Laplacian::run

template<int Dim, int Order, int RDim, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Laplacian<Dim, Order, RDim, Entity>::solve( sparse_matrix_ptrtype& D,
                                            element_type& u,
                                            vector_ptrtype& F )
{


    vector_ptrtype U( backend->newVector( u.functionSpace() ) );
    backend->solve( D, D, U, F );
    u = *U;
} // Laplacian::solve


template<int Dim, int Order, int RDim, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Laplacian<Dim, Order, RDim, Entity>::exportResults( element_type& U, element_type& v )
{
    if ( exporter->doExport() )
        {
            Log() << "exportResults starts\n";

            exporter->step(0)->setMesh( U.functionSpace()->mesh() );

            exporter->step(0)->add( "pid",
                           regionProcess( boost::shared_ptr<p0_space_type>( new p0_space_type( U.functionSpace()->mesh() ) ) ) );
            exporter->step(0)->add( "u", U );
            exporter->step(0)->add( "exact", v );

            exporter->save();
        }
} // Laplacian::export


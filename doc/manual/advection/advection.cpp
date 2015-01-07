/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-11-23

  Copyright (C) 2006-2011 University Joseph Fourier

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
#include <feel/feel.hpp>


/**
 * Create the ring geometry
 */
std::pair<std::string,std::string> createRing( int Dim, double meshSize );


inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description advectionoptions( "Advection options" );
    advectionoptions.add_options()
    ( "penal", Feel::po::value<double>()->default_value( 0.5 ), "penalisation parameter" )
    ( "f", Feel::po::value<double>()->default_value( 0 ), "forcing term" )
    //        ("g", Feel::po::value<double>()->default_value( 0 ), "boundary term")
    ( "bx", Feel::po::value<double>()->default_value( 1.0 ), "convection X component" )
    ( "by", Feel::po::value<double>()->default_value( 1.0 ), "convection Y component" )
    ( "mu", Feel::po::value<double>()->default_value( 1.0 ), "reaction coefficient component" )
    ( "geomap", Feel::po::value<int>()->default_value( 0 ), "type of geomap for integrals" )
    ( "stiff", Feel::po::value<double>()->default_value( 1.0 ), "stiffness parameter of solution" )
    ( "ring", Feel::po::value<bool>()->default_value( 0 ), "0 = square computational domain, 1 = quarter of a ring as computational domain" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.5 ), "first h value to start convergence" )
    ( "bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
    ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
    ( "export-matlab", "export matrices and vectors in matlab format" )

    ;
    return advectionoptions.add( Feel::feel_options() ) ;
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "advection" ,
                           "advection" ,
                           "0.3",
                           "nD(n=1,2,3)Advection equation on simplices or simplex products using cG and dG",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2006-2011 University Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Benjamin Stamm", "developer", "benjamin.stamm@epfl.ch", "" );
    return about;

}


namespace Feel
{
gmsh_ptrtype
createRing( int Dim, int Order, double meshSize, std::string const& convex )
{
    std::ostringstream ostr;
    std::ostringstream nameStr;
    gmsh_ptrtype gmshp( new Gmsh( Dim, Order ) );
    gmshp->setCharacteristicLength( meshSize );
    ostr << gmshp->preamble() << "\n";

    switch ( Dim )
    {
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
        nameStr << "ring-" << convex;
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

    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( ostr.str() );
    return gmshp;
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

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Entity<Dim,Order,Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    template<typename Conti = Cont, template<uint16_type D> class FieldType = Scalar>
    struct space
    {
        typedef bases<Lagrange<Order, FieldType,Conti> > basis_type;
        typedef FunctionSpace<mesh_type, basis_type> type;
    };

    /* export */
    typedef Exporter<mesh_type,Order> export_type;

    Advection()
        :
        super(),
        backend( backend_type::build( soption("backend") ) ),
        meshSize( doption("hsize") ),
        bcCoeff(  doption("bccoeff") ),
        geomap( ( GeomapStrategyType )ioption("geomap") ),
        exporter( export_type::New( this->vm(), "advection" ) )
    {
        LOG(INFO) << "[Advection] hsize = " << meshSize << "\n";
        LOG(INFO) << "[Advection] bccoeff = " << bcCoeff << "\n";
    }

    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    template<typename f1_type, typename f2_type, typename f3_type>
    void exportResults( f1_type& u,
                        f2_type& v,
                        f3_type& e );

private:

    backend_ptrtype backend;
    mesh_ptrtype mesh;

    double meshSize;
    double bcCoeff;

    GeomapStrategyType geomap;

    boost::shared_ptr<export_type> exporter;
}; // Advection

template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Advection<Dim, Order, Cont, Entity>::run()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    using namespace Feel::vf;

    this->changeRepository( boost::format( "%1%/%2%/P%3%/h_%4%/" )
                            % this->about().appName()
                            % entity_type::name()
                            % Order
                            % meshSize
                          );
    value_type penalisation = this->vm()["penal"].template as<value_type>();
    int bctype =       ioption("bctype");
    double beta_x =    this->vm()["bx"].template as<value_type>();
    double beta_y =    this->vm()["by"].template as<value_type>();
    value_type mu =    this->vm()["mu"].template as<value_type>();
    value_type stiff = this->vm()["stiff"].template as<value_type>();
    bool ring =        boption("ring");

    /*
     * First we create the mesh
     */

    if ( ring )
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                               _desc=createRing( Dim,Order,meshSize,entity_type::name() ) );

    else
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                               _desc=domain( _name=( boost::format( "hypercube-%1%" ) % Dim ).str() ,
                                             _order=Order,
                                             _shape="hypercube",
                                             _dim=Dim,
                                             _h=meshSize ) );


    /*
     * The function space and some associate elements are then defined
     */
    auto Xh = space<Cont>::type::New( mesh );
    auto u = Xh->element();
    auto v = Xh->element();


    auto r = sqrt( Px()*Px()+Py()*Py() );
    //auto r = sqrt(trans(P())*P());
    auto beta = beta_x*( ring*Py()/r + !ring )*oneX()+beta_y*( -ring*Px()/r + !ring )*oneY();
    //auto beta = (ones<Dim,1>());
    auto beta_N = ( trans( N() )*beta );
    auto beta_abs = abs( beta_N );
    auto beta_minus = constant( 0.5 )*( beta_abs-beta_N );
    auto beta_plus = constant( 0.5 )*( beta_abs+beta_N );
    auto g = ( ( !ring*exp( -mu*Px() )*atan( ( Py()-	0.5 )/stiff ) + ring*exp( -mu*r*acos( Py()/r ) )*atan( ( r-0.5 )/stiff ) ) );
    auto f = ( constant( 0.0 ) );

    auto F = backend->newVector( Xh );
    form1( _test=Xh, _vector=F, _init=true )  =
        integrate( elements( mesh ), trans( f )*id( v ) );

    if ( bctype == 1 || !Cont::is_continuous )
    {
        form1( _test=Xh, _vector=F ) +=
            integrate( _range=boundaryfaces( mesh ), _expr=trans( beta_minus*g )*id( v ),_geomap=geomap );
    }

    auto D = backend->newMatrix( Xh, Xh );
    //size_type pattern = Pattern::COUPLED|Pattern::EXTENDED;
    size_type pattern = Pattern::COUPLED;
    form2( _test=Xh, _trial=Xh, _matrix=D, _init=true, _pattern=pattern ) =
        integrate( _range=elements( mesh ), _quad=_Q<2*Order>(),
                   // -(u,beta*grad(v))+(mu*u,v)-(u,div(beta)*v)
                   _expr=-trans( idt( u ) )*( grad( v )*beta ) + mu*trans( idt( u ) )*id( v ),
                   //- idt(u)*id(v)*(dx(beta_x)+dy(beta_y))
                   _geomap=geomap
                 );

    if ( !Cont::is_continuous )
    {
        form2( _test=Xh, _trial=Xh, _matrix=D ) +=
            integrate( _range=internalfaces( mesh ), _quad=_Q<2*Order>(),
                       // {beta u} . [v]
                       //( trans(averaget(trans(beta)*idt(u))) * jump(trans(id(v))) )
                       _expr=( averaget( trans( beta )*idt( u ) ) * jump( id( v ) ) )
                             // penal*[u] . [v]
                             + penalisation*beta_abs*( trans( jumpt( trans( idt( u ) ) ) )*jump( trans( id( v ) ) ) ),
                       _geomap=geomap );
    }

    else // continuous case: stabilization by interior penalty
    {
        form2( _test=Xh, _trial=Xh, _matrix=D ) +=
            integrate( internalfaces( mesh ),
                       // penal*[grad(u)] . [grad(v)]
                       + penalisation*beta_abs*hFace()*hFace()*( trans( jumpt( gradt( u ) ) )*jump( grad( v ) ) ) );
    }

    if ( bctype == 1 || !Cont::is_continuous )
    {
        form2( _test=Xh, _trial=Xh, _matrix=D ) +=
            integrate( _range=boundaryfaces( mesh ), _expr=beta_plus*trans( idt( u ) )*id( v ),_geomap=geomap );

        D->close();
        F->close();
    }

    else if ( bctype == 0 )
    {
        D->close();
        F->close();
        form2( _test=Xh, _trial=Xh, _matrix=D ) +=
            on( boundaryfaces( mesh ), u, F, g );
    }

    if ( this->vm().count( "export-matlab" ) )
    {
        F->printMatlab( "F.m" );
        D->printMatlab( "D.m" );
    }

    backend->solve( _matrix=D, _solution=u, _rhs=F );

    double c = integrate( internalfaces( mesh ), trans( jumpv( idv( u ) ) )*jumpv( idv( u ) )  ).evaluate()( 0, 0 );
    double error = integrate( elements( mesh ), trans( idv( u )-g )*( idv( u )-g ) ).evaluate()( 0, 0 );

    LOG(INFO) << "||error||_0 =" << error << "\n";
    std::cout << "c =" << c << "\n";
    std::cout << "||error||_0 =" << error << "\n";

    this->exportResults( u, g, beta );
} // Advection::run

template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity>
template<typename f1_type, typename f2_type, typename f3_type>
void
Advection<Dim, Order, Cont, Entity>::exportResults( f1_type& U,
        f2_type& E,
        f3_type& beta )
{
    auto Xch = space<Continuous>::type::New( mesh );
    auto uEx = Xch->element();
    auto uC = Xch->element();
    auto Xvch = space<Continuous,Vectorial>::type::New( mesh );
    auto betaC = Xvch->element();

    auto L2Proj = projector( Xch, Xch );
    uEx = L2Proj->project( E );
    uC = L2Proj->project( vf::idv( U ) );
    auto L2Projv = projector( Xvch, Xvch );
    betaC = L2Projv->project( ( beta ) );

    exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );
    exporter->step( 0 )->add( "u", U );
    exporter->step( 0 )->add( "uC", uC );
    exporter->step( 0 )->add( "exact", uEx );
    exporter->step( 0 )->add( "beta", betaC );
    exporter->save();

} // Advection::export
} // Feel




int
main( int argc, char** argv )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,_desc=makeOptions(),_about=makeAbout() );

    /* change parameters below */
    const int nDim = 2;
    const int nOrder = 2;
    //typedef Continuous MyContinuity;
    typedef Discontinuous MyContinuity;
    typedef Feel::Advection<nDim, nOrder, MyContinuity, Simplex> advection_type;

    /* define and run application */
    advection_type advection;
    advection.run();
}





/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
 */
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>
/** use Feel namespace */
using namespace Feel;
using namespace Feel::vf;


inline
po::options_description
makeOptions()
{
    po::options_description gridoptions( "Grid options" );
    gridoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
    ( "kappa", po::value<double>()->default_value( 1 ), "coefficient" )
    ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)" )
    ( "nu", po::value<double>()->default_value( 1 ), "grad.grad coefficient" )
    ;
    return gridoptions.add( Feel::feel_options() );
}

template<int Dim>
class Grid
    :
public Simget
{
    typedef Simget super;
public:
    static const uint16_type Order = 2;
    typedef double value_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::vector_type vector_type;
    typedef Simplex<Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef bases<Lagrange<Order,Scalar> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef Exporter<mesh_type> export_type;

    /**
     * Constructor
     */
    Grid()
        :
        super(),
        M_backend( backend_type::build( soption("backend") ) ),
        meshSize( doption("hsize") ),
        shape( soption("shape") ),
        eigen( SolverEigen<value_type>::build( this->vm() ) )
    {
    }

    void run();
private:

    backend_ptrtype M_backend;
    double meshSize;
    std::string shape;
    std::vector<int> flags;
    boost::shared_ptr<SolverEigen<value_type> > eigen;
}; // Grid

template<int Dim> const uint16_type Grid<Dim>::Order;

template<int Dim>
void
Grid<Dim>::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute Grid<" << Dim << ">\n";
    Environment::changeRepository( boost::format( "doc/manual/%1%/%2%-%3%/P%4%/h_%5%/" )
                                   % this->about().appName()
                                   % shape
                                   % Dim
                                   % Order
                                   % meshSize );

#if 0
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 1,1 );
    GeoTool::Rectangle Omega( meshSize,"Omega",x1,x2 );
    Omega.setMarker( _type="line",_name="Paroi",_markerAll=true );
    Omega.setMarker( _type="surface",_name="Omega",_markerAll=true );

    auto mesh = Omega.createMesh<mesh_type>( "omega_"+ mesh_type::shape_type::name() );
#endif

    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                        _addmidpoint=false,
                                        _usenames=false,
                                        _shape=shape,
                                        _dim=Dim,
                                        _h=meshSize,
                                        _xmin=0.,
                                        _xmax=1.,
                                        _ymin=0.,
                                        _ymax=1.,
                                        _zmin=0.,
                                        _zmax=1. ) );


    auto Xh = space_type::New( mesh );
    auto u = Xh->element();
    auto v = Xh->element();


    if ( Dim == 1 )
    {
        using namespace boost::assign;
        flags += 1,3;
    }

    else if ( Dim == 2 )
    {
        using namespace boost::assign;
        flags += 1,2,3,4;
    }

    else if ( Dim == 3 )
    {
        using namespace boost::assign;
        flags += 6,15,19,23,27,28;
    }

    value_type kappa = doption("kappa");
    value_type nu = doption("nu");

    auto A = M_backend->newMatrix( Xh, Xh ) ;
    form2( _test=Xh, _trial=Xh, _matrix=A ) =
        integrate( elements( mesh ), kappa*gradt( u )*trans( grad( v ) ) + nu*idt( u )*id( v ) );

    auto B = M_backend->newMatrix( Xh, Xh ) ;
    form2( _test=Xh, _trial=Xh, _matrix=B );
    BOOST_FOREACH( int marker, flags )
    {
        form2( Xh, Xh, _matrix=B ) +=
            integrate( markedfaces( mesh,marker ), kappa*idt( u )*id( v ) );
    }

    int maxit = ioption("solvereigen-maxiter");
    int tol = doption("solvereigen-tol");

    int nev = ioption("solvereigen-nev");

    int ncv = ioption("solvereigen-ncv");;

    double eigen_real, eigen_imag;

    SolverEigen<double>::eigenmodes_type modes;

    std::cout << "nev= " << nev <<std::endl;
    std::cout << "ncv= " << ncv <<std::endl;

    modes=
        eigs( _matrixA=A,
              _matrixB=B,
              _nev=nev,
              _ncv=ncv,
              _transform=SINVERT,
              _spectrum=SMALLEST_MAGNITUDE,
              _verbose = true );

    auto femodes = std::vector<decltype( Xh->element() )>( modes.size(), Xh->element() );

    if ( !modes.empty() )
    {
        LOG(INFO) << "eigenvalue " << 0 << " = (" << modes.begin()->second.get<0>() << "," <<  modes.begin()->second.get<1>() << ")\n";

        int i = 0;
        BOOST_FOREACH( auto mode, modes )
        {
            std::cout << " -- eigenvalue " << i << " = (" << mode.second.get<0>() << "," <<  mode.second.get<1>() << ")\n";
            femodes[i++] = *mode.second.get<2>();
        }
    }

    auto exporter =  export_type::New( this->vm(),
                                       ( boost::format( "%1%-%2%-%3%" )
                                         % this->about().appName()
                                         % shape
                                         % Dim ).str() ) ;

    if ( exporter->doExport() )
    {
        LOG(INFO) << "exportResults starts\n";

        exporter->step( 0 )->setMesh( mesh );

        int i = 0;
        BOOST_FOREACH( auto mode, femodes )
        {
            exporter->step( 0 )->add( ( boost::format( "mode-%1%" ) % i++ ).str(), mode );
        }

        exporter->save();
        LOG(INFO) << "exportResults done\n";
    }
}

int
main( int argc, char** argv )
{
    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="mortar",
                                  _author="Abdoulaye Samake",
                                  _email="samakeablo@gmail.com") );

    Application app;

    app.add( new Grid<2>() );
    app.run();
}






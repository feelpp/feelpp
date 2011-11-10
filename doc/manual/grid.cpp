/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
 */
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelfilters/gmsh.hpp>
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
    po::options_description gridoptions("Grid options");
    gridoptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.5 ), "mesh size")
        ("kappa", po::value<double>()->default_value( 1 ), "coefficient")
        ("shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)")
        ("nu", po::value<double>()->default_value( 1 ), "grad.grad coefficient")
        ;
    return gridoptions.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "grid" ,
                     "grid" ,
                     "0.2",
                     "nD(n=1,2,3) Grid on simplices or simplex products",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2008-2009 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

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
    Grid( po::variables_map const& vm, AboutData const& about )
        :
        super( vm, about ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        eigen( SolverEigen<value_type>::build( this->vm() ) ),
        shape( this->vm()["shape"].template as<std::string>() )
    {
    }

    void run();

    void run( const double* X, unsigned long P, double* Y, unsigned long N );

private:

    backend_ptrtype M_backend;
    double meshSize;
    std::string shape;
    boost::shared_ptr<SolverEigen<value_type> > eigen;
}; // Grid

template<int Dim> const uint16_type Grid<Dim>::Order;

template<int Dim>
void
Grid<Dim>::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute Grid<" << Dim << ">\n";
    std::vector<double> X( 2 );
    X[0] = meshSize;
    if ( shape == "hypercube" )
        X[1] = 1;
    else // default is simplex
        X[1] = 0;
    std::vector<double> Y( 3 );
    run( X.data(), X.size(), Y.data(), Y.size() );
}
template<int Dim>
void
Grid<Dim>::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    if ( X[1] == 0 ) shape = "simplex";
    if ( X[1] == 1 ) shape = "hypercube";

    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "doc/manual/%1%/%2%-%3%/P%4%/h_%5%/" )
                                       % this->about().appName()
                                       % shape
                                       % Dim
                                       % Order
                                       % meshSize );

    GeoTool::Node x1(0,0);
    GeoTool::Node x2(1,1);
    GeoTool::Rectangle Omega(meshSize,"Omega",x1,x2);
    Omega.setMarker(_type="line",_name="Paroi",_markerAll=true);
    Omega.setMarker(_type="surface",_name="Omega",_markerAll=true);

    auto mesh = Omega.createMesh<mesh_type>( "omega_"+ mesh_type::shape_type::name() );
#if 0
    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _desc=domain( _name=(boost::format( "%1%-%2%" ) % shape % Dim).str() ,
                                              _usenames=true,
                                              _shape=shape,
                                              _dim=Dim,
                                              _h=X[0],
                                              _xmin=0.,
                                              _xmax=1.,
                                              _ymin=0.,
                                              _ymax=1. ) );
#endif
    auto meshb = createSubmesh(mesh,boundaryelements(mesh));
    auto Xh = space_type::New( mesh );
    auto u = Xh->element();
    auto v = Xh->element();

    auto Xhb = space_type::New( meshb );
    auto ub = Xhb->element();
    auto vb = Xhb->element();

    value_type kappa = this->vm()["kappa"].template as<double>();
    value_type nu = this->vm()["nu"].template as<double>();

    using namespace Feel::vf;

    auto A = M_backend->newMatrix( Xhb, Xhb ) ;
    form2( Xhb, Xhb, A, _init=true ) =
        integrate( elements(meshb), kappa*gradt(u)*trans(grad(v)) + nu*idt(u)*id(v) );

    auto B = M_backend->newMatrix( Xhb, Xhb ) ;
    form2( Xhb, Xhb, B, _init=true ) =
        integrate( markedfaces(mesh,"Paroi"), kappa*idt(u)*id(v) );


    int maxit = this->vm()["solvereigen-maxiter"].template as<int>();
    int tol = this->vm()["solvereigen-tol"].template as<double>();

    int nev = this->vm()["solvereigen-nev"].template as<int>();

    int ncv = this->vm()["solvereigen-ncv"].template as<int>();;

    double eigen_real, eigen_imag;

    SolverEigen<double>::eigenmodes_type modes;


    modes=
        eigs( _matrixA=A,
              _matrixB=B,
              _nev=nev,
              _ncv=ncv,
              _spectrum=SMALLEST_MAGNITUDE );
    std::cout <<"pass"<< std::endl;

    auto mode = Xhb->element() ;

    if ( !modes.empty() )
    {
        Log() << "eigenvalue " << 0 << " = (" << modes.begin()->second.get<0>() << "," <<  modes.begin()->second.get<1>() << ")\n";
        std::cout << "eigenvalue " << 0 << " = (" << modes.begin()->second.get<0>()
                  << "," <<  modes.begin()->second.get<1>() << ")\n";

        mode = *modes.begin()->second.get<2>();
    }

    auto exporter =  export_type::New( this->vm(),
                                       (boost::format( "%1%-%2%-%3%" )
                                        % this->about().appName()
                                        % shape
                                        % Dim).str() ) ;

    if ( exporter->doExport() )
    {
        Log() << "exportResults starts\n";

        exporter->step(0)->setMesh( meshb );

        //exporter->step(0)->add( "u", u );
        exporter->step(0)->add( "mode", mode );

        exporter->save();
        Log() << "exportResults done\n";
    }
}


int
main( int argc, char** argv )
{
    Environment env( argc, argv );

    Application app( argc, argv, makeAbout(), makeOptions() );
    if ( app.vm().count( "help" ) )
    {
        std::cout << app.optionsDescription() << "\n";
        return 0;
    }

    app.add( new Grid<2>( app.vm(), app.about() ) );
    app.run();
}






/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-07

  Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file Lshape.cpp
   \author Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   \date 2012-05-07
 */
#include <feel/feel.hpp>
#include <feel/feeldiscr/elementdiv.hpp>
#include <feel/feelmesh/meshadaptation.hpp>

/** use Feel namespace */
using namespace Feel;

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
inline
po::options_description
makeOptions()
{
    po::options_description lShapeoptions( "LShape options" );
    lShapeoptions.add_options()
        ( "nDim", po::value<int>()->default_value( 2 ), "dimension" )
        ( "hsize", po::value<double>()->default_value( 0.5 ), "mesh size" )
        ( "Lx", po::value<double>()->default_value( 2.0 ), "length (x) of Lshape" )
        ( "Ly", po::value<double>()->default_value( 2.0 ), "length (y) of Lshape" )
        ( "Lz", po::value<double>()->default_value( 1.0 ), "length (z) of Lshape" )
        ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)" )
        ( "weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
        ( "penaldir", Feel::po::value<double>()->default_value( 60 ),
          "penalisation parameter for the weak boundary Dirichlet formulation" )
        ("meshadapt_type", Feel::po::value<std::string>()->default_value( "isotropic" ), "type of mesh adaptation (isotropic, anisotropic)" )
        ("tolerance", Feel::po::value<double>()->default_value( 0.5 ), "tolerance parameter mesh adaptation criterion")
        ("geo_tolerance", Feel::po::value<double>()->default_value( 0.5 ), "geometrical tolerance parameter for mesh adaptation")
        ("max_iterations", Feel::po::value<int>()->default_value( 10 ), "max of mesh adaptation iterations")
        ;
    return lShapeoptions.add( Feel::feel_options() );
}

/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
inline
AboutData
makeAbout()
{
    AboutData about( "lShape" ,
                     "lShape" ,
                     "0.2",
                     "nD(n=1,2,3) LShape on simplices or simplex products",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2008-2009 Universite Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}


/**
 * \class LShape
 *
 * LShape Solver using continuous approximation spaces
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 *
 * \tparam Dim the geometric dimension of the problem (e.g. Dim=1, 2 or 3)
 */
template<int Dim>
class LShape
    :
public Simget
{
    typedef Simget super;
public:

    //! Polynomial order \f$P_2\f$
    static const uint16_type Order = 2;

    //! Geometrical order
    static const uint16_type GOrder = 1;

    //! numerical type is double
    typedef double value_type;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order 1
    typedef Simplex<Dim> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //! the basis type of our approximation space
    typedef bases<Lagrange<Order,Scalar> > basis_type;
    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<space_type> space_ptrtype;
    //! an element type of the approximation function space
    typedef typename space_type::element_type element_type;

    //! Scalar P0 space
    typedef bases<Lagrange<0, Scalar, Discontinuous > > p0_basis_type;
    typedef FunctionSpace<mesh_type, p0_basis_type> p0_space_type;
    typedef boost::shared_ptr<p0_space_type> p0_space_ptrtype;
    typedef typename p0_space_type::element_type p0_element_type;

    //! Vectorial P1 space
    typedef bases<Lagrange<1, Vectorial> > p1vec_basis_type;
    typedef FunctionSpace<mesh_type, p1vec_basis_type> p1vec_space_type;
    typedef boost::shared_ptr<p1vec_space_type> p1vec_space_ptrtype;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /**
     * Constructor
     */
    LShape()
        :
        super(),
        M_backend( backend_type::build( soption("backend") ) ),
        meshSize( doption("hsize") ),
        shape( soption("shape") )
    {
    }

    typedef MeshAdaptation<Dim, Order, GOrder> MeshAdapt;

    void run();
    void run( const double* X, unsigned long P, double* Y, unsigned long N );
    gmsh_ptrtype createLShapeGeo( double meshSize, double Lx, double Ly, double Lz);
    boost::tuple<double, double, p0_element_type> zz_estimator(const element_type& U, const mesh_ptrtype& mesh);


private:

    //! backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    double meshSize;

    //! shape of the domain
    std::string shape;
}; // LShape

template<int Dim> const uint16_type LShape<Dim>::Order;

template<int Dim>
gmsh_ptrtype
LShape<Dim>::createLShapeGeo( double meshSize, double Lx, double Ly, double Lz)
{
    std::ostringstream costr;
    costr << "hs=" << meshSize << ";\n"
          << "Lx=" << Lx << ";\n"
          << "Ly=" << Ly << ";\n";

    if (Dim == 3)
        costr << "Lz=" << Lz << ";\n";

    costr << "Point(1) = {0, 0, 0, hs};\n"
          << "Point(2) = {Lx, 0, 0, hs};\n"
          << "Point(3) = {Lx/2, Ly/2, 0, hs};\n"
          << "Point(4) = {Lx, Ly/2, 0, hs};\n"
          << "Point(5) = {0, Ly, 0, hs};\n"
          << "Point(6) = {Lx/2, Ly, 0, hs};\n"
          << "Line(1) = {1, 2};\n"
          << "Line(2) = {2, 4};\n"
          << "Line(3) = {4, 3};\n"
          << "Line(4) = {3, 6};\n"
          << "Line(5) = {6, 5};\n"
          << "Line(6) = {5, 1};\n"
          << "Line Loop(8) = {4, 5, 6, 1, 2, 3};\n"
          << "Plane Surface(8) = {8};\n";

    if (Dim == 2)
        {
            costr << "Physical Line(\"boundary\") = {4, 5, 6, 1, 2, 3};\n"
                  << "Physical Surface(\"omega\") = {8};\n";
        }
    else if (Dim == 3 )
        {
            // Extrude Lshape2D -> Lshape 3D
            costr << "Extrude{0,0,Lz}{ Surface{8}; }\n"
                  << "Physical Surface(\"boundary\") = {19, 23, 27, 31, 35, 39};\n"
                  << "Physical Volume(\"omega\") = {1};\n";
        }

    std::ostringstream nameStr;
    nameStr << "geoLShape";
    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( costr.str() );
    return gmshp;
}


template<int Dim>
boost::tuple<double, double, typename LShape<Dim>::p0_element_type>
LShape<Dim>::zz_estimator(const element_type& U, const mesh_ptrtype& mesh)
{
    using namespace Feel::vf;

    p0_space_ptrtype P0h = p0_space_type::New( mesh );
    p1vec_space_ptrtype P1hvec = p1vec_space_type::New( mesh );

    auto GhUh = div( vf::sum( P1hvec, trans(gradv(U))*vf::meas()), vf::sum( P1hvec, vf::meas()*vf::one()) );
    auto eta_k_U = integrate(elements(mesh), trans(idv(GhUh) - trans(gradv(U)))*(idv(GhUh) - trans(gradv(U))), _Q<10>() ).broken(P0h).sqrt();

    auto h=vf::project(P0h, elements(mesh), vf::h() );

    auto estimatorH1_U = math::sqrt( (eta_k_U.pow(2)).sum() );
    auto estimatorL2_U = math::sqrt( (element_product(eta_k_U,h).pow(2)).sum() );

    return boost::make_tuple( estimatorL2_U, estimatorH1_U, eta_k_U);
}


template<int Dim>
void
LShape<Dim>::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute LShape<" << Dim << ">\n";
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
LShape<Dim>::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    if ( X[1] == 0 ) shape = "simplex";
    if ( X[1] == 1 ) shape = "hypercube";

    Environment::changeRepository( boost::format( "doc/manual/adapt/%1%/%2%-%3%/P%4%/h_%5%/" )
                                   % this->about().appName()
                                   % shape
                                   % Dim
                                   % Order
                                   % meshSize );

    //! Set dimensions of Lshape geometry
    double Lx = doption("Lx");
    double Ly = doption("Ly");
    double Lz = doption("Lz");

    mesh_ptrtype mesh = createGMSHMesh ( _mesh = new mesh_type,
                                         _desc = createLShapeGeo( meshSize, Lx, Ly, Lz ),
                                         _update=MESH_UPDATE_FACES | MESH_UPDATE_EDGES );

    element_type u; //solution
    element_type v; //test function

    // Loop for calculation and mesh adaptation
    int adapt_iter = 0; //current iteration
    double mesh_eps = doption("geo_tolerance"); //geometrical tolerance
    double tol = doption("tolerance"); //tolerance for stop criterion
    int max_iter = ioption("max_iterations");
    std::string meshadapt_type = soption("meshadapt_type");

    // Mesh adaptation stop criterion
    boost::tuple<double, double, p0_element_type> estimator_U;

    double criterion_U; // Relative error
    bool criterion = true;

    MeshAdapt mesh_adaptation;
    do{
        space_ptrtype Xh = space_type::New( mesh );
        u = Xh->element();
        v = Xh->element();

        // print some information (number of local/global dof in logfile)
        Xh->printInfo();

        value_type pi = M_PI;

        //! Resolve -\Delta u = 1 on \Omega, with u = 0 on \delta \Omega
        auto f = cst(1.);
        auto g = cst(0.);

        bool weak_dirichlet = ioption("weakdir");
        value_type penaldir = doption("penaldir");

        //! Define right hand side :
        // \int_{\Omega} fv
        auto F = M_backend->newVector( Xh );
        form1( _test=Xh, _vector=F, _init=true ) =
            integrate( _range=elements( mesh ), _expr=f*id(v));

        //! Weak Dirichlet boundary conditions treatment :
        // \int_{d\Omega} \frac{\gamma}{h} gv - \int_{d\Omega} (\nabla v \cdot n) g
        if ( weak_dirichlet )
            {
                form1( _test=Xh, _vector=F ) +=
                    integrate( _range=markedfaces( mesh,"boundary" ),
                               _expr=g*( -grad( v )*vf::N()+penaldir*id( v )/hFace() ) );
            }

        F->close();

        //! Define left hand side
        auto D = M_backend->newMatrix( _test=Xh, _trial=Xh  );

        // \int_{\Omega} \nabla u \nabla v
        form2( _test=Xh, _trial=Xh, _matrix=D ) =
            integrate( _range=elements( mesh ), _expr=gradt( u )*trans( grad( v ) ) );

        //! Weak Dirichlet boundary conditions treatment :
        // - \int_{d\Omega} (\nabla u \cdot n) v + \int_{d\Omega} \frac{\gamma}{h} uv - \int_{d\Omega} (\nabla v \cdot n) u
        if ( weak_dirichlet )
            {
                form2( _test=Xh, _trial=Xh, _matrix=D ) +=
                    integrate( _range=markedfaces( mesh,"boundary" ),
                               _expr= ( -( gradt( u )*vf::N() )*id( v )
                                        -( grad( v )*vf::N() )*idt( u )
                                        +penaldir*id( v )*idt( u )/hFace() ) );
                D->close();

            }
        else
            {
                // Set u = 0 on \delta \Omega
                D->close();
                form2( _test=Xh, _trial=Xh, _matrix=D ) +=
                    on( _range=markedfaces( mesh, "boundary" ),
                        _element=u, _rhs=F, _expr=g );
            }

        //! Solve the system
        backend_type::build(soption("backend"))->solve( _matrix=D, _solution=u, _rhs=F );

        //! Exportation of results
        export_ptrtype exporter( export_type::New( this->vm(),
                                                   ( boost::format( "%1%-%2%-%3%" )
                                                     % this->about().appName()
                                                     % shape
                                                     % Dim ).str() ) );

        if ( exporter->doExport() )
            {
                LOG(INFO) << "exportResults starts\n";
                exporter->step( 0 )->setMesh( mesh );
                exporter->step( 0 )->add( "u", u );
                exporter->save();
                LOG(INFO) << "exportResults done\n";
            }

        // ****** Mesh adaptation for Potential and Temperature ****** ///

        auto norm2_U = integrate( elements(mesh), idv(u)*idv(u) ).evaluate()(0,0);
        auto norm_U = math::sqrt( norm2_U );

        estimator_U = zz_estimator(u, mesh);
        criterion_U = math::abs( estimator_U.template get<0>() - norm_U  )/norm_U;

        std::cout << "criterion = " << criterion_U << std::endl;
        std::cout << "ZZ estimator = " << estimator_U.template get<0>() << std::endl;

        LOG(INFO) << "criterion = " << criterion_U;
        LOG(INFO) << "ZZ estimator = " << estimator_U.template get<0>();

        criterion = criterion && (criterion_U > (1+tol) || criterion_U < (1-tol) );

        //// Use interface of mesh adaptation
        std::list<std::pair<element_type, std::string> > var_list;
        std::pair<element_type, std::string> solution_pair = std::make_pair( u, "solution");
        var_list.push_back(solution_pair);

        std::string geofile_path = "./geoLShape.geo";

        if (criterion)
            mesh = mesh_adaptation.adaptMesh(_initMesh=mesh, _geofile=geofile_path, _adaptType=meshadapt_type,
                                             _var=var_list, _tol=mesh_eps);

        /// Make a pause between two adaptations : testing purpose
        // cout << "Press any key to continue: " << endl;
        // char c;
        // cin >> c;

        adapt_iter++;

        /// update mesh_eps
        mesh_eps = mesh_eps/5.0;

    }while(criterion && adapt_iter < max_iter);
    // End mesh adaptation loop

} // LShape::run

/**
 * main function: entry point of the program
 */
int
main( int argc, char** argv )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );

    /**
     * create an application
     */
    Application app;

    /**
     * register the simgets
     */
    switch (app.vm()["nDim"].as<int>()) {
    case(2) : {
        app.add( new LShape<2>() );
        break;
    }
    case(3) : {
        app.add( new LShape<3>() );
        break;
    }
    default: {
        std::cerr << "wrong pb dimension - should be either 2 or 3\n";
        return 1;
    }
    }
    //app.add( new LShape<3>() );

    /**
     * run the application
     */
    app.run();
}






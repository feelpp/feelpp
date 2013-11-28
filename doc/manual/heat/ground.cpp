/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2011-07-25

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file ground.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2011-07-25
 */

#include <feel/feel.hpp>

namespace Feel
{
//# marker1 #
inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description groundoptions("ground options");
    groundoptions.add_options()
        // mesh parameters
        ("hsize", Feel::po::value<double>()->default_value( 0.1 ),
         "first h value to start convergence")
        ("D", Feel::po::value<double>()->default_value( 1 ),
         "depth")
        ("W", Feel::po::value<double>()->default_value( 0.5 ),
         "width")

        ("soil.0.k", Feel::po::value<double>()->default_value( 0.2 ),"conductivity")
        ("soil.0.c", Feel::po::value<double>()->default_value( 1 ),"capacity")
        ("soil.0.rho", Feel::po::value<double>()->default_value( 1 ),"density")
        ("soil.1.k", Feel::po::value<double>()->default_value( 0.2 ),"conductivity")
        ("soil.1.c", Feel::po::value<double>()->default_value( 1 ),"capacity")
        ("soil.1.rho", Feel::po::value<double>()->default_value( 1 ),"density")

        ("TA", Feel::po::value<double>()->default_value( 1 ),"TA")
        ("TR", Feel::po::value<double>()->default_value( 0 ),"TR")
        ("frequency", Feel::po::value<double>()->default_value( 2*M_PI ),"w")

        ("gamma", Feel::po::value<double>()->default_value( 20 ),"gamma");


        return groundoptions.add( Feel::feel_options() );
}
//# endmarker1 #

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "ground" ,
                           "ground" ,
                           "0.1",
                           "nD(n=1,2,3) heating and cooling of the ground due to surface temperature",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2011 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;
}

struct Soil
{
    double k,rho,c;
};

/**
 * \class Ground
 * \brief Heating and Cooling of the ground due to surface temperature
 *
 * @author Christophe Prud'homme
 * @see
 */
template<int Dim, int Order = 1>
class Ground
    :
    public Application
{
    typedef Application super;
public:

    // -- TYPEDEFS --
    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef backend_type::sparse_matrix_type sparse_matrix_type;
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_type vector_type;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Simplex<Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpace<mesh_type, bases<Lagrange<0, Scalar> >, Discontinuous > p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

	/* BDF discretization */
	typedef Bdf<space_type>  bdf_type;
	typedef boost::shared_ptr<bdf_type> bdf_ptrtype;

    /* export */
    typedef Exporter<mesh_type> export_type;

	/* constructor */
	Ground();

    /* run the simulation */
    void run();

private:

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( double time, element_type& u );

private:

    backend_ptrtype M_backend;

	/* mesh parameters */
    double meshSize;
    double W,D;
    double k0,c0,rho0;
    double k1,c1,rho1;
    double TR,TA,frequency;
    double gamma;

    mesh_ptrtype mesh;
    space_ptrtype Xh;

    sparse_matrix_ptrtype M;
    vector_ptrtype F;

	bdf_ptrtype M_bdf;

    boost::shared_ptr<export_type> exporter;

}; // Ground class

/* Constructor */
template<int Dim, int Order>
Ground<Dim,Order>::Ground()
    :
    super(),
    M_backend( backend_type::build( this->vm() ) ),
    meshSize( this->vm()["hsize"].template as<double>() ),
    W( this->vm()["W"].template as<double>() ),
    D( this->vm()["D"].template as<double>() ),
    k0( this-> vm()["soil.0.k"].template as<double>() ),
    c0( this-> vm()["soil.0.c"].template as<double>() ),
    rho0( this-> vm()["soil.0.rho"].template as<double>() ),
    k1( this-> vm()["soil.1.k"].template as<double>() ),
    c1( this-> vm()["soil.1.c"].template as<double>() ),
    rho1( this-> vm()["soil.1.rho"].template as<double>() ),
    TR( this-> vm()["TR"].template as<double>() ),
    TA( this-> vm()["TA"].template as<double>() ),
    frequency( this-> vm()["frequency"].template as <double>() ),
    gamma( this-> vm()["gamma"].template as <double>() ),
    exporter( export_type::New( this->vm(), this->about().appName() ) )
{
    this->changeRepository( boost::format( "%1%/%2%/%3%/" )
                            % this->about().appName()
                            % entity_type::name()
                            % this->vm()["hsize"].template as<double>()
        );
    using namespace Feel::vf;


    LOG(INFO) << "meshSize = " << meshSize << "\n"
              << "D = "<<D<<"\n"
              << "W = " << W << "\n"
              << "TR = " << TR << "\n"
              << "TA = " << TA << "\n"
              << "rho0 = " << rho0 << "\n"
              << "k0 = " << k0 << "\n"
              << "c0 = " << c0 << "\n"
              << "rho1 = " << rho1 << "\n"
              << "k1 = " << k1 << "\n"
              << "c1 = " << c1 << "\n"
              << "frequency = " << frequency << "\n";

    /*
     * First we create the mesh
     */
#if 0
    std::cout << "Omega0\n";
    GeoTool::Rectangle Omega0( meshSize,"Omega0",GeoTool::Node(-W/2,-D),GeoTool::Node(W/2,0));
    Omega0.setMarker( _type="surface", _name="soil.0", _markerAll=true );
    Omega0.setMarker( _type="line", _name="ground", _marker3=true );
    std::cout << "Omega1\n";
    GeoTool::Rectangle Omega1( meshSize,"Omega1",GeoTool::Node(-W/4,-D/2),GeoTool::Node(W/4,-D/2+D/4));
    Omega1.setMarker( _type="surface", _name="soil.1", _markerAll=true );

    auto Omega = (Omega0+Omega1).createMesh<mesh_type>( "Ground" );
#else
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=geo( _filename="ground.geo",
                                      _dim=2,
                                      _order=1,
                                      _h=meshSize ) );
#endif
    /*
     * The function space associated to the mesh
     */
    Xh = space_type::New( mesh );
    M_bdf = bdf( _space=Xh, _vm=this->vm(), _name="Temperature" );

    /*
     * Right hand side
     */
    F = M_backend->newVector( Xh );

    /*
     * Left hand side
     */
    M = M_backend->newMatrix( Xh, Xh );

}


template<int Dim, int Order>
void
Ground<Dim, Order>::run()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    using namespace Feel::vf;


    /*
     * T is the unknown, v the test function
     */
    auto T = Xh->element();
    auto v = Xh->element();



    /*
     * Left hand side construction (steady state)
     */
    form2( Xh, Xh, M, _init=true ) = integrate( _range= markedelements(mesh,"soil.0"), _expr= k0*gradt(T)*trans(grad(v)) );
    form2( Xh, Xh, M) += integrate( _range= markedelements(mesh,"soil.1"), _expr= k1*gradt(T)*trans(grad(v)) );
    form2( Xh, Xh, M) +=
        integrate( _range= markedfaces(mesh, "ground"),
                   _expr=k0*(-gradt(T)*N()*id(v)-grad(T)*N()*idt(v)+gamma*idt(T)*id(v)/hFace() ) );
    form2(Xh, Xh, M) +=
        integrate( _range=markedelements(mesh, "soil.0"), _expr=rho0*c0*idt(T)*id(v)*M_bdf->polyDerivCoefficient(0) )
        + integrate( _range=markedelements(mesh, "soil.1"), _expr=rho1*c1*idt(T)*id(v)*M_bdf->polyDerivCoefficient(0) );

    M->close();

    if ( this->vm().count( "export-matlab" ) )
    {
        M->printMatlab( "D.m" );
    }

    /*
     * Left and right hand sides construction (non-steady state) with BDF
     */
    T = vf::project( _space=Xh, _expr=cst(TR) );

    M_bdf->initialize(T);

    std::cout << "The step is : " << M_bdf->timeStep() << "\n"
              << "The initial time is : " << M_bdf->timeInitial() << "\n"
              << "The final time is  : " << M_bdf->timeFinal() << "\n";

    for ( M_bdf->start(); M_bdf->isFinished()==false; M_bdf->next(T) )
    {
        std::cout << "============================================================\n";
        std::cout << "T =" << M_bdf->time() << "s\n";
        auto T0=TR+TA*sin(frequency*cst(M_bdf->time()) );
        // update right hand side with time dependent terms
        auto bdf_poly = M_bdf->polyDeriv();
        form1( _test=Xh, _vector=F) =
            integrate( _range= markedfaces(mesh, "ground"), _expr=T0*k0*(-grad(v)*N()+gamma*id(v)/hFace() ) )+
            integrate( _range=markedelements(mesh, "soil.0"), _expr=rho0*c0*idv(bdf_poly)*id(v)) +
            integrate( _range=markedelements(mesh, "soil.1"), _expr=rho1*c1*idv(bdf_poly)*id(v) );

		M_backend->solve( _matrix=M, _solution=T, _rhs=F );

		this->exportResults( M_bdf->time(), T );

     }

    std::cout << "Resolution ended, export done \n";

} // Ground::run


template<int Dim, int Order>
void
Ground<Dim, Order>::exportResults( double time, element_type& U )
{
    exporter->step(time)->setMesh( U.functionSpace()->mesh() );
    exporter->step(time)->add( "Temperature", U );
    exporter->save();
} // Ground::exportResults

}

int
main( int argc, char** argv )
{
    try {
        using namespace Feel;
        Environment env( _argc=argc, _argv=argv,
                         _desc=makeOptions(),
                         _about=makeAbout() );
        /* Parameters to be changed */
        const int nDim = 2;
        const int nOrder = 1;

        Ground<nDim, nOrder> ground;
        ground.run();
    }
    catch( std::exception const& e )
    {
        std::cout << "Exception caught: " << e.what() << "\n";
    }

}

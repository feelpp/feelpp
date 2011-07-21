/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-06-11

  Copyright (C) 2007-2011 Université Joseph Fourier (Grenoble I)

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
   \file heatsink.cpp
   \author Baptiste Morin <baptistemorin@gmail.com>
   \date 2011-06-28
 */
#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feelalg/backendgmm.hpp>
#include <feel/feelalg/backendpetsc.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelpoly/polynomialset.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/bdf2.hpp>

Feel::gmsh_ptrtype makefin( double hsize, double deep );

//# marker1 #
inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description heatsinkoptions("heatsink options");
    heatsinkoptions.add_options()
    // mesh parameters
    ("hsize", Feel::po::value<double>()->default_value( 0.5 ),
     "first h value to start convergence")

	// 3D parameter
	("deep", Feel::po::value<double>()->default_value( 0 ),
	 "depth of the mesh, only in 3D simulation")

    // thermal conductivities parameters
    ("lambda_s", Feel::po::value<double>()->default_value( 386 ),
     "thermal conductivity of the base spreader")
    ("lambda_f", Feel::po::value<double>()->default_value( 386 ),
     "thermal conductivity of the fin")

	// density parameter
	("rho_s", Feel::po::value<int>()->default_value( 8940 ),
	 "density of the spreader's material in SI unit kg.m^{-3}")
	("rho_f", Feel::po::value<int>()->default_value( 8940 ),
	 "density of the fin's material in SI unit kg.m^{-3}")

	// physical coeff
    ("Bimin", Feel::po::value<double>()->default_value( 0.01 ),
     "minimum value of Biot number")
    ("Bimax", Feel::po::value<double>()->default_value( 1 ),
     "maximum value of Biot number")
    ("therm_coeff", Feel::po::value<double>()->default_value(60000),
     "thermal coefficient for the biot number")

    // export
    ("export", "export results(ensight, data file(1D)")
    ("export-matlab", "export matrix and vectors in matlab" );

    return heatsinkoptions.add( Feel::feel_options() );
}
//# endmarker1 #

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "heatsink" ,
                            "heatsink" ,
                            "0.1",
                            "nD(n=1,2,3) Heat sink thermal fin on simplices or simplex products",
                            Feel::AboutData::License_GPL,
                            "Copyright (c) 2006-2011 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    about.addAuthor("Baptiste Morin", "junior developer", "baptistemorin@gmail.com","");
    return about;
}


namespace Feel
{
/**
 * Heat sink application
 */
template<int Dim, int Order>
class HeatSink
    :
    public Application
{
    typedef Application super;
public:

#define Entity Simplex

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
    typedef Entity<Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpace<mesh_type, fusion::vector<Lagrange<0, Scalar> >, Discontinuous > p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;

    /*basis*/
    typedef fusion::vector<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

	/* BDF discretization */
	//typedef Bdf<space_type>  bdf_type;
	//typedef boost::shared_ptr<bdf_type> bdf_ptrtype;

    /* export */
    typedef Exporter<mesh_type> export_type;

    /**
     * create the mesh
     */
    mesh_ptrtype createMesh();

	/* constructor */
	HeatSink( int argc, char** argv, AboutData const& ad, po::options_description const& od );

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
    double depth;

	/* thermal conductivities */
	double lambda_s;
	double lambda_f;

	/* density of the material */
	int rho_s;
	int rho_f;

	/* Biot number */
	double Bi;
    double therm_coeff;

    /* characteristical length of the body */
	double charact_length;

    /* messh, pointers and spaces*/
    mesh_ptrtype mesh;
    space_ptrtype Xh;

    sparse_matrix_ptrtype D;
    vector_ptrtype F;

	//bdf_ptrtype M_bdf;

    boost::shared_ptr<export_type> exporter;

}; // HeatSink class

/* Constructor */
template<int Dim, int Order>
HeatSink<Dim,Order>::HeatSink( int argc, char** argv, AboutData const& ad, po::options_description const& od )
		:
		super( argc, argv, ad, od ),
		M_backend( backend_type::build( this->vm() ) ),
		meshSize( this->vm()["hsize"].template as<double>() ),
		depth( this->vm()["deep"].template as<double>() ),
		lambda_s( this-> vm()["lambda_s"].template as<double>() ),
		lambda_f( this-> vm()["lambda_f"].template as<double>() ),
		rho_s( this-> vm()["rho_s"].template as<int>() ),
		rho_f( this-> vm()["rho_f"].template as<int>() ),
        therm_coeff( this-> vm()["therm_coeff"].template as <double>() ),
		exporter( export_type::New( this->vm(), this->about().appName() ) )
		{
			this->changeRepository( boost::format( "%1%/%2%/%3%/" )
								   % this->about().appName()
								   % entity_type::name()
								   % this->vm()["hsize"].template as<double>()
								  );
		using namespace Feel::vf;

	    /*
		 * First we create the mesh
		 */
        mesh = createMesh();

		/*
		 * Calculate the characterisical length of the mesh = length of the base
		 */
		charact_length = integrate( _range= markedelements(mesh, "gamma4"), _expr= cst(1.) ).evaluate()(0,0);

		/*
		 * calculate the biot number
		 */
		Bi = therm_coeff * charact_length / lambda_f;

        /*
         * The function space associated to the mesh
         */
        Xh = space_type::New( mesh );
        //M_bdf = bdf_ptrtype( new bdf_type( Xh ) );

        /*
         * Right hand side
         */
		F = M_backend->newVector( Xh );

		/*
		 * Left hand side
		 */
		D = M_backend->newMatrix( Xh, Xh );

    }

template<int Dim, int Order>
typename HeatSink<Dim,Order>::mesh_ptrtype
HeatSink<Dim,Order>::createMesh( )
{
	mesh_ptrtype mesh = createGMSHMesh ( _mesh = new mesh_type,
								         _desc = makefin( meshSize, depth ),
                                         _h = meshSize,
                                         _update=MESH_UPDATE_FACES | MESH_UPDATE_EDGES );
    return mesh;
} // HeatSink::createMesh


template<int Dim, int Order>
void
HeatSink<Dim, Order>::run()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    using namespace Feel::vf;

    Log() << "Bi = " << Bi << "\n"
          << "meshSize = " << meshSize << "\n"
          << "depth = " << depth << "\n"
          << "lambda_spreader = " << lambda_s << "\n"
          << "lambda_fin = " << lambda_f << "\n"
          << "rho_spreader = " << rho_s << "\n"
          << "rho_fin = " << rho_f << "\n"
          << "thermal coefficient = " << therm_coeff << "\n";

    /*
     * T is the unknown, v the test function
     */
    element_type T( Xh, "T" );
    element_type v( Xh, "v" );

    /*
     * Right hand side construction (steady state)
     */
    form1( Xh, F, _init=true ) = integrate( _range= markedfaces(mesh,"gamma4"), _expr= id(v) );

    /*
     * Left hand side construction (steady state)
     */
    form2( Xh, Xh, D, _init=true ) = integrate( _range= markedelements(mesh,"spreader_mesh"), _expr= lambda_s*gradt(T)*trans(grad(v)) );
    form2( Xh, Xh, D) += integrate( _range= markedelements(mesh,"fin_mesh"), _expr= lambda_f*gradt(T)*trans(grad(v)) );
    form2 (Xh, Xh, D) += integrate( _range= markedfaces(mesh, "gamma1"), _expr= lambda_f*Bi*idt(T)*id(v));


    //form2(Xh, Xh, D) +=
    //    integrate( _range=markedelements(mesh, "spreader_mesh"), _expr=rho_s*idt(T)*id(v)*M_bdf->polyDerivCoefficient(0) )
    //    + integrate( _range=markedelements(mesh, "fin_mesh"), _expr=rho_f*idt(T)*id(v)*M_bdf->polyDerivCoefficient(0) );

    D->close();

    if ( this->vm().count( "export-matlab" ) )
    {
        D->printMatlab( "D.m" );
    }

    /*
     * Left and right hand sides construction (non-steady state) with BDF
     *//*
    for ( M_bdf->start(); M_bdf->isFinished(); M_bdf->next() )
    {
        auto Ft = M_backend->newVector( Xh );
        auto bdf_poly = M_bdf->polyDeriv();
        form1( _test=Xh, _vector=Ft) =
            integrate( _range=markedelements(mesh, "spreader_mesh"), _expr=rho_s*idv(bdf_poly)*id(v) )+
            integrate( _range=markedelements(mesh, "fin_mesh"), _expr=rho_f*idv(bdf_poly)*id(v) );


        std::cout << "Begin of the resolution \n";
		M_backend->solve( _matrix=D, _solution=T, _rhs=F );
        std::cout << "Resolution ended \n";

		this->exportResults( M_bdf->time(), T );
     }*/

} // HeatSink::run


template<int Dim, int Order>
void
HeatSink<Dim, Order>::exportResults( double time, element_type& U )
{
    if ( this->vm().count( "export" ) )
        {
            exporter->step(1.)->setMesh( U.functionSpace()->mesh() );
            exporter->step(1.)->add( "Temperature", U );
            //exporter->step(time)->add( "Temperature", U);
            exporter->save();

        }
} // HeatSink::exportResults

} // Feel


int
main( int argc, char** argv )
{
    using namespace Feel;

	/* Parameters to be changed */
	const int nDim = 2;
	const int nOrder = 1;

	/* define application */
	typedef Feel::HeatSink<nDim, nOrder> heat_sink_type;

	/* instanciate */
	heat_sink_type heatsink( argc, argv, makeAbout(), makeOptions() );

	/* run */
    heatsink.run();
}






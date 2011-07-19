
/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-06-11

  Copyright (C) 2007-2008 Université Joseph Fourier (Grenoble I)

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
	"depth of the fin, only in 3D simulation")

        // thermal conductivities parameters
        ("spreader", Feel::po::value<double>()->default_value( 386 ),
         "thermal conductivity of the base spreader")
        ("fin", Feel::po::value<double>()->default_value( 386 ),
         "thermal conductivity of the fin")
	
	// density parameter
	("rho_s", Feel::po::value<int>()->default_value( 8940 ),
	"density of the spreader's material in SI unit kg.m^{-3}")
	("rho_f", Feel::po::value<int>()->default_value( 8940 ),
	 "density of the fin's material in SI unit kg.m^{-3}")
		 
		// physical coeff
        ("Bimin", Feel::po::value<double>()->default_value( 0.01 ), "minimum value of Biot number")
        ("Bimax", Feel::po::value<double>()->default_value( 1 ), "maximum value of Biot number")

        ("N", Feel::po::value<int>()->default_value( 1 ), "number of samples withing parameter space")

        // export
        ("export", "export results(ensight, data file(1D)")
        ("export-matlab", "export matrix and vectors in matlab" )
        ;
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
    static const uint16_type imOrder = Dim*Order;

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
	typedef Bdf<space_type>  bdf_type;
	typedef boost::shared_ptr<bdf_type> bdf_ptrtype;

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
     * solve the system
     */
    void solve( sparse_matrix_ptrtype& D, element_type& u, vector_ptrtype& F );


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u );

private:

    backend_ptrtype M_backend;

	/* mesh parameters */
    double meshSize;
    double depth;

	/* thermal conductivities */
	double lambda_spreader;
	double lambda_fin;
	double ratio;
	
	/* density of the material */
	int density_s;
	int density_f;
	
	/* Biot number */
	double Bi;
	
	double charact_length;
	
    mesh_ptrtype mesh;
    space_ptrtype Xh;

    sparse_matrix_ptrtype Dcst;
    vector_ptrtype Fcst;

	bdf_ptrtype M_bdf;
	
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
		lambda_spreader( this-> vm()["spreader"].template as<double>() ),
		lambda_fin( this-> vm()["fin"].template as<double>() ),
		density_s( this-> vm()["rho_s"].template as<int>() ),
		density_f( this-> vm()["rho_f"].template as<int>() ),
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
		 * Calculate the characterisical length of the mesh = surface / perimeter
		 */
		double surface = integrate( markedelements(mesh,"spreader_mesh"), cst(1.) ).evaluate()(0,0) +
			integrate( markedelements(mesh,"fin_mesh"), cst(1.) ).evaluate()(0,0);
	
			
		double perimeter = integrate( markedfaces(mesh,"gamma4"), cst(1.) ).evaluate()(0,0) +
						integrate( markedfaces(mesh,"gamma7"), cst(1.) ).evaluate()(0,0) +
						integrate( markedfaces(mesh,"gamma5"), cst(1.) ).evaluate()(0,0) +
						integrate( markedfaces(mesh,"gamma3"), cst(1.) ).evaluate()(0,0) +
						integrate( markedfaces(mesh,"gamma6"), cst(1.) ).evaluate()(0,0) +
						integrate( markedfaces(mesh,"gamma1"), cst(1.) ).evaluate()(0,0) +
						integrate( markedfaces(mesh,"gamma2"), cst(1.) ).evaluate()(0,0);
			
		charact_length = surface/perimeter;
			
		/*
		 * N.B : you can also calculate the surface with a p0 space.
		 * To do it, you have to name a p0 space ptrtype (such as p0h)
		 * and then P0h = p0_space_type::New( mesh );
		 * The value is obtain in the field 'Value' of the string return of
		 *        auto surface2 = vf::sum( P0h, meas() );
		 */
			
			
		/*
		 * calculate the biot number
		 */
		Bi = cst(60000.) * charact_length / lambda_fin;

		/*
		 * calculate the ratio for the equation
		 */
		ratio = lambda_spreader/lambda_fin;

        /*
         * The function space and some associate elements are then defined
         */
        Xh = space_type::New( mesh );
        element_type u( Xh, "u" );
        element_type v( Xh, "v" );
			
		/*
		 * Construction of the right hand side (the linear form)
		 */	
        Fcst = M_backend->newVector( Xh );
		form1( Xh, Fcst, _init=true ) = integrate( markedfaces(mesh,"gamma4"), cst(1.) );

        if ( this->vm().count( "export-matlab" ) )
            Fcst->printMatlab( "F.m" );

        /*
         * Construction of the left hand side
         */
        Dcst = M_backend->newMatrix( Xh, Xh );
		form2( Xh, Xh, Dcst, _init=true ) = integrate( markedelements(mesh,"spreader_mesh"), _Q<imOrder>(), ( ratio*gradt(T)*trans(grad(v))) )
										+ integrate( markedelements(mesh,"fin_mesh"), _Q<imOrder>(), ( ratio*gradt(T)*trans(grad(v))) )
										+ integrate( markedfaces(mesh, "gamma1"), _Q<imOrder>(), Bi*idv(T));

		// At this step, all the elements have been introduced to the equation EXCEPT the discretization terms (right and left ones)
        // this part is included in the run() method

		Dcst->close();
        if ( this->vm().count( "export-matlab" ) )
            Dcst->printMatlab( "Dcst" );
			
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

		element_type u( Xh, "u" );
		element_type v( Xh, "v" );

		/*
		 * Construction of the left hand side
		 */
		//sparse_matrix_ptrtype D( M_backend->newMatrix( Xh, Xh ) );

		Log() << "Bi = " << Bi << "\n"
			  << "meshSize = " << meshSize << "\n"
			  << "depth = " << depth << "\n"
			  << "lambda_spreader = " << lambda_spreader << "\n"
			  << "lambda_fin = " << lambda_fin << "\n"
			  << "density_s = " << density_s << "\n"
			  << "density_f = " << density_f << "\n";

		std::cout << " DEBUT DE LA BOUCLE POUR BDF \n";

		/* discretization with BDF */
		for ( M_bdf->start(); M_bdf->isFinished(); M_bdf->next() ) {

			auto bdf_poly = M_bdf->polyDeriv();
			form2(Xh, Xh, Dcst) += integrate( markedelements(mesh), density_s*trans(idv(T))*id(v)*M_bdf->polyDerivCoefficient(0) );
			form1( Xh, Fcst) += integrate( elements(mesh), density_s*idv(bdf_poly)*idv(v) );

		}

		std::cout << " FIN DE LA BOUCLE \n";

		Dcst->close();

		if ( this->vm().count( "export-matlab" ) )
			Dcst->printMatlab( "D" );

		this->solve( Dcst, u, Fcst );

		double moy_u = ( integrate( markedfaces(mesh,1), _Q<imOrder>(), idv(u) ).evaluate()(0,0) /
							integrate( markedfaces(mesh,1), _Q<imOrder>(), constant(1.0) ).evaluate()(0,0) );

		std::cout.precision( 5 );
		std::cout << std::setw( 5 ) << Bi << " "
				  << std::setw( 10 ) << moy_u << "\n";
		this->exportResults( u );

} // HeatSink::run
	
template<int Dim, int Order>
void
HeatSink<Dim, Order>::solve( sparse_matrix_ptrtype& D, element_type& u, vector_ptrtype& F  )
{
    vector_ptrtype U( M_backend->newVector( u.functionSpace()->map() ) );
    M_backend->solve( D, D, U, F );
    u = *U;
} // HeatSink::solve


template<int Dim, int Order>
void
HeatSink<Dim, Order>::exportResults( element_type& U )
{
    if ( this->vm().count( "export" ) )
        {
            exporter->step(1.)->setMesh( U.functionSpace()->mesh() );
            exporter->step(1.)->add( "ProcessId",
                           regionProcess( boost::shared_ptr<p0_space_type>( new p0_space_type( U.functionSpace()->mesh() ) ) ) );
            exporter->step(1.)->add( "Temperature", U );
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






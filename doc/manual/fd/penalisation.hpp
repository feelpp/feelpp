#ifndef PENALISATION_HPP
#define PENALISATION_HPP

#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feelalg/backend.hpp>

#include <iostream>
#include <fstream>

#define VELOCITY_ORDER 3
#define PRESSURE_ORDER 2
#define CARAC_ORDER 2


namespace Feel
{

inline po::options_description makeOptions()
{
	po::options_description mesOptions( "Mes options pour l'appli" );
	mesOptions.add_options()
		( "OutFolder", po::value<std::string>()->default_value( "" ), "Dossier Sortie" )
		( "hsize", po::value<double>()->default_value( 0.1 ), "h size" )
		( "x0", po::value<double>()->default_value( 0 ), "x0" )
		( "y0", po::value<double>()->default_value( 0.333 ), "y0" )
		( "z0", po::value<double>()->default_value( 0 ), "z0" )
		( "nu", po::value<double>()->default_value( 0.001002 ),"viscosity" )
		( "radius", po::value<double>()->default_value( 0.35 ), "Radius" )
		( "Q", po::value<double>()->default_value( 2. ), "In Flow rate" )
		( "RQ", po::value<double>()->default_value( 0.5 ), "Flow rate ratio Q1/Q0" )
		( "epsilonpress", po::value<double>()->default_value( 0.00001 ), "penalisation term : div u = epsilonpress p  " )
		( "epsilon", po::value<double>()->default_value( 0.8 ), "penalisation term for particle" )
		( "Tfinal", po::value<double>()->default_value( 3 ), "Final time" )
		( "Ylimit", po::value<double>()->default_value( 0 ), "Stop the simulation if y reaches Ylimit (never stop if 0)" )
		( "DT", po::value<double>()->default_value( 0.01 ), "time step" )
		( "test", po::value<int>()->default_value( 1 ), "test application" );
	return mesOptions
		.add( Feel::feel_options() )
		.add( backend_options( "stokes_backend" ) )
		;
}//makeoptions()

inline AboutData makeAbout()
{
	AboutData about( "Penalisation", "penalisation", "0.1", "desc" );
	return about;
}

template <int Dim>
class Penalisation : public Application
{
	// mesh
	typedef Simplex<Dim,1,Dim> entity_type;

	typedef Mesh<entity_type> mesh_type;
	typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

	// bases
	typedef Lagrange<VELOCITY_ORDER, Vectorial> basis_veloc_type;
	typedef Lagrange<PRESSURE_ORDER, Scalar> basis_pressure_type;
	typedef Lagrange<0, Scalar> basis_lag_type;
	// typedef bases<Lagrange<CARAC_ORDER,Scalar,Discontinuous> > basis_carac_type;
	typedef bases<Lagrange<CARAC_ORDER,Scalar> > basis_carac_type;
	typedef bases<Lagrange<0,Scalar,Discontinuous> > basis_p0_type;
	typedef bases<basis_veloc_type, basis_pressure_type, basis_lag_type> basis_type;
	//typedef bases<basis_veloc_type, basis_pressure_type> basis_type;


	// space
	typedef FunctionSpace<mesh_type, basis_type> space_type;
	typedef boost::shared_ptr<space_type> space_ptrtype;
	typedef FunctionSpace<mesh_type, basis_carac_type,Discontinuous> space_carac_type;
	typedef boost::shared_ptr<space_carac_type> space_carac_ptrtype;
	typedef FunctionSpace<mesh_type, basis_p0_type,Discontinuous> space_p0_type;
	typedef boost::shared_ptr<space_p0_type> space_p0_ptrtype;

	// elements
	typedef typename space_type::element_type element_type;
	typedef boost::shared_ptr<element_type> element_ptrtype;
	typedef typename element_type::template sub_element<0>::type element_veloc_type;
	typedef typename element_type::template sub_element<1>::type element_pressure_type;
	typedef typename element_type::template sub_element<2>::type element_lag_type;
	typedef typename space_carac_type::element_type element_carac_type;
	typedef boost::shared_ptr<element_carac_type> element_carac_ptrtype;
	typedef typename space_p0_type::element_type element_p0_type;
	typedef boost::shared_ptr<element_p0_type> element_p0_ptrtype;

	// backend
	typedef Backend<double> backend_type;
	typedef boost::shared_ptr<backend_type> backend_ptrtype;

	//algebra
	typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
	typedef typename backend_type::vector_type vector_type;
	typedef boost::shared_ptr<sparse_matrix_type> sparse_matrix_ptrtype;
	typedef boost::shared_ptr<vector_type> vector_ptrtype;

	//export
	typedef Exporter<mesh_type> export_type;
	typedef boost::shared_ptr<export_type> export_ptrtype;

public :
	Penalisation();
	void run();

private :
	void exportResults();
	void updatePosition();
	void updateChi();
	void stokes();
	void addCL();
	void initStokes();

	const double nu;
	const double radius;

	const double H1;
	const double H2;
	const double L1;
	const double E;

	const double Q;
	const double Qtop;
	const double Qbot;

	space_ptrtype Xh;
	space_carac_ptrtype Yh;
	space_p0_ptrtype P0h;

	mesh_ptrtype mesh,mesh_visu;
	export_ptrtype exporter;
	backend_ptrtype M_backend;

	sparse_matrix_ptrtype C;
	sparse_matrix_ptrtype D;
	vector_ptrtype F;

	element_ptrtype U;
	element_carac_type carac;
	element_p0_type p0;

	double xp,yp,zp;
	double Vx, Vy, Vz;
	const double epsilonpress, epsilon;
	double Tfinal,t;
	const double dt;
	const double Ylimit;
	int iter;

	boost::timer chrono;

}; //Penalisation

}//namespace Feel

#endif //PENALISATION_HPP

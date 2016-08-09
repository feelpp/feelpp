#ifndef _MIXEDPOISSONELASTICITY_HPP
#define _MIXEDPOISSONELASTICITY_HPP

#include <mixedpoisson2.hpp>
#include <mixedelasticity.hpp>

namespace Feel
{

namespace FeelModels
{

inline
po::options_description
makeMixedPoissonElasticityOptions( std::string prefix = "mixedpoissonelasticity" )
{
	po::options_description mpOptions( "Mixed Poisson Elasticity HDG options");
	mpOptions.add ( makeMixedPoissonOptions("mixedpoisson") );
	mpOptions.add ( makeMixedElasticityOptions("mixedelasticity") );

	return mpOptions;
}

inline po::options_description
makeMixedPoissonElasticityLibOptions( std::string prefix = "mixedpoissonelasticity" )
{
    po::options_description mpLibOptions( "Mixed Poisson Elasticity HDG Lib options");
    // if ( !prefix.empty() )
    //     mpLibOptions.add( backend_options( prefix ) );
    return mpLibOptions;
}


template <int Dim, int Order>
class MixedPoissonElasticity
{
 public:
	typedef double 						value_type;
	typedef MixedPoisson<Dim,Order,1> 	mp_type;
	typedef MixedElasticity<Dim,Order> 	me_type;
	typedef boost::shared_ptr<mp_type>  mp_ptrtype;
	typedef boost::shared_ptr<me_type>  me_ptrtype;
	typedef MixedPoissonElasticity<Dim,Order> 	self_type;
	typedef boost::shared_ptr<self_type> 		self_ptrtype;

    using sparse_matrix_type = backend_type::sparse_matrix_type;
    using sparse_matrix_ptrtype = backend_type::sparse_matrix_ptrtype;
    using vector_type = backend_type::vector_type;
    using vector_ptrtype = backend_type::vector_ptrtype;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order G_order
    typedef Simplex<Dim,1> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    // The Lagrange multiplier lives in R^n-1
    typedef Simplex<Dim-1,1,Dim> face_convex_type;
    typedef Mesh<face_convex_type> face_mesh_type;
    typedef boost::shared_ptr<face_mesh_type> face_mesh_ptrtype;
    
   	/*
// ---- //
    using Vh_t =  Pdhms_type<mesh_type,Order>;
    using Vh_ptr_t =  Pdhms_ptrtype<mesh_type,Order>;
    using Vh_element_t = typename Vh_t::element_type;
    using Vh_element_ptr_t = typename Vh_t::element_ptrtype;
// ---- //
    using Wh_t =  Pdhv_type<mesh_type,Order>;
    using Wh_ptr_t =  Pdhv_ptrtype<mesh_type,Order>;
    using Wh_element_t = typename Wh_t::element_type;
    using Wh_element_ptr_t = typename Wh_t::element_ptrtype;
// ---- //
    using Mh_t =  Pdhv_type<face_mesh_type,Order>;
    using Mh_ptr_t =  Pdhv_ptrtype<face_mesh_type,Order>;
    using Mh_element_t = typename Mh_t::element_type;
    using Mh_element_ptr_t = typename Mh_t::element_ptrtype;
// ---- //
    using M0h_t =  Pdh_type<face_mesh_type,0>;
    using M0h_ptr_t =  Pdh_ptrtype<face_mesh_type,0>;
    using M0h_element_t = typename M0h_t::element_type;
    using M0h_element_ptr_t = typename M0h_t::element_ptrtype;
    
	using op_interp_ptrtype = boost::shared_ptr<OperatorInterpolation<Wh_t, Pdhv_type<mesh_type,Order>>>;
    using opv_interp_ptrtype = boost::shared_ptr<OperatorInterpolation<Vh_t, Pdhms_type<mesh_type,Order>>>;
 	*/
    
	//! Model properties type
    using model_prop_type = ModelProperties;
    using model_prop_ptrtype = std::shared_ptr<model_prop_type>;

    // using product_space_std = ProductSpaces<Vh_ptr_t,Wh_ptr_t,Mh_ptr_t>;
    // using bilinear_block_std = BlockBilinearForm<product_space_std>;

    typedef Exporter<mesh_type,1> exporter_type;
    typedef boost::shared_ptr <exporter_type> exporter_ptrtype;
    
    // typedef Newmark<space_mixedelasticity_type>  newmark_type;
    // typedef Newmark <Wh_t> newmark_type;
    // typedef boost::shared_ptr<newmark_type> newmark_ptrtype;
 
 private:
	mesh_ptrtype M_mesh;
	
	mp_ptrtype M_PoissonModel;
	me_ptrtype M_ElasticityModel;


 public:
	MixedPoissonElasticity(mesh_ptrtype _mesh) : M_mesh(_mesh) {
		M_PoissonModel = mp_type::New("mixedpoisson");
		M_PoissonModel->init( _mesh );

		M_ElasticityModel = me_type::New("mixedelasticity");
		M_ElasticityModel->init( _mesh );
	}

	
	void assembleF_Poisson();  				// this is the assembleF of MixedPoisson
	
	void assembleF_Elasticity();			// this is the assembleF of MixedElasticity

	void solvePoisson() { this->assembleF_Poisson(); M_PoissonModel->solve(); }
	void solveElasticity() { this->assembleF_Elasticity(); M_ElasticityModel->solve(); }
 	
	void run();
	
}; // end class declaration


template <int Dim, int Order>
void
MixedPoissonElasticity<Dim,Order>::assembleF_Poisson()
{
	auto ps = M_PoissonModel->getPS();
	auto F = M_PoissonModel->getF();

	auto blf = blockform1 (*ps, F);
	auto w = M_PoissonModel->potentialSpace()->element();
	auto dt = M_PoissonModel->timeStepBDF()->timeStep();

	// - <d/dt div(u),w> 
	blf(1_c) += integrate( _range=elements(M_mesh), 
	 					   _expr= -( div(M_ElasticityModel-> potentialField())-div(M_ElasticityModel->timeStepNM()->previousUnknown()) ) * id(w) / dt );

}

template <int Dim, int Order>
void
MixedPoissonElasticity<Dim,Order>::assembleF_Elasticity()
{
	 
	auto ps = product( M_ElasticityModel->fluxSpace(), M_ElasticityModel->potentialSpace(), M_ElasticityModel->traceSpace() );

	auto blf = blockform1 ( ps, M_ElasticityModel->getF() );
	auto v = M_ElasticityModel->fluxSpace()->element( "v" );

	// Get data from Poisson model
	auto pressure = M_PoissonModel->potentialField();

	// Feel::cout << "Pressure: " << pressure << std::endl;
	// - < pI , v>
	blf( 0_c ) += integrate( _range=elements(M_mesh), 
	 						 _expr= - inner( idv(M_PoissonModel->potentialField())*eye<Dim,Dim>(),  id(v)) );

}

template <int Dim, int Order>
void MixedPoissonElasticity<Dim,Order>::run()
{
	auto t_init = M_PoissonModel->timeStepBDF()->timeInitial();
	auto dt = M_PoissonModel->timeStepBDF()->timeStep();
	auto t_fin = M_PoissonModel->timeStepBDF()->timeFinal();	
	
	for (; !M_PoissonModel->timeStepBase()->isFinished() && !M_ElasticityModel->timeStepBase()->isFinished() ; M_PoissonModel->updateTimeStep() )
	{
		Feel::cout << "===============================================" << std::endl;
		Feel::cout << "time simulation: " << M_PoissonModel->time() << "s \n";	
		Feel::cout << "===============================================" << std::endl;

		// Elasticity problem
		this->assembleF_Elasticity();
		M_ElasticityModel->solve();	
		M_ElasticityModel->exportResults(M_mesh);

		// Poisson problem
		this->assembleF_Poisson();
		M_PoissonModel->solve();	
		M_PoissonModel->exportResults(M_mesh);

		// update
		M_ElasticityModel->updateTimeStep();
	}

}

} // end namespace FeelModels

} // end namespace Feel

#endif


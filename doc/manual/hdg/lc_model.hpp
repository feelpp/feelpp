/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel++ library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
             Lorenzo Sala <sala@unistra.fr>
       Date: 2016-02-10

  Copyright (C) 2016 Feel++ Consortium

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

#ifndef LC_MODEL_HPP
#define	LC_MODEL_HPP

#include <feel/feelmodels/hdg/mixedpoisson.hpp>
#include <boost/numeric/odeint.hpp>

namespace Feel {

namespace FeelModels {

using value_type = double;
typedef boost::array <value_type,3> state_type;


struct ode_model {

    // Inverse matrix of capacitance
    const boost::numeric::ublas::matrix<value_type> Cinv ;
    // matrix of resistance
    const boost::numeric::ublas::matrix<value_type> A ;
    // rhs of the initial equation
    const boost::numeric::ublas::vector<value_type> g ;

    ode_model( boost::numeric::ublas::matrix<value_type> &Cinv,
               boost::numeric::ublas::matrix<value_type> &A,
               boost::numeric::ublas::vector<value_type> &g ) : Cinv(Cinv),A(A),g(g) {}

    void operator()( const state_type &x , state_type &dxdt , double t ) const {
		/*
		Feel::cout << "Cinv: [" << Cinv(0,0) << " , " << Cinv(1,1) << " , " << Cinv(2,2) << " ]" << std::endl;

		Feel::cout << "A: " << std::endl << A(0,0) << " , " << A(0,1) << " , " << A(0,2)  << std::endl;
		Feel::cout << A(1,0) << " , " << A(1,1) << " , " << A(1,2)  << std::endl;
        Feel::cout << A(2,0) << " , " << A(2,1) << " , " << A(2,2)  << std::endl << std::endl;

        Feel::cout << "g: [" << g(0) << " , " << g(1) << " , " << g(2) << " ]" << std::endl;
		*/
		// Feel::cout << "Initial solution second step: " << x[0] << std::endl;

		// solution
        dxdt[0] = Cinv(0,0) * ( A(0,0)*x[0] + A(0,1)*x[1] + A(0,2)*x[2] ) + Cinv(0,0) * g(0) ;
        dxdt[1] = Cinv(1,1) * ( A(1,0)*x[0] + A(1,1)*x[1] + A(1,2)*x[2] ) + Cinv(1,1) * g(1) ;
        dxdt[2] = Cinv(2,2) * ( A(2,0)*x[0] + A(2,1)*x[1] + A(2,2)*x[2] ) + Cinv(2,2) * g(2) ;
    }
};

template<int Dim, int Order, int G_Order = 1, int E_Order = 4>
class LaminaCribrosa : public MixedPoisson<Dim,Order,G_Order,E_Order>
{

public:
    typedef MixedPoisson<Dim,Order,G_Order,E_Order> super_type;


    typedef LaminaCribrosa<Dim,Order,G_Order,E_Order> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

   	static const uint16_type expr_order = Order + E_Order;

    using convex_type = typename super_type::convex_type;
    using mesh_type = typename super_type::mesh_type;
    using mesh_ptrtype = typename super_type::mesh_ptrtype;
    using op_interp_ptrtype = typename super_type::op_interp_ptrtype;
    using opv_interp_ptrtype = typename super_type::opv_interp_ptrtype;

    // Ch
    using Ch_t = typename super_type::Ch_t;
    using Ch_ptr_t = typename super_type::Ch_ptr_t;
    using Ch_element_t = typename super_type::Ch_element_t;
    using Ch_element_ptr_t = typename super_type::Ch_element_ptr_t;
    using Ch_element_vector_type = typename super_type::Ch_element_vector_type;

    typedef Bdf <Ch_t> statevar_bdf_type;
    typedef std::shared_ptr<statevar_bdf_type> statevar_bdf_ptrtype;

    using product2_space_type = typename super_type::product2_space_type;
    using integral_boundary_list_type = typename super_type::integral_boundary_list_type;

	typedef boost::numeric::ublas::vector<value_type> state_vector_type;
	typedef std::shared_ptr<state_vector_type> state_vector_ptrtype;

private:

    Ch_element_ptr_t M_Y;
    // Ch_element_vector_type  M_Y;

    statevar_bdf_ptrtype M_bdf_statevariable;
    state_type M_statevar_solution;
    boost::numeric::ublas::matrix<value_type> M_A0d ;
    boost::numeric::ublas::matrix<value_type> M_Cinv ;
    // boost::numeric::ublas::vector<value_type> M_g ;
    state_vector_type M_g ;

    int M_0dCondition;
    integral_boundary_list_type M_0dList;


public:

	LaminaCribrosa( std::string const& prefix = "hdg.poisson",
                    worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                    std::string const& subPrefix = "",
                    ModelBaseRepository const& modelRep = ModelBaseRepository() )
		: super_type( prefix, MixedPoissonPhysics::None, worldComm, subPrefix, modelRep ),
          ModelBase( prefix, worldComm, subPrefix, modelRep )
        {}
    LaminaCribrosa( self_type const& LC ) = default;


    static self_ptrtype New( std::string const& prefix = "hdg.poisson",
                             worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                             std::string const& subPrefix = "",
                             ModelBaseRepository const& modelRep = ModelBaseRepository() );


    virtual void assembleCstPart();
    virtual void assembleNonCstPart();
    void assemble0d( int i );
    void assembleRhs0d( int i );

    virtual void initModel();
    virtual void initSpaces();

    void exportResults( double time, mesh_ptrtype mesh = nullptr, op_interp_ptrtype Idh = nullptr, opv_interp_ptrtype Idhv = nullptr );
    void exportResults(mesh_ptrtype mesh = nullptr, op_interp_ptrtype Idh = nullptr, opv_interp_ptrtype Idhv = nullptr)
	{
	   this->exportResults (this->currentTime(), mesh, Idh, Idhv );
	   this->exporterMP()->save();
	}

    virtual void solve();


    // time step scheme
    virtual void createTimeDiscretization() ;
    virtual void updateTimeStepBDF();
    virtual void initTimeStep();
    statevar_bdf_ptrtype timeStepBDF_statevar() { return M_bdf_statevariable; }
    statevar_bdf_ptrtype const& timeStepBDF_statevar() const { return M_bdf_statevariable; }
    std::shared_ptr<TSBase> timeStepBase_statevar() { return this->timeStepBDF_statevar(); }
    std::shared_ptr<TSBase> timeStepBase_statevar() const { return this->timeStepBDF_statevar(); }

    // For the second step
	void odeForceTermEvaluation( double time );
    void second_step();
    // void ode_model(const state_type x , state_type &dxdt , const value_type t );

};


template<int Dim, int Order, int G_Order, int E_Order>
typename LaminaCribrosa<Dim,Order, G_Order, E_Order>::self_ptrtype
LaminaCribrosa<Dim,Order,G_Order, E_Order>::New( std::string const& prefix,
                                                 worldcomm_ptr_t const& worldComm,
                                                 std::string const& subPrefix,
                                                 ModelBaseRepository const& modelRep )
{
    return std::make_shared<self_type>( prefix, worldComm, subPrefix, modelRep );
}




template<int Dim, int Order, int G_Order, int E_Order>
void
LaminaCribrosa<Dim, Order, G_Order, E_Order>::initTimeStep()
{
    super_type::initTimeStep();
    // start or restart time step scheme
    if (!this->doRestart())
    {
        // start time step
        M_bdf_statevariable -> start( *M_Y );

        // up current time
        // this->updateTime( M_bdf_statevariable -> time() );
    }
    else
    {
        // start time step
        M_bdf_statevariable->restart();
        // load a previous solution as current solution
        *M_Y = M_bdf_statevariable->unknown(0);
        // up initial time
        this->setTimeInitial( M_bdf_statevariable->timeInitial() );
        // up current time
        // this->updateTime( M_bdf_statevariable->time() );

        this->log("LaminaCribrosa","initTimeStep", "restart bdf/exporter done" );
    }
}


template<int Dim, int Order, int G_Order, int E_Order>
void
LaminaCribrosa<Dim, Order, G_Order, E_Order>::updateTimeStepBDF()
{
    super_type::updateTimeStepBDF();

    this->log("LaminaCribrosa","updateTimeStepBDF", "start" );
    this->timerTool("TimeStepping").setAdditionalParameter("time_Y",this->currentTime());
    this->timerTool("TimeStepping").start();

    int previousTimeOrder_statevar = this->timeStepBDF_statevar()->timeOrder();

    M_bdf_statevariable->next( *M_Y );

    int currentTimeOrder_statevar = this->timeStepBDF_statevar()->timeOrder();

    // this->updateTime( M_bdf_statevariable->time() );
    this->timerTool("TimeStepping").stop("updateBdf");
    if ( this->scalabilitySave() ) this->timerTool("TimeStepping").save();
    this->log("LaminaCribrosa","updateTimeStepBDF", "finish" );
}


template<int Dim, int Order, int G_Order, int E_Order>
void
LaminaCribrosa<Dim,Order, G_Order, E_Order>::createTimeDiscretization()
{
    super_type::createTimeDiscretization();

    this->log("LaminaCribrosa","createTimeDiscretization", "start" );
    this->timerTool("Constructor").start();

    std::string myFileFormat = soption(_name="ts.file-format");// without prefix
    std::string suffixName = "";
    if ( myFileFormat == "binary" )
         suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();

   M_bdf_statevariable = bdf( _vm=Environment::vm(), _space=this->M_Ch ,
                       _name=prefixvm(this->prefix(),prefixvm(this->subPrefix(),"Y"+suffixName)) ,
                       _prefix="",
                       _initial_time=this->timeInitial(),
                       _final_time=this->timeFinal(),
                       _time_step=this->timeStep(),
                       _restart=this->doRestart(),
                       _restart_path=this->restartPath(),
                       _restart_at_last_save=this->restartAtLastSave(),
                       _save=this->tsSaveInFile(), _freq=this->tsSaveFreq() );
    M_bdf_statevariable->setfileFormat( myFileFormat );
    M_bdf_statevariable->setPathSave( (fs::path(this->rootRepository()) /
                               fs::path( prefixvm(this->prefix(), (boost::format("bdfY_o_%1%_dt_%2%")%M_bdf_statevariable->bdfOrder()%this->timeStep() ).str() ) ) ).string() );

    double tElapsed = this->timerTool("Constructor").stop("createTimeDiscr");
    this->log("LaminaCribrosa","createTimeDiscretization", (boost::format("finish in %1% s") %tElapsed).str() );

}


// Overriding of init Model
template<int Dim, int Order, int G_Order, int E_Order>
void
LaminaCribrosa<Dim, Order, G_Order, E_Order>::initModel(){

    super_type::initModel();

    M_0dList.clear();
    auto itField = this->modelProperties().boundaryConditions().find( "flux");
    if ( itField != this->modelProperties().boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Integral_coupled_with_0d" );
        if ( itType != mapField.end() )
        {
            Feel::cout << "Integral coupled with 0d equation:";
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if ( this->mesh()->hasFaceMarker(marker) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
                this->M_IBCList.push_back(exAtMarker);
				M_0dList.push_back(exAtMarker);
            }
            Feel::cout << std::endl;
        }
    }

    if ( this->M_IBCList.empty() )
        this->M_integralCondition = 0;
    else
        this->M_integralCondition = this->M_IBCList.size();

    if ( M_0dList.empty() )
        M_0dCondition = 0;
    else
        M_0dCondition = M_0dList.size();

    // Initialization for second step: the 0d equation
    M_A0d.resize(3,3);
    M_Cinv.resize(3,3);
    M_g.resize(3);


    for( auto const& pairMat : this->modelProperties().materials() )
    {
        auto material = pairMat.second;
	    auto C1 = material.getDouble("C1");
	    auto C2 = material.getDouble("C2");
	    auto C3 = material.getDouble("C3");
    	auto R12 = material.getDouble("R12");
	    auto R23 = material.getDouble("R23");
		auto Rout = material.getDouble("Rout");

    	// Initialize matrices and vector of the ODE
    	M_A0d(0,0) = -1/R12;
    	M_A0d(0,1) = 1/R12;
    	M_A0d(0,2) = 0;
    	M_A0d(1,0) = 1/R12;
    	M_A0d(1,1) = -1/R12 - 1/R23;
    	M_A0d(1,2) = +1/R23;
		M_A0d(2,0) = 0;
    	M_A0d(2,1) = 1/R23;
    	M_A0d(2,2) = -1/R23 -1/Rout;

    	M_Cinv(0,0) = 1/C1;
    	M_Cinv(1,1) = 1/C2;
    	M_Cinv(2,2) = 1/C3;

    }

}

// Overriding of init Spaces

template<int Dim, int Order, int G_Order, int E_Order>
void
LaminaCribrosa<Dim, Order, G_Order, E_Order>::initSpaces(){


    // for( int i = 0; i < M_0dCondition; i++)
    //    this->M_IBCList.push_back(M_0dList[i]);

    super_type::initSpaces();

    // for( int i = 0; i < M_0dCondition; i++)
    //    this->M_IBCList.pop_back();

    for( int i = 0; i < M_0dCondition; i++ )
        (this->M_mup).push_back(this->M_Ch->element("mup"));


	if(M_0dCondition)
	{
		Feel::cout << "Number of 0d: " << M_0dCondition << std::endl;
        solve::strategy s = this->M_useSC ? solve::strategy::static_condensation : solve::strategy::monolithic;
    	auto ibcSpaces = std::make_shared<ProductSpace<Ch_ptr_t,true> >( this->integralCondition() + M_0dCondition, this->M_Ch);
    	this->M_ps = std::make_shared<product2_space_type>(product2(ibcSpaces,this->M_Vh,this->M_Wh,this->M_Mh));

        this->M_A_cst = makeSharedMatrixCondensed<value_type>(s, csrGraphBlocks(*this->M_ps, (s>=solve::strategy::static_condensation)?Pattern::ZERO:Pattern::COUPLED), *this->M_backend );
#ifndef USE_SAME_MAT
        this->M_A = makeSharedMatrixCondensed<value_type>(s,  csrGraphBlocks(*this->M_ps, (s>=solve::strategy::static_condensation)?Pattern::ZERO:Pattern::COUPLED), *this->M_backend );
#endif
        this->M_F = makeSharedVectorCondensed<value_type>(s, blockVector(*this->M_ps), *this->M_backend, false);
	}

    // Init for second step
    M_Y = this->constantSpace()->elementPtr( "yy" );

	// Initialization variables second step
    M_statevar_solution.fill(0);
	if(M_0dCondition)
	{
    	auto itField = this->modelProperties().boundaryConditions().find( "InitialCondition_CircuitModel");
    	if ( itField != this->modelProperties().boundaryConditions().end() )
    	{
        	auto mapField = (*itField).second;
        	auto itType = mapField.find( "Pi1" );
        	if ( itType != mapField.end() )
        	{
            	for ( auto const& exAtMarker : (*itType).second )
            	{
					auto expression =  expr(exAtMarker.expression());
					expression.setParameterValues( { {"t", this->time() } } );
					auto P1 = mean( _range = markedfaces(this->mesh(),exAtMarker.marker()), _expr = expression )(0,0);
					M_statevar_solution[0] = P1;
				}
			}

			itType = mapField.find( "Pi2" );
        	if ( itType != mapField.end() )
        	{
            	for ( auto const& exAtMarker : (*itType).second )
            	{
					auto expression =  expr(exAtMarker.expression());
					expression.setParameterValues( { {"t", this->time() } } );
					auto P2 = mean( _range = markedfaces(this->mesh(),exAtMarker.marker()), _expr = expression )(0,0);
					M_statevar_solution[1] = P2;
				}
			}

			itType = mapField.find( "Pi3" );
        	if ( itType != mapField.end() )
        	{
            	for ( auto const& exAtMarker : (*itType).second )
            	{
					auto expression =  expr(exAtMarker.expression());
					expression.setParameterValues( { {"t", this->time() } } );
					auto P3 = mean( _range = markedfaces(this->mesh(),exAtMarker.marker()), _expr = expression )(0,0);
					M_statevar_solution[2] = P3;
				}
			}
		}
	}
}



template<int Dim, int Order, int G_Order, int E_Order>
void
LaminaCribrosa<Dim, Order, G_Order, E_Order>::assembleCstPart( )
{
    for( int i = 0; i < M_0dCondition; i++)
        this->M_IBCList.push_back(M_0dList[i]);

	super_type::assembleCstPart();

    for( int i = 0; i < M_0dCondition; i++)
        this->M_IBCList.pop_back();

    for ( int i = 0; i < M_0dList.size(); i++ )
        this->assemble0d( i );

}

template<int Dim, int Order, int G_Order, int E_Order>
void
LaminaCribrosa<Dim, Order, G_Order, E_Order>::assembleNonCstPart( )
{
    super_type::assembleNonCstPart();

    for ( int i = 0; i < M_0dList.size(); i++ )
        this->assembleRhs0d( i );

}


template<int Dim, int Order, int G_Order, int E_Order>
void
LaminaCribrosa<Dim, Order, G_Order, E_Order>::assemble0d( int i )
{

    this->log("LaminaCribrosa","assembleMatrix0d", "start" );
    this->timerTool("Constructor").start();

    auto bbf = blockform2( *(this->getPS()), this->M_A_cst);

    auto u = this->fluxSpace()->element( "u" );
    auto p = this->potentialSpace()->element( "p" );
     auto w = this->potentialSpace()->element( "w" );
    auto nu = this->constantSpace()->element("nu");
    auto uI = this->constantSpace()->element( "uI" );
    auto yy = this->constantSpace()->element( "yy" );

    auto H = this->traceSpaceOrder0()->element( "H" );
    if ( ioption(prefixvm(this->prefix(), "hface") ) == 0 )
        H.on( _range=elements(this->traceSpaceOrder0()->mesh()), _expr=cst(this->fluxSpace()->mesh()->hMax()) );
    else if ( ioption(prefixvm(this->prefix(), "hface") ) == 1 )
        H.on( _range=elements(this->traceSpaceOrder0()->mesh()), _expr=cst(this->fluxSpace()->mesh()->hMin()) );
    else if ( ioption(prefixvm(this->prefix(), "hface") ) == 2 )
        H.on( _range=elements(this->traceSpaceOrder0()->mesh()), _expr=cst(this->fluxSpace()->mesh()->hAverage()) );
    else
        H.on( _range=elements(this->traceSpaceOrder0()->mesh()), _expr=h() );
    // stabilisation parameter
    auto tau_constant = cst(doption(prefixvm(this->prefix(), "tau_constant")));


    auto marker = M_0dList[i].marker();
    int j = this->integralCondition()+i; // index where to start

	// Feel::cout << "j: " << j << "\t i: " << i << std::endl;

    // 0D EQUATION
    double meas = integrate( _range=markedfaces(this->mesh(),marker),_expr=cst(1.0)).evaluate()(0,0);

    // - <j.n, mu3>_Gamma_I
    bbf ( 3_c, 0_c, j, 0) += integrate( _range=markedfaces(this->mesh(),marker), _expr= -trans(idt(u))*N()*id(nu) );

    // - <tau p, mu3>_Gamma_I
    bbf ( 3_c, 1_c, j, 1) += integrate( _range=markedfaces(this->mesh(),marker), _expr= -tau_constant*idt(p)*id(nu) );

    // + <tau u_I, mu3>_Gamma_I
    bbf ( 3_c, 3_c, j, j-1) += integrate( _range=markedfaces(this->mesh(),marker), _expr= tau_constant*idt(uI)*id(nu) );

	for( auto const& pairMat : this->modelProperties().materials() )
    {
       	auto material = pairMat.second;
        auto RR = material.getScalar("RR");       // Resistance of the buffer
        auto CC = material.getScalar("CC");       // Capacitance of the buffer

	    // -1/(R |Gamma_I|) <u_I, mu2>_Gamma_I
        bbf ( 3_c, 3_c, j-1, j-1) += integrate( _range=markedfaces(this->mesh(),marker), _expr = - idt(uI)*id(nu)/RR/meas ) ;

        // +1/(R |Gamma_I|) <Y,mu2>_Gamma_I
        bbf ( 3_c, 3_c, j-1, j) += integrate( _range=markedfaces(this->mesh(),marker), _expr = idt(yy)*id(nu)/RR/meas ) ;

      	// < C/|Gamma_I| Y/dt, mu3>_Gamma_I
        bbf ( 3_c, 3_c, j, j) += integrate( _range=markedfaces(this->mesh(),marker), _expr= CC*this->timeStepBDF_statevar()->polyDerivCoefficient(0) * idt(yy)*id(nu)/meas );
    }

    double tElapsed = this->timerTool("Constructor").stop("assembleMatrix0d");
    this->log("LaminaCribrosa","assembleMatrix0d", (boost::format("finish in %1% s") %tElapsed).str() );
}



template<int Dim, int Order, int G_Order, int E_Order>
void
LaminaCribrosa<Dim, Order, G_Order, E_Order>::assembleRhs0d( int i )
{

    this->log("LaminaCribrosa","assembleRhs0d", "start" );
    this->timerTool("Constructor").start();

    auto blf = blockform1( *(this->getPS()), this->getF() );
    auto nu = this->constantSpace()->element( "nu" );


    int j = this->integralCondition() + i; // index where to start


    auto exAtMarker = M_0dList[i];
    auto marker = exAtMarker.marker();
    auto g = expr(exAtMarker.expression());

    double meas = integrate( _range=markedfaces(this->mesh(),marker),_expr=cst(1.0)).evaluate()(0,0);;

	// 0d PART
    for( auto const& pairMat : this->modelProperties().materials() )
    {
        auto material = pairMat.second;
        auto CC = material.getScalar("CC");       // Capacitance of the buffer

        auto bdf_poly = M_bdf_statevariable->polyDeriv();
		// < C/|Gamma_I| Yold/dt, mu3>
        blf( 3_c, j ) += integrate( _range = markedfaces(this->mesh(),marker), _expr = CC*idv(bdf_poly) * id(nu)/meas );
    }

	double tElapsed = this->timerTool("Constructor").stop("assembleRhs0d");
    this->log("LaminaCribrosa","assembleRhs0d", (boost::format("finish in %1% s") %tElapsed).str() );
}



template <int Dim, int Order, int G_Order, int E_Order>
void
LaminaCribrosa<Dim,Order, G_Order, E_Order>::exportResults( double time, mesh_ptrtype mesh, op_interp_ptrtype Idh, opv_interp_ptrtype Idhv)
{
    super_type::exportResults( time, mesh, Idh, Idhv );
    this->log("LaminaCribrosa","exportResults", "start");

    // Export computed solutions
    for ( auto const& field : this->modelProperties().postProcess().exports().fields() )
    {
        if ( field == "state variable" )
        {
            this->exporterMP()->step( time )->add(prefixvm(this->prefix(), "Pi_1"), M_statevar_solution[0] );
            this->exporterMP()->step( time )->add(prefixvm(this->prefix(), "Pi_2"), M_statevar_solution[1] );
            this->exporterMP()->step( time )->add(prefixvm(this->prefix(), "Pi_3"), M_statevar_solution[2] );

            auto itField = this->modelProperties().boundaryConditions().find( "CircuitModel");
            if ( itField != this->modelProperties().boundaryConditions().end() )
            {
                auto mapField = (*itField).second;
                auto itType = mapField.find( "Pi1_exact" );
                if ( itType != mapField.end() )
                {
                    for ( auto const& exAtMarker : (*itType).second )
                    {
                        auto exprP1_exact =  expr(exAtMarker.expression());
                        exprP1_exact.setParameterValues( { {"t", this->time() } } );
                        auto P1_exact = exprP1_exact.evaluate()(0,0);
                        double mean_ex = std::abs(P1_exact);
                        if (mean_ex < 1e-10)
                            mean_ex = 1;
                        Feel::cout << "||P1-P1_ex|=\t" << std::abs(P1_exact - M_statevar_solution[0])/mean_ex << std::endl;
                        this->exporterMP() -> step( time )->add(prefixvm(this->prefix(), "P1_error"),  std::abs(P1_exact - M_statevar_solution[0])/mean_ex );
                    }
                }

                itType = mapField.find( "Pi2_exact" );
                if ( itType != mapField.end() )
                {
                    for ( auto const& exAtMarker : (*itType).second )
                    {
                        auto exprP2_exact =  expr(exAtMarker.expression());
                        exprP2_exact.setParameterValues( { {"t", this->time() } } );
                        auto P2_exact = exprP2_exact.evaluate()(0,0);
                        double mean_ex = std::abs(P2_exact);
                        if (mean_ex < 1e-10)
                            mean_ex = 1;
                        Feel::cout << "||P2-P2_ex|=\t" << std::abs(P2_exact - M_statevar_solution[1])/mean_ex << std::endl;
                        this->exporterMP() -> step( time )->add(prefixvm(this->prefix(), "P2_error"),  std::abs(P2_exact - M_statevar_solution[1])/mean_ex );
                    }
                }

                itType = mapField.find( "Pi3_exact" );
                if ( itType != mapField.end() )
                {
                    for ( auto const& exAtMarker : (*itType).second )
                    {
                        auto exprP3_exact =  expr(exAtMarker.expression());
                        exprP3_exact.setParameterValues( { {"t", this->time() } } );
                        auto P3_exact = exprP3_exact.evaluate()(0,0);
                        double mean_ex = std::abs(P3_exact);
                        if (mean_ex < 1e-10)
                            mean_ex = 1;
                        Feel::cout << "||P3-P3_ex|=\t" << std::abs(P3_exact - M_statevar_solution[2])/mean_ex << std::endl;
                        this->exporterMP() -> step( time )->add(prefixvm(this->prefix(), "P3_error"),  std::abs(P3_exact - M_statevar_solution[2])/mean_ex );
                    }
                }
                Feel::cout << "---------------------------" << std::endl;
            }
        }
    }
    this->log("LaminaCribrosa","exportResults", "finish");
}


template<int Dim, int Order, int G_Order, int E_Order>
void
LaminaCribrosa<Dim, Order, G_Order, E_Order>::solve()
{

#ifdef USE_SAME_MAT
    auto bbf = blockform2(*(this->getPS()), this->M_A_cst);
#else
    auto bbf = blockform2(*(this->getPS()), this->M_A);
#endif

    auto blf = blockform1(*(this->getPS()), this->M_F);

    auto U = this->getPS()->element();


    tic();
    bbf.solve(_solution=U, _rhs=blf, _condense=boption(prefixvm(this->prefix(), "use-sc")), _name=this->prefix());
    toc("LaminaCribrosa : solve");

    this->M_up = U(0_c);
    this->M_pp = U(1_c);

    for( int i = 0; i < this->integralCondition(); i++ )
        (this->M_mup)[i] = U(3_c,i);

	if (M_0dCondition)
	    *M_Y = U(3_c,this->integralCondition());


}



template<int Dim, int Order, int G_Order, int E_Order>
void
LaminaCribrosa<Dim, Order, G_Order, E_Order>::odeForceTermEvaluation( double time )
{
	value_type Piout = 0;

    auto itField = this->modelProperties().boundaryConditions().find( "CircuitModel");
    if ( itField != this->modelProperties().boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Piout" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
				auto exprPiout = expr(exAtMarker.expression());
                exprPiout.setParameterValues( { {"t", time } } );
				Piout = exprPiout.evaluate()(0,0);
				// Piout = mean( _range = markedfaces(this->mesh(),exAtMarker.marker()), _expr = exprPiout )(0,0);
			}
		}
	}

	for( auto const& pairMat : this->modelProperties().materials() )
	{
        auto material = pairMat.second;
		auto Rout = material.getDouble("Rout");
		M_g(0) = 0;
		M_g(1) = 0;
		M_g(2) = Piout/Rout;
	}
}



template<int Dim, int Order, int G_Order, int E_Order>
void
LaminaCribrosa<Dim, Order, G_Order, E_Order>::second_step()
{
	this->log("LaminaCribrosa","0D model", "start");
    tic();

	using namespace boost::numeric::odeint;
    using namespace boost::numeric::ublas;

    if (M_0dCondition)
    {
    	// Update the initial solution for Pi1 (for step 2)
    	// M_statevar_solution[0] = mean( _quad=_Q<expr_order>(), _range= elements(this->mesh()), _expr=idv(*M_Y) )(0,0) ;
    	M_statevar_solution[0] = (*M_Y)[0];

		auto u_i = mean( _quad=_Q<expr_order>(), _range= elements(this->mesh()), _expr=idv((this->M_mup)[0]) )(0,0) ;
    	Feel::cout << "Value of U_I after first step : \t " << u_i << std::endl;


    	Feel::cout << "Value of P1 after first step : \t " << (*M_Y)[0] << std::endl;

        auto marker = M_0dList[0].marker();
        double j_integral = integrate( _quad=_Q<expr_order>(), _range=markedfaces(this->mesh(),marker),_expr=trans(idv(this->M_up))*N()).evaluate()(0,0);


    	Feel::cout << "Integral value of the flow: \t " << j_integral << std::endl;

	    // solve the problem
		this->odeForceTermEvaluation( this->time() );
		runge_kutta4< state_type > stepper;
		// euler < state_type > stepper;
		auto initial_time = this->time()-this->timeStep();
		/*
		Feel::cout << "Initial time second step: " << initial_time << std::endl;
		Feel::cout << "Final time second step: " << this->time() << std::endl;
		Feel::cout << "Time step second step: " << this->timeStep() << std::endl;
		*/
	    boost::numeric::odeint::integrate_const(stepper, ode_model(M_Cinv,M_A0d,M_g), M_statevar_solution,
		    	this->time(),		 				// initial time
		    	this->time()+this->timeStep(), 		// final time
		    	this->timeStep()/50					// time step
				);

	    Feel::cout << "Pi1: \t" << M_statevar_solution[0] << std::endl;
	    Feel::cout << "Pi2: \t" << M_statevar_solution[1] << std::endl;
    	Feel::cout << "Pi3: \t" << M_statevar_solution[2] << std::endl;

		/*
       	auto help = expr("2*t:t");
        help.setParameterValues( { {"t", this->time() } } );
		M_statevar_solution[0] = help.evaluate()(0,0);
		*/

		// Update the initial solution for Pi1 (for step 1)
	    *M_Y = project ( _space = this->M_Ch, _expr = cst(M_statevar_solution[0]) );
	    M_bdf_statevariable -> setUnknown(0,*M_Y);

	    Feel::cout << "Value of P1 after second step : \t " << (*M_Y)[0] << std::endl;
    }
    else
    {
        Feel::cout << "No need to make the second step: no integral boundary condition coupled with a 0d circuit." << std::endl;
    }

    this->log("LaminaCribrosa","0D model", "finish");
	toc("0D model");

}


} // FeelModels

} // Feel


#endif	/* LC_MODEL_HPP */


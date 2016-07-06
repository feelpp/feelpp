/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel++ library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
             Daniele Prada <daniele.prada85@gmail.com>
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

#include "mixedpoisson.hpp"
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
	// solution
        dxdt[0] = -Cinv(0,0) * ( A(0,0)*x[0] + A(0,1)*x[1] + A(0,2)*x[2] );
        dxdt[1] = -Cinv(1,1) * ( A(1,0)*x[0] + A(1,1)*x[1] + A(1,2)*x[2] );
        dxdt[2] = Cinv(2,2) * g(2) - Cinv(2,2) * ( A(2,0)*x[0] + A(2,1)*x[1] + A(2,2)*x[2] );
    }
};

/*
inline
po::options_description
makeLCHDGOptions()
{
    po::options_description options ( "Lamina cribrosa options");
    options.add( makeMixedPoissonOptions("LC"));
    return options;
}

inline
po::options_description
makeLCHDGLibOptions()
{
    po::options_description options ( "Lamina cribrosa lib options");
    options.add( makeMixedPoissonLibOptions("LC"));
    return options.add( feel_options() );
}
*/


 
template<int Dim, int Order>
class LaminaCribrosa : public MixedPoisson<Dim,Order,1>
{
public:
    typedef MixedPoisson<Dim,Order,1> super_type;
    
    typedef LaminaCribrosa<Dim,Order> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;
    
    
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
    
    typedef Bdf <Ch_t> statevar_bdf_type;
    typedef boost::shared_ptr<statevar_bdf_type> statevar_bdf_ptrtype;


private:
    
    Ch_element_ptr_t M_Y;
    statevar_bdf_ptrtype M_bdf_statevariable;
    state_type M_statevar_solution;
    boost::numeric::ublas::matrix<value_type> M_A0d ;
    boost::numeric::ublas::matrix<value_type> M_Cinv ;
    boost::numeric::ublas::vector<value_type> M_g ;

public:
   
    LaminaCribrosa() : super_type() {  }


    virtual void initModel();
    virtual void initSpaces();
    virtual void initGraphs(int extraRow, int extraCol);
    virtual void initExporter(mesh_ptrtype meshVisu = nullptr);
    virtual void assembleA();
    virtual void assembleF();
    virtual void exportResults( double time, mesh_ptrtype mesh = nullptr, op_interp_ptrtype Idh = nullptr, opv_interp_ptrtype Idhv = nullptr );
    void exportResults(mesh_ptrtype mesh = nullptr, op_interp_ptrtype Idh = nullptr, opv_interp_ptrtype Idhv = nullptr) 
	{ 
	   this->exportResults (this->currentTime(), mesh, Idh, Idhv ); 
	   this->exporterMP()->save(); 
	}
    
       
    // time step scheme
    virtual void createTimeDiscretization() ;
    virtual void updateTimeStepBDF();
    virtual void initTimeStep();    
    statevar_bdf_ptrtype timeStepBDF_statevar() { return M_bdf_statevariable; }
    statevar_bdf_ptrtype const& timeStepBDF_statevar() const { return M_bdf_statevariable; }
    boost::shared_ptr<TSBase> timeStepBase_statevar() { return this->timeStepBDF_statevar(); }
    boost::shared_ptr<TSBase> timeStepBase_statevar() const { return this->timeStepBDF_statevar(); }
   
    // For the second step 
    void second_step();
    // void ode_model(const state_type x , state_type &dxdt , const value_type t );

};

template<int Dim, int Order> 
void 
LaminaCribrosa<Dim, Order>::initTimeStep()
{
    super_type::initTimeStep();
    // start or restart time step scheme
    if (!this->doRestart())
    {
        // start time step
        M_bdf_statevariable -> start( *M_Y );

        // up current time
        this->updateTime( M_bdf_statevariable -> time() );
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
        this->updateTime( M_bdf_statevariable->time() );

        this->log("LaminaCribrosa","initTimeStep", "restart bdf/exporter done" );
    }
    
}


template<int Dim, int Order>
void 
LaminaCribrosa<Dim, Order>::updateTimeStepBDF()
{
    super_type::updateTimeStepBDF();

    this->log("LaminaCribrosa","updateTimeStepBDF", "start" );
    this->timerTool("TimeStepping").setAdditionalParameter("time_Y",this->currentTime());
    this->timerTool("TimeStepping").start();

    int previousTimeOrder_statevar = this->timeStepBDF_statevar()->timeOrder();
    
     
    M_bdf_statevariable->next( *M_Y );

    int currentTimeOrder_statevar = this->timeStepBDF_statevar()->timeOrder();

    this->updateTime( M_bdf_statevariable->time() );

    this->timerTool("TimeStepping").stop("updateBdf");
    if ( this->scalabilitySave() ) this->timerTool("TimeStepping").save();
    this->log("LaminaCribrosa","updateTimeStepBDF", "finish" );
}


template<int Dim, int Order>
void
LaminaCribrosa<Dim,Order>::createTimeDiscretization()
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
template<int Dim, int Order>
void
LaminaCribrosa<Dim, Order>::initModel(){
    super_type::initModel();
    if (!this->integralCondition()){
        Feel::cout << std::endl << "ERROR Lamina Cribrosa: no integral conditions found" << std::endl << std::endl;
    }
    
    M_A0d.resize(3,3);
    M_Cinv.resize(3,3);
    M_g.resize(3);

    for( auto const& pairMat : this->M_modelProperties->materials() )
    {
        auto material = pairMat.second;
        auto Piout = material.getDouble("Piout");  
	auto C1 = material.getDouble("C1");
	auto C2 = material.getDouble("C2");
	auto C3 = material.getDouble("C3");
    	auto R12 = material.getDouble("R12");
	auto R23 = material.getDouble("R23");
	auto Rout = material.getDouble("Rout");
    
    	// Initialize matrices and vector of the ODE
    	M_A0d(0,0) = 1/R12;
    	M_A0d(0,1) = -1/R12;
    	M_A0d(1,0) = -1/R12;
    	M_A0d(1,1) = 1/R12 + 1/R23;
    	M_A0d(1,2) = -1/R23;
    	M_A0d(2,1) = -1/R23;
    	M_A0d(2,2) = 1/R23 + 1/Rout;

    	M_Cinv(0,0) = 1/C1; 
    	M_Cinv(1,1) = 1/C2; 
    	M_Cinv(2,2) = 1/C3; 

    	M_g(2) = Piout/Rout;
    }

}

// Overriding of init Spaces
template<int Dim, int Order>
void
LaminaCribrosa<Dim, Order>::initSpaces(){
    super_type::initSpaces();

    M_Y = this->constantSpace()->elementPtr( "yy" );
    M_statevar_solution.fill(0); // initializtion 
}


// Overriding of init Exporter
template<int Dim, int Order>
void
LaminaCribrosa<Dim, Order>::initExporter(mesh_ptrtype meshVisu)  {  
    super_type::initExporter(meshVisu);
    
}

// Overriding of the graph method
template<int Dim, int Order>
void
LaminaCribrosa<Dim, Order>::initGraphs(int extraRow, int extraCol)
{
    super_type::initGraphs(extraRow,extraCol);

 
    auto Vh = this->fluxSpace();
    auto Wh = this->potentialSpace();
    auto Mh = this->traceSpace();
    auto Ch = this->constantSpace();

    this->M_hdg_graph(4,0) = stencil( _test=Ch, _trial=Vh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    this->M_hdg_graph(4,1) = stencil( _test=Ch, _trial=Wh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    this->M_hdg_graph(4,2) = stencil( _test=Ch, _trial=Mh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    this->M_hdg_graph(4,3) = stencil( _test=Ch, _trial=Ch, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();

    this->M_hdg_graph(0,4) = stencil( _test=Vh, _trial=Ch, _diag_is_nonzero=false, _close=false)->graph();
    this->M_hdg_graph(1,4) = stencil( _test=Wh, _trial=Ch, _diag_is_nonzero=false, _close=false)->graph();
    this->M_hdg_graph(2,4) = stencil( _test=Mh, _trial=Ch, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    this->M_hdg_graph(3,4) = stencil( _test=Ch, _trial=Ch, _diag_is_nonzero=false, _close=false)->graph();
    this->M_hdg_graph(4,4) = stencil( _test=Ch, _trial=Ch, _diag_is_nonzero=false, _close=false)->graph();

    this->M_hdg_vec(4,0) = this->get_backend()->newVector( Ch );
    this->M_hdg_sol(4,0) = M_Y;
    
}

template<int Dim, int Order>
void
LaminaCribrosa<Dim, Order>::assembleA( ) 
{
    super_type::assembleA();
    
    auto u = this->fluxSpace()->element( "u" );
    auto p = this->potentialSpace()->element( "p" );

    auto mu2 = this->constantSpace()->element("mu2");
    auto mu3 = this->constantSpace()->element( "mu3" );
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

    // Add extra equation for the coupling 
    auto a15 = form2(_trial=this->constantSpace(), _test=this->fluxSpace(),_matrix=this->M_A_cst,
                     _rowstart=0,
                     _colstart=4);
    auto a25 = form2(_trial=this->constantSpace(), _test=this->potentialSpace(),_matrix=this->M_A_cst,
                     _rowstart=1,
                     _colstart=4);    
    auto a35 = form2(_trial=this->constantSpace(), _test=this->traceSpace(),_matrix=this->M_A_cst,
                     _rowstart=2,
                     _colstart=4);
    auto a44 = form2(_trial=this->constantSpace(), _test=this->constantSpace(),_matrix=this->M_A_cst,
                     _rowstart=3,
                     _colstart=3);
    auto a45 = form2(_trial=this->constantSpace(), _test=this->constantSpace(),_matrix=this->M_A_cst,
                     _rowstart=3,
                     _colstart=4);  
    auto a51 = form2(_trial=this->fluxSpace(), _test=this->constantSpace(),_matrix=this->M_A_cst,
                     _rowstart=4,
                     _colstart=0);
    auto a52 = form2(_trial=this->potentialSpace(), _test=this->constantSpace(),_matrix=this->M_A_cst,
                     _rowstart=4,
                     _colstart=1);
    auto a53 = form2(_trial=this->traceSpace(), _test=this->constantSpace(),_matrix=this->M_A_cst,
                     _rowstart=4,
                     _colstart=2);
    auto a54 = form2(_trial=this->constantSpace(), _test=this->constantSpace(),_matrix=this->M_A_cst,
                     _rowstart=4,
                     _colstart=3);  
    auto a55 = form2(_trial=this->constantSpace(), _test=this->constantSpace(),_matrix=this->M_A_cst,
                     _rowstart=4,
                     _colstart=4);

    double meas = 0;
    for( auto marker : this->M_integralMarkersList)
    {
	meas += integrate(_range=markedfaces(this->mesh(),marker),_expr=cst(1.0)).evaluate()(0,0);
    }
 
    auto itField = this->M_modelProperties->boundaryConditions().find( "flux");
    if ( itField != this->M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Integral" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {

            std::string marker = exAtMarker.marker();
	    	    
            // - <j.n, mu3>_Gamma_I
            a51 += integrate( _range=markedfaces(this->mesh(),marker), _expr= -trans(idt(u))*N()*id(mu3) );

            // - <tau p, mu3>_Gamma_I
            a52 += integrate( _range=markedfaces(this->mesh(),marker), _expr= -tau_constant*( pow(idv(H),this->tau_order())*idt(p) )*id(mu3) );
            
            // + <tau u_I, mu3>_Gamma_I
            a54 += integrate( _range=markedfaces(this->mesh(),marker), _expr= tau_constant * (pow(idv(H),this->tau_order())*idt(uI)*id(mu3)) );
            
            for( auto const& pairMat : this->M_modelProperties->materials() )
            {
            	auto material = pairMat.second;
        	auto RR = material.getScalar("RR");       // Resistence of the buffer
        	auto CC = material.getScalar("CC");       // Capacitance of the buffer
	        
		// -1/(R |Gamma_I|) <u_I, mu2>_Gamma_I
                a44 += integrate( _range=markedfaces(this->mesh(),marker), _expr = - idt(uI)*id(mu2)/RR/meas ) ;
		// +1/(R |Gamma_I|) <Y,mu2>_Gamma_I
		a45 += integrate( _range=markedfaces(this->mesh(),marker), _expr = idt(yy)*id(mu2)/RR/meas ) ;

            	// < C/|Gamma_I| Y/dt, mu3>_Gamma_I
            	a55 += integrate(_range=markedfaces(this->mesh(),marker), 
				 _expr= CC*this->timeStepBDF_statevar()->polyDerivCoefficient(0) * idt(yy)*id(mu3)/meas );

	    }
            }
	}    
    }    
}


template<int Dim, int Order>
void
LaminaCribrosa<Dim, Order>::assembleF()
{
    super_type::assembleF();

    auto mu3 = this->constantSpace()->element( "mu3" );
    int RowStart = this->integralCondition() ? 4 : 3;
 
    auto rhs5 = form1( _test=this->constantSpace(), _vector=this->M_F, 
						   _rowstart=RowStart );

    double meas = 0;
    for( auto marker : this->M_integralMarkersList)
    {
	meas += integrate(_range=markedfaces(this->mesh(),marker),_expr=cst(1.0)).evaluate()(0,0);
    }
    
    for( auto const& pairMat : this->M_modelProperties->materials() )
    {
        auto material = pairMat.second;
        auto CC = material.getScalar("CC");       // Capacitance of the buffer
	for (auto marker : this->M_integralMarkersList)
	{
             // < C/|Gamma_I| Yold/dt, mu3>
	     rhs5 += integrate( _range = markedfaces(this->mesh(),marker),
                               _expr = CC*idv(this->timeStepBDF_statevar()->polyDeriv()) * id(mu3)/meas);
	}
    }		
}

template <int Dim, int Order>
void
LaminaCribrosa<Dim,Order>::exportResults( double time, mesh_ptrtype mesh, op_interp_ptrtype Idh, opv_interp_ptrtype Idhv)
{
    super_type::exportResults( time, mesh, Idh, Idhv );
    this->log("LaminaCribrosa","exportResults", "start");

    
    // Export computed solutions
    auto postProcess = this->modelProperties().postProcess();
    auto itField = postProcess.find( "Fields");
    if ( itField != postProcess.end() )
    {
        for ( auto const& field : (*itField).second )
        {
            if ( field == "state variable" )
	    {
		this->exporterMP()->step( time )->add(prefixvm(this->prefix(), "Pi_1"), M_statevar_solution[0] );
		this->exporterMP()->step( time )->add(prefixvm(this->prefix(), "Pi_2"), M_statevar_solution[1] );
		this->exporterMP()->step( time )->add(prefixvm(this->prefix(), "Pi_3"), M_statevar_solution[2] );
		
		for( auto const& pairMat : this->M_modelProperties->materials() )
                {
                   auto material = pairMat.second;
                   auto P1_exact = material.getDouble( "P1_exact" ); 
                   auto P2_exact = material.getDouble( "P2_exact" );
                   auto P3_exact = material.getDouble( "P3_exact" );
		   						
    		   Feel::cout << "||P1-P1_ex|=\t" << std::abs(P1_exact - M_statevar_solution[0]) << std::endl;
    		   Feel::cout << "||P2-P2_ex|=\t" << std::abs(P2_exact - M_statevar_solution[1]) << std::endl;
    		   Feel::cout << "||P3-P3_ex|=\t" << std::abs(P3_exact - M_statevar_solution[2]) << std::endl;
		   Feel::cout << "---------------------------" << std::endl;
                }
	    }
        }
    }
   
    this->log("LaminaCribrosa","exportResults", "finish");
}

template<int Dim, int Order>
void
LaminaCribrosa<Dim, Order>::second_step(){

	this->log("LaminaCribrosa","0D model", "start");
	tic();

	using namespace boost::numeric::odeint;
	using namespace boost::numeric::ublas;

	// Update the initial solution for Pi1 (for step 2)  
	M_statevar_solution[0] = mean( _range= elements(this->mesh()), _expr=idv(*M_Y) )(0,0) ;
	M_statevar_solution[0] = (*M_Y)[0];
        	
	Feel::cout << "Value of P1 after first step : \t " << (*M_Y)[0] << std::endl;

        double j_integral = 0;
        for( auto marker : this->M_integralMarkersList)
        {
            j_integral += integrate(_range=markedfaces(this->mesh(),marker),_expr=trans(idv(this->M_up))*N()).evaluate()(0,0);
        }
	
	Feel::cout << "Integral value of the flow: \t " << j_integral << std::endl;

	// solve the problem
	boost::numeric::odeint::integrate(ode_model(M_Cinv,M_A0d,M_g), M_statevar_solution, 
			M_bdf_statevariable->time(),                      		// initial time
			M_bdf_statevariable->time()+M_bdf_statevariable->timeStep(), 	// final time
			M_bdf_statevariable->timeStep()/100 ); 				// time step
	Feel::cout << "Pi1: \t" << M_statevar_solution[0] << std::endl;
	Feel::cout << "Pi2: \t" << M_statevar_solution[1] << std::endl;
	Feel::cout << "Pi3: \t" << M_statevar_solution[2] << std::endl;

	// Update the initial solution for Pi1 (for step 1)
	*M_Y = project ( _space = this->M_Ch, _expr = cst(M_statevar_solution[0]) );
	M_bdf_statevariable -> setUnknown(0,*M_Y);
  
	Feel::cout << "Value of P1 after second step : \t " << (*M_Y)[0] << std::endl;

	this->log("LaminaCribrosa","0D model", "finish");
	toc("0D model");
}


} // FeelModels

} // Feel


#endif	/* LC_MODEL_HPP */


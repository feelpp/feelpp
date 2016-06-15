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

namespace Feel {

namespace FeelModels {
    
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
    


public:
    /*
    LaminaCribrosa() : M_LcModel("mixedpoisson") {
        M_mesh = loadMesh( new mesh_type );
        this->init(M_mesh);
    }*/ 
    
    LaminaCribrosa() : super_type() {
	// with one extra row and 1 extra column
        this -> init(this->mesh(),1,1);
    }


    virtual void initModel();
    virtual void initSpaces();
    virtual void initGraphs(int extraRow, int extraCol);
    virtual void initExporter();
    virtual void assembleA();
    virtual void assembleF();
    virtual void exportResults(double time);
    void exportResults() { this->exportResults (this->currentTime() ); this->exporterMP()->save(); }
    
    // void updateLcAssembly();
    void run();
    // void updateLcAssembly( sparse_matrix_ptrtype& A, vector_ptrtype& F) const;
        
    // time step scheme
    virtual void createTimeDiscretization() ;
    virtual void updateTimeStepBDF();
    virtual void initTimeStep();    
    statevar_bdf_ptrtype timeStepBDF_statevar() { return M_bdf_statevariable; }
    statevar_bdf_ptrtype const& timeStepBDF_statevar() const { return M_bdf_statevariable; }
    boost::shared_ptr<TSBase> timeStepBase_statevar() { return this->timeStepBDF(); }
    boost::shared_ptr<TSBase> timeStepBase_statevar() const { return this->timeStepBDF(); }
    
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

    // int previousTimeOrder = this->timeStepBDF()->timeOrder();
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
}

// Overriding of init Spaces
template<int Dim, int Order>
void
LaminaCribrosa<Dim, Order>::initSpaces(){
    super_type::initSpaces();

    M_Y = this->constantSpace()->elementPtr( "yy" );
}


// Overriding of init Exporter
template<int Dim, int Order>
void
LaminaCribrosa<Dim, Order>::initExporter()  {  
    super_type::initExporter();
    
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
    this->M_hdg_graph(4,3) = stencil( _test=Ch, _trial=Ch, _diag_is_nonzero=false, _close=false)->graph();
    
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

    Feel::cout << __LINE__ << std::endl;

    
    auto itField = this->M_modelProperties->boundaryConditions().find( "flux");
    if ( itField != this->M_modelProperties->boundaryConditions().end() )
    {
        Feel::cout << __LINE__ << std::endl; 
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Integral" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
		Feel::cout << __LINE__ << std::endl;

            std::string marker = exAtMarker.marker();
            Feel::cout << __LINE__ << std::endl;
	    double meas = integrate( _range=markedfaces(this->mesh(),marker), _expr=cst(1.0)).evaluate()(0,0);
    	    // -<tau uI, mu2>_Gamma_I
            a44 += integrate( _range=markedfaces(this->mesh(),marker), _expr=-pow(idv(H),this->tau_order())*idt(uI)*id(mu2) );
        
            // - <j.n, mu3>_Gamma_I
            a51 += integrate( _range=markedfaces(this->mesh(),marker), _expr= -trans(idt(u))*N()*id(mu3) );

	    Feel::cout << __LINE__ << std::endl;

            // - <tau p, mu3>_Gamma_I
            a52 += integrate( _range=markedfaces(this->mesh(),marker), _expr= -tau_constant*( pow(idv(H),this->tau_order())*idt(p) )*id(mu3) );
            
            if (this->integralCondition()){
                // + <tau u_I, mu3>_Gamma_I
                a54 += integrate( _range=markedfaces(this->mesh(),marker), _expr= pow(idv(H),this->tau_order())*id(mu3)*idt(uI) );
            }

		Feel::cout << __LINE__ << std::endl;

            for( auto const& pairMat : this->M_modelProperties->materials() )
            {
            	auto material = pairMat.second;
        	auto RR = material.getScalar("RR");       // Resistence of the buffer
        	auto CC = material.getScalar("CC");       // Capacitance of the buffer
	        // -1/(R |Gamma_I|) <u_I, mu2>_Gamma_I
                a44 += integrate( _range=markedfaces(this->mesh(),marker), _expr = - idt(uI)*id(mu2)/(RR*meas) ) ;
		// +1/(R |Gamma_I|) <Y,mu2>_Gamma_I
		a45 += integrate( _range=markedfaces(this->mesh(),marker), _expr = idt(yy)*id(mu2)/(RR*meas) ) ;
            	// < C/|Gamma_I| Y/dt, mu3>_Gamma_I
            	a55 += integrate( _range=markedfaces(this->mesh(),marker), _expr= CC/meas*this->timeStepBDF_statevar()->polyDerivCoefficient(0)*idt(yy)*id(mu3) );
	    }
Feel::cout << __LINE__ << std::endl;             
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
						   _rowstart=RowStart);
						   
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
                double meas = integrate( _range=markedfaces(this->mesh(),marker), _expr=cst(1.0)).evaluate()(0,0);
		for( auto const& pairMat : this->M_modelProperties->materials() )
    		{
        	    auto material = pairMat.second;
                    auto CC = material.getScalar("CC");       // Capacitance of the buffer
                    // < C/|Gamma_I| Yold/dt, mu3>
		    rhs5 += integrate( _range = markedelements(this->mesh(),marker),
                                       _expr = CC/meas*idv(this->timeStepBDF_statevar()->polyDeriv()) * id(mu3));
		}
            }
        }
    }
}

template <int Dim, int Order>
void
LaminaCribrosa<Dim,Order>::exportResults( double time )
{
    super_type::exportResults( time );
    this->log("LaminaCribrosa","exportResults", "start");
    this->timerTool("PostProcessing").start();
   
    // Export computed solutions
    auto postProcess = this->modelProperties().postProcess();
    auto itField = postProcess.find( "Fields");
    if ( itField != postProcess.end() )
    {
        for ( auto const& field : (*itField).second )
        {
            if ( field == "state variable" )
                this->exporterMP()->step( time )->add(prefixvm(this->prefix(), "state variable"), *M_Y);
        }
    }
   
    this->log("LaminaCribrosa","exportResults", "finish");
}

template<int Dim, int Order>
void
LaminaCribrosa<Dim, Order>::run()
{
    auto ModelProp = this->modelProperties();
    
    
    for ( ; !( this->timeStepBase()->isFinished() && this->timeStepBase_statevar()->isFinished() ) ; this -> updateTimeStep() ) { // start time cycle
   
        Feel::cout << "============================================================" << std::endl;
        Feel::cout << "time simulation: \t" << this->time() << "s " << std::endl;
        Feel::cout << "============================================================" << std::endl;
    
        // this->M_updateAssembly = boost::bind( &LaminaCribrosa<Dim, Order>::updateLcAssembly, this, _1, _2 );
        

        this->solve();
        //M_J = M_LcModel.fluxField();
        //M_u = M_LcModel.potentialField();
        this->exportResults(); 
        
    } // end time cycle
}


} // FeelModels

} // Feel


#endif	/* LC_MODEL_HPP */


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

#ifndef COUPLING_HPP
#define COUPLING_HPP

#include <feel/feelfmi/fmu.hpp>
#include <feel/feelmodels/hdg/mixedpoisson.hpp>

namespace Feel
{

namespace FeelModels
{

using value_type = double;

inline po::options_description
makeCoupledMixedPoissonOptions( std::string const& _prefix = "", std::string const& _toolbox_prefix = "hdg.poisson" )
{
    std::string prefix = _toolbox_prefix.empty() ? "hdg.poisson" : _toolbox_prefix;
    if ( !_prefix.empty() )
        prefix = prefixvm( prefix, _prefix );
    po::options_description cmpOptions( "Coupled Mixed Poisson HDG options" );
    cmpOptions.add_options()( "gmsh.submesh", po::value<std::string>()->default_value( "" ), "submesh extraction" )
        ( prefixvm( prefix, "tau_constant" ).c_str(), po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( prefixvm( prefix, "tau_order" ).c_str(), po::value<int>()->default_value( 0 ), "order of the stabilization function on the selected edges" ) // -1, 0, 1 ==> h^-1, h^0, h^1
        ( prefixvm( prefix, "hface" ).c_str(), po::value<int>()->default_value( 0 ), "hface" )
        ( prefixvm( prefix, "conductivity_json" ).c_str(), po::value<std::string>()->default_value( "cond" ), "key for conductivity in json" )
        ( prefixvm( prefix, "conductivityNL_json" ).c_str(), po::value<std::string>()->default_value( "condNL" ), "key for non linear conductivity in json (depends on potential p)" )
        ( prefixvm( prefix, "use-sc" ).c_str(), po::value<bool>()->default_value( true ), "use static condensation" )
        ( prefixvm( prefix, "error-quadrature" ).c_str(), po::value<int>()->default_value( 10 ), "quadrature to compute errors" )
        ( prefixvm( prefix, "set-zero-by-init" ).c_str(), po::value<bool>()->default_value( true ), "reinit matrix and vector when setting to zero" )
        ( prefixvm( "coupling", "Cbuffer_name" ).c_str(), po::value<std::string>()->default_value( "" ), "Name of the C buffer in the circuit" )
        ( prefixvm( "coupling", "var_buffer" ).c_str(), po::value<std::string>()->default_value( "" ), "Name of the buffer variable in the circuit" )
        ( prefixvm( "coupling", "Rbuffer_name" ).c_str(), po::value<std::string>()->default_value( "" ), "Name of the R buffer in the circuit" );
    cmpOptions.add( modelnumerical_options( prefix ) );
    cmpOptions.add( backend_options( prefix + ".sc" ) );
    return cmpOptions;
}

inline po::options_description
makeCoupledMixedPoissonLibOptions( std::string const& prefix = "", std::string const& _toolbox_prefix = "hdg.poisson" )
{
    po::options_description cmpLibOptions( "Coupled Mixed Poisson HDG Lib options" );
    // if ( !prefix.empty() )
    //     mpLibOptions.add( backend_options( prefix ) );
    return cmpLibOptions;
}

template <int Dim, int Order, int G_Order = 1, int E_Order = 4>
class CoupledMixedPoisson : public MixedPoisson<Dim, Order, G_Order, E_Order>
{

  public:
    typedef MixedPoisson<Dim, Order, G_Order, E_Order> super_type;

    typedef CoupledMixedPoisson<Dim, Order, G_Order, E_Order> self_type;
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

    typedef Bdf<Ch_t> buffer_bdf_type;
    typedef std::shared_ptr<buffer_bdf_type> buffer_bdf_ptrtype;
    using Bdf_element_vector_t = typename std::vector<buffer_bdf_ptrtype>;

    using product2_space_type = typename super_type::product2_space_type;
    using integral_boundary_list_type = typename super_type::integral_boundary_list_type;

  private:
    // Ch_element_ptr_t M_Y;
    Ch_element_vector_type M_Y;

    //buffer_bdf_ptrtype M_bdf_buffer;
    Bdf_element_vector_t M_bdf_buffer;

    int M_0dCondition;
    integral_boundary_list_type M_0dList;
    std::shared_ptr<FMU> M_circuit;

  public:
    CoupledMixedPoisson( std::string const& prefix = "hdg.poisson",
                         worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                         std::string const& subPrefix = "",
                         ModelBaseRepository const& modelRep = ModelBaseRepository() )
        : super_type( prefix, MixedPoissonPhysics::None, worldComm, subPrefix, modelRep ),
          ModelBase( prefix, MixedPoissonPhysicsMap[MixedPoissonPhysics::None]["keyword"], worldComm, subPrefix, modelRep )
        {}

    static self_ptrtype New( std::string const& prefix = "hdg.poisson",
                             worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                             std::string const& subPrefix = "" );

    std::shared_ptr<FMU> circuit() { return M_circuit; }

    virtual void assembleCstPart();
    virtual void assembleNonCstPart();
    void assemble0d( int i );
    void assembleRhs0d( int i );

    virtual void initModel();
    virtual void initSpaces();

    void exportResults( double time, mesh_ptrtype mesh = nullptr, op_interp_ptrtype Idh = nullptr, opv_interp_ptrtype Idhv = nullptr );
    void exportResults( mesh_ptrtype mesh = nullptr, op_interp_ptrtype Idh = nullptr, opv_interp_ptrtype Idhv = nullptr )
    {
        this->exportResults( this->currentTime(), mesh, Idh, Idhv );
        this->exporterMP()->save();
    }

    virtual void solve();

    // time step scheme
    virtual void createTimeDiscretization();
    virtual void updateTimeStepBDF();
    virtual void initTimeStep();

    Bdf_element_vector_t timeStepBDF_buffer() { return M_bdf_buffer; }
    Bdf_element_vector_t const& timeStepBDF_buffer() const { return M_bdf_buffer; }
    // boost::shared_ptr<TSBase> timeStepBase_buffer() { return this->timeStepBDF_buffer(); }
    // boost::shared_ptr<TSBase> timeStepBase_buffer() const { return this->timeStepBDF_buffer(); }

    // For the second step
    void second_step( double final_time, int i );

    // Run simulation
    void run( op_interp_ptrtype Idh_poi = nullptr, opv_interp_ptrtype Idhv_poi = nullptr );
};

template <int Dim, int Order, int G_Order, int E_Order>
typename CoupledMixedPoisson<Dim, Order, G_Order, E_Order>::self_ptrtype
CoupledMixedPoisson<Dim, Order, G_Order, E_Order>::New( std::string const& prefix,
                                                        worldcomm_ptr_t const& worldComm,
                                                        std::string const& subPrefix )
{
    return std::make_shared<self_type>( prefix, worldComm, subPrefix );
}

template <int Dim, int Order, int G_Order, int E_Order>
void CoupledMixedPoisson<Dim, Order, G_Order, E_Order>::initTimeStep()
{
    super_type::initTimeStep();
    // start or restart time step scheme
    if ( !this->doRestart() )
    {
        auto itField = this->modelProperties().boundaryConditions().find( "buffer" );
        if ( itField != this->modelProperties().boundaryConditions().end() )
        {
            for ( int i = 0; i < M_0dList.size(); i++ )
            {
                auto mapField = ( *itField ).second;
                auto itType = mapField.find( "InitialSolution" );
                if ( itType != mapField.end() )
                {
                    for ( auto const& exAtMarker : ( *itType ).second )
                    {
                        auto marker = M_0dList[i].marker();
                        marker += "-buffer";
                        auto buffer = expr( exAtMarker.expression() );
                        if ( exAtMarker.marker() == marker )
                        {
                            M_Y[i] = project( _space = this->M_Ch, _expr = buffer );
                            Feel::cout << "Initial buffer condition value of potential on " << marker << " : \t " << buffer << std::endl;

                            for ( auto time : M_bdf_buffer[i]->priorTimes() )
                            {
                                buffer.setParameterValues( {{"t", time.second}} );
                                auto buffer_e = project( _space = this->M_Ch, _expr = buffer );
                                M_bdf_buffer[i]->setUnknown( time.first, buffer_e );
                            }
                        }
                    }
                }
            }
        }

        // start time step
        for ( int i = 0; i < M_0dCondition; i++ )
            M_bdf_buffer[i]->start();

        // Feel::cout << "Pressure on buffer before 2nd step: \t " << *M_Y << std::endl;
        // up current time
        this->updateTime( M_bdf_buffer[0]->time() );
    }
    else
    {
        for ( int i = 0; i < M_0dCondition; i++ )
        {
            // start time step
            M_bdf_buffer[i]->restart();
            // load a previous solution as current solution
            M_Y[i] = M_bdf_buffer[i]->unknown( 0 );
            // up initial time
            this->setTimeInitial( M_bdf_buffer[i]->timeInitial() );
            // up current time
            this->updateTime( M_bdf_buffer[0]->time() );
        }
        this->log( "CoupledMixedPoisson", "initTimeStep", "restart bdf/exporter done" );
    }
}

template <int Dim, int Order, int G_Order, int E_Order>
void CoupledMixedPoisson<Dim, Order, G_Order, E_Order>::updateTimeStepBDF()
{
    super_type::updateTimeStepBDF();

    this->log( "CoupledMixedPoisson", "updateTimeStepBDF", "start" );
    this->timerTool( "TimeStepping" ).setAdditionalParameter( "time_Y", this->currentTime() );
    this->timerTool( "TimeStepping" ).start();

    int previousTimeOrder_buffer = this->timeStepBDF_buffer()[0]->timeOrder();

    for ( int i = 0; i < M_0dCondition; i++ )
        M_bdf_buffer[i]->next( M_Y[i] );

    int currentTimeOrder_buffer = this->timeStepBDF_buffer()[0]->timeOrder();

    this->updateTime( M_bdf_buffer[0]->time() );

    this->timerTool( "TimeStepping" ).stop( "updateBdf" );
    if ( this->scalabilitySave() ) this->timerTool( "TimeStepping" ).save();
    this->log( "CoupledMixedPoisson", "updateTimeStepBDF", "finish" );
}

template <int Dim, int Order, int G_Order, int E_Order>
void CoupledMixedPoisson<Dim, Order, G_Order, E_Order>::createTimeDiscretization()
{
    super_type::createTimeDiscretization();

    this->log( "CoupledMixedPoisson", "createTimeDiscretization", "start" );
    this->timerTool( "Constructor" ).start();

    std::string myFileFormat = soption( _name = "ts.file-format" ); // without prefix
    std::string suffixName = "";
    if ( myFileFormat == "binary" )
        suffixName = ( boost::format( "_rank%1%_%2%" ) % this->worldComm().rank() % this->worldComm().size() ).str();

    for ( int i = 0; i < M_0dCondition; i++ )
    {
        M_bdf_buffer[i] = bdf( _vm = Environment::vm(), _space = this->M_Ch,
                               _name = prefixvm( this->prefix(), prefixvm( this->subPrefix(), "Y" + suffixName ) ),
                               _prefix = "",
                               _initial_time = this->timeInitial(),
                               _final_time = this->timeFinal(),
                               _time_step = this->timeStep(),
                               _restart = this->doRestart(),
                               _restart_path = this->restartPath(),
                               _restart_at_last_save = this->restartAtLastSave(),
                               _save = this->tsSaveInFile(), _freq = this->tsSaveFreq() );
        M_bdf_buffer[i]->setfileFormat( myFileFormat );
        M_bdf_buffer[i]->setPathSave( ( fs::path( this->rootRepository() ) /
                                        fs::path( prefixvm( this->prefix(), ( boost::format( "bdfY_o_%1%_dt_%2%" ) % M_bdf_buffer[i]->bdfOrder() % this->timeStep() ).str() ) ) )
                                          .string() );
    }

    double tElapsed = this->timerTool( "Constructor" ).stop( "createTimeDiscr" );
    this->log( "CoupledMixedPoisson", "createTimeDiscretization", ( boost::format( "finish in %1% s" ) % tElapsed ).str() );
}

// Overriding of init Model
template <int Dim, int Order, int G_Order, int E_Order>
void CoupledMixedPoisson<Dim, Order, G_Order, E_Order>::initModel()
{

    super_type::initModel();

    M_0dList.clear();
    auto itField = this->modelProperties().boundaryConditions().find( "flux" );
    if ( itField != this->modelProperties().boundaryConditions().end() )
    {
        auto mapField = ( *itField ).second;
        auto itType = mapField.find( "Integral_coupled_with_0d" );
        if ( itType != mapField.end() )
        {
            Feel::cout << "Integral coupled with 0d equation:";
            for ( auto const& exAtMarker : ( *itType ).second )
            {
                std::string marker = exAtMarker.marker();
                if ( this->mesh()->hasFaceMarker( marker ) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl
                               << "WARNING!! marker " << marker << "does not exist!" << std::endl;
                this->M_IBCList.push_back( exAtMarker );
                M_0dList.push_back( exAtMarker );
            }
            Feel::cout << std::endl;
        }
    }

    /*
	if ( this->M_IBCList.empty() )
    	this->M_integralCondition = 0;
   	else
       	this->M_integralCondition = this->M_IBCList.size();	
	*/

    if ( M_0dList.empty() )
        M_0dCondition = 0;
    else
    {
        M_0dCondition = M_0dList.size();
        this->M_integralCondition = this->M_IBCList.size();
    }
}

// Overriding of init Spaces
template <int Dim, int Order, int G_Order, int E_Order>
void CoupledMixedPoisson<Dim, Order, G_Order, E_Order>::initSpaces()
{

    // for (int i = 0; i < M_0dCondition; i++)
    // 	this->M_IBCList.push_back(M_0dList[i]);

    super_type::initSpaces();

    // for (int i = 0; i < M_0dCondition; i++)
    //	this->M_IBCList.pop_back();

    if ( M_0dCondition )
    {
        Feel::cout << "Number of 0d: " << M_0dCondition << std::endl;
        /*/ Mh only on the faces whitout 0d condition
    	auto complement_integral_bdy = complement(faces(M_mesh),[this]( auto const& ewrap ) {
        	auto const& e = unwrap_ref( ewrap );
            for( auto exAtMarker : this->M_0dList)
            {
                if ( e.hasMarker() && e.marker().value() == this->M_mesh->markerName( exAtMarker.marker() ) )
                    return true;
            }
            return false; });
		/
    	auto face_mesh = createSubmesh( _mesh=M_mesh, _range=complement_integral_bdy, _update=0 );

    	this->M_Mh = Pdh<Order>( face_mesh, true );
    	this->M_M0h = Pdh<0>( face_mesh );

	    std::vector<std::string> 0d_markers(M_0dCondition);
    	for( int i = 0; i < M_0dCondition; i++)
     	{
         	0d_markers.push_back(M_0dList[i].marker());
     	}

     	auto 0d_mesh = createSubmesh( _mesh=this->mesh(), _range=markedfaces(this->mesh(), 0d_markers), _update=0 );
    	this->M_Ch = Pch<0>( 0d_mesh, true );
		*/
        // auto zeroSpaces = boost::make_shared<ProductSpace<Ch_ptr_t,true> >( 2*M_0dCondition, this->M_Ch);
        auto zeroSpaces = std::make_shared<ProductSpace<Ch_ptr_t, true>>( this->integralCondition() + M_0dCondition, this->M_Ch );
        this->M_ps = std::make_shared<product2_space_type>( product2( zeroSpaces, this->M_Vh, this->M_Wh, this->M_Mh ) );

        solve::strategy s = boption( prefixvm( this->prefix(), "use-sc" ) ) ? solve::strategy::static_condensation : solve::strategy::monolithic;

        this->M_A_cst = makeSharedMatrixCondensed<value_type>( s, csrGraphBlocks( *( this->getPS() ) ), *( this->get_backend() ) );
#ifndef USE_SAME_MAT
        this->M_A = makeSharedMatrixCondensed<value_type>( s, csrGraphBlocks( *( this->getPS() ) ), *( this->get_backend() ) );
#endif
        this->M_F = makeSharedVectorCondensed<value_type>( s, blockVector( *( this->getPS() ) ), *( this->get_backend() ), false );

        M_bdf_buffer.resize( M_0dCondition );
    }

    // Initialization of the buffer
    for ( int i = 0; i < M_0dCondition; i++ )
        M_Y.push_back( this->constantSpace()->element( "yy" ) );

    // Initialization for second step: the 0d equation
    if ( M_0dCondition )
    {
        // Environment::setOptionValue<std::string>("fmu.filename");
        M_circuit = std::make_shared<FMU>();

        M_circuit->load();
        M_circuit->initialize( this->timeInitial(), this->timeFinal() );
    }
}

template <int Dim, int Order, int G_Order, int E_Order>
void CoupledMixedPoisson<Dim, Order, G_Order, E_Order>::assembleCstPart()
{

    for ( int i = 0; i < M_0dCondition; i++ )
        this->M_IBCList.push_back( M_0dList[i] );

    super_type::assembleCstPart();

    for ( int i = 0; i < M_0dCondition; i++ )
        this->M_IBCList.pop_back();

    for ( int i = 0; i < M_0dList.size(); i++ )
        this->assemble0d( i );
}

template <int Dim, int Order, int G_Order, int E_Order>
void CoupledMixedPoisson<Dim, Order, G_Order, E_Order>::assembleNonCstPart()
{
    super_type::assembleNonCstPart();

    for ( int i = 0; i < M_0dList.size(); i++ )
        this->assembleRhs0d( i );
}

template <int Dim, int Order, int G_Order, int E_Order>
void CoupledMixedPoisson<Dim, Order, G_Order, E_Order>::assemble0d( int i )
{

    this->log( "CoupledMixedPoisson", "assembleMatrix0d", "start" );
    this->timerTool( "Constructor" ).start();

    auto bbf = blockform2( *( this->getPS() ), this->M_A_cst );

    auto u = this->fluxSpace()->element( "u" );
    auto p = this->potentialSpace()->element( "p" );
    auto w = this->potentialSpace()->element( "w" );
    auto nu = this->constantSpace()->element( "nu" );
    auto uI = this->constantSpace()->element( "uI" );
    auto yy = this->constantSpace()->element( "yy" );

    auto H = this->traceSpaceOrder0()->element( "H" );
    if ( ioption( prefixvm( this->prefix(), "hface" ) ) == 0 )
        H.on( _range = elements( this->traceSpaceOrder0()->mesh() ), _expr = cst( this->fluxSpace()->mesh()->hMax() ) );
    else if ( ioption( prefixvm( this->prefix(), "hface" ) ) == 1 )
        H.on( _range = elements( this->traceSpaceOrder0()->mesh() ), _expr = cst( this->fluxSpace()->mesh()->hMin() ) );
    else if ( ioption( prefixvm( this->prefix(), "hface" ) ) == 2 )
        H.on( _range = elements( this->traceSpaceOrder0()->mesh() ), _expr = cst( this->fluxSpace()->mesh()->hAverage() ) );
    else
        H.on( _range = elements( this->traceSpaceOrder0()->mesh() ), _expr = h() );
    // stabilisation parameter
    auto tau_constant = cst( doption( prefixvm( this->prefix(), "tau_constant" ) ) );

    auto marker = M_0dList[i].marker();
    int j = i + 1;
    // this->assembleIBC(i,marker);

    // Feel::cout << "j: " << j << "\t i: " << i << std::endl;

    // 0D EQUATION
    double meas = integrate( _range = markedfaces( this->mesh(), marker ), _expr = cst( 1.0 ) ).evaluate()( 0, 0 );

    std::string Cbuffer_str = soption( "coupling.Cbuffer_name" );
    Cbuffer_str += ".C";
    std::string Rbuffer_str = soption( "coupling.Rbuffer_name" );
    Rbuffer_str += ".R";

    double Cbuffer = M_circuit->getValue<double>( Cbuffer_str );
    ;
    double Rbuffer = M_circuit->getValue<double>( Rbuffer_str );

    Feel::cout << "Resistance of the buffer " << Rbuffer_str << ": " << Rbuffer << std::endl;
    Feel::cout << "Capacitance of the buffer " << Cbuffer_str << ": " << Cbuffer << std::endl;

    // Part on the IBC line
    // -1/(R |Gamma_I|) <u_I, mu2>_Gamma_I
    bbf( 3_c, 3_c, i, i ) += integrate( _range = markedfaces( this->mesh(), marker ), _expr = -idt( uI ) * id( nu ) / Rbuffer / meas );

    // 1/(R |Gamma_I|) <Y,mu2>_Gamma_I
    bbf( 3_c, 3_c, i, j ) += integrate( _range = markedfaces( this->mesh(), marker ), _expr = idt( yy ) * id( nu ) / Rbuffer / meas );

    // Part on the 0d line

#if 1
    // <j.n, m>_Gamma_I
    bbf( 3_c, 0_c, j, 0 ) += integrate( _range = markedfaces( this->mesh(), marker ), _expr = ( trans( idt( u ) ) * N() ) * id( nu ) );

    // <tau p, m>_Gamma_I
    bbf( 3_c, 1_c, j, 1 ) += integrate( _range = markedfaces( this->mesh(), marker ), _expr = tau_constant * idt( p ) * id( nu ) * ( pow( idv( H ), this->tauOrder() ) ) );

    // -<tau lambda2, m>_Gamma_I
    bbf( 3_c, 3_c, j, i ) += integrate( _range = markedfaces( this->mesh(), marker ), _expr = -tau_constant * id( nu ) * idt( uI ) * ( pow( idv( H ), this->tauOrder() ) ) );
#else
    // 1/(R |Gamma_I|) <u_I, mu2>_Gamma_I
    bbf( 3_c, 3_c, j, i ) += integrate( _range = markedfaces( this->mesh(), marker ), _expr = idt( uI ) * id( nu ) / Rbuffer / meas );

    // -1/(R |Gamma_I|) <Y,mu2>_Gamma_I
    bbf( 3_c, 3_c, j, j ) += integrate( _range = markedfaces( this->mesh(), marker ), _expr = -idt( yy ) * id( nu ) / Rbuffer / meas );
#endif

    // -< C/|Gamma_I| Y/dt, mu3>_Gamma_I
    bbf( 3_c, 3_c, j, j ) += integrate( _range = markedfaces( this->mesh(), marker ), _expr = -Cbuffer * M_bdf_buffer[i]->polyDerivCoefficient( 0 ) * idt( yy ) * id( nu ) / meas );

    double tElapsed = this->timerTool( "Constructor" ).stop( "assembleMatrix0d" );
    this->log( "CoupledMixedPoisson", "assembleMatrix0d", ( boost::format( "finish in %1% s" ) % tElapsed ).str() );
}

template <int Dim, int Order, int G_Order, int E_Order>
void CoupledMixedPoisson<Dim, Order, G_Order, E_Order>::assembleRhs0d( int i )
{

    this->log( "CoupledMixedPoisson", "assembleRhs0d", "start" );
    this->timerTool( "Constructor" ).start();

    auto blf = blockform1( *( this->getPS() ), this->getF() );
    auto nu = this->constantSpace()->element( "nu" );

    // int j = this->integralCondition() + i; // index where to start
    int j = i + 1;

    auto exAtMarker = M_0dList[i];
    auto marker = exAtMarker.marker();
    // double g = 0.0; // exAtMarker.expression();

    // this->assembleRhsIBC(i, marker, g);

    double meas = integrate( _range = markedfaces( this->mesh(), marker ), _expr = cst( 1.0 ) ).evaluate()( 0, 0 );
    ;

    // 0d PART
    std::string Cbuffer_str = soption( "coupling.Cbuffer_name" );
    Cbuffer_str += ".C";
    double Cbuffer = M_circuit->getValue<double>( Cbuffer_str );

    auto bdf_poly = M_bdf_buffer[i]->polyDeriv();
    Feel::cout << "Value rhs 0d: " << mean( _range = markedfaces( this->mesh(), marker ), _expr = idv( bdf_poly ) )( 0, 0 ) << std::endl;
    // Feel::cout << "Value measure on 0d: " << meas << std::endl;

    // -< C/|Gamma_I| Yold/dt, mu3>
    blf( 3_c, j ) += integrate( _range = markedfaces( this->mesh(), marker ), _expr = -Cbuffer * idv( bdf_poly ) * id( nu ) / meas );

    double tElapsed = this->timerTool( "Constructor" ).stop( "assembleRhs0d" );
    this->log( "CoupledMixedPoisson", "assembleRhs0d", ( boost::format( "finish in %1% s" ) % tElapsed ).str() );
}

template <int Dim, int Order, int G_Order, int E_Order>
void CoupledMixedPoisson<Dim, Order, G_Order, E_Order>::exportResults( double time, mesh_ptrtype mesh, op_interp_ptrtype Idh, opv_interp_ptrtype Idhv )
{
    super_type::exportResults( time, mesh, Idh, Idhv );
    this->log( "CoupledMixedPoisson", "exportResults", "start" );

    // Export computed solutions
    for ( auto const& field : this->modelProperties().postProcess().exports().fields() )
    {
        if ( field == "buffer_variable" )
        {
            LOG( INFO ) << "exporting buffer variable at time " << time;
            // this->exporterMP() -> step( time )->add(prefixvm(this->prefix(), "buffer_variable"), M_Y );
        }
    }

    this->log( "CoupledMixedPoisson", "exportResults", "finish" );
}

template <int Dim, int Order, int G_Order, int E_Order>
void CoupledMixedPoisson<Dim, Order, G_Order, E_Order>::solve()
{

#ifdef USE_SAME_MAT
    auto bbf = blockform2( *( this->getPS() ), this->M_A_cst );
#else
    auto bbf = blockform2( *( this->getPS() ), this->M_A );
    //this->M_A->printMatlab( "A-" + std::to_string( this->currentTime() ) + ".m" );
#endif

    auto blf = blockform1( *( this->getPS() ), this->M_F );

    auto U = this->getPS()->element();

    //this->M_F->printMatlab( "F-" + std::to_string( this->currentTime() ) + ".m" );

    tic();
    bbf.solve( _solution = U, _rhs = blf, _condense = boption( prefixvm( this->prefix(), "use-sc" ) ), _name = this->prefix() );
    toc( "CoupledMixedPoisson : solve", this->verbose() || FLAGS_v > 0 );

    this->M_up = U( 0_c );
    this->M_pp = U( 1_c );

    for ( int i = 0; i < this->integralCondition(); i++ )
        this->M_mup[i] = U( 3_c, i );

    for ( int i = 0; i < M_0dCondition; i++ )
    {
        M_Y[i] = U( 3_c, i + 1 );
        this->second_step( this->currentTime() + this->timeStep(), i );
    }
}

// Second step solved with OpenModelica
template <int Dim, int Order, int G_Order, int E_Order>
void CoupledMixedPoisson<Dim, Order, G_Order, E_Order>::second_step( double final_time, int i )
{
    tic();
    this->log( "CoupledMixedPoisson", "0D model", "start" );

    // Update the buffer variable after first step
    std::string buffer_str = soption( "coupling.var_buffer" );

    Feel::cout << "Buffer variable: \t" << buffer_str << std::endl;
    // double Pi_1 = M_Mmean(_range=markedelements(this->mesh(),M_0dList[i].marker()),_expr= idv(M_Y[i]))(0,0);

    M_circuit->setValue<double>( buffer_str, ( M_Y[i] ).max() );
    // M_circuit->setValue<double>( buffer_str, Pi_1 );
    Feel::cout << "Pressure on buffer before 2nd step: \t " << M_circuit->getValue<double>( buffer_str ) << std::endl;

    // double j_integral = integrate( _quad=_Q<expr_order>(), _range=markedfaces(this->mesh(),marker),_expr=trans(idv(this->M_up))*N()).evaluate()(0,0);
    // Feel::cout << "Integral value of the flow: \t " << j_integral << std::endl;

    // solve the problem with OpenModelica
    Feel::cout << "2nd step: from " << this->currentTime() << " to " << final_time << std::endl;
    M_circuit->doSteps( final_time );

    // Update the buffer variable after second step
    double buffer = M_circuit->getValue<double>( buffer_str );
    Feel::cout << "Pressure on buffer after 2nd step: \t " << buffer << std::endl;

    M_Y[i] = project( _space = this->M_Ch, _expr = cst( buffer ) );
    M_bdf_buffer[i]->setUnknown( 0, M_Y[i] );

    this->log( "CoupledMixedPoisson", "0D model", "finish" );
    toc( "0D model", this->verbose() || FLAGS_v > 0 );
}

template <int Dim, int Order, int G_Order, int E_Order>
void CoupledMixedPoisson<Dim, Order, G_Order, E_Order>::run( op_interp_ptrtype Idh_poi, opv_interp_ptrtype Idhv_poi )
{

    if ( this->isStationary() )
    {
        Feel::cout << std::endl
                   << "ERROR: this model has to be unsteady." << std::endl
                   << std::endl;
        return;
    }

    for ( ; !this->timeStepBase()->isFinished(); this->updateTimeStep() )
    {
        Feel::cout << "===============================================" << std::endl;
        Feel::cout << "time simulation: " << this->time() << "s \n";
        Feel::cout << "===============================================" << std::endl;
        this->setMatricesAndVectorToZero();
        this->assembleAll();
        this->solve();
        this->exportResults( this->mesh(), Idh_poi, Idhv_poi );
    }
}

} // namespace FeelModels

} // namespace Feel

#endif /* COUPLING_HPP */

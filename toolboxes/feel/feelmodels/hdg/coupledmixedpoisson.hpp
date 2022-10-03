#ifndef FEELPP_COUPLEDMIXEDPOISSON_HPP
#define FEELPP_COUPLEDMIXEDPOISSON_HPP

#include <feel/feelfmi/fmi4cpp.hpp>
#include <feel/feelmodels/hdg/mixedpoisson.hpp>

namespace Feel
{
#if 0
class ModelBoundaryConditionCoupling : public ModelBoundaryCondition
{
public:
    using super_type = ModelBoundaryCondition;
    ModelBoundaryConditionCoupling( ModelBoundaryCondition const& bc ) : super_type(bc) {
        M_circuit = M_pt.get<std::string>("circuit");
        M_capacitor = M_pt.get<std::string>("capacitor");
        M_resistor = M_pt.get<std::string>("resistor");
        M_buffer = M_pt.get<std::string>("buffer");
    }
    std::string const& circuit() const { return M_circuit; }
    std::string const& capacitor() const { return M_capacitor; }
    std::string const& resistor() const { return M_resistor; }
    std::string const& buffer() const { return M_buffer; }

private:
    std::string M_circuit;
    std::string M_capacitor;
    std::string M_resistor;
    std::string M_buffer;
}; // class ModelBoundaryConditionCoupling
#endif
namespace FeelModels
{

/**
 * Toolbox CoupledMixedPoisson
 * @ingroup Toolboxes
 */
template<typename ConvexType, int Order>
class CoupledMixedPoisson : public MixedPoisson<ConvexType, Order>// ,
                            // public std::enable_shared_from_this<CoupledMixedPoisson<ConvexType, Order> >
{
public:
    using super_type = MixedPoisson<ConvexType, Order>;
    using super2_type = ModelNumerical;
    using self_type = CoupledMixedPoisson<ConvexType, Order>;
    using self_ptrtype = std::shared_ptr<self_type>;
    using value_type = double;

    using space_constant_type = typename super_type::space_traceibc_type;
    using space_constant_ptrtype = std::shared_ptr<space_constant_type>;
    using element_constant_type = typename super_type::element_traceibc_type;
    using element_constant_ptrtype = typename super_type::element_traceibc_ptrtype;
    using element_flux_type = typename super_type::element_flux_type;
    using element_potential_type = typename super_type::element_potential_type;
    using product_space_ptrtype = typename super_type::product_space_ptrtype;

    using bdf_buffer_type = Bdf<space_constant_type>;
    using bdf_buffer_ptrtype = std::shared_ptr<bdf_buffer_type>;
    using bdfs_buffer_type = std::vector<bdf_buffer_ptrtype>;

    using fmu_ptrtype = std::unique_ptr<fmi2::cs_slave>;
    using fmus_type = std::vector<fmu_ptrtype>;

    CoupledMixedPoisson( std::string const& prefix = "hdg.poisson",
                         MixedPoissonPhysics const& physic = MixedPoissonPhysics::None,
                         worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
                         std::string const& subPrefix = "",
                         ModelBaseRepository const& modelRep = ModelBaseRepository()
                         )
        : super_type(prefix, physic, _worldComm, subPrefix, modelRep),
          ModelBase(prefix, MixedPoissonPhysicsMap[physic]["keyword"], _worldComm, subPrefix, modelRep ) {
    }
    self_ptrtype shared_from_this()
    { return std::static_pointer_cast<self_type>(super_type::shared_from_this()); }
    std::shared_ptr<const self_type> shared_from_this() const
    { return std::static_pointer_cast<self_type>(super_type::shared_from_this()); }

    void initBoundaryConditions() override {
        super_type::initBoundaryConditions();
        for ( auto const& [bcName,bcData] : this->M_boundaryConditions->couplingODEs() )
        {
            auto fmu = fmi2::fmu( Environment::expand( bcData->circuit() ) ).as_cs_fmu()->new_instance();
            fmu->setup_experiment();
            fmu->enter_initialization_mode();
            fmu->exit_initialization_mode();
            this->M_circuits.push_back( std::move(fmu) );
        }
    }

    void initTimeStep() override {
        super_type::initTimeStep();

        std::string myFileFormat = soption(_name="ts.file-format");// without prefix

        int bdfOrder = 1;
        if ( this->M_timeStepping == "BDF" )
            bdfOrder = ioption(_prefix=this->prefix(),_name="bdf.order");
        int nConsecutiveSave = std::max( 3, bdfOrder ); // at least 3 is required when restart with theta scheme

        int i = 0;
        int nbIbc = this->nbIbc();
        for ( auto const& [bcName,bcData] : this->M_boundaryConditions->couplingODEs() )
        {
            M_bdfsBuffer.push_back(this->createBdf( this->spaceTraceIbc(), bcName, bdfOrder, nConsecutiveSave, myFileFormat) );
            if( this->doRestart() )
            {
                double tir = M_bdfsBuffer.back()->restart();
                *(this->M_mup[nbIbc+2*i+1]) = M_bdfsBuffer.back()->unknown(0);
            }
            ++i;
        }

        this->updateParameterValues();
        this->updateInitialConditions(this->symbolsExpr());
    }

    void startTimeStep() override {
        super_type::startTimeStep();
        if(!this->doRestart())
        {
            for(auto const& thebdf : M_bdfsBuffer )
                thebdf->start( thebdf->unknowns() );
        }
    }

    void updateTimeStep() override {
        super_type::updateTimeStep();

        int i = 0;
        int nbIbc = this->nbIbc();
        for(auto const& thebdf : M_bdfsBuffer )
        {
            bool rebuildCstAssembly = false;
            if ( this->M_timeStepping == "BDF" )
            {
                int previousTimeOrder = thebdf->timeOrder();
                thebdf->next( *(this->M_mup[nbIbc+2*i+1]) );
                int currentTimeOrder = thebdf->timeOrder();
                rebuildCstAssembly = previousTimeOrder != currentTimeOrder && thebdf->strategy() == TS_STRATEGY_DT_CONSTANT;
                this->updateTime( thebdf->time() );
            }
            else if ( this->M_timeStepping == "Theta" )
            {
                thebdf->next( *(this->M_mup[nbIbc+2*i+1]) );
                this->updateTime( thebdf->time() );
            }

            if ( rebuildCstAssembly )
                this->setNeedToRebuildCstPart(true);
            ++i;
        }
        this->updateParameterValues();
    }

    template<typename SymbolsExprType>
    void updateInitialConditions( SymbolsExprType const& se ) {
        super_type::updateInitialConditions( se );
        if( !this->doRestart() )
        {
            int i = 0;
            int nbIbc = this->nbIbc();
            for ( auto const& [bcName,bcData] : this->M_boundaryConditions->couplingODEs() )
            {
                std::vector<element_constant_ptrtype> icConstantFields;
                if( this->isStationary() )
                    icConstantFields = { this->M_mup[nbIbc+2*i+1] };
                else
                    icConstantFields = M_bdfsBuffer[i]->unknowns();

                super2_type::updateInitialConditions( bcName, markedfaces(this->mesh(), bcData->markers()), se, icConstantFields );

                if( !this->isStationary() )
                    *(this->M_mup[nbIbc+2*i+1]) = M_bdfsBuffer[i]->unknown(0);
                ++i;
            }
        }
    }

    void initPostProcess() override {
        super_type::initPostProcess();

        int i = 0;
        for ( auto const& [bcName,bcData] : this->M_boundaryConditions->couplingODEs() )
        {
            auto fields = this->modelProperties().postProcess().exports( bcName ).fields();
            M_export_list.insert(M_export_list.end(),fields.begin(),fields.end());
        }

        std::string export_dir = Environment::expand(soption(_name="fmu.export-directory"));
        if( !export_dir.empty() )
            M_export_path = fs::path(export_dir)/"fmu_values.csv";
        else
            M_export_path = fs::path("fmu_values.csv");

        // need to add option with specific prefix for each circuit
        // M_export_list = vsoption(_name="fmu.exported-variables");

        std::ofstream data;
        data.open( M_export_path.string() , std::ios::trunc );
        data << "t";
        for( auto const& var : M_export_list )
            data << "," << var;
        data << std::endl;
        data.close();
    }

    using super_type::exportResults;
    void exportResults( double time ) override {
        super_type::exportResults(time);

        std::ofstream data;
        data.open( M_export_path.string() , std::ios::app );
        data << time;
        int i = 0;
        for ( auto const& [bcName,bcData] : this->M_boundaryConditions->couplingODEs() )
        {
            auto fields = this->modelProperties().postProcess().exports( bcName ).fields();
            auto md = M_circuits[i]->get_model_description();
            std::vector<fmi2ValueReference> vr_buffer;
            std::transform(fields.begin(), fields.end(), std::back_inserter(vr_buffer),
                           [&md](std::string const& f) { return md->get_variable_by_name(f).value_reference; });
            std::vector<double> buffer(vr_buffer.size());
            M_circuits[i]->read_real( vr_buffer, buffer );
            for( auto const& v : buffer )
                data << "," << v;
            ++i;
        }

        data << std::endl;
        data.close();
    }


    int constantSpacesSize() const override { return this->M_boundaryConditions->integral().size() + this->M_boundaryConditions->couplingODEs().size(); }
    int nbIbc() const { return this->constantSpacesSize() - 2*this->M_boundaryConditions->couplingODEs().size(); }

    void setSpaceProperties(product_space_ptrtype const& ibcSpaces) override {
        std::vector<std::string> props(this->M_boundaryConditions->integral().size(), "Ibc");
        for ( auto const& [bcName,bcData] : this->M_boundaryConditions->couplingODEs() )
        {
            props.push_back("Ibc");
            props.push_back("Ode");
        }
        ibcSpaces->setProperties( props );
    }

    void updateLinearPDE( ModelAlgebraic::DataUpdateLinear& data ) const override {
        super_type::updateLinearPDE( data );

        auto A = std::dynamic_pointer_cast<condensed_matrix_t<typename super_type::value_type>>(data.matrix());
        auto F = std::dynamic_pointer_cast<condensed_vector_t<typename super_type::value_type>>(data.rhs());
        bool buildCstPart = data.buildCstPart();
        bool buildNonCstPart = !buildCstPart;
        if( buildNonCstPart )
            return;

        auto ps = this->spaceProduct();
        auto bbf = blockform2( ps, A );
        auto blf = blockform1( ps, F );
        auto u = this->fieldFlux();
        auto p = this->fieldPotential();
        auto l = this->M_Ch->element();
        auto y = this->M_Ch->element();
        auto tau_constant = cst(this->M_tauCst);

        int i = 0;
        int nbIbc = this->nbIbc();

        for ( auto const& [bcName,bcData] : this->M_boundaryConditions->couplingODEs() )
        {
            auto bcRangeFaces = markedfaces(support(this->M_Wh), bcData->markers());
            // Feel::cout << "M_mup["<<nbIbc+2*i<<"] = " << this->M_mup[nbIbc+2*i]->max() << std::endl;
            // Feel::cout << "M_mup["<<nbIbc+2*i+1<<"] = " << this->M_mup[nbIbc+2*i+1]->max() << std::endl;
            double meas = integrate( _range=bcRangeFaces,
                                     _expr=cst(1.)).evaluate()(0,0);
            auto md = M_circuits[i]->get_model_description();

            std::vector<fmi2ValueReference> vr_buffer = { md->get_variable_by_name( bcData->capacitor() ).value_reference,
                                                          md->get_variable_by_name( bcData->resistor() ).value_reference };
            std::vector<double> buffer(2);
            M_circuits[i]->read_real( vr_buffer, buffer );
            double Cbuffer=buffer[0], Rbuffer=buffer[1];
            // Feel::cout << "Cbuffer = " << Cbuffer << "\tRbuffer = " << Rbuffer << std::endl;

            // <lambda, v.n>_Gamma_I
            bbf( 0_c, 3_c, 0, nbIbc+2*i ) += integrate( _range=bcRangeFaces,
                                                        _expr= idt(l) * normal(u) );

            // -<lambda, tau w>_Gamma_I
            bbf( 1_c, 3_c, 0, nbIbc+2*i ) += integrate( _range=bcRangeFaces,
                                                        _expr=-tau_constant*inner(idt(l), id(p)) );

            // <j.n, m>_Gamma_I
            bbf( 3_c, 0_c, nbIbc+2*i, 0 ) += integrate( _range=bcRangeFaces,
                                                        _expr=normalt(u) * id(l) );

            // <tau p, m>_Gamma_I
            bbf( 3_c, 1_c, nbIbc+2*i, 0 ) += integrate( _range=bcRangeFaces,
                                                        _expr=tau_constant*inner(idt(p), id(l)) );

            // -<lambda2, m>_Gamma_I
            bbf( 3_c, 3_c, nbIbc+2*i, nbIbc+2*i ) += integrate( _range=bcRangeFaces,
                                                                _expr=-tau_constant*inner(id(l), idt(l)) );

            bbf( 3_c, 3_c, nbIbc+2*i, nbIbc+2*i ) += integrate( _range=bcRangeFaces,
                                                                _expr=-inner(idt(l), id(l))/meas/Rbuffer );
            bbf( 3_c, 3_c, nbIbc+2*i, nbIbc+2*i+1 ) += integrate( _range=bcRangeFaces,
                                                                  _expr=inner(idt(y),id(l))/meas/Rbuffer );


            bbf( 3_c, 3_c, nbIbc+2*i+1, nbIbc+2*i ) += integrate( _range=bcRangeFaces,
                                                                  _expr=-inner(idt(l), id(y))/meas );
            bbf( 3_c, 3_c, nbIbc+2*i+1, nbIbc+2*i+1 ) += integrate( _range=bcRangeFaces,
                                                                    _expr=inner(idt(y), id(y))/meas );

            bbf( 3_c, 3_c, nbIbc+2*i+1, nbIbc+2*i+1 ) += integrate( _range=bcRangeFaces,
                                                                    _expr=Cbuffer*M_bdfsBuffer[i]->polyDerivCoefficient(0)*inner(idt(y), id(y))/meas );
            blf( 3_c, nbIbc+2*i+1 ) += integrate( _range=bcRangeFaces,
                                                  _expr=Cbuffer*inner(idv(this->M_bdfsBuffer[i]->polyDeriv()), id(y))/meas );
            ++i;
        }
    }

    void solve() override {
        solve::strategy s = this->M_useSC ? solve::strategy::static_condensation : solve::strategy::monolithic;
        this->M_A = makeSharedMatrixCondensed<value_type>(s, csrGraphBlocks(*this->M_ps, (s>=solve::strategy::static_condensation)?Pattern::ZERO:Pattern::COUPLED), *this->backend(), (s>=solve::strategy::static_condensation)?false:true );
        this->M_F = makeSharedVectorCondensed<value_type>(s, blockVector(*this->M_ps), *this->backend(), false);
        // this->M_A->zero();
        // this->M_F->zero();
        auto A = std::dynamic_pointer_cast<typename super_type::super_type::backend_type::sparse_matrix_type>(this->M_A);
        auto F = std::dynamic_pointer_cast<typename super_type::super_type::backend_type::vector_type>(this->M_F);
        CHECK(A) << "Dynamic cast for M_A not ok !!!";
        CHECK(F) << "Dynamic cast for M_F not ok !!!";
        auto U = this->M_ps->element();
        U.buildVector(this->backend());
        this->algebraicFactory()->applyAssemblyLinear(U.vectorMonolithic(), A, F);

        auto bbf = blockform2( *this->M_ps, this->M_A);
        auto blf = blockform1( *this->M_ps, this->M_F);

        bbf.solve(_solution=U, _rhs=blf, _condense=this->M_useSC, _name=this->prefix());

        this->M_up = std::make_shared<element_flux_type>( U(0_c) );
        this->M_pp = std::make_shared<element_potential_type>( U(1_c));
        for ( int i = 0; i < this->constantSpacesSize(); i++ )
            this->M_mup[i] = std::make_shared<element_constant_type>( U( 3_c, i ) );

        // Feel::cout << "M_mup["<<0<<"] = " << this->M_mup[0]->max() << std::endl;
        // Feel::cout << "M_mup["<<1<<"] = " << this->M_mup[1]->max() << std::endl;

        int i = 0;
        int nbIbc = this->nbIbc();
        for ( auto const& [bcName,bcData] : this->M_boundaryConditions->couplingODEs() )
        {
            auto md = M_circuits[i]->get_model_description();
            fmi2ValueReference vr_buffer = md->get_variable_by_name( bcData->buffer() ).value_reference;
            fmi2Real v_buffer = this->M_mup[nbIbc+2*i+1]->max();
            M_circuits[i]->write_real( vr_buffer, v_buffer );
            // Feel::cout << "M_mup["<<i<<"] = " << v_buffer << std::endl;
            LOG(INFO) << fmt::format("bc {} time pde system: {}, time ode system: {}, time step: {}", 
                                     i, this->currentTime(), M_circuits[i]->get_simulation_time(), this->timeStepBase()->timeStep() );
            M_circuits[i]->step( this->timeStepBase()->timeStep() );
            M_circuits[i]->read_real( vr_buffer, v_buffer );
            // Feel::cout << "M_cir["<<i<<"] = " << v_buffer << std::endl;
            *this->M_mup[nbIbc+2*i+1] = project( _space=this->M_Ch, _expr=cst(v_buffer) );
            M_bdfsBuffer[i]->setUnknown( 0, *this->M_mup[nbIbc+2*i+1] );
            ++i;
        }

#if 0
        auto functions = this->modelProperties().parameters();
        functions.setParameterValues({{"t",this->currentTime()}});

        auto p = this->M_mup[0]->max();
        auto p_exp = functions["p"].expression();
        auto p_err = normL2(_range=elements(this->mesh()), _expr=idv(this->M_pp)-p_exp);
        auto p_ex = mean(_range=markedfaces(this->mesh(), "top"), _expr=p_exp)(0,0);
        Feel::cout << "p\t exact = " << p_ex << "\tfeel = " << p
                   << "\terr = " << std::abs(p-p_ex) << std::endl;
        Feel::cout << "p_err = " << p_err << std::endl;

        auto j_exp = functions["j"].template expression<3>();
        auto j_err = normL2(_range=elements(this->mesh()), _expr=idv(this->M_up)-j_exp);
        auto j_ex = mean(_range=markedfaces(this->mesh(), "top"), _expr=j_exp)(0,0);
        auto j = mean(_range=markedfaces(this->mesh(), "top"), _expr=idv(this->M_up))(0,0);
        Feel::cout << "j\t exact = " << j_ex << "\tfeel = " << j
                   << "\terr = " << std::abs(j-j_ex) << std::endl;
        Feel::cout << "j_err = " << j_err << std::endl;

        if( M_circuits.size() )
        {
            auto md = M_circuits[0]->get_model_description();
            std::vector<std::string> vs = { "Pi_1.phi", "Pi_2.phi", "Pi_out.p.v" };
            std::vector<fmi2ValueReference> vr_buffer;
            std::vector<double> fs_exact;
            std::transform(vs.begin(), vs.end(), std::back_inserter(vr_buffer),
                           [&md](std::string const& f) { return md->get_variable_by_name(f).value_reference; });
            std::transform(vs.begin(), vs.end(), std::back_inserter(fs_exact),
                           [&functions](std::string const& f) { return functions[f].expression().evaluate()(0,0); });
            std::vector<double> fs_buffer(vr_buffer.size());
            M_circuits[0]->read_real( vr_buffer, fs_buffer );

            for(int i = 0; i < vs.size(); ++i )
                Feel::cout << vs[i] << " exact = " << fs_exact[i] << "\tfmu = " << fs_buffer[i]
                           << "\terr = " << std::abs(fs_exact[i]-fs_buffer[i]) << std::endl;
        }
#endif
    }

private:
    bdfs_buffer_type M_bdfsBuffer;
    fmus_type M_circuits;
    fs::path M_export_path;
    std::vector<std::string> M_export_list;
}; // class CoupledMixedPoisson

} // namespace FeelModels
} // namespace Feel

#endif

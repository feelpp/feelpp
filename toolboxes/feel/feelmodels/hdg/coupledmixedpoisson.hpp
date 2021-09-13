#ifndef FEELPP_COUPLEDMIXEDPOISSON_HPP
#define FEELPP_COUPLEDMIXEDPOISSON_HPP

#include <feel/feelfmi/fmu.hpp>
#include <feel/feelmodels/hdg/mixedpoisson.hpp>

namespace Feel
{

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

namespace FeelModels
{

/**
 * Toolbox CopledMixedPoisson
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

    using fmu_ptrtype = std::shared_ptr<FMU>;
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
        for( auto const& [name, bc] : this->modelProperties().boundaryConditions2().byFieldType( this->M_fluxKey, "Coupling") )
        {
            this->M_bcIntegralMarkerManagement.addMarkerIntegralBC(name, bc.markers() );
            this->M_coupledBC.push_back(ModelBoundaryConditionCoupling(bc));

            auto thefmu = std::make_shared<FMU>();
            thefmu->load(Environment::expand(this->M_coupledBC.back().circuit()));
            thefmu->initialize( this->timeInitial(), this->timeFinal() );
            this->M_circuits.push_back(thefmu);
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
        for( auto const& [name, bc] : this->modelProperties().boundaryConditions2().byFieldType( this->M_fluxKey, "Coupling") )
        {
            M_bdfsBuffer.push_back(this->createBdf( this->spaceTraceIbc(), bc.name(), bdfOrder, nConsecutiveSave, myFileFormat) );
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
            for( auto const& [name, bc] : this->modelProperties().boundaryConditions2().byFieldType( this->M_fluxKey, "Coupling") )
            {
                std::vector<element_constant_ptrtype> icConstantFields;
                if( this->isStationary() )
                    icConstantFields = { this->M_mup[nbIbc+2*i+1] };
                else
                    icConstantFields = M_bdfsBuffer[i]->unknowns();

                super2_type::updateInitialConditions( bc.name(), markedfaces(this->mesh(), bc.markers()), se, icConstantFields );

                if( !this->isStationary() )
                    *(this->M_mup[nbIbc+2*i+1]) = M_bdfsBuffer[i]->unknown(0);
                ++i;
            }
        }
    }

    int constantSpacesSize() const override { return this->M_bcIntegralMarkerManagement.markerIntegralBC().size() + this->M_coupledBC.size(); }
    int nbIbc() const { return this->constantSpacesSize() - 2*this->M_coupledBC.size(); }

    void setSpaceProperties(product_space_ptrtype const& ibcSpaces) override {
        std::vector<std::string> props(this->M_bcIntegralMarkerManagement.markerIntegralBC().size(), "Ibc");
        for( auto const& [name, bc] : this->modelProperties().boundaryConditions2().byFieldType( this->M_fluxKey, "Coupling") )
        {
            props.push_back("Ibc");
            props.push_back("Ode");
        }
        ibcSpaces->setProperties( props );
    }

    void updateLinearPDE( ModelAlgebraic::DataUpdateHDG& data ) const override {
        super_type::updateLinearPDE( data );

        condensed_matrix_ptr_t<double>& A = data.matrix();
        condensed_vector_ptr_t<double>& F = data.rhs();
        bool buildCstPart = data.buildCstPart();
        bool buildNonCstPart = !buildCstPart;

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

        for( auto const& bc : this->M_coupledBC )
        {
            double meas = integrate( _range=markedfaces(support(this->M_Wh), bc.markers()),
                                     _expr=cst(1.)).evaluate()(0,0);
            double Cbuffer = M_circuits[i]->template getValue<double>( bc.capacitor() );
            double Rbuffer = M_circuits[i]->template getValue<double>( bc.resistor() );

            // <lambda, v.n>_Gamma_I
            bbf( 0_c, 3_c, 0, nbIbc+2*i ) += integrate( _range=markedfaces(support(this->M_Wh), bc.markers()),
                                                        _expr= idt(l) * normal(u) );

            // -<lambda, tau w>_Gamma_I
            bbf( 1_c, 3_c, 0, nbIbc+2*i ) += integrate( _range=markedfaces(support(this->M_Wh), bc.markers()),
                                                        _expr=-tau_constant*inner(idt(l), id(p)) );

            // <j.n, m>_Gamma_I
            bbf( 3_c, 0_c, nbIbc+2*i, 0 ) += integrate( _range=markedfaces(support(this->M_Wh), bc.markers()),
                                                        _expr=normalt(u) * id(l) );

            // <tau p, m>_Gamma_I
            bbf( 3_c, 1_c, nbIbc+2*i, 0 ) += integrate( _range=markedfaces(support(this->M_Wh), bc.markers()),
                                                        _expr=tau_constant*inner(idt(p), id(l)) );

            // -<lambda2, m>_Gamma_I
            bbf( 3_c, 3_c, nbIbc+2*i, nbIbc+2*i ) += integrate( _range=markedfaces(support(this->M_Wh), bc.markers()),
                                                                _expr=-tau_constant*inner(id(l), idt(l)) );

            bbf( 3_c, 3_c, nbIbc+2*i, nbIbc+2*i ) += integrate( _range=markedfaces(support(this->M_Wh), bc.markers()),
                                                                _expr=-inner(idt(l), id(l))/meas/Rbuffer );
            bbf( 3_c, 3_c, nbIbc+2*i, nbIbc+2*i+1 ) += integrate( _range=markedfaces(support(this->M_Wh), bc.markers()),
                                                                  _expr=inner(idt(y),id(l))/meas/Rbuffer );


            bbf( 3_c, 3_c, nbIbc+2*i+1, nbIbc+2*i ) += integrate( _range=markedfaces(support(this->M_Wh), bc.markers()),
                                                                  _expr=-inner(idt(l), id(y))/meas );
            bbf( 3_c, 3_c, nbIbc+2*i+1, nbIbc+2*i+1 ) += integrate( _range=markedfaces(support(this->M_Wh), bc.markers()),
                                                                    _expr=inner(idt(y), id(y))/meas );

            bbf( 3_c, 3_c, nbIbc+2*i+1, nbIbc+2*i+1 ) += integrate( _range=markedfaces(support(this->M_Wh), bc.markers()),
                                                                    _expr=Cbuffer*M_bdfsBuffer[i]->polyDerivCoefficient(0)*inner(idt(y), id(y))/meas );
            blf( 3_c, nbIbc+2*i+1 ) += integrate( _range=markedfaces(support(this->M_Wh), bc.markers()),
                                                  _expr=Cbuffer*inner(idv(this->M_bdfsBuffer[i]->polyDeriv()), id(y))/meas );
            ++i;
        }
    }

    void solve() override {
        this->M_A->zero();
        this->M_F->zero();
        ModelAlgebraic::DataUpdateHDG dataHDGCst(this->M_A, this->M_F, true);
        this->updateLinearPDE(dataHDGCst);

        auto bbf = blockform2( *this->M_ps, this->M_A);
        auto blf = blockform1( *this->M_ps, this->M_F);
        auto U = this->M_ps->element();

        bbf.solve(_solution=U, _rhs=blf, _condense=this->M_useSC, _name=this->prefix());

        this->M_up = std::make_shared<element_flux_type>( U(0_c) );
        this->M_pp = std::make_shared<element_potential_type>( U(1_c));
        for ( int i = 0; i < this->constantSpacesSize(); i++ )
            this->M_mup[i] = std::make_shared<element_constant_type>( U( 3_c, i ) );

        int i = 0;
        int nbIbc = this->nbIbc();
        for( auto const& bc : M_coupledBC )
        {
            M_circuits[i]->setValue( bc.buffer(), this->M_mup[nbIbc+2*i+1]->max() );
            M_circuits[i]->doSteps(this->currentTime() );
            double buffer = M_circuits[i]->template getValue<double>( bc.buffer() );
            *this->M_mup[nbIbc+2*i+1] = project( _space=this->M_Ch, _expr=cst(buffer) );
            M_bdfsBuffer[i]->setUnknown( 0, *this->M_mup[nbIbc+2*i+1] );
            ++i;
        }

#if 0
        auto functions = this->modelProperties().functions();
        functions.setParameterValues({{"t",this->currentTime()}});
        auto p = this->M_mup[0]->max();
        auto p_exp = functions["p"].expressionScalar();
        auto p_err = normL2(_range=elements(this->mesh()), _expr=idv(this->M_pp)-p_exp);
        auto p_ex = mean(_range=markedfaces(this->mesh(), "top"), _expr=p_exp)(0,0);
        auto j_ex = functions["j"].expressionVectorial3();
        auto j_err = normL2(_range=elements(this->mesh()), _expr=idv(this->M_up)-j_ex);
        Feel::cout << "get Pi_1.phi" << std::endl;
        auto pi1 = M_circuits[0]->template getValue<double>( "Pi_1.phi");
        auto pi1_ex = functions["Pi_1.phi"].expressionScalar().evaluate()(0,0);
        Feel::cout << "get Pi_2.phi" << std::endl;
        auto pi2 = M_circuits[0]->template getValue<double>( "Pi_2.phi");
        auto pi2_ex = functions["Pi_2.phi"].expressionScalar().evaluate()(0,0);
        Feel::cout << "get Pi_out.p.v" << std::endl;
        auto piOut = M_circuits[0]->template getValue<double>( "Pi_out.p.v");
        auto piOut_ex = functions["Pi_out.p.v"].expressionScalar().evaluate()(0,0);
        Feel::cout << "p\t exact = " << p_ex << "\tfeel = " << p
                   << "\terr = " << std::abs(p-p_ex) << std::endl;
        Feel::cout << "p_err = " << p_err << "\tj_err = " << j_err << std::endl;
        Feel::cout << "Pi1\t exact = " << pi1_ex << "\tfmu = " << pi1
                   << "\terr = " << std::abs(pi1_ex-pi1) << std::endl;
        Feel::cout << "Pi2\t exact = " << pi2_ex << "\tfmu = " << pi2
                   << "\terr = " << std::abs(pi2_ex-pi2) << std::endl;
        Feel::cout << "PiOut\t exact = " << piOut_ex << "\tfmu = " << piOut
                   << "\terr = " << std::abs(piOut_ex-piOut) << std::endl;
#endif
    }

private:
    std::vector<ModelBoundaryConditionCoupling> M_coupledBC;
    bdfs_buffer_type M_bdfsBuffer;
    fmus_type M_circuits;
}; // class CoupledMixedPoisson

} // namespace FeelModels
} // namespace Feel

#endif

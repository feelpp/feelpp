/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_VF_FLUIDMEC_CONVECTIVETERM_H
#define FEELPP_TOOLBOXES_VF_FLUIDMEC_CONVECTIVETERM_H

#include <feel/feelmodels/modelvf/exprtensorbase.hpp>
#include <feel/feelmodels/modelvf/exprevaluatefieldoperators.hpp>

namespace Feel {
namespace FeelModels {

enum class FluidMecConvectiveTermFormulation { Convective=0, SkewSymmetric, Conservative, Rotational, EMAC };

template<typename FluidMecConvectiveTermBaseType,FluidMecConvectiveTermFormulation ConvTermType,bool ApplyLinearization>
struct FluidMecConvectiveTermImpl;

template<typename ExprEvaluateFieldOperatorsVelocityType,typename FiniteElementVelocityType,
         typename ModelPhysicFluidType, typename ExprEvaluateFieldOperatorsVelocityLinearizationType, ExprApplyType ExprApplied>
struct FluidMecConvectiveTermBase
{
    using this_type = FluidMecConvectiveTermBase<ExprEvaluateFieldOperatorsVelocityType,FiniteElementVelocityType,ModelPhysicFluidType,ExprEvaluateFieldOperatorsVelocityLinearizationType,ExprApplied>;
    using model_physic_fluid_type = ModelPhysicFluidType;
    using expr_evaluate_velocity_opertors_type = ExprEvaluateFieldOperatorsVelocityType;
    using expr_evaluate_velocity_linearization_opertors_type = ExprEvaluateFieldOperatorsVelocityLinearizationType;
    using value_type = typename expr_evaluate_velocity_opertors_type::value_type;

    static constexpr bool is_applied_as_eval = (ExprApplied == ExprApplyType::EVAL);
    static constexpr bool is_applied_as_jacobian = (ExprApplied == ExprApplyType::JACOBIAN);

    static constexpr uint16_type nDim = model_physic_fluid_type::nDim;
    static constexpr uint16_type nRealDim = nDim;
    typedef Shape<nDim, Vectorial, false, false> shape_type;

    static const bool is_terminal = true;

    static const size_type context = vm::KB|vm::GRAD;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = false;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = is_applied_as_jacobian && std::is_same_v<Func,FiniteElementVelocityType>;
    };
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    static constexpr auto tuple_specialization_type_object = hana::to_tuple(hana::tuple_t<
                                                                            FluidMecConvectiveTermImpl<this_type,FluidMecConvectiveTermFormulation::Convective,false>,
                                                                            FluidMecConvectiveTermImpl<this_type,FluidMecConvectiveTermFormulation::Convective,true>,
                                                                            FluidMecConvectiveTermImpl<this_type,FluidMecConvectiveTermFormulation::SkewSymmetric,false>,
                                                                            FluidMecConvectiveTermImpl<this_type,FluidMecConvectiveTermFormulation::SkewSymmetric,true>,
                                                                            FluidMecConvectiveTermImpl<this_type,FluidMecConvectiveTermFormulation::Conservative,false>,
                                                                            FluidMecConvectiveTermImpl<this_type,FluidMecConvectiveTermFormulation::Conservative,true>,
                                                                            FluidMecConvectiveTermImpl<this_type,FluidMecConvectiveTermFormulation::Rotational,false>,
                                                                            FluidMecConvectiveTermImpl<this_type,FluidMecConvectiveTermFormulation::Rotational,true>,
                                                                            FluidMecConvectiveTermImpl<this_type,FluidMecConvectiveTermFormulation::EMAC,false>,
                                                                            FluidMecConvectiveTermImpl<this_type,FluidMecConvectiveTermFormulation::EMAC,true>
                                                                            >);


    FluidMecConvectiveTermBase( std::shared_ptr<expr_evaluate_velocity_opertors_type> exprEvaluateVelocityOperators,
                                std::shared_ptr<model_physic_fluid_type> const& physicFluidData,
                                std::shared_ptr<expr_evaluate_velocity_linearization_opertors_type> exprEvaluateVelocityLinearizationOperators )
        :
        M_exprEvaluateVelocityOperators( exprEvaluateVelocityOperators ),
        M_physicFluidData( physicFluidData ),
        M_exprEvaluateVelocityLinearizationOperators( exprEvaluateVelocityLinearizationOperators )
        {}

    FluidMecConvectiveTermBase( FluidMecConvectiveTermBase const& ) = default;
    FluidMecConvectiveTermBase( FluidMecConvectiveTermBase && ) = default;
    virtual ~FluidMecConvectiveTermBase() {}


    static std::shared_ptr<this_type> New( std::shared_ptr<expr_evaluate_velocity_opertors_type> exprEvaluateVelocityOperators,
                                           std::shared_ptr<model_physic_fluid_type> const& physicFluidData,
                                           std::shared_ptr<expr_evaluate_velocity_linearization_opertors_type> exprEvaluateVelocityLinearizationOperators,
                                           bool applyLinearization );


    virtual bool applyLinearization() const = 0;

    size_type dynamicContext() const { return 0; }

    uint16_type polynomialOrder() const
        {
            if ( this->applyLinearization() )
                return expr_evaluate_velocity_linearization_opertors_type::polynomialOrderId()+expr_evaluate_velocity_opertors_type::polynomialOrderGrad();
            else
                return expr_evaluate_velocity_opertors_type::polynomialOrderId()+expr_evaluate_velocity_opertors_type::polynomialOrderGrad();
        }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    using tensor_base_type = tensorBase<Geo_t,Basis_i_t,Basis_j_t,shape_type,value_type>;

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor : public this_type::template tensor_base_type<Geo_t,Basis_i_t,Basis_j_t>
    {
    protected :
        using super_type = typename this_type::template tensor_base_type<Geo_t,Basis_i_t,Basis_j_t>;
    public :

        typedef typename super_type::gmc_type gmc_type;
        typedef typename super_type::gm_type gm_type;
        typedef typename super_type::value_type value_type;

        using shape = typename super_type::shape_type;
        using matrix_shape_type = typename super_type::matrix_shape_type;

        using array_shape_type = typename super_type::new_array_shape_type;
        using ret_type = typename super_type::ret_type;

        template<typename... TheArgsType>
        tensor( this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            super_type( geom, theInitArgs... )
            {}

        void update( Geo_t const& geom ) override { CHECK( false ) << "should be override"; };
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu ) override { CHECK( false ) << "should be override"; }
        void update( Geo_t const& geom, Basis_i_t const& fev ) override { CHECK( false ) << "should be override"; }

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs );


        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const override
            {
                if constexpr ( this_type::is_applied_as_eval )
                    return this->evalq( c1,c2,q );
                else
                {
                    CHECK( false) << "not allow";
                    return value_type(0);
                }
            }
        ret_type
        evaliq( uint16_type i, uint16_type q ) const override
            {
                if constexpr ( this_type::is_applied_as_eval )
                    return this->evalq( q );
                else
                {
                    CHECK( false) << "not allow";
                    return ret_type(this->locMatrixShape().data());
               }
            }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const override
            {
                return M_localEval[q](c1,c2);
            }
        ret_type
        evalq( uint16_type q ) const override
            {
                return ret_type(M_localEval[q].data());
            }

    protected :
        std::vector<matrix_shape_type> M_localEval;
    };

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename... TheArgsType>
    std::shared_ptr<tensor<Geo_t,Basis_i_t,Basis_j_t> >
    evaluator( const TheArgsType&... theInitArgs ) const;
#if 0
    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename TheExprExpandedType,typename TupleTensorSymbolsExprType,typename... TheArgsType>
    std::shared_ptr<tensor<Geo_t,Basis_i_t,Basis_j_t> >
    evaluator( std::true_type/**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse, const TheArgsType&... theInitArgs ) const;
#endif

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    friend struct tensor;

protected :
    std::shared_ptr<expr_evaluate_velocity_opertors_type> M_exprEvaluateVelocityOperators;
    std::shared_ptr<model_physic_fluid_type> M_physicFluidData;
    std::shared_ptr<expr_evaluate_velocity_linearization_opertors_type> M_exprEvaluateVelocityLinearizationOperators;
};



template<typename FluidMecConvectiveTermBaseType,FluidMecConvectiveTermFormulation ConvTermFormulation,bool ApplyLinearization>
struct FluidMecConvectiveTermImpl : public FluidMecConvectiveTermBaseType
{
    using super_type = FluidMecConvectiveTermBaseType;
public :
    using this_type = FluidMecConvectiveTermImpl<FluidMecConvectiveTermBaseType,ConvTermFormulation,ApplyLinearization>;
    static constexpr bool apply_linearization = ApplyLinearization;
    static constexpr FluidMecConvectiveTermFormulation formulation = ConvTermFormulation;

    FluidMecConvectiveTermImpl( std::shared_ptr<typename super_type::expr_evaluate_velocity_opertors_type> exprEvaluateVelocityOperators,
                                std::shared_ptr<typename super_type::model_physic_fluid_type> const& physicFluidData,
                                std::shared_ptr<typename super_type::expr_evaluate_velocity_linearization_opertors_type> exprEvaluateVelocityLinearizationOperators )
        :
        super_type( exprEvaluateVelocityOperators, physicFluidData,exprEvaluateVelocityLinearizationOperators )
        {
            if constexpr ( !apply_linearization )
            {
                this->M_exprEvaluateVelocityOperators->setEnableId( true );
                this->M_exprEvaluateVelocityOperators->setEnableGrad( true );
            }
            else
            {
                this->M_exprEvaluateVelocityLinearizationOperators->setEnableId( true );
                if constexpr ( this_type::is_applied_as_eval )
                {
                     this->M_exprEvaluateVelocityOperators->setEnableGrad( true );
                }
            }
        }

    FluidMecConvectiveTermImpl( FluidMecConvectiveTermImpl const& ) = default;
    FluidMecConvectiveTermImpl( FluidMecConvectiveTermImpl && ) = default;

    bool applyLinearization() const override { return apply_linearization; }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor : public this_type::super_type::template tensor<Geo_t,Basis_i_t,Basis_j_t>
    {
        using super_type = typename this_type::super_type::template tensor<Geo_t,Basis_i_t,Basis_j_t>;
        using tensor_expr_evaluate_velocity_opertors_type = typename this_type::expr_evaluate_velocity_opertors_type::template tensor<Geo_t,Basis_i_t,Basis_j_t>;
        using tensor_expr_evaluate_velocity_linearization_opertors_type = typename this_type::expr_evaluate_velocity_linearization_opertors_type::template tensor<Geo_t,Basis_i_t,Basis_j_t>;

    public :
        typedef typename super_type::gmc_type gmc_type;
        typedef typename super_type::gm_type gm_type;
        typedef typename super_type::value_type value_type;

        //using shape = typename super_type::shape_type;
        using matrix_shape_type = typename super_type::matrix_shape_type;
        using array_shape_type = typename super_type::new_array_shape_type;
        using ret_type = typename super_type::ret_type;

        template<typename... TheArgsType>
        tensor( this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            super_type( expr, geom, theInitArgs... ),
            M_expr( expr ),
            M_sameTensorExprWithLinearization( false )
            {
                if constexpr ( this_type::apply_linearization )
                {
                    M_tensorExprEvaluateVelocityLinearizationOperators = std::make_shared<tensor_expr_evaluate_velocity_linearization_opertors_type>( *(expr.M_exprEvaluateVelocityLinearizationOperators), geom, theInitArgs... );
                    if constexpr ( this_type::is_applied_as_eval )
                    {
                        if ( expr.M_exprEvaluateVelocityOperators != expr.M_exprEvaluateVelocityLinearizationOperators )
                            M_tensorExprEvaluateVelocityOperators = std::make_shared<tensor_expr_evaluate_velocity_opertors_type>( *(expr.M_exprEvaluateVelocityOperators), geom, theInitArgs... );
                        else if constexpr ( std::is_same_v<tensor_expr_evaluate_velocity_opertors_type,tensor_expr_evaluate_velocity_linearization_opertors_type> )
                        {
                            M_tensorExprEvaluateVelocityOperators = M_tensorExprEvaluateVelocityLinearizationOperators;
                            M_sameTensorExprWithLinearization= true;
                        }
                   }
                }
                else
                    M_tensorExprEvaluateVelocityOperators = std::make_shared<tensor_expr_evaluate_velocity_opertors_type>( *(expr.M_exprEvaluateVelocityOperators), geom, theInitArgs... );
            }

        void update( Geo_t const& geom ) override
            {
                if constexpr ( !apply_linearization )// this_type::is_applied_as_eval || this_type::is_applied_as_jacobian ) // no linearization
                {
                    this->M_tensorExprEvaluateVelocityOperators->update( geom );
                }
                else
                {
                    this->M_tensorExprEvaluateVelocityLinearizationOperators->update( geom );
                    if constexpr ( this_type::is_applied_as_eval )
                    {
                        if ( !M_sameTensorExprWithLinearization )
                            this->M_tensorExprEvaluateVelocityOperators->update( geom );
                    }
                }
                this->updateImpl();
            }

        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const override
            {
                return evalijqImpl<this_type>( i,j,q );
            }
    private :
        void updateImpl()
        {
            if constexpr ( !this_type::is_applied_as_eval )
                 return;

            uint16_type nPoint = this->gmc()->nPoints();
            //this->M_localEval.resize( nPoint, matrix_shape_type::Zero() );
            this->M_localEval.assign( nPoint, matrix_shape_type::Zero() );
            for ( uint16_type q = 0; q < nPoint; ++q )
            {
                matrix_shape_type & theLocEval = this->M_localEval[q];
                auto const& gradVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalGrad( q );

                // term : ( u \dot \nabla ) u
                if constexpr ( this_type::formulation == FluidMecConvectiveTermFormulation::Convective ||
                               this_type::formulation == FluidMecConvectiveTermFormulation::SkewSymmetric ||
                               this_type::formulation == FluidMecConvectiveTermFormulation::Conservative )
                {
                    if constexpr ( apply_linearization )
                    {
                        auto const& idVelocityEval = M_tensorExprEvaluateVelocityLinearizationOperators->localEvalId( q );
                        //auto const& gradVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalGrad( q );
                        theLocEval += gradVelocityEval*idVelocityEval;
                    }
                    else
                    {
                        auto const& idVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalId( q );
                        //auto const& gradVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalGrad( q );
                        theLocEval += gradVelocityEval*idVelocityEval;
                    }
                }
                // term : alpha* div( u ) u
                if constexpr ( this_type::formulation == FluidMecConvectiveTermFormulation::SkewSymmetric ||
                               this_type::formulation == FluidMecConvectiveTermFormulation::Conservative ||
                               this_type::formulation == FluidMecConvectiveTermFormulation::EMAC )
                {
                    double coeff = this_type::formulation == FluidMecConvectiveTermFormulation::SkewSymmetric? 0.5 : 1.0;
                    if constexpr ( apply_linearization )
                    {
                        auto const& idVelocityEval = M_tensorExprEvaluateVelocityLinearizationOperators->localEvalId( q );
                        //auto const& gradVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalGrad( q );
                        theLocEval += coeff*gradVelocityEval.trace()*idVelocityEval;
                    }
                    else
                    {
                        auto const& idVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalId( q );
                        //auto const& gradVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalGrad( q );
                        theLocEval += coeff*gradVelocityEval.trace()*idVelocityEval;
                    }
                }

                // term : ( \nabla x u ) x u
                if constexpr ( this_type::formulation == FluidMecConvectiveTermFormulation::Rotational )
                {
                    if constexpr ( this_type::nRealDim == 2 )
                    {
                        value_type _curlv = -gradVelocityEval(0,1) + gradVelocityEval(1,0);
                        if constexpr ( apply_linearization )
                        {
                            auto const& idVelocityEval = M_tensorExprEvaluateVelocityLinearizationOperators->localEvalId( q );
                            theLocEval(0) += -_curlv*idVelocityEval(1);
                            theLocEval(1) += _curlv*idVelocityEval(0);
                        }
                        else
                        {
                            auto const& idVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalId( q );
                            theLocEval(0) += -_curlv*idVelocityEval(1);
                            theLocEval(1) += _curlv*idVelocityEval(0);
                        }
                    }
                    else
                    {
                        CHECK( false ) << "TODO";
                    }
                }

                // term : 2 D(u) u
                if constexpr ( this_type::formulation == FluidMecConvectiveTermFormulation::EMAC )
                {
                    //auto const& gradVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalGrad( q );
                    if constexpr ( this_type::nRealDim == 2 )
                    {
                        value_type d00 = 2*gradVelocityEval(0,0);
                        value_type d01 = gradVelocityEval(0,1)+gradVelocityEval(1,0);
                        value_type d11 = 2*gradVelocityEval(1,1);
                        if constexpr ( apply_linearization )
                        {
                            auto const& idVelocityEval = M_tensorExprEvaluateVelocityLinearizationOperators->localEvalId( q );
                            theLocEval(0) += d00*idVelocityEval(0) + d01*idVelocityEval(1);
                            theLocEval(1) += d01*idVelocityEval(0) + d11*idVelocityEval(1);
                        }
                        else
                        {
                            auto const& idVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalId( q );
                            theLocEval(0) += d00*idVelocityEval(0) + d01*idVelocityEval(1);
                            theLocEval(1) += d01*idVelocityEval(0) + d11*idVelocityEval(1);
                        }
                    }
                    else // 3d
                    {
                        CHECK( false ) << "TODO";
                    }
                }

            }
        }

        template <typename TT,std::enable_if_t< TT::is_applied_as_eval, bool> = true>
        ret_type
        evalijqImpl( uint16_type i, uint16_type j, uint16_type q ) const
            {
                CHECK( false ) << "not allow";
                return ret_type(this->locMatrixShape().data());
            }

        template <typename TT,std::enable_if_t< TT::is_applied_as_jacobian, bool> = true>
        ret_type
        evalijqImpl( uint16_type i, uint16_type j, uint16_type q ) const
        {
            auto const& idTrial = this->fecTrial()->id( j, q );
            auto const& gradTrial = this->fecTrial()->grad( j, q );

            auto const& gradVelocityEval = this->M_tensorExprEvaluateVelocityOperators->localEvalGrad( q );
            Eigen::Map< const Eigen::Matrix<typename super_type::value_type,this_type::nRealDim,this_type::nRealDim/*,Eigen::ColMajor*/ > >
                gradTrial_matrix( gradTrial.data() );

            typename super_type::matrix_shape_type & thelocMat = this->locMatrixShape();
            thelocMat.setZero();
            if constexpr ( apply_linearization )
            {
                auto const& idVelocityEval = M_tensorExprEvaluateVelocityLinearizationOperators->localEvalId( q );
                // term : ( u \dot \nabla ) u
                if constexpr ( this_type::formulation == FluidMecConvectiveTermFormulation::Convective ||
                               this_type::formulation == FluidMecConvectiveTermFormulation::SkewSymmetric ||
                               this_type::formulation == FluidMecConvectiveTermFormulation::Conservative )
                {
                    thelocMat += gradTrial_matrix*idVelocityEval;
                }
                // term : alpha* div( u ) u
                if constexpr ( this_type::formulation == FluidMecConvectiveTermFormulation::SkewSymmetric ||
                               this_type::formulation == FluidMecConvectiveTermFormulation::Conservative ||
                               this_type::formulation == FluidMecConvectiveTermFormulation::EMAC )
                {
                    double coeff = this_type::formulation == FluidMecConvectiveTermFormulation::SkewSymmetric? 0.5 : 1.0;
                    thelocMat += coeff*(gradTrial_matrix.trace()*idVelocityEval);
                }
                // term : ( \nabla x u ) x u
                if constexpr ( this_type::formulation == FluidMecConvectiveTermFormulation::Rotational )
                {
                    if constexpr ( this_type::nRealDim == 2 )
                    {
                        value_type _curlt = -gradTrial_matrix(0,1) + gradTrial_matrix(1,0);
                        thelocMat(0) += -_curlt*idVelocityEval(1);
                        thelocMat(1) +=  _curlt*idVelocityEval(0);
                    }
                    else
                    {
                        CHECK( false ) << "TODO";
                    }
                }

                // term : 2 D(u) u
                if constexpr ( this_type::formulation == FluidMecConvectiveTermFormulation::EMAC )
                {
                    if constexpr ( this_type::nRealDim == 2 )
                    {
                        value_type d00t = 2*gradTrial_matrix(0,0);
                        value_type d01t = gradTrial_matrix(0,1)+gradTrial_matrix(1,0);
                        value_type d11t = 2*gradTrial_matrix(1,1);
                        thelocMat(0) += d00t*idVelocityEval(0) + d01t*idVelocityEval(1);
                        thelocMat(1) += d01t*idVelocityEval(0) + d11t*idVelocityEval(1);
                    }
                    else // 3d
                    {
                        CHECK( false ) << "TODO";
                    }
                }

            }
            else
            {
                auto const& idVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalId( q );
                Eigen::Map< const Eigen::Matrix<typename super_type::value_type,this_type::nRealDim,1/*,Eigen::ColMajor*/ > > idTrial_matrix( idTrial.data() );
                // term : ( u \dot \nabla ) u
                if constexpr ( this_type::formulation == FluidMecConvectiveTermFormulation::Convective ||
                               this_type::formulation == FluidMecConvectiveTermFormulation::SkewSymmetric ||
                               this_type::formulation == FluidMecConvectiveTermFormulation::Conservative )
                {
                    thelocMat += gradTrial_matrix*idVelocityEval + gradVelocityEval*idTrial_matrix;
                }
                // term : alpha* div( u ) u
                if constexpr ( this_type::formulation == FluidMecConvectiveTermFormulation::SkewSymmetric ||
                               this_type::formulation == FluidMecConvectiveTermFormulation::Conservative ||
                               this_type::formulation == FluidMecConvectiveTermFormulation::EMAC )
                {
                    double coeff = this_type::formulation == FluidMecConvectiveTermFormulation::SkewSymmetric? 0.5 : 1.0;
                    thelocMat += coeff*( gradTrial_matrix.trace()*idVelocityEval + gradVelocityEval.trace()*idTrial_matrix);
                }
                // term : ( \nabla x u ) x u
                if constexpr ( this_type::formulation == FluidMecConvectiveTermFormulation::Rotational )
                {
                    if constexpr ( this_type::nRealDim == 2 )
                    {
                        value_type _curlv = -gradVelocityEval(0,1) + gradVelocityEval(1,0);
                        value_type _curlt = -gradTrial_matrix(0,1) + gradTrial_matrix(1,0);
                        thelocMat(0) += -_curlv*idTrial_matrix(1) - _curlt*idVelocityEval(1);
                        thelocMat(1) +=  _curlv*idTrial_matrix(0) + _curlt*idVelocityEval(0);
                    }
                    else
                    {
                        CHECK( false ) << "TODO";
                    }
                }

                // term : 2 D(u) u
                if constexpr ( this_type::formulation == FluidMecConvectiveTermFormulation::EMAC )
                {
                    if constexpr ( this_type::nRealDim == 2 )
                    {
                        value_type d00 = 2*gradVelocityEval(0,0);
                        value_type d01 = gradVelocityEval(0,1)+gradVelocityEval(1,0);
                        value_type d11 = 2*gradVelocityEval(1,1);
                        value_type d00t = 2*gradTrial_matrix(0,0);
                        value_type d01t = gradTrial_matrix(0,1)+gradTrial_matrix(1,0);
                        value_type d11t = 2*gradTrial_matrix(1,1);
                        thelocMat(0) += d00*idTrial_matrix(0) + d00t*idVelocityEval(0) + d01*idTrial_matrix(1) + d01t*idVelocityEval(1);
                        thelocMat(1) += d01*idTrial_matrix(0) + d01t*idVelocityEval(0) + d11*idTrial_matrix(1) + d11t*idVelocityEval(1);
                    }
                    else // 3d
                    {
                        CHECK( false ) << "TODO";
                    }
                }

            }

            return ret_type(thelocMat.data());
        }

    private :
        this_type const& M_expr;
        std::shared_ptr<tensor_expr_evaluate_velocity_opertors_type> M_tensorExprEvaluateVelocityOperators;
        std::shared_ptr<tensor_expr_evaluate_velocity_linearization_opertors_type> M_tensorExprEvaluateVelocityLinearizationOperators;
        bool M_sameTensorExprWithLinearization;
    };

};



template<typename ExprEvaluateFieldOperatorsVelocityType,typename FiniteElementVelocityType,
         typename ModelPhysicFluidType, typename ExprEvaluateFieldOperatorsVelocityLinearizationType, ExprApplyType ExprApplied>
std::shared_ptr<FluidMecConvectiveTermBase<ExprEvaluateFieldOperatorsVelocityType,FiniteElementVelocityType,ModelPhysicFluidType,ExprEvaluateFieldOperatorsVelocityLinearizationType,ExprApplied>>
FluidMecConvectiveTermBase<ExprEvaluateFieldOperatorsVelocityType,FiniteElementVelocityType,ModelPhysicFluidType,
                           ExprEvaluateFieldOperatorsVelocityLinearizationType,ExprApplied>::New( std::shared_ptr<expr_evaluate_velocity_opertors_type> exprEvaluateVelocityOperators,
                                                                                                  std::shared_ptr<model_physic_fluid_type> const& physicFluidData,
                                                                                                  std::shared_ptr<expr_evaluate_velocity_linearization_opertors_type> exprEvaluateVelocityLinearizationOperators,
                                                                                                  bool applyLinearization )

{
    switch ( physicFluidData->navierStokesFormulation() )
    {
    case model_physic_fluid_type::NavierStokesFormulation::Convective:
        if ( applyLinearization )
            return std::make_shared<FluidMecConvectiveTermImpl<this_type,FluidMecConvectiveTermFormulation::Convective,true>>( exprEvaluateVelocityOperators,physicFluidData,exprEvaluateVelocityLinearizationOperators );
        else
            return std::make_shared<FluidMecConvectiveTermImpl<this_type,FluidMecConvectiveTermFormulation::Convective,false>>( exprEvaluateVelocityOperators,physicFluidData,exprEvaluateVelocityLinearizationOperators );
    case model_physic_fluid_type::NavierStokesFormulation::SkewSymmetric:
        if ( applyLinearization )
            return std::make_shared<FluidMecConvectiveTermImpl<this_type,FluidMecConvectiveTermFormulation::SkewSymmetric,true>>( exprEvaluateVelocityOperators,physicFluidData,exprEvaluateVelocityLinearizationOperators );
        else
            return std::make_shared<FluidMecConvectiveTermImpl<this_type,FluidMecConvectiveTermFormulation::SkewSymmetric,false>>( exprEvaluateVelocityOperators,physicFluidData,exprEvaluateVelocityLinearizationOperators );
    case model_physic_fluid_type::NavierStokesFormulation::Conservative:
        if ( applyLinearization )
            return std::make_shared<FluidMecConvectiveTermImpl<this_type,FluidMecConvectiveTermFormulation::Conservative,true>>( exprEvaluateVelocityOperators,physicFluidData,exprEvaluateVelocityLinearizationOperators );
        else
            return std::make_shared<FluidMecConvectiveTermImpl<this_type,FluidMecConvectiveTermFormulation::Conservative,false>>( exprEvaluateVelocityOperators,physicFluidData,exprEvaluateVelocityLinearizationOperators );
    case model_physic_fluid_type::NavierStokesFormulation::Rotational:
        if ( applyLinearization )
            return std::make_shared<FluidMecConvectiveTermImpl<this_type,FluidMecConvectiveTermFormulation::Rotational,true>>( exprEvaluateVelocityOperators,physicFluidData,exprEvaluateVelocityLinearizationOperators );
        else
            return std::make_shared<FluidMecConvectiveTermImpl<this_type,FluidMecConvectiveTermFormulation::Rotational,false>>( exprEvaluateVelocityOperators,physicFluidData,exprEvaluateVelocityLinearizationOperators );
    case model_physic_fluid_type::NavierStokesFormulation::EMAC:
        if ( applyLinearization )
            return std::make_shared<FluidMecConvectiveTermImpl<this_type,FluidMecConvectiveTermFormulation::EMAC,true>>( exprEvaluateVelocityOperators,physicFluidData,exprEvaluateVelocityLinearizationOperators );
        else
            return std::make_shared<FluidMecConvectiveTermImpl<this_type,FluidMecConvectiveTermFormulation::EMAC,false>>( exprEvaluateVelocityOperators,physicFluidData,exprEvaluateVelocityLinearizationOperators );
    }

    return std::shared_ptr<this_type>{};
}



template<typename ExprEvaluateFieldOperatorsVelocityType,typename FiniteElementVelocityType,
         typename ModelPhysicFluidType, typename ExprEvaluateFieldOperatorsVelocityLinearizationType, ExprApplyType ExprApplied>
template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename... TheArgsType>
std::shared_ptr<typename FluidMecConvectiveTermBase<ExprEvaluateFieldOperatorsVelocityType,FiniteElementVelocityType,ModelPhysicFluidType,ExprEvaluateFieldOperatorsVelocityLinearizationType,ExprApplied>::template tensor<Geo_t,Basis_i_t,Basis_j_t> >
FluidMecConvectiveTermBase<ExprEvaluateFieldOperatorsVelocityType,FiniteElementVelocityType,ModelPhysicFluidType,
                           ExprEvaluateFieldOperatorsVelocityLinearizationType,ExprApplied>::evaluator( const TheArgsType&... theInitArgs ) const
{
    std::shared_ptr<tensor<Geo_t,Basis_i_t,Basis_j_t>> ret;
    hana::for_each( tuple_specialization_type_object, [this,&ret,&theInitArgs...]( auto e ) {
                                                          if ( ret )
                                                              return;
                                                          using the_expr_type = typename std::decay_t<decltype( e ) >::type;
                                                          if ( auto objptr = dynamic_cast<the_expr_type const*>(this) )
                                                              ret = std::make_shared< typename the_expr_type::template tensor<Geo_t,Basis_i_t,Basis_j_t> >( *objptr, theInitArgs... );
                                                      });
    CHECK( ret ) << "cast failed : something wrong";
    return ret;
}



template <typename ExprType>
struct PolymorphicExpr
{
public:

    using this_type = PolymorphicExpr<ExprType>;

    using expr_type = ExprType;
    static constexpr size_type context = expr_type::context;
    using value_type = typename expr_type::value_type;
    static constexpr bool is_terminal = expr_type::is_terminal;

    template<typename Func>
    struct HasTestFunction { static const bool result = expr_type::template HasTestFunction<Func>; };
    template<typename Func>
    struct HasTrialFunction { static const bool result = expr_type::template HasTrialFunction<Func>; };
    using test_basis = typename expr_type::test_basis;
    using trial_basis = typename expr_type::trial_basis;

    PolymorphicExpr( std::shared_ptr<expr_type> expr )
        :
        M_expr( expr )
        {}
    PolymorphicExpr( PolymorphicExpr const& ) = default;
    PolymorphicExpr( PolymorphicExpr && ) = default;

    size_type dynamicContext() const { return M_expr->dynamicContext(); }

    //! polynomial order
    uint16_type polynomialOrder() const { return M_expr->polynomialOrder(); }

    //! expression is polynomial?
    bool isPolynomial() const { return true; }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
        {
            return M_expr->hasSymbolDependency( symb,se );
        }

    template <typename TheSymbolExprType>
    void dependentSymbols( std::string const& symb, std::map<std::string,std::set<std::string>> & res, TheSymbolExprType const& se ) const
        {
            CHECK( false )  << "TODO";
        }
#if 0
    template <typename OtherSymbolsExprType>
    auto applySymbolsExpr( OtherSymbolsExprType const& se ) const
        {
            CHECK( false )  << "TODO";
            return *this;
        }
#endif
    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
        {
            CHECK( false ) << "not implemented";
            return *this;
        }

    std::shared_ptr<expr_type> exprPtr() const { return M_expr; }
    expr_type const& expr() const { return *M_expr; }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        using tensor_expr_type = typename expr_type::template tensor<Geo_t,Basis_i_t,Basis_j_t>;
        using value_type = typename tensor_expr_type::value_type;
        using shape = typename tensor_expr_type::shape;
        struct is_zero { static const bool value = false; };
        using ret_type = typename tensor_expr_type::ret_type;

        template<typename ... TheArgsType>
        tensor( this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            M_expr( expr )
            {
                this->initTensor( expr,geom, theInitArgs... );
            }
#if 0
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            M_expr( expr )
            {
                this->initTensor( std::true_type{}, exprExpanded, ttse, expr, geom, theInitArgs... );
            }
#endif

        tensor( tensor const& t ) = default;

        void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        {
            this->update( geom );
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )
        {
            this->update( geom );
        }
        void update( Geo_t const& geom )
        {
            M_tensor->update( geom );
        }

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
            {
                M_tensor->update(  std::true_type{}, *(exprExpanded.exprFirstPiolaKirchhoffBase()), ttse, geom, theUpdateArgs... );
            }

        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const
        {
            return M_tensor->evalijq( i,j,q );
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensor->evalijq( i,j,c1,c2,q );
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensor->evaliq( i,c1,c2,q );
        }
        ret_type
        evaliq( uint16_type i, uint16_type q ) const
        {
            return M_tensor->evaliq( i, q );
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensor->evalq( c1,c2,q );
        }
        ret_type
        evalq( uint16_type q ) const
        {
            return M_tensor->evalq( q );
        }
    private :
        template<typename... TheArgsType>
        void initTensor( this_type const& expr, const TheArgsType&... theInitArgs )
            {
                M_tensor = expr.expr().template evaluator<Geo_t, Basis_i_t, Basis_j_t>( theInitArgs... );
            }
#if 0
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void initTensor( std::true_type /**/,TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                         this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            {
                M_tensor = expr.expr().template evaluator<Geo_t, Basis_i_t, Basis_j_t>( std::true_type{}, *(exprExpanded.exprFirstPiolaKirchhoffBase()), ttse, geom, theInitArgs... );
            }
#endif
    private:
        this_type const& M_expr;
        std::shared_ptr<tensor_expr_type> M_tensor;
    };

private :
    std::shared_ptr<expr_type> M_expr;
};



template<ExprApplyType ExprApplied,typename ElementVelocityType,typename ModelPhysicFluidType,typename ElementVelocityLinearizationType>
inline
auto
fluidMecConvectiveTermImpl( ElementVelocityType const& u, std::shared_ptr<ModelPhysicFluidType> const& physicFluidData,
                            ElementVelocityLinearizationType const& beta_u, bool applyLinearization )
{
    using fe_velocity_type = typename unwrap_ptr_t<ElementVelocityType>::functionspace_type::reference_element_type;
    using expr_evaluate_velocity_opertors_type = ExprEvaluateFieldOperators<ElementVelocityType>;
    auto exprEvaluateVelocityOperators = std::make_shared<expr_evaluate_velocity_opertors_type>( u );

    using expr_evaluate_velocity_linearization_opertors_type = ExprEvaluateFieldOperators<ElementVelocityLinearizationType>;
    std::shared_ptr<expr_evaluate_velocity_linearization_opertors_type> exprEvaluateVelocityLinearizationOperators;
    if ( applyLinearization )
    {
        if ( static_cast<void const*>( std::addressof(u) ) != static_cast<void const*>( std::addressof( beta_u ) ) )
            exprEvaluateVelocityLinearizationOperators = std::make_shared<expr_evaluate_velocity_linearization_opertors_type>( beta_u );
        else if constexpr( std::is_same_v<expr_evaluate_velocity_opertors_type,expr_evaluate_velocity_linearization_opertors_type> )
            exprEvaluateVelocityLinearizationOperators = exprEvaluateVelocityOperators;
        CHECK( exprEvaluateVelocityLinearizationOperators ) << "something wrong";
    }

    using expr_convectionterm_base_type = FluidMecConvectiveTermBase<expr_evaluate_velocity_opertors_type,fe_velocity_type,ModelPhysicFluidType,expr_evaluate_velocity_linearization_opertors_type,ExprApplied>;
    using expr_convectionterm_type = PolymorphicExpr<expr_convectionterm_base_type>;

    return Expr<expr_convectionterm_type>( expr_convectionterm_type{ expr_convectionterm_base_type::New( exprEvaluateVelocityOperators, physicFluidData, exprEvaluateVelocityLinearizationOperators, applyLinearization ) } );
}


template<typename ElementVelocityType,typename ModelPhysicFluidType,typename ElementVelocityLinearizationType>
inline
auto
fluidMecConvectiveTerm( ElementVelocityType const& u, std::shared_ptr<ModelPhysicFluidType> const& physicFluidData,
                        ElementVelocityLinearizationType const& beta_u, bool applyLinearization )
{
    return fluidMecConvectiveTermImpl<ExprApplyType::EVAL>( u, physicFluidData, beta_u, applyLinearization );
}

template<class ElementVelocityType,typename ModelPhysicFluidType>
inline
auto
fluidMecConvectiveTerm( ElementVelocityType const& u, std::shared_ptr<ModelPhysicFluidType> const& physicFluidData )
{
    return fluidMecConvectiveTerm( u, physicFluidData, u, false );
}

template<class ElementVelocityType,typename ModelPhysicFluidType,typename ElementVelocityLinearizationType>
inline
auto
fluidMecConvectiveTermJacobian( ElementVelocityType const& u, std::shared_ptr<ModelPhysicFluidType> const& physicFluidData,
                                ElementVelocityLinearizationType const& beta_u, bool applyLinearization )
{
    return fluidMecConvectiveTermImpl<ExprApplyType::JACOBIAN>( u, physicFluidData, beta_u, applyLinearization );
}

template<class ElementVelocityType,typename ModelPhysicFluidType>
inline
auto
fluidMecConvectiveTermJacobian( ElementVelocityType const& u, std::shared_ptr<ModelPhysicFluidType> const& physicFluidData, bool applyLinearization = false )
{
    return fluidMecConvectiveTermJacobian( u,physicFluidData,u,applyLinearization );
}

} // namespace FeelModels
} // namespace Feel

#endif

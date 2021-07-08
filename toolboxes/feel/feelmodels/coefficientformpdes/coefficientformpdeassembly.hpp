/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#ifndef FEELPP_TOOLBOXES_COEFFICIENTFORMPDE_ASSEMBLY_HPP
#define FEELPP_TOOLBOXES_COEFFICIENTFORMPDE_ASSEMBLY_HPP 1

#include <feel/feelmodels/coefficientformpdes/coefficientformpde.hpp>
#include <feel/feelmodels/modelcore/diffsymbolicexpr.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisUnknownType>
template <typename ModelContextType>
void
CoefficientFormPDE<ConvexType,BasisUnknownType>::updateLinearPDE( ModelAlgebraic::DataUpdateLinear & data, ModelContextType const& mctx ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("CoefficientFormPDE","updateLinearPDE", "start"+sc);
    this->timerTool("Solve").start();

    double timeSteppingScaling = 1.;
    if ( !this->isStationary() )
    {
        if ( this->timeStepping() == "Theta" )
            timeSteppingScaling = this->M_timeStepThetaValue;
        data.addDoubleInfo( prefixvm(this->prefix(),"time-stepping.scaling"), timeSteppingScaling );
    }

    auto const& se = mctx.symbolsExpr();

    auto mesh = this->mesh();
    auto Xh = this->spaceUnknown();
    auto const& u = mctx.field( FieldTag::unknown(this), this->unknownName() );
    auto const& v = this->fieldUnknown();

    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );
    auto linearForm = form1( _test=Xh, _vector=F,
                             _rowstart=this->rowStartInVector() );

    for ( std::string const& matName : this->materialsProperties()->physicToMaterials( this->physicDefault() ) )
    {
        auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );

        // Convection
        if ( this->materialsProperties()->hasProperty( matName, this->convectionCoefficientName() ) )
        {
            auto const& coeff_beta = this->materialsProperties()->materialProperty( matName, this->convectionCoefficientName() );
            auto coeff_beta_expr = expr( coeff_beta.template expr<nDim,1>(), se );
            bool build_ConvectionTerm = coeff_beta_expr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
            if ( build_ConvectionTerm )
            {
                bilinearForm +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*inner( gradt(u)*coeff_beta_expr, id(v) ),
                               _geomap=this->geomap() );
            }
        }

        // Diffusion
        if ( this->materialsProperties()->hasProperty( matName, this->diffusionCoefficientName() ) )
        {
            auto const& coeff_c = this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() );
            if ( coeff_c.template hasExpr<nDim,nDim>() )
            {
                auto coeff_c_expr = expr( coeff_c.template expr<nDim,nDim>(), se );
                bool build_DiffusionTerm = coeff_c_expr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
                if ( build_DiffusionTerm )
                {
                    if constexpr ( unknown_is_vectorial )
                        {
                            bilinearForm +=
                                integrate( _range=range,
                                           _expr= timeSteppingScaling*inner(coeff_c_expr*gradt(u),grad(v)),
                                           _geomap=this->geomap() );
                        }
                    else
                    {
                            bilinearForm +=
                                integrate( _range=range,
                                           _expr= timeSteppingScaling*grad(v)*(coeff_c_expr*trans(gradt(u))),
                                           _geomap=this->geomap() );
                    }
                }
            }
            else
            {
                auto coeff_c_expr = expr( coeff_c.expr(), se );
                bool build_DiffusionTerm = coeff_c_expr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
                if ( build_DiffusionTerm )
                {
                    bilinearForm +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*coeff_c_expr*inner(gradt(u),grad(v)),
                                   _geomap=this->geomap() );
                }
            }
        }

        // conservative flux source
        if ( this->materialsProperties()->hasProperty( matName, this->conservativeFluxSourceCoefficientName() ) )
        {
            auto const& coeff_gamma = this->materialsProperties()->materialProperty( matName, this->conservativeFluxSourceCoefficientName() );
            auto getexpr = [&]( auto const& c ) {
                if constexpr ( unknown_is_scalar )
                {
                    auto coeff_gamma_expr = expr( c.template expr<nDim, 1>(), se );
                    return std::pair{
                        trans( coeff_gamma_expr ), coeff_gamma_expr.expression().isNumericExpression()};
                }
                else if constexpr ( unknown_is_vectorial )
                {
                    auto coeff_gamma_expr = expr( c.template expr<nDim, nDim>(), se );
                    return std::pair{ coeff_gamma_expr, coeff_gamma_expr.expression().isNumericExpression() };
                }
            };
            auto [coeff_gamma_expr,is_numeric] = getexpr( coeff_gamma );
            bool build_conservativeFluxSourceTerm = is_numeric? buildCstPart : buildNonCstPart;
            if ( build_conservativeFluxSourceTerm )
            {
                linearForm +=
                    integrate( _range=range,
                                _expr= timeSteppingScaling*inner(coeff_gamma_expr,grad(v)),
                                _geomap=this->geomap() );
            }
        }

        // conservative flux convection
        if ( this->materialsProperties()->hasProperty( matName, this->conservativeFluxConvectionCoefficientName() ) )
        {
            auto const& coeff_alpha = this->materialsProperties()->materialProperty( matName, this->conservativeFluxConvectionCoefficientName() );
            auto coeff_alpha_expr = expr( coeff_alpha.template expr<nDim,1>(), se );
            bool build_conservativeFluxConvectionTerm = coeff_alpha_expr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
            if ( build_conservativeFluxConvectionTerm )
            {
                if constexpr ( unknown_is_scalar )
                {
                    bilinearForm +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*idt(u)*inner(trans(coeff_alpha_expr),grad(v)),
                                   _geomap=this->geomap() );
                }
                else if constexpr ( unknown_is_vectorial )
                {
                    bilinearForm +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*inner(idt(u)*trans(coeff_alpha_expr),grad(v)),
                                   _geomap=this->geomap() );
                }
            }
        }

        // Reaction
        if ( this->materialsProperties()->hasProperty( matName, this->reactionCoefficientName() ) )
        {
            auto const& coeff_a = this->materialsProperties()->materialProperty( matName, this->reactionCoefficientName() );
            auto coeff_a_expr = expr( coeff_a.expr(), se );
            bool build_ReactionTerm = coeff_a_expr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
            if ( build_ReactionTerm )
            {
                bilinearForm +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*coeff_a_expr*inner(idt(u), id(v)),
                               _geomap=this->geomap() );
            }
        }

        // Transient terms
        if ( !this->isStationary() && this->materialsProperties()->hasProperty( matName, this->firstTimeDerivativeCoefficientName() ) )
        {
            auto const& coeff_d = this->materialsProperties()->materialProperty( matName, this->firstTimeDerivativeCoefficientName() );
            auto coeff_d_expr = expr( coeff_d.expr(), se );

            bool build_Form2TransientTerm = buildNonCstPart;
            bool build_Form1TransientTerm = buildNonCstPart;
            if ( coeff_d_expr.expression().isNumericExpression() && this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
            {
                build_Form2TransientTerm = buildCstPart;
            }

            if (build_Form2TransientTerm)
            {
                bilinearForm +=
                    integrate( _range=range,
                               _expr=coeff_d_expr*M_bdfUnknown->polyDerivCoefficient(0)*inner(idt(u),id(v)),
                               _geomap=this->geomap() );
            }
            if (build_Form1TransientTerm)
            {
                linearForm +=
                    integrate( _range=range,
                               _expr=coeff_d_expr*inner(idv(M_bdfUnknown->polyDeriv()),id(v)),
                               _geomap=this->geomap() );
            }
        }

        // Source
        if ( this->materialsProperties()->hasProperty( matName, this->sourceCoefficientName() ) )
        {
            auto const& coeff_f = this->materialsProperties()->materialProperty( matName, this->sourceCoefficientName() );
            auto coeff_f_expr = hana::eval_if( hana::bool_c<unknown_is_scalar>,
                                               [&coeff_f,&se] { return expr( coeff_f.expr(), se ); },
                                               [&coeff_f,&se] { return expr( coeff_f.template expr<nDim,1>(), se ); } );
            //auto const& coeff_f_expr = expr( coeff_f.expr(), se );
            bool buildSourceTerm = coeff_f_expr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
            if ( buildSourceTerm )
            {
                linearForm +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*inner(coeff_f_expr,id(v)),
                               _geomap=this->geomap() );
            }
        }

        // Stab
        if constexpr ( is_lagrange_polynomialset_v<typename space_unknown_type::fe_type> && space_unknown_type::fe_type::isContinuous )
        {
            if ( this->M_applyStabilization && buildNonCstPart )
            {
                this->updateLinearPDEStabilizationGLS( data, mctx, matName, range );
            }
        }

        if constexpr ( is_lagrange_polynomialset_v<typename space_unknown_type::fe_type> && !space_unknown_type::fe_type::isContinuous )
        {

            auto bilinearFormPatternExtended = form2( _test=Xh,_trial=Xh,_matrix=A,
                                                      _pattern=size_type(Pattern::EXTENDED),
                                                      _rowstart=this->rowStartInMatrix(),
                                                      _colstart=this->colStartInMatrix() );
            if ( this->materialsProperties()->hasProperty( matName, this->convectionCoefficientName() ) )
            {
                auto const& coeff_beta = this->materialsProperties()->materialProperty( matName, this->convectionCoefficientName() );
                auto coeff_beta_expr = expr( coeff_beta.template expr<nDim,1>(), se );
                bool build_ConvectionTerm = coeff_beta_expr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
                if ( build_ConvectionTerm )
                {
                    // see ERN, A. and GUERMOND, J. L. Discontinuous Galerkin methods for Friedrichs’ symmetric systems. I. General theory. SIAM J. Numer. Anal.
                    double alpha=0.5;
                    bilinearFormPatternExtended +=
                        integrate( _range=internalfaces(mesh), // TODO here
                                   _expr= timeSteppingScaling*(alpha*abs(inner(N(),coeff_beta_expr))-0.5*inner(leftface(N())+rightface(N()),coeff_beta_expr) )*inner( jumpt( idt(u) ), jump( id(v) ) ),
                                   _geomap=this->geomap() );

                    bilinearForm +=
                        integrate( _range=boundaryfaces(mesh), // TODO here : only on material + precompute before the range
                                   _expr= -timeSteppingScaling*( inner(N(),coeff_beta_expr)<0 )*inner(N(),coeff_beta_expr)*inner(idt(u),id(v)),
                                   _geomap=this->geomap() );
                }
            }
        }


    } // for each material


    // update weak bc
    if ( buildNonCstPart )
    {
        // k \nabla u  n = g
        for( auto const& d : this->M_bcNeumann )
        {
            auto theExpr = hana::eval_if( hana::bool_c<unknown_is_scalar>,
                                          [&d,&se] { return expression( d,se ); },
                                          [&d,&se] { return expression( d,se )*N(); } );
            linearForm +=
                integrate( _range=markedfaces(this->mesh(),M_bcNeumannMarkerManagement.markerNeumannBC(MarkerManagementNeumannBC::NeumannBCShape::SCALAR,name(d)) ),
                           _expr= timeSteppingScaling*inner(theExpr,id(v)),
                           _geomap=this->geomap() );
        }

        // k \nabla u  n + g u = h  -> - k \nabla u  n = gu - h
        for( auto const& d : this->M_bcRobin )
        {
            if constexpr( unknown_is_scalar )
            {
                auto theExpr1 = expression1( d,se );
                bilinearForm +=
                    integrate( _range=markedfaces(mesh,M_bcRobinMarkerManagement.markerRobinBC( name(d) ) ),
                               _expr= timeSteppingScaling*theExpr1*idt(v)*id(v),
                               _geomap=this->geomap() );
                auto theExpr2 = expression2( d,se );
                linearForm +=
                    integrate( _range=markedfaces(mesh,M_bcRobinMarkerManagement.markerRobinBC( name(d) ) ),
                               _expr= timeSteppingScaling*theExpr2*id(v),
                               _geomap=this->geomap() );
            }
            else
                CHECK( false ) << "robin bc with vectorial unknown can not be implemented for now";
        }
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("CoefficientFormPDE","updateLinearPDE",
              "finish in "+(boost::format("%1% s") % timeElapsed).str() );

}

namespace detail
{

// TODO : create struct for DistributeMarker
template <ElementsType ET>
constexpr int indexInDistributeMarker()
{
    if constexpr ( ET == MESH_FACES )
        return 0;
    else if constexpr ( ET == MESH_EDGES )
        return 1;
    else if constexpr ( ET == MESH_POINTS )
        return 2;
    else //if constexpr ( ET == MESH_ELEMENTS )
        return 3;
}

template <ElementsType ET, typename MeshType, typename MarkerType>
auto rangeOfMarkedEntity( MeshType const& mesh, MarkerType const& markers )
{
    if constexpr ( ET == MESH_FACES )
        return markedfaces(mesh,markers );
    else if constexpr ( ET == MESH_EDGES )
        return markededges(mesh,markers );
    else if constexpr ( ET == MESH_POINTS )
        return markedpoints(mesh,markers );
    else //if constexpr ( ET == MESH_ELEMENTS )
        return markedelements(mesh,markers );
}

template <ElementsType ET, typename BfType, typename MeshType, typename EltType, typename SymbolsExprType, typename BcType, typename BcCompType>
void
applyDofEliminationLinear( BfType& bilinearForm, vector_ptrtype& F, MeshType const& mesh, EltType const& u, SymbolsExprType const& se,
                           std::map<ComponentType, std::map<std::string, std::tuple< std::set<std::string>,std::set<std::string>,std::set<std::string>,std::set<std::string> > > > const& M_meshMarkersDofEliminationUnknown,
                           BcType const& M_bcDirichlet, BcCompType const& M_bcDirichletComponents )
{
    static const bool unknown_is_vectorial = EltType::functionspace_type::is_vectorial;
    static const int indexDistrib = indexInDistributeMarker<ET>();
    for ( auto const& [comp,mapMarkerBCToEntitiesMeshMarker] : M_meshMarkersDofEliminationUnknown )
    {
        if ( comp == ComponentType::NO_COMPONENT )
        {
            for( auto const& d : M_bcDirichlet )
            {
                auto itFindMarker = mapMarkerBCToEntitiesMeshMarker.find( name(d) );
                if ( itFindMarker == mapMarkerBCToEntitiesMeshMarker.end() )
                    continue;
                auto const& listMarkedEntities = std::get</*0*/indexDistrib >( itFindMarker->second );
                if ( listMarkedEntities.empty() )
                    continue;
                auto theExpr = expression(d,se);
                bilinearForm +=
                    on( _range=rangeOfMarkedEntity<ET>(mesh,listMarkedEntities),
                        _element=u,_rhs=F,_expr=theExpr );
            }
        }
        else if constexpr ( unknown_is_vectorial )
        {
            auto itFindComp = M_bcDirichletComponents.find( comp );
            if ( itFindComp == M_bcDirichletComponents.end() )
                continue;
            for( auto const& d : itFindComp->second )
            {
                auto itFindMarker = mapMarkerBCToEntitiesMeshMarker.find( name(d) );
                if ( itFindMarker == mapMarkerBCToEntitiesMeshMarker.end() )
                    continue;
                auto const& listMarkedEntities = std::get</*0*/indexDistrib>( itFindMarker->second );
                if ( listMarkedEntities.empty() )
                    continue;
                auto theExpr = expression(d,se);
                bilinearForm +=
                    on( _range=rangeOfMarkedEntity<ET>(mesh, listMarkedEntities),
                        _element=u.comp( comp ),_rhs=F,_expr=theExpr );
            }
        }
    }
}

template <ElementsType ET, typename MeshType, typename EltType, typename SymbolsExprType, typename BcType, typename BcCompType>
void
applyNewtonInitialGuess( MeshType const& mesh, EltType & u, SymbolsExprType const& se,
                         std::map<ComponentType, std::map<std::string, std::tuple< std::set<std::string>,std::set<std::string>,std::set<std::string>,std::set<std::string> > > > const& M_meshMarkersDofEliminationUnknown,
                         BcType const& M_bcDirichlet, BcCompType const& M_bcDirichletComponents )
{
    static const bool unknown_is_vectorial = EltType::functionspace_type::is_vectorial;
    static const int indexDistrib = indexInDistributeMarker<ET>();
    for ( auto const& [comp,mapMarkerBCToEntitiesMeshMarker] : M_meshMarkersDofEliminationUnknown )
    {
        if ( comp == ComponentType::NO_COMPONENT )
        {
            for( auto const& d : M_bcDirichlet )
            {
                auto itFindMarker = mapMarkerBCToEntitiesMeshMarker.find( name(d) );
                if ( itFindMarker == mapMarkerBCToEntitiesMeshMarker.end() )
                    continue;
                auto const& listMarkedEntities = std::get</*0*/indexDistrib>( itFindMarker->second );
                if ( listMarkedEntities.empty() )
                    continue;
                auto theExpr = expression(d,se);
                u.on(_range=rangeOfMarkedEntity<ET>(mesh, listMarkedEntities),
                     _expr=theExpr );
            }
        }
        else if constexpr ( unknown_is_vectorial )
        {
            auto itFindComp = M_bcDirichletComponents.find( comp );
            if ( itFindComp == M_bcDirichletComponents.end() )
                continue;
            for( auto const& d : itFindComp->second )
            {
                auto itFindMarker = mapMarkerBCToEntitiesMeshMarker.find( name(d) );
                if ( itFindMarker == mapMarkerBCToEntitiesMeshMarker.end() )
                    continue;
                auto const& listMarkedEntities = std::get</*0*/indexDistrib>( itFindMarker->second );
                if ( listMarkedEntities.empty() )
                    continue;
                auto theExpr = expression(d,se);
                u.comp( comp ).on(_range=rangeOfMarkedEntity<ET>(mesh, listMarkedEntities),
                                  _expr=theExpr );
            }
        }
    }
}

} // namespace detail

template< typename ConvexType, typename BasisUnknownType>
template <typename ModelContextType>
void
CoefficientFormPDE<ConvexType,BasisUnknownType>::updateLinearPDEDofElimination( ModelAlgebraic::DataUpdateLinear & data, ModelContextType const& mctx ) const
{
    if ( !M_bcDirichletMarkerManagement.hasMarkerDirichletBCelimination() ) return;

    this->log("CoefficientFormPDE","updateLinearPDEDofElimination","start" );

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();

    auto const& se = mctx.symbolsExpr();
    auto mesh = this->mesh();
    auto Xh = this->spaceUnknown();
    auto const& u = this->fieldUnknown();

    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );


    Feel::FeelModels::detail::applyDofEliminationLinear<MESH_ELEMENTS>( bilinearForm, F, mesh, u, se, M_meshMarkersDofEliminationUnknown, M_bcDirichlet, M_bcDirichletComponents );
    Feel::FeelModels::detail::applyDofEliminationLinear<MESH_FACES>( bilinearForm, F, mesh, u, se, M_meshMarkersDofEliminationUnknown, M_bcDirichlet, M_bcDirichletComponents );
    if constexpr ( !is_hcurl_conforming_v<typename space_unknown_type::fe_type> )
    {
    if constexpr ( nDim == 3 )
         Feel::FeelModels::detail::applyDofEliminationLinear<MESH_EDGES>( bilinearForm, F, mesh, u, se, M_meshMarkersDofEliminationUnknown, M_bcDirichlet, M_bcDirichletComponents );
    Feel::FeelModels::detail::applyDofEliminationLinear<MESH_POINTS>( bilinearForm, F, mesh, u, se, M_meshMarkersDofEliminationUnknown, M_bcDirichlet, M_bcDirichletComponents );
    }
    this->log("CoefficientFormPDE","updateLinearPDEDofElimination","finish" );
}


template< typename ConvexType, typename BasisUnknownType>
template <typename ModelContextType>
void
CoefficientFormPDE<ConvexType,BasisUnknownType>::updateNewtonInitialGuess( ModelAlgebraic::DataNewtonInitialGuess & data, ModelContextType const& mctx ) const
{
    if ( !M_bcDirichletMarkerManagement.hasMarkerDirichletBCelimination() ) return;

    this->log("CoefficientFormPDE","updateNewtonInitialGuess","start" );

    vector_ptrtype& U = data.initialGuess();
    auto mesh = this->mesh();
    size_type startBlockIndexUnknown = this->startSubBlockSpaceIndex( this->unknownName() );
    auto u = this->spaceUnknown()->element( U, this->rowStartInVector()+startBlockIndexUnknown );
    auto const& se = mctx.symbolsExpr();

    Feel::FeelModels::detail::applyNewtonInitialGuess<MESH_ELEMENTS>( mesh, u, se, M_meshMarkersDofEliminationUnknown, M_bcDirichlet, M_bcDirichletComponents );
    Feel::FeelModels::detail::applyNewtonInitialGuess<MESH_FACES>( mesh, u, se, M_meshMarkersDofEliminationUnknown, M_bcDirichlet, M_bcDirichletComponents );
    if constexpr ( !is_hcurl_conforming_v<typename space_unknown_type::fe_type> )
    {
    if constexpr ( nDim == 3 )
        Feel::FeelModels::detail::applyNewtonInitialGuess<MESH_EDGES>( mesh, u, se, M_meshMarkersDofEliminationUnknown, M_bcDirichlet, M_bcDirichletComponents );
    Feel::FeelModels::detail::applyNewtonInitialGuess<MESH_POINTS>( mesh, u, se, M_meshMarkersDofEliminationUnknown, M_bcDirichlet, M_bcDirichletComponents );
    }
    // update info for synchronization
    this->updateDofEliminationIds( this->unknownName(), data );

    this->log("CoefficientFormPDE","updateNewtonInitialGuess","finish" );

}

template< typename ConvexType, typename BasisUnknownType>
template <typename ModelContextType>
void
CoefficientFormPDE<ConvexType,BasisUnknownType>::updateJacobian( ModelAlgebraic::DataUpdateJacobian & data, ModelContextType const& mctx ) const
{
    sparse_matrix_ptrtype& J = data.jacobian();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("CoefficientFormPDE","updateJacobian", "start"+sc);

    double timeSteppingScaling = 1.;
    if ( !this->isStationary() )
    {
        if ( this->timeStepping() == "Theta" )
            timeSteppingScaling = this->M_timeStepThetaValue;
        data.addDoubleInfo( prefixvm(this->prefix(),"time-stepping.scaling"), timeSteppingScaling );
    }

    auto const& se = mctx.symbolsExpr();
    auto const& tse = mctx.trialSymbolsExpr();
    auto trialSymbolNames = tse.names();

    auto mesh = this->mesh();
    auto Xh = this->spaceUnknown();
    auto const& u = mctx.field( FieldTag::unknown(this), this->unknownName() );
    auto const& v = this->fieldUnknown();

    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=J,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );

    for ( std::string const& matName : this->materialsProperties()->physicToMaterials( this->physicDefault() ) )
    {
        auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );

        // Convection
        if ( this->materialsProperties()->hasProperty( matName, this->convectionCoefficientName() ) )
        {
            auto const& coeff_beta = this->materialsProperties()->materialProperty( matName, this->convectionCoefficientName() );
            auto coeff_beta_expr = expr( coeff_beta.template expr<nDim,1>(), se );
            bool build_ConvectionTerm = coeff_beta_expr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
            if ( build_ConvectionTerm )
            {
                bilinearForm +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*inner( gradt(u)*coeff_beta_expr, id(v) ),
                               _geomap=this->geomap() );
            }
            bool coeffConvectionDependOnUnknown = coeff_beta.hasSymbolDependency( trialSymbolNames, se );
            if ( coeffConvectionDependOnUnknown && buildNonCstPart )
            {
                hana::for_each( tse.map(), [this,&coeff_beta_expr,&u,&v,&J,&range,&Xh,&timeSteppingScaling]( auto const& e )
                {
                    // NOTE : a strange compilation error related to boost fusion if we use [trialXh,trialBlockIndex] in the loop for
                    for ( auto const& trialSpacePair /*[trialXh,trialBlockIndex]*/ : hana::second(e).blockSpaceIndex() )
                    {
                        auto trialXh = trialSpacePair.first;
                        auto trialBlockIndex = trialSpacePair.second;

                        auto coeff_beta_diff_expr = diffSymbolicExpr( coeff_beta_expr, hana::second(e), trialXh, trialBlockIndex, this->worldComm(), this->repository().expr() );

                        if ( !coeff_beta_diff_expr.expression().hasExpr() )
                            continue;

                        form2( _test=Xh,_trial=trialXh,_matrix=J,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=trialBlockIndex ) +=
                            integrate( _range=range,
                                       _expr= timeSteppingScaling*inner( gradv(u)*coeff_beta_diff_expr, id(v) ),
                                       _geomap=this->geomap() );
                    }
                });
            }
        }

        // Diffusion
        if ( this->materialsProperties()->hasProperty( matName, this->diffusionCoefficientName() ) )
        {
            auto const& coeff_c = this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() );
            bool coeffDiffusionDependOnUnknown = coeff_c.hasSymbolDependency( trialSymbolNames, se );
            if ( coeff_c.template hasExpr<nDim,nDim>() )
            {
                auto coeff_c_expr = expr( coeff_c.template expr<nDim,nDim>(), se );
                bool build_DiffusionTerm = coeff_c_expr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
                if ( build_DiffusionTerm )
                {
                    if constexpr ( unknown_is_vectorial )
                        {
                            bilinearForm +=
                                integrate( _range=range,
                                           _expr= timeSteppingScaling*inner(coeff_c_expr*gradt(u),grad(v)),
                                           _geomap=this->geomap() );
                        }
                    else
                    {
                            bilinearForm +=
                                integrate( _range=range,
                                           _expr= timeSteppingScaling*grad(v)*(coeff_c_expr*trans(gradt(u))),
                                           _geomap=this->geomap() );
                    }
                }
                if ( coeffDiffusionDependOnUnknown && buildNonCstPart )
                {
                    hana::for_each( tse.map(), [this,&coeff_c_expr,&u,&v,&J,&range,&Xh,&timeSteppingScaling]( auto const& e )
                    {
                        // NOTE : a strange compilation error related to boost fusion if we use [trialXh,trialBlockIndex] in the loop for
                        for ( auto const& trialSpacePair /*[trialXh,trialBlockIndex]*/ : hana::second(e).blockSpaceIndex() )
                        {
                            auto trialXh = trialSpacePair.first;
                            auto trialBlockIndex = trialSpacePair.second;

                            auto coeff_c_diff_expr = diffSymbolicExpr( coeff_c_expr, hana::second(e), trialXh, trialBlockIndex, this->worldComm(), this->repository().expr() );

                            if ( !coeff_c_diff_expr.expression().hasExpr() )
                                continue;

                            if constexpr ( unknown_is_vectorial )
                            {
                                form2( _test=Xh,_trial=trialXh,_matrix=J,
                                       _pattern=size_type(Pattern::COUPLED),
                                       _rowstart=this->rowStartInMatrix(),
                                       _colstart=trialBlockIndex ) +=
                                    integrate( _range=range,
                                               _expr= timeSteppingScaling*inner(coeff_c_diff_expr*gradv(u),grad(v)),
                                               _geomap=this->geomap() );
                            }
                            else
                            {
                                form2( _test=Xh,_trial=trialXh,_matrix=J,
                                       _pattern=size_type(Pattern::COUPLED),
                                       _rowstart=this->rowStartInMatrix(),
                                       _colstart=trialBlockIndex ) +=
                                    integrate( _range=range,
                                               _expr= timeSteppingScaling*grad(v)*(coeff_c_diff_expr*trans(gradv(u))),
                                               _geomap=this->geomap() );
                            }
                        }
                    });
                }
            }
            else
            {
                auto coeff_c_expr = expr( coeff_c.expr(), se );
                bool build_DiffusionTerm = coeff_c_expr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
                if ( build_DiffusionTerm )
                {
                    bilinearForm +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*coeff_c_expr*inner(gradt(u),grad(v)),
                                   _geomap=this->geomap() );
                }
                if ( coeffDiffusionDependOnUnknown && buildNonCstPart )
                {
                    hana::for_each( tse.map(), [this,&coeff_c_expr,&u,&v,&J,&range,&Xh,&timeSteppingScaling]( auto const& e )
                    {
                        // NOTE : a strange compilation error related to boost fusion if we use [trialXh,trialBlockIndex] in the loop for
                        for ( auto const& trialSpacePair /*[trialXh,trialBlockIndex]*/ : hana::second(e).blockSpaceIndex() )
                        {
                            auto trialXh = trialSpacePair.first;
                            auto trialBlockIndex = trialSpacePair.second;

                            auto coeff_c_diff_expr = diffSymbolicExpr( coeff_c_expr, hana::second(e), trialXh, trialBlockIndex, this->worldComm(), this->repository().expr() );

                            if ( !coeff_c_diff_expr.expression().hasExpr() )
                                continue;

                            form2( _test=Xh,_trial=trialXh,_matrix=J,
                                   _pattern=size_type(Pattern::COUPLED),
                                   _rowstart=this->rowStartInMatrix(),
                                   _colstart=trialBlockIndex ) +=
                                integrate( _range=range,
                                           _expr= timeSteppingScaling*coeff_c_diff_expr*inner(gradv(u),grad(v)),
                                           _geomap=this->geomap() );
                        }
                    });
                }
            }
        }
        // conservative flux source
        if ( buildNonCstPart && this->materialsProperties()->hasProperty( matName, this->conservativeFluxSourceCoefficientName() ) )
        {
            auto const& coeff_gamma = this->materialsProperties()->materialProperty( matName, this->conservativeFluxSourceCoefficientName() );
            bool coeffConservativeFluxSourceDependOnUnknown = coeff_gamma.hasSymbolDependency( trialSymbolNames, se );
            if ( coeffConservativeFluxSourceDependOnUnknown )
            {
                auto getexpr = [&se]( auto const& c ) {
                    if constexpr ( unknown_is_scalar )
                    {
                        auto coeff_gamma_expr = expr( c.template expr<nDim, 1>(), se );
                        return trans( coeff_gamma_expr );
                    }
                    else if constexpr ( unknown_is_vectorial )
                    {
                        auto coeff_gamma_expr = expr( c.template expr<nDim, nDim>(), se );
                        return coeff_gamma_expr;
                    }
                };
                auto coeff_gamma_expr = getexpr( coeff_gamma );
                hana::for_each( tse.map(), [this, &coeff_gamma_expr, &v, &J, &range, &Xh, &timeSteppingScaling]( auto const& e ) {
                    for ( auto const& trialSpacePair /*[trialXh,trialBlockIndex]*/ : hana::second( e ).blockSpaceIndex() )
                    {
                        auto trialXh = trialSpacePair.first;
                        auto trialBlockIndex = trialSpacePair.second;
                        auto coeff_gamma_diff_expr = diffSymbolicExpr( coeff_gamma_expr, hana::second( e ), trialXh, trialBlockIndex, this->worldComm(), this->repository().expr() );
                        if ( !coeff_gamma_diff_expr.expression().hasExpr() )
                            continue;

                        form2( _test = Xh, _trial = trialXh, _matrix = J,
                                _pattern = size_type( Pattern::COUPLED ),
                                _rowstart = this->rowStartInMatrix(),
                                _colstart = trialBlockIndex ) +=
                            integrate( _range = range,
                                        _expr = -timeSteppingScaling * inner( coeff_gamma_diff_expr, grad( v ) ),
                                        _geomap = this->geomap() );
                    }
                } );
            }
        }

        // conservative flux convection
        if ( this->materialsProperties()->hasProperty( matName, this->conservativeFluxConvectionCoefficientName() ) )
        {
            auto const& coeff_alpha = this->materialsProperties()->materialProperty( matName, this->conservativeFluxConvectionCoefficientName() );
            auto coeff_alpha_expr = expr( coeff_alpha.template expr<nDim,1>(), se );
            bool build_conservativeFluxConvectionTerm = coeff_alpha_expr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;

            auto getConservativeFluxConvectionAssemblyExpr = [&timeSteppingScaling,&v]( auto const& _id_u_expr, auto const& _coeff_alpha_expr ) {
                if constexpr ( unknown_is_scalar )
                    return timeSteppingScaling*_id_u_expr*inner(trans(_coeff_alpha_expr),grad(v));
                else
                    return timeSteppingScaling*inner( _id_u_expr*trans(_coeff_alpha_expr) , grad(v) ); // outer product inside the inner product
            };
            if ( build_conservativeFluxConvectionTerm )
            {
                bilinearForm +=
                    integrate( _range=range,
                               _expr=getConservativeFluxConvectionAssemblyExpr(idt(u),coeff_alpha_expr),//   timeSteppingScaling*idt(u)*inner(trans(coeff_alpha_expr),grad(v)),
                               _geomap=this->geomap() );
            }
            bool coeffConservativeFluxConvectionDependOnUnknown = coeff_alpha.hasSymbolDependency( trialSymbolNames, se );
            if ( coeffConservativeFluxConvectionDependOnUnknown && buildNonCstPart )
            {
                hana::for_each( tse.map(), [this,&coeff_alpha_expr,&u,&v,&J,&range,&Xh,&getConservativeFluxConvectionAssemblyExpr]( auto const& e )
                {
                    for ( auto const& trialSpacePair /*[trialXh,trialBlockIndex]*/ : hana::second(e).blockSpaceIndex() )
                    {
                        auto trialXh = trialSpacePair.first;
                        auto trialBlockIndex = trialSpacePair.second;
                        auto coeff_alpha_diff_expr = diffSymbolicExpr( coeff_alpha_expr, hana::second(e), trialXh, trialBlockIndex, this->worldComm(), this->repository().expr() );
                        if ( !coeff_alpha_diff_expr.expression().hasExpr() )
                            continue;

                        form2( _test=Xh,_trial=trialXh,_matrix=J,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=trialBlockIndex ) +=
                            integrate( _range=range,
                                       _expr=getConservativeFluxConvectionAssemblyExpr(idv(u),coeff_alpha_diff_expr),// timeSteppingScaling*idv(u)*inner(trans(coeff_alpha_diff_expr),grad(v)),
                                       _geomap=this->geomap() );
                    }
                });
            }
        }


        // Reaction
        if ( this->materialsProperties()->hasProperty( matName, this->reactionCoefficientName() ) )
        {
            auto const& coeff_a = this->materialsProperties()->materialProperty( matName, this->reactionCoefficientName() );
            bool coeffReactionDependOnUnknown = coeff_a.hasSymbolDependency( trialSymbolNames, se );
            auto coeff_a_expr = expr( coeff_a.expr(), se );
            bool build_ReactionTerm = coeff_a_expr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
            if ( build_ReactionTerm )
            {
                bilinearForm +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*coeff_a_expr*inner(idt(u), id(v)),
                               _geomap=this->geomap() );
            }
            if ( coeffReactionDependOnUnknown && buildNonCstPart )
            {
                hana::for_each( tse.map(), [this,&coeff_a_expr,&u,&v,&J,&range,&Xh,&timeSteppingScaling]( auto const& e )
                {
                    // NOTE : a strange compilation error related to boost fusion if we use [trialXh,trialBlockIndex] in the loop for
                    for ( auto const& trialSpacePair /*[trialXh,trialBlockIndex]*/ : hana::second(e).blockSpaceIndex() )
                    {
                        auto trialXh = trialSpacePair.first;
                        auto trialBlockIndex = trialSpacePair.second;

                        auto coeff_a_diff_expr = diffSymbolicExpr( coeff_a_expr, hana::second(e), trialXh, trialBlockIndex, this->worldComm(), this->repository().expr() );

                        if ( !coeff_a_diff_expr.expression().hasExpr() )
                            continue;

                        form2( _test=Xh,_trial=trialXh,_matrix=J,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=trialBlockIndex ) +=
                            integrate( _range=range,
                                       _expr= timeSteppingScaling*coeff_a_diff_expr*inner(idv(u), id(v)),
                                       _geomap=this->geomap() );
                    }
                });
            }
        }

        // Transient terms
        if ( !this->isStationary() && this->materialsProperties()->hasProperty( matName, this->firstTimeDerivativeCoefficientName() ) )
        {
            auto const& coeff_d = this->materialsProperties()->materialProperty( matName, this->firstTimeDerivativeCoefficientName() );
            auto coeff_d_expr = expr( coeff_d.expr(), se );

            bool build_Form2TransientTerm = buildNonCstPart;
            if ( coeff_d_expr.expression().isNumericExpression() && this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
            {
                build_Form2TransientTerm = buildCstPart;
            }

            if (build_Form2TransientTerm)
            {
                bilinearForm +=
                    integrate( _range=range,
                               _expr=coeff_d_expr*M_bdfUnknown->polyDerivCoefficient(0)*inner(idt(u),id(v)),
                               _geomap=this->geomap() );
            }
        }

        // Source
        if ( this->materialsProperties()->hasProperty( matName, this->sourceCoefficientName() ) )
        {
            auto const& coeff_f = this->materialsProperties()->materialProperty( matName, this->sourceCoefficientName() );

            bool sourceDependOnUnknown = coeff_f.hasSymbolDependency( trialSymbolNames, se );
            if ( sourceDependOnUnknown && buildNonCstPart )
            {
                auto coeff_f_expr = hana::eval_if( hana::bool_c<unknown_is_scalar>,
                                                   [&coeff_f,&se] { return expr( coeff_f.expr(), se ); },
                                                   [&coeff_f,&se] { return expr( coeff_f.template expr<nDim,1>(), se ); } );

                hana::for_each( tse.map(), [this,&coeff_f_expr,&v,&J,&range,&Xh,&timeSteppingScaling]( auto const& e )
                {
                    // NOTE : a strange compilation error related to boost fusion if we use [trialXh,trialBlockIndex] in the loop for
                    for ( auto const& trialSpacePair /*[trialXh,trialBlockIndex]*/ : hana::second(e).blockSpaceIndex() )
                    {
                        auto trialXh = trialSpacePair.first;
                        auto trialBlockIndex = trialSpacePair.second;

                        auto coeff_f_diff_expr = diffSymbolicExpr( coeff_f_expr, hana::second(e), trialXh, trialBlockIndex, this->worldComm(), this->repository().expr() );
                        if ( !coeff_f_diff_expr.expression().hasExpr() )
                            continue;

                        form2( _test=Xh,_trial=trialXh,_matrix=J,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=trialBlockIndex ) +=
                            integrate( _range=range,
                                       _expr= -timeSteppingScaling*inner(coeff_f_diff_expr,id(v)),
                                       _geomap=this->geomap() );
                    }
                });
            }
        }

        // Stab
        if constexpr ( is_lagrange_polynomialset_v<typename space_unknown_type::fe_type> && space_unknown_type::fe_type::isContinuous )
        {
            if ( this->M_applyStabilization && buildNonCstPart )
            {
                this->updateJacobianStabilizationGLS( data, mctx, matName, range );
            }
        }

        if constexpr ( is_lagrange_polynomialset_v<typename space_unknown_type::fe_type> && !space_unknown_type::fe_type::isContinuous )
        {

            auto bilinearFormPatternExtended = form2( _test=Xh,_trial=Xh,_matrix=J,
                                                      _pattern=size_type(Pattern::EXTENDED),
                                                      _rowstart=this->rowStartInMatrix(),
                                                      _colstart=this->colStartInMatrix() );
            if ( this->materialsProperties()->hasProperty( matName, this->convectionCoefficientName() ) )
            {
                auto const& coeff_beta = this->materialsProperties()->materialProperty( matName, this->convectionCoefficientName() );
                auto coeff_beta_expr = expr( coeff_beta.template expr<nDim,1>(), se );
                bool build_ConvectionTerm = coeff_beta_expr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
                if ( build_ConvectionTerm )
                {
                    // see ERN, A. and GUERMOND, J. L. Discontinuous Galerkin methods for Friedrichs’ symmetric systems. I. General theory. SIAM J. Numer. Anal.
                    double alpha=0.5;
                    bilinearFormPatternExtended +=
                        integrate( _range=internalfaces(mesh), // TODO here
                                   _expr= timeSteppingScaling*(alpha*abs(inner(N(),coeff_beta_expr))-0.5*inner(leftface(N())+rightface(N()),coeff_beta_expr) )*inner( jumpt( idt(u) ), jump( id(v) ) ),
                                   _geomap=this->geomap() );

                    bilinearForm +=
                        integrate( _range=boundaryfaces(mesh), // TODO here : only on material + precompute before the range
                                   _expr= -timeSteppingScaling*( inner(N(),coeff_beta_expr)<0 )*inner(N(),coeff_beta_expr)*inner(idt(u),id(v)),
                                   _geomap=this->geomap() );
                }
            }
        }


    } // for each material


    // Neumann bc :  k \nabla u  n = g
    if ( buildNonCstPart )
    {
        for( auto const& d : M_bcNeumann )
        {
            auto neumannExprBase = expression( d );
            bool neumannnBcDependOnUnknown = neumannExprBase.hasSymbolDependency( trialSymbolNames, se );
            if ( neumannnBcDependOnUnknown )
            {
                auto neumannExpr = expr( neumannExprBase, se );
                hana::for_each( tse.map(), [this,&d,&neumannExpr,&u,&v,&J,&Xh,&timeSteppingScaling]( auto const& e )
                {
                    for ( auto const& trialSpacePair : hana::second(e).blockSpaceIndex() )
                    {
                        auto trialXh = trialSpacePair.first;
                        auto trialBlockIndex = trialSpacePair.second;

                        auto neumannDiffExpr = diffSymbolicExpr( neumannExpr, hana::second(e), trialXh, trialBlockIndex, this->worldComm(), this->repository().expr() );

                        if ( !neumannDiffExpr.expression().hasExpr() )
                            continue;

                        auto neumannDiffExprUsed = hana::eval_if( hana::bool_c<unknown_is_scalar>,
                                                                  [&neumannDiffExpr] { return neumannDiffExpr; },
                                                                  [&neumannDiffExpr] { return neumannDiffExpr*N(); } );

                        form2( _test=Xh,_trial=trialXh,_matrix=J,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=trialBlockIndex ) +=
                            integrate( _range=markedfaces(this->mesh(),M_bcNeumannMarkerManagement.markerNeumannBC(MarkerManagementNeumannBC::NeumannBCShape::SCALAR,name(d)) ),
                                       _expr= -timeSteppingScaling*inner(neumannDiffExprUsed, id(v)),
                                       _geomap=this->geomap() );
                    }
                });
            }
        }
    }

    // Robin bc
    for( auto const& d : this->M_bcRobin )
    {
        if constexpr( unknown_is_scalar )
            {
                auto theExpr1 = expression1( d,se );
                bool build_BcRobinTermLhs = theExpr1.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
                if ( !build_BcRobinTermLhs )
                    continue;

                bilinearForm +=
                    integrate( _range=markedfaces(mesh,M_bcRobinMarkerManagement.markerRobinBC( name(d) ) ),
                               _expr= timeSteppingScaling*theExpr1*idt(u)*id(v),
                               _geomap=this->geomap() );

                if ( theExpr1.hasSymbolDependency( trialSymbolNames ) )
                    CHECK( false ) << "TODO";
                auto theExpr2base = expression2( d );
                if ( theExpr2base.hasSymbolDependency( trialSymbolNames, se ) )
                    CHECK( false ) << "TODO";
            }
        else
            CHECK( false ) << "robin bc with vectorial unknown can not be implemented for now";
    }


}

template< typename ConvexType, typename BasisUnknownType>
template <typename ModelContextType>
void
CoefficientFormPDE<ConvexType,BasisUnknownType>::updateResidual( ModelAlgebraic::DataUpdateResidual & data, ModelContextType const& mctx ) const
{
    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("CoefficientFormPDE","updateResidual", "start"+sc);

    bool timeSteppingEvaluateResidualWithoutTimeDerivative = data.hasInfo( "time-stepping.evaluate-residual-without-time-derivative" );
    if ( timeSteppingEvaluateResidualWithoutTimeDerivative )
    {
        if ( this->isStationary() )
            return;
        if ( this->timeStepping() != "Theta" )
            return;
    }

    double timeSteppingScaling = 1.;
    if ( !this->isStationary() )
    {
        if ( this->timeStepping() == "Theta" )
        {
            if ( timeSteppingEvaluateResidualWithoutTimeDerivative )
                timeSteppingScaling = 1. - this->M_timeStepThetaValue;
            else
                timeSteppingScaling = this->M_timeStepThetaValue;
        }
        data.addDoubleInfo( prefixvm(this->prefix(),"time-stepping.scaling"), timeSteppingScaling );
    }

    auto const& se = mctx.symbolsExpr();
    auto const& tse = mctx.trialSymbolsExpr();
    auto trialSymbolNames = tse.names();
    auto mesh = this->mesh();
    auto Xh = this->spaceUnknown();
    auto const& u = mctx.field( FieldTag::unknown(this), this->unknownName() );
    auto const& v = this->fieldUnknown();

    auto linearForm = form1( _test=Xh, _vector=R,
                             _rowstart=this->rowStartInVector() );

    for ( std::string const& matName : this->materialsProperties()->physicToMaterials( this->physicDefault() ) )
    {
        auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );

        // Convection
        if ( buildNonCstPart && this->materialsProperties()->hasProperty( matName, this->convectionCoefficientName() ) )
        {
            auto const& coeff_beta = this->materialsProperties()->materialProperty( matName, this->convectionCoefficientName() );
            auto coeff_beta_expr = expr( coeff_beta.template expr<nDim,1>(), se );
            bool build_ConvectionTerm = coeff_beta_expr.expression().isNumericExpression()? buildNonCstPart && !UseJacobianLinearTerms : buildNonCstPart;
            if ( build_ConvectionTerm )
            {
                linearForm +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*inner( gradv(u)*coeff_beta_expr, id(v) ),
                               _geomap=this->geomap() );
            }
        }

        // Diffusion
        if ( buildNonCstPart && this->materialsProperties()->hasProperty( matName, this->diffusionCoefficientName() ) )
        {
            auto const& coeff_c = this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() );
            if ( coeff_c.template hasExpr<nDim,nDim>() )
            {
                auto coeff_c_expr = expr( coeff_c.template expr<nDim,nDim>(), se );
                bool build_DiffusionTerm = coeff_c_expr.expression().isNumericExpression()? buildNonCstPart && !UseJacobianLinearTerms : buildNonCstPart;
                if ( build_DiffusionTerm )
                {
                    if constexpr ( unknown_is_vectorial )
                        {
                            linearForm +=
                                integrate( _range=range,
                                           _expr= timeSteppingScaling*inner(coeff_c_expr*gradv(u),grad(v)),
                                           _geomap=this->geomap() );
                        }
                    else
                    {
                            linearForm +=
                                integrate( _range=range,
                                           _expr= timeSteppingScaling*grad(v)*(coeff_c_expr*trans(gradv(u))),
                                           _geomap=this->geomap() );
                    }
                }
            }
            else
            {
                auto coeff_c_expr = expr( coeff_c.expr(), se );
                bool build_DiffusionTerm = coeff_c_expr.expression().isNumericExpression()?  buildNonCstPart && !UseJacobianLinearTerms : buildNonCstPart;
                if ( build_DiffusionTerm )
                {
                    linearForm +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*coeff_c_expr*inner(gradv(u),grad(v)),
                                   _geomap=this->geomap() );
                }
            }
        }
        // conservative flux source
        if ( this->materialsProperties()->hasProperty( matName, this->conservativeFluxSourceCoefficientName() ) )
        {
            auto const& coeff_gamma = this->materialsProperties()->materialProperty( matName, this->conservativeFluxSourceCoefficientName() );
            auto getexpr = [&se]( auto const& c ){
                if constexpr ( unknown_is_scalar )
                {
                    auto coeff_gamma_expr = expr( c.template expr<nDim, 1>(), se );
                    return std::pair{
                        trans( coeff_gamma_expr ), coeff_gamma_expr.expression().isNumericExpression()};
                }
                else if constexpr ( unknown_is_vectorial )
                {
                    auto coeff_gamma_expr = expr( c.template expr<nDim, nDim>(), se );
                    return std::pair{ coeff_gamma_expr, coeff_gamma_expr.expression().isNumericExpression() };
                }
            };
            auto [coeff_gamma_expr,is_numeric] = getexpr( coeff_gamma );
            bool build_conservativeFluxSourceTerm =  is_numeric ? buildNonCstPart && !UseJacobianLinearTerms : buildNonCstPart;
            if ( build_conservativeFluxSourceTerm )
            {
                linearForm +=
                    integrate( _range = range,
                               _expr = -timeSteppingScaling * inner( coeff_gamma_expr, grad( v ) ),
                               _geomap = this->geomap() );
            }
        }

        // conservative flux convection
        if ( this->materialsProperties()->hasProperty( matName, this->conservativeFluxConvectionCoefficientName() ) )
        {
            auto const& coeff_alpha = this->materialsProperties()->materialProperty( matName, this->conservativeFluxConvectionCoefficientName() );
            auto coeff_alpha_expr = expr( coeff_alpha.template expr<nDim,1>(), se );
            bool build_conservativeFluxConvectionTerm = coeff_alpha_expr.expression().isNumericExpression()? buildNonCstPart && !UseJacobianLinearTerms : buildNonCstPart;
            if ( build_conservativeFluxConvectionTerm )
            {
                if constexpr ( unknown_is_scalar )
                {
                    linearForm +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*idv(u)*inner(trans(coeff_alpha_expr),grad(v)),
                                   _geomap=this->geomap() );
                }
                else if constexpr ( unknown_is_vectorial )
                {
                    linearForm +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*inner(idv(u)*trans(coeff_alpha_expr),grad(v)),
                                   _geomap=this->geomap() );
                }
            }
        }

        // Reaction
        if ( buildNonCstPart && this->materialsProperties()->hasProperty( matName, this->reactionCoefficientName() ) )
        {
            auto const& coeff_a = this->materialsProperties()->materialProperty( matName, this->reactionCoefficientName() );
            auto coeff_a_expr = expr( coeff_a.expr(), se );
            bool build_ReactionTerm = coeff_a_expr.expression().isNumericExpression()?  buildNonCstPart && !UseJacobianLinearTerms : buildNonCstPart;
            if ( build_ReactionTerm )
            {
                linearForm +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*coeff_a_expr*inner(idv(u), id(v)),
                               _geomap=this->geomap() );
            }
        }

        // Transient terms
        if ( !this->isStationary() && !timeSteppingEvaluateResidualWithoutTimeDerivative && this->materialsProperties()->hasProperty( matName, this->firstTimeDerivativeCoefficientName() ) )
        {
            auto const& coeff_d = this->materialsProperties()->materialProperty( matName, this->firstTimeDerivativeCoefficientName() );
            auto coeff_d_expr = expr( coeff_d.expr(), se );

            bool build_TransientTermLhs = buildNonCstPart;
            bool build_Form1TransientRhs = buildCstPart; // TODO : if coeff_d_expr depends on unkwnon, then it's buildNonCstPart
            if ( coeff_d_expr.expression().isNumericExpression() && this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
            {
                build_TransientTermLhs = buildNonCstPart && !UseJacobianLinearTerms;
            }

            if ( build_TransientTermLhs )
            {
                linearForm +=
                    integrate( _range=range,
                               _expr=coeff_d_expr*M_bdfUnknown->polyDerivCoefficient(0)*inner(idv(u),id(v)),
                               _geomap=this->geomap() );
            }
            if (build_Form1TransientRhs)
            {
                linearForm +=
                    integrate( _range=range,
                               _expr=-coeff_d_expr*inner(idv(M_bdfUnknown->polyDeriv()),id(v)),
                               _geomap=this->geomap() );
            }
        }

        // Source
        if ( this->materialsProperties()->hasProperty( matName, this->sourceCoefficientName() ) )
        {
            auto const& coeff_f = this->materialsProperties()->materialProperty( matName, this->sourceCoefficientName() );

            bool sourceDependOnUnknown = coeff_f.hasSymbolDependency( trialSymbolNames, se );

            auto coeff_f_expr = hana::eval_if( hana::bool_c<unknown_is_scalar>,
                                               [&coeff_f,&se] { return expr( coeff_f.expr(), se ); },
                                               [&coeff_f,&se] { return expr( coeff_f.template expr<nDim,1>(), se ); } );
            bool buildSourceTerm = sourceDependOnUnknown ? buildNonCstPart : buildCstPart;

            if ( buildSourceTerm )
            {
                linearForm +=
                    integrate( _range=range,
                               _expr= -timeSteppingScaling*inner(coeff_f_expr,id(v)),
                               _geomap=this->geomap() );
            }
        }

        // Stab
        if constexpr ( is_lagrange_polynomialset_v<typename space_unknown_type::fe_type> && space_unknown_type::fe_type::isContinuous )
        {
            if ( this->M_applyStabilization && buildNonCstPart )
            {
                this->updateResidualStabilizationGLS( data, mctx, matName, range );
            }
        }

        if constexpr ( is_lagrange_polynomialset_v<typename space_unknown_type::fe_type> && !space_unknown_type::fe_type::isContinuous )
        {
            if ( buildNonCstPart && this->materialsProperties()->hasProperty( matName, this->convectionCoefficientName() ) )
            {
                auto const& coeff_beta = this->materialsProperties()->materialProperty( matName, this->convectionCoefficientName() );
                auto coeff_beta_expr = expr( coeff_beta.template expr<nDim,1>(), se );
                bool build_ConvectionTerm = coeff_beta_expr.expression().isNumericExpression()? buildNonCstPart && !UseJacobianLinearTerms : buildNonCstPart;
                if ( build_ConvectionTerm )
                {
                    // see ERN, A. and GUERMOND, J. L. Discontinuous Galerkin methods for Friedrichs’ symmetric systems. I. General theory. SIAM J. Numer. Anal.
                    double alpha=0.5;
                    linearForm +=
                        integrate( _range=internalfaces(mesh), // TODO here
                                   _expr= timeSteppingScaling*(alpha*abs(inner(N(),coeff_beta_expr))-0.5*inner(leftface(N())+rightface(N()),coeff_beta_expr) )*inner( jumpv( idv(u) ), jump( id(v) ) ),
                                   _geomap=this->geomap() );

                    linearForm +=
                        integrate( _range=boundaryfaces(mesh), // TODO here : only on material + precompute before the range
                                   _expr= -timeSteppingScaling*( inner(N(),coeff_beta_expr)<0 )*inner(N(),coeff_beta_expr)*inner(idv(u),id(v)),
                                   _geomap=this->geomap() );
                }
            }
        }


    } // for each material


    // Neumann bc
    for( auto const& d : this->M_bcNeumann )
    {
        auto neumannExprBase = expression( d );
        bool neumannnBcDependOnUnknown = neumannExprBase.hasSymbolDependency( trialSymbolNames, se );
        bool assembleNeumannBcTerm = neumannnBcDependOnUnknown? buildNonCstPart : buildCstPart;
        if ( assembleNeumannBcTerm )
        {
            auto theExpr = hana::eval_if( hana::bool_c<unknown_is_scalar>,
                                          [&neumannExprBase,&se] { return expr( neumannExprBase, se ); },
                                          [&neumannExprBase,&se] { return expr( neumannExprBase, se )*N(); } );
            linearForm +=
                integrate( _range=markedfaces(this->mesh(),M_bcNeumannMarkerManagement.markerNeumannBC(MarkerManagementNeumannBC::NeumannBCShape::SCALAR,name(d)) ),
                           _expr= -timeSteppingScaling*inner(theExpr,id(v)),
                           _geomap=this->geomap() );
        }
    }

    // Robin bc
    for( auto const& d : this->M_bcRobin )
    {
        if constexpr( unknown_is_scalar )
            {
                auto theExpr1 = expression1( d,se );
                bool build_BcRobinTermLhs = theExpr1.expression().isNumericExpression()? buildNonCstPart && !UseJacobianLinearTerms : buildNonCstPart;
                if ( build_BcRobinTermLhs )
                {
                    linearForm +=
                        integrate( _range=markedfaces(mesh,M_bcRobinMarkerManagement.markerRobinBC( name(d) ) ),
                                   _expr= timeSteppingScaling*theExpr1*idv(u)*id(v),
                                   _geomap=this->geomap() );
                }
                auto theExpr2 = expression2( d,se );
                bool build_BcRobinTermRhs = buildCstPart; // TODO : theExpr2 depends on unkwnon, then it's buildNonCstPart
                if ( build_BcRobinTermRhs )
                {
                    linearForm +=
                        integrate( _range=markedfaces(mesh,M_bcRobinMarkerManagement.markerRobinBC( name(d) ) ),
                                   _expr= -timeSteppingScaling*theExpr2*id(v),
                                   _geomap=this->geomap() );
                }
            }
        else
            CHECK( false ) << "robin bc with vectorial unknown can not be implemented for now";
    }

}

} // namespace Feel
} // namespace FeelModels

#endif

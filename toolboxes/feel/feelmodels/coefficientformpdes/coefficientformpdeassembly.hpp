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
        if ( M_timeStepping == "Theta" )
            timeSteppingScaling = M_timeStepThetaValue;
        data.addDoubleInfo( prefixvm(this->prefix(),"time-stepping.scaling"), timeSteppingScaling );
    }

    auto const& se = mctx.symbolsExpr();

    auto mesh = this->mesh();
    auto Xh = this->spaceUnknown();
    auto const& u = this->fieldUnknown();
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
        if ( this->M_applyStabilization && buildNonCstPart )
        {
            this->updateLinearPDEStabilizationGLS( data, mctx, matName, range );
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

    // store markers for each entities in order to apply strong bc with priority (points erase edges erace faces)
    std::map<std::string, std::tuple< std::set<std::string>,std::set<std::string>,std::set<std::string>,std::set<std::string> > > mapMarkerBCToEntitiesMeshMarker;
    for( auto const& d : M_bcDirichlet )
    {
        mapMarkerBCToEntitiesMeshMarker[name(d)] =
            detail::distributeMarkerListOnSubEntity(mesh,M_bcDirichletMarkerManagement.markerDirichletBCByNameId( "elimination",name(d) ) );
    }

    for( auto const& d : M_bcDirichlet )
    {
        auto itFindMarker = mapMarkerBCToEntitiesMeshMarker.find( name(d) );
        if ( itFindMarker == mapMarkerBCToEntitiesMeshMarker.end() )
            continue;
        auto const& listMarkedFaces = std::get<0>( itFindMarker->second );
        if ( listMarkedFaces.empty() )
            continue;
        auto theExpr = expression(d,se);
        bilinearForm +=
            on( _range=markedfaces(mesh, listMarkedFaces ),
                _element=u,_rhs=F,_expr=theExpr );
    }
    if constexpr ( nDim == 3 )
    {
        for( auto const& d : M_bcDirichlet )
        {
            auto itFindMarker = mapMarkerBCToEntitiesMeshMarker.find( name(d) );
            if ( itFindMarker == mapMarkerBCToEntitiesMeshMarker.end() )
                continue;
            auto const& listMarkedEdges = std::get<1>( itFindMarker->second );
            if ( listMarkedEdges.empty() )
                continue;
            auto theExpr = expression(d,se);
            bilinearForm +=
                on( _range=markededges(mesh, listMarkedEdges ),
                    _element=u,_rhs=F,_expr=theExpr );
        }
    }
    for( auto const& d : M_bcDirichlet )
    {
        auto itFindMarker = mapMarkerBCToEntitiesMeshMarker.find( name(d) );
        if ( itFindMarker == mapMarkerBCToEntitiesMeshMarker.end() )
            continue;
        auto const& listMarkedPoints = std::get<2>( itFindMarker->second );
        if ( listMarkedPoints.empty() )
            continue;
        auto theExpr = expression(d,se);
        bilinearForm +=
            on( _range=markedpoints(mesh, listMarkedPoints ),
                _element=u,_rhs=F,_expr=theExpr );
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

    // store markers for each entities in order to apply strong bc with priority (points erase edges erace faces)
    std::map<std::string, std::tuple< std::set<std::string>,std::set<std::string>,std::set<std::string>,std::set<std::string> > > mapMarkerBCToEntitiesMeshMarker;
    for( auto const& d : M_bcDirichlet )
    {
        mapMarkerBCToEntitiesMeshMarker[name(d)] =
            detail::distributeMarkerListOnSubEntity(mesh,M_bcDirichletMarkerManagement.markerDirichletBCByNameId( "elimination",name(d) ) );
    }

    for( auto const& d : M_bcDirichlet )
    {
        auto itFindMarker = mapMarkerBCToEntitiesMeshMarker.find( name(d) );
        if ( itFindMarker == mapMarkerBCToEntitiesMeshMarker.end() )
            continue;
        auto const& listMarkedFaces = std::get<0>( itFindMarker->second );
        if ( listMarkedFaces.empty() )
            continue;
        auto theExpr = expression(d,se);
        u.on(_range=markedfaces(mesh, listMarkedFaces ),
             _expr=theExpr );
    }
    if constexpr ( nDim == 3 )
    {
        for( auto const& d : M_bcDirichlet )
        {
            auto itFindMarker = mapMarkerBCToEntitiesMeshMarker.find( name(d) );
            if ( itFindMarker == mapMarkerBCToEntitiesMeshMarker.end() )
                continue;
            auto const& listMarkedEdges = std::get<1>( itFindMarker->second );
            if ( listMarkedEdges.empty() )
                continue;
            auto theExpr = expression(d,se);
            u.on(_range=markededges(mesh, listMarkedEdges ),
                 _expr=theExpr );
        }
    }
    for( auto const& d : M_bcDirichlet )
    {
        auto itFindMarker = mapMarkerBCToEntitiesMeshMarker.find( name(d) );
        if ( itFindMarker == mapMarkerBCToEntitiesMeshMarker.end() )
            continue;
        auto const& listMarkedPoints = std::get<2>( itFindMarker->second );
        if ( listMarkedPoints.empty() )
            continue;
        auto theExpr = expression(d,se);
        u.on(_range=markedpoints(mesh, listMarkedPoints ),
             _expr=theExpr );
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
        if ( M_timeStepping == "Theta" )
            timeSteppingScaling = M_timeStepThetaValue;
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

                hana::for_each( tse.map(), [this,&coeff_f_expr,&v,&J,&range,&Xh]( auto const& e )
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
                                       _expr= -inner(coeff_f_diff_expr,id(v)),
                                       _geomap=this->geomap() );
                    }
                });
            }
        }

        // Stab
        if ( this->M_applyStabilization && buildNonCstPart )
        {
            this->updateJacobianStabilizationGLS( data, mctx, matName, range );
        }


    } // for each material

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

    double timeSteppingScaling = 1.;
    bool timeSteppingEvaluateResidualWithoutTimeDerivative = false;
    if ( !this->isStationary() )
    {
        timeSteppingEvaluateResidualWithoutTimeDerivative = data.hasInfo( prefixvm( this->prefix(),"time-stepping.evaluate-residual-without-time-derivative") );
        if ( M_timeStepping == "Theta" )
        {
            if ( timeSteppingEvaluateResidualWithoutTimeDerivative )
                timeSteppingScaling = 1. - M_timeStepThetaValue;
            else
                timeSteppingScaling = M_timeStepThetaValue;
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

            bool sourceDependOnUnknown = coeff_f.hasSymbolDependency( trialSymbolNames );

            auto coeff_f_expr = hana::eval_if( hana::bool_c<unknown_is_scalar>,
                                               [&coeff_f,&se] { return expr( coeff_f.expr(), se ); },
                                               [&coeff_f,&se] { return expr( coeff_f.template expr<nDim,1>(), se ); } );
            bool buildSourceTerm = sourceDependOnUnknown ? buildNonCstPart : buildCstPart;

            if ( buildSourceTerm )
            {
                linearForm +=
                    integrate( _range=range,
                               _expr= -inner(coeff_f_expr,id(v)),
                               _geomap=this->geomap() );
            }
        }

        // Stab
        if ( this->M_applyStabilization && buildNonCstPart )
        {
            this->updateResidualStabilizationGLS( data, mctx, matName, range );
        }

    } // for each material


    // Neumann bc
    for( auto const& d : this->M_bcNeumann )
    {
        auto theExpr = hana::eval_if( hana::bool_c<unknown_is_scalar>,
                                      [&d,&se] { return expression( d,se ); },
                                      [&d,&se] { return expression( d,se )*N(); } );
        bool build_BcNeumann = buildCstPart; // TODO : if theExpr depends on unkwnon, then it's buildNonCstPart
        if ( build_BcNeumann )
        {
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

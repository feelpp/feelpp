#ifndef FEELPP_TOOLBOXES_ELECTRIC_ASSEMBLY_HPP
#define FEELPP_TOOLBOXES_ELECTRIC_ASSEMBLY_HPP 1

#include <feel/feelmodels/modelcore/diffsymbolicexpr.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisPotentialType>
template <typename ModelContextType>
void
Electric<ConvexType,BasisPotentialType>::updateLinearPDE( DataUpdateLinear & data, ModelContextType const& mctx ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("Electric","updateLinearPDE", "start"+sc);
    boost::mpi::timer thetimer;

    auto mesh = this->mesh();
    auto XhV = this->spaceElectricPotential();
    auto const& v = this->fieldElectricPotential();
    auto const& symbolsExpr = mctx.symbolsExpr();

    auto bilinearForm_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix() ,
                                              _colstart=this->colStartInMatrix() );
    auto myLinearForm = form1( _test=XhV, _vector=F,
                               _rowstart=this->rowStartInVector() );

    //--------------------------------------------------------------------------------------------------//

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
            auto const& electricConductivity = this->materialsProperties()->electricConductivity( matName );

            auto sigma = expr( electricConductivity.expr(), symbolsExpr);
            //if ( sigma.expression().hasSymbol( "heat_T" ) )
            //    continue;
            //auto sigma = idv(M_electricProperties->fieldElectricConductivity());
            bool buildDiffusion = sigma.expression().isConstant()? buildCstPart : buildNonCstPart;
            if ( buildDiffusion )
            {
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= sigma*inner(gradt(v),grad(v)),
                               _geomap=this->geomap() );
            }
        }
    }

    // update source term
    if ( buildNonCstPart )
    {
#if 0 //VINCENT
        for( auto const& d : this->M_volumicForcesProperties )
        {
            auto rangeEltUsed = (markers(d).empty())? M_rangeMeshElements : markedelements(this->mesh(),markers(d));
            myLinearForm +=
                integrate( _range=rangeEltUsed,
                           _expr= expression(d,symbolsExpr)*id(v),
                           _geomap=this->geomap() );
        }
#endif
    }

    // update weak bc
    if ( buildNonCstPart )
    {
        for ( auto const& [bcName,bcData] : M_boundaryConditions->surfaceChargeDensity() )
        {
            auto theExpr = bcData->expr( symbolsExpr );
            myLinearForm +=
                integrate( _range=markedfaces(mesh,bcData->markers()),
                           _expr=-theExpr*id(v),
                           _geomap=this->geomap() );
        }
    }
}


template< typename ConvexType, typename BasisPotentialType>
template <typename ModelContextType>
void
Electric<ConvexType,BasisPotentialType>::updateLinearPDEDofElimination( DataUpdateLinear & data, ModelContextType const& mctx ) const
{
    if ( !M_boundaryConditions->hasTypeDofElimination() )
        return;

    this->log("Electric","updateLinearPDEDofElimination","start" );

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    auto const& se = mctx.symbolsExpr();
    auto XhV = this->spaceElectricPotential();
    auto const& u = this->fieldElectricPotential();
    auto mesh = XhV->mesh();

    auto bilinearForm = form2( _test=XhV,_trial=XhV,_matrix=A,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );

    M_boundaryConditions->applyDofEliminationLinear( bilinearForm, F, mesh, u, se );

    this->log("Electric","updateLinearPDEDofElimination","finish" );
}

template< typename ConvexType, typename BasisPotentialType>
template <typename ModelContextType>
void
Electric<ConvexType,BasisPotentialType>::updateNewtonInitialGuess( DataNewtonInitialGuess & data, ModelContextType const& mctx ) const
{
    if ( !M_boundaryConditions->hasTypeDofElimination() )
        return;
    this->log("Electric","updateNewtonInitialGuess","start" );

    vector_ptrtype& U = data.initialGuess();
    auto mesh = this->mesh();
    size_type startBlockIndexElectricPotential = this->startSubBlockSpaceIndex( "electric-potential" );
    auto u = this->spaceElectricPotential()->element( U, this->rowStartInVector()+startBlockIndexElectricPotential );
    auto const& se = mctx.symbolsExpr();

    M_boundaryConditions->applyNewtonInitialGuess( mesh, u, se );

    // update info for synchronization
    this->updateDofEliminationIds( "electric-potential", data );

    this->log("Electric","updateNewtonInitialGuess","finish" );
}


template< typename ConvexType, typename BasisPotentialType>
template <typename ModelContextType>
void
Electric<ConvexType,BasisPotentialType>::updateJacobian( DataUpdateJacobian & data, ModelContextType const& mctx ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    vector_ptrtype& RBis = data.vectorUsedInStrongDirichlet();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("Electric","updateJacobian", "start"+sc);
    size_type startBlockIndexElectricPotential = this->startSubBlockSpaceIndex( "electric-potential" );

    auto mesh = this->mesh();
    auto XhV = this->spaceElectricPotential();
    auto const& v = mctx.field( FieldTag::potential(this), "electric-potential" );
    auto const& se = mctx.symbolsExpr();
    auto const& tse = mctx.trialSymbolsExpr();
    auto trialSymbolNames = tse.names();

    auto bilinearForm_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix()+startBlockIndexElectricPotential,
                                              _colstart=this->colStartInMatrix()+startBlockIndexElectricPotential );

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
            auto const& electricConductivity =  this->materialsProperties()->electricConductivity( matName );
            bool electricConductivityDependOnTrialSymbol = electricConductivity.hasSymbolDependency( trialSymbolNames,se );
            auto sigmaExpr = expr(electricConductivity.expr(),se);
            bool buildDiffusion = sigmaExpr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
            if ( buildDiffusion )
            {
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= sigmaExpr*inner(gradt(v),grad(v)),
                               _geomap=this->geomap() );
            }
            if ( electricConductivityDependOnTrialSymbol && buildNonCstPart )
            {
                hana::for_each( tse.map(), [this,&sigmaExpr,&v,&J,&range,&XhV]( auto const& e )
                    {
                        // NOTE : a strange compilation error related to boost fusion if we use [trialXh,trialBlockIndex] in the loop for
                        for ( auto const& trialSpacePair /*[trialXh,trialBlockIndex]*/ : hana::second(e).blockSpaceIndex() )
                        {
                            auto trialXh = trialSpacePair.second;
                            auto trialBlockIndex = trialSpacePair.first;

                            auto sigmaDiffExpr = diffSymbolicExpr( sigmaExpr, hana::second(e), trialXh, trialBlockIndex, this->worldComm(), this->repository().expr() );

                            if ( !sigmaDiffExpr.expression().hasExpr() )
                                continue;

                            form2( _test=XhV,_trial=trialXh,_matrix=J,
                                   _pattern=size_type(Pattern::COUPLED),
                                   _rowstart=this->rowStartInMatrix(),
                                   _colstart=trialBlockIndex ) +=
                                integrate( _range=range,
                                           _expr= sigmaDiffExpr*inner(gradv(v),grad(v)),
                                           _geomap=this->geomap() );
                        }
                    });
            }
        }
    }

}

template< typename ConvexType, typename BasisPotentialType>
template <typename ModelContextType>
void
Electric<ConvexType,BasisPotentialType>::updateResidual( DataUpdateResidual & data, ModelContextType const& mctx ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("Electric","updateResidual", "start"+sc);

    size_type startBlockIndexElectricPotential = this->startSubBlockSpaceIndex( "electric-potential" );

    auto mesh = this->mesh();
    auto XhV = this->spaceElectricPotential();
    // auto const& v = this->fieldElectricPotential();
    //auto const v = XhV->element(XVec, this->rowStartInVector()+startBlockIndexElectricPotential );
    auto const& v = mctx.field( FieldTag::potential(this), "electric-potential" );
    auto const& symbolsExpr = mctx.symbolsExpr();

    auto myLinearForm = form1( _test=XhV, _vector=R,
                               _rowstart=this->rowStartInVector() + startBlockIndexElectricPotential );

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
            auto const& electricConductivity =  this->materialsProperties()->electricConductivity( matName );

            auto sigmaExpr = expr(electricConductivity.expr(),symbolsExpr);
            bool buildDiffusion = sigmaExpr.expression().isNumericExpression()? buildNonCstPart && !UseJacobianLinearTerms : buildNonCstPart;
            if ( buildDiffusion )
            {
                myLinearForm +=
                    integrate( _range=range,
                               _expr= sigmaExpr*inner(gradv(v),grad(v)),
                               _geomap=this->geomap() );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------//
    // source term
    if ( buildCstPart )
    {
#if 0 //VINCENT
        for( auto const& d : this->M_volumicForcesProperties )
        {
            auto rangeEltUsed = (markers(d).empty())? M_rangeMeshElements : markedelements(this->mesh(),markers(d));
            myLinearForm +=
                integrate( _range=rangeEltUsed,
                           _expr= -expression(d,symbolsExpr)*id(v),
                           _geomap=this->geomap() );
        }
#endif
    }

    //--------------------------------------------------------------------------------------------------//
    // weak bc
    if ( buildCstPart )
    {
        for ( auto const& [bcName,bcData] : M_boundaryConditions->surfaceChargeDensity() )
        {
            auto theExpr = bcData->expr( symbolsExpr );
            myLinearForm +=
                integrate( _range=markedfaces(mesh,bcData->markers()),
                           _expr= theExpr*id(v),
                           _geomap=this->geomap() );
        }
    }

    //--------------------------------------------------------------------------------------------------//
    this->log("Electric","updateResidual", "finish"+sc);
}


} // namespace Feel
} // namespace FeelModels

#endif

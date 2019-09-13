#ifndef FEELPP_TOOLBOXES_ELECTRIC_ASSEMBLY_HPP
#define FEELPP_TOOLBOXES_ELECTRIC_ASSEMBLY_HPP 1

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisPotentialType>
template <typename SymbolsExpr>
void
Electric<ConvexType,BasisPotentialType>::updateLinearPDE( DataUpdateLinear & data, SymbolsExpr const& symbolsExpr ) const
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

    auto bilinearForm_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix() ,
                                              _colstart=this->colStartInMatrix() );
    auto myLinearForm = form1( _test=XhV, _vector=F,
                               _rowstart=this->rowStartInVector() );

    //--------------------------------------------------------------------------------------------------//

    for ( auto const& rangeData : M_electricProperties->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& electricConductivity = M_electricProperties->electricConductivity( matName );

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

    // update source term
    if ( !buildCstPart )
    {
        for( auto const& d : this->M_volumicForcesProperties )
        {
            auto rangeEltUsed = (markers(d).empty())? M_rangeMeshElements : markedelements(this->mesh(),markers(d));
            myLinearForm +=
                integrate( _range=rangeEltUsed,
                           _expr= expression(d)*id(v),
                           _geomap=this->geomap() );
        }
    }

    // update bc
    this->updateLinearPDEWeakBC(A,F,buildCstPart);
}


template< typename ConvexType, typename BasisPotentialType>
template <typename SymbolsExpr>
void
Electric<ConvexType,BasisPotentialType>::updateJacobian( DataUpdateJacobian & data, SymbolsExpr const& symbolsExpr ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    vector_ptrtype& RBis = data.vectorUsedInStrongDirichlet();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("Electric","updateJacobian", "start"+sc);
    size_type startBlockIndexElectricPotential = this->startSubBlockSpaceIndex( "potential-electric" );

    auto mesh = this->mesh();
    auto XhV = this->spaceElectricPotential();
    // auto const& v = this->fieldElectricPotential();
    auto const v = XhV->element(XVec, this->rowStartInVector()+startBlockIndexElectricPotential );

    auto bilinearForm_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix()+startBlockIndexElectricPotential,
                                              _colstart=this->colStartInMatrix()+startBlockIndexElectricPotential );

    for ( auto const& rangeData : M_electricProperties->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& electricConductivity = M_electricProperties->electricConductivity( matName );
        bool buildDiffusion = ( electricConductivity.isConstant() )? buildCstPart : buildNonCstPart;
        if ( buildDiffusion )
        {
            auto sigmaExpr = expr(electricConductivity.expr(),symbolsExpr);
            bilinearForm_PatternCoupled +=
                integrate( _range=range,
                           _expr= sigmaExpr*inner(gradt(v),grad(v)),
                           _geomap=this->geomap() );
        }
    }

    this->updateJacobianWeakBC( v,J,buildCstPart );

}

template< typename ConvexType, typename BasisPotentialType>
template <typename SymbolsExpr>
void
Electric<ConvexType,BasisPotentialType>::updateResidual( DataUpdateResidual & data, SymbolsExpr const& symbolsExpr ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool _BuildCstPart = data.buildCstPart();
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();

    bool buildNonCstPart = !_BuildCstPart;
    bool buildCstPart = _BuildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("Electric","updateResidual", "start"+sc);

    size_type startBlockIndexElectricPotential = this->startSubBlockSpaceIndex( "potential-electric" );

    auto mesh = this->mesh();
    auto XhV = this->spaceElectricPotential();
    // auto const& v = this->fieldElectricPotential();
    auto const v = XhV->element(XVec, this->rowStartInVector()+startBlockIndexElectricPotential );


    auto myLinearForm = form1( _test=XhV, _vector=R,
                               _rowstart=this->rowStartInVector() + startBlockIndexElectricPotential );


    for ( auto const& rangeData : M_electricProperties->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& electricConductivity = M_electricProperties->electricConductivity( matName );

        bool buildDiffusion = ( electricConductivity.isConstant() )? buildNonCstPart && !UseJacobianLinearTerms : buildNonCstPart;
        if ( buildDiffusion )
        {
            auto sigmaExpr = expr(electricConductivity.expr(),symbolsExpr);
            myLinearForm +=
                integrate( _range=range,
                           _expr= sigmaExpr*inner(gradv(v),grad(v)),
                           _geomap=this->geomap() );
        }
    }
    // source term
    if ( buildCstPart )
    {
        for( auto const& d : this->M_volumicForcesProperties )
        {
            auto rangeEltUsed = (markers(d).empty())? M_rangeMeshElements : markedelements(this->mesh(),markers(d));
            myLinearForm +=
                integrate( _range=rangeEltUsed,
                           _expr= -expression(d,symbolsExpr)*id(v),
                           _geomap=this->geomap() );
        }
    }

    // weak bc
    this->updateResidualWeakBC( v,R,buildCstPart );

    this->log("Electric","updateResidual", "finish"+sc);
}


} // namespace Feel
} // namespace FeelModels

#endif

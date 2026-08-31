//!

#ifndef FEELPP_TOOLBOXES_ELECTROMAGNETIC_ASSEMBLY_HPP
#define FEELPP_TOOLBOXES_ELECTROMAGNETIC_ASSEMBLY_HPP

namespace Feel::FeelModels
{

template< typename ElectricType, typename MagneticType>
template <typename ModelContextType>
void
Electromagnetic<ElectricType,MagneticType>::updateLinearPDE( DataUpdateLinear & data, ModelContextType const& mctx ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    bool doAssemblyRhs = !data.hasInfo( "ignore-assembly.rhs" );
    bool doAssemblyLhs = !data.hasInfo( "ignore-assembly.lhs" );

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("Electromagnetic","updateLinearPDE", "start"+sc);
    this->timerTool("Solve").start();

    M_electricModel->updateLinearPDE( data,mctx );
    M_magneticModel->updateLinearPDE( data,mctx );

    if ( buildNonCstPart && doAssemblyLhs )
    {
        auto const& se = mctx.symbolsExpr();
        size_type blockIndexElectricPotential = this->startSubBlockSpaceIndex( "electric" ) + 0; // TODO
        size_type blockIndexMagneticVectorPotential = this->startSubBlockSpaceIndex( "magnetic" ) +
            M_magneticModel->startSubBlockSpaceIndex( magnetic_model_type::FieldTag::vectorPotential(M_magneticModel.get()).identifier() );

        auto const& electricV = mctx.field( electric_model_type::FieldTag::potential(this->electricModel().get()),
                                            "electric-potential" );//electric_model_type::FieldTag::potential(this->electricModel().get()).identifier() );
        auto const& magneticA = mctx.field( magnetic_model_type::FieldTag::vectorPotential(M_magneticModel.get()),
                                            magnetic_model_type::FieldTag::vectorPotential(M_magneticModel.get()).identifier() );

       auto mybfAV = form2( _test=M_magneticModel->spaceVectorPotential(),
                            _trial=M_electricModel->spaceElectricPotential(),
                            _matrix=A,
                            _pattern=size_type(Pattern::COUPLED),
                            _rowstart=this->rowStartInMatrix()+blockIndexMagneticVectorPotential ,
                            _colstart=this->colStartInMatrix()+blockIndexElectricPotential );

        for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        {
            for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
            {
                auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(), matName );
                auto const& electricConductivity = this->materialsProperties()->electricConductivity( matName );
                auto sigmaExpr = expr( electricConductivity.expr(), se );
                mybfAV +=
                    integrate( _range=range,
                               _expr= sigmaExpr*inner(trans(gradt(electricV)),id(magneticA)),
                               _geomap=this->geomap() );
            }
        }
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("Electromagnetic","updateLinearPDE",
              "finish in "+(boost::format("%1% s") % timeElapsed).str() );
}

template< typename ElectricType, typename MagneticType>
template <typename ModelContextType>
void
Electromagnetic<ElectricType,MagneticType>::updateLinearPDEDofElimination( DataUpdateLinear & data, ModelContextType const& mctx ) const
{
    M_electricModel->updateLinearPDEDofElimination( data, mctx );
    M_magneticModel->updateLinearPDEDofElimination( data, mctx );
}

}

#endif

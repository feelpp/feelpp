/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/modelcore/modelgenericpde.hpp>


namespace Feel
{
namespace FeelModels
{

template <uint16_type Dim>
ModelGenericPDE<Dim>::ModelGenericPDE( infos_type const& infos )
    :
    super_type( "GenericPDE" ),
    ModelBase(""),
    M_infos( infos )
{
    this->setupGenericPDE();
}

template <uint16_type Dim>
ModelGenericPDE<Dim>::Infos::Infos( std::string const& name, pt::ptree const& eqPTree )
{
    M_equationName = name;

    if ( auto nameEqOpt =  eqPTree.template get_optional<std::string>( "name" ) )
        M_equationName = *nameEqOpt;

    if ( auto unknownPTree = eqPTree.get_child_optional("unknown") )
    {
        if ( auto unknownNameOpt =  unknownPTree->template get_optional<std::string>( "name" ) )
            M_unknownName = *unknownNameOpt;
        else
            CHECK( false ) << "require to define unknown.name";
        if ( auto unknownSymbolOpt =  unknownPTree->template get_optional<std::string>( "symbol" ) )
            M_unknownSymbol = *unknownSymbolOpt;
        else
            M_unknownSymbol = M_unknownName;
        if ( auto unknownBasisOpt =  unknownPTree->template get_optional<std::string>( "basis" ) )
            M_unknownBasis = *unknownBasisOpt;
        else
            M_unknownBasis = "Pch1";
    }

    std::string unknownShape;
    if ( M_unknownBasis == "Pch1" || M_unknownBasis == "Pch2" || M_unknownBasis == "Pdh1" )
        unknownShape = "scalar";
    else if ( M_unknownBasis == "Pchv1" || M_unknownBasis == "Pchv2" || M_unknownBasis == "Ned1h0" )
        unknownShape = "vectorial";
    else
        CHECK( false ) << "invalid unknown.basis : " << M_unknownBasis;
}

template <uint16_type Dim>
void
ModelGenericPDE<Dim>::setupGenericPDE()
{
    this->M_physicDefault = M_infos.equationName();

    auto mphysic = std::make_shared<ModelPhysic<Dim>>( this->physicType(),  this->physicDefault() );

    std::string unknownShape;
    if ( this->unknownBasis() == "Pch1" ||  this->unknownBasis() == "Pch2" || this->unknownBasis() == "Pdh1" )
        unknownShape = "scalar";
    else if ( this->unknownBasis() == "Pchv1" || this->unknownBasis() == "Pchv2" || this->unknownBasis() == "Ned1h0" )
        unknownShape = "vectorial";
    else
        CHECK( false ) << "invalid unknown.basis : " << this->unknownBasis();

    material_property_shape_dim_type scalarShape = std::make_pair(1,1);
    material_property_shape_dim_type vectorialShape = std::make_pair(nDim,1);
    material_property_shape_dim_type matrixShape = std::make_pair(nDim,nDim);

    mphysic->addMaterialPropertyDescription( this->convectionCoefficientName(), this->convectionCoefficientName(), { vectorialShape } );
    mphysic->addMaterialPropertyDescription( this->diffusionCoefficientName(), this->diffusionCoefficientName(), { scalarShape, matrixShape } );
    mphysic->addMaterialPropertyDescription( this->reactionCoefficientName(), this->reactionCoefficientName(), { scalarShape } );
    mphysic->addMaterialPropertyDescription( this->firstTimeDerivativeCoefficientName(), this->firstTimeDerivativeCoefficientName(), { scalarShape } );
    mphysic->addMaterialPropertyDescription( this->secondTimeDerivativeCoefficientName(), this->secondTimeDerivativeCoefficientName(), { scalarShape } );

    mphysic->addMaterialPropertyDescription( this->conservativeFluxConvectionCoefficientName(), this->conservativeFluxConvectionCoefficientName(), { vectorialShape } );

    if ( unknownShape == "scalar" )
    {
        mphysic->addMaterialPropertyDescription( this->sourceCoefficientName(), this->sourceCoefficientName(), { scalarShape } );
        mphysic->addMaterialPropertyDescription( this->conservativeFluxSourceCoefficientName(), this->conservativeFluxSourceCoefficientName(), { vectorialShape } );
    }
    else if ( unknownShape == "vectorial" )
    {
        mphysic->addMaterialPropertyDescription( this->sourceCoefficientName(), this->sourceCoefficientName(), { vectorialShape } );
        mphysic->addMaterialPropertyDescription( this->conservativeFluxSourceCoefficientName(), this->conservativeFluxSourceCoefficientName(), { matrixShape } );
        mphysic->addMaterialPropertyDescription( this->curlCurlCoefficientName(), this->curlCurlCoefficientName(), { scalarShape } );
    }

    this->M_physics.emplace( mphysic->name(), mphysic );
}

template <uint16_type Dim>
ModelGenericPDEs<Dim>::ModelGenericPDEs( /*std::string const& physic*/ )
    :
    super_type( "GenericPDEs" ),
    ModelBase("")
{}

template <uint16_type Dim>
void
ModelGenericPDEs<Dim>::setupGenericPDEs( pt::ptree const& modelPTree )
{
    if ( auto equationsOpt = modelPTree.get_child_optional("equations") )
    {
        if ( equationsOpt->empty() )
        {
            std::string equationName = equationsOpt->get_value<std::string>();
            CHECK( false ) << "TODO";
        }
        else
        {
            for ( auto const& itemEq : *equationsOpt )
            {
                CHECK( itemEq.first.empty() ) << "should be an array, not a subtree";
                std::string nameEqDefault = (boost::format("equation%1%")%M_pdes.size()).str();
                typename ModelGenericPDE<nDim>::infos_type infos( nameEqDefault, itemEq.second );
                M_pdes.push_back( std::make_tuple( std::move( infos ), std::shared_ptr<ModelGenericPDE<nDim>>{} ) );
            }
        }
    }
}

template <uint16_type Dim>
void
ModelGenericPDEs<Dim>::initGenericPDEs( std::string const& name )
{
    this->M_physicDefault = name;
    auto mphysic = std::make_shared<ModelPhysic<Dim>>( this->physicType(), name );
    this->M_physics.emplace( name, mphysic );
}

template <uint16_type Dim>
void
ModelGenericPDEs<Dim>::addGenericPDE( typename ModelGenericPDE<nDim>::infos_type const& infos )
{
    M_pdes.push_back( std::make_tuple( infos, std::shared_ptr<ModelGenericPDE<nDim>>{} ) );
}

template <uint16_type Dim>
void
ModelGenericPDEs<Dim>::updateForUseGenericPDEs()
{
    auto & mphysic = this->M_physics[this->physicDefault()];
    for ( auto const& [infos,pde] : M_pdes )
    {
        CHECK( pde ) <<"pde not defined";
        this->M_physics.insert( pde->physics().begin(), pde->physics().end() );
        for ( auto const& subPhysic : pde->physics() ) // normally only one
            mphysic->addSubphysic( subPhysic.second );
    }
}

template class ModelGenericPDE<2>;
template class ModelGenericPDE<3>;
template class ModelGenericPDEs<2>;
template class ModelGenericPDEs<3>;

} // namespace FeelModels
} // namespace Feel

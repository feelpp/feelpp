/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels/modelcore/modelgenericpde.hpp>


namespace Feel
{
namespace FeelModels
{

template <uint16_type Dim>
ModelGenericPDE<Dim>::ModelGenericPDE( /*std::string const& physic*/ )
{}

template <uint16_type Dim>
ModelGenericPDE<Dim>::ModelGenericPDE( std::string const& name, pt::ptree const& p )
{
    this->setupGenericPDE( name, p );
}

template <uint16_type Dim>
void
ModelGenericPDE<Dim>::setupGenericPDE( std::string const& name, pt::ptree const& eqPTree )
{
    if ( auto nameEqOpt =  eqPTree.template get_optional<std::string>( "name" ) )
         this->M_physic = *nameEqOpt;
    else
        this->M_physic  = name;

    this->M_physics = { this->M_physic };
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
        //CHECK( M_unknownShape == "scalar" || M_unknownShape == "vectorial" ) << "invalid unknown.shape : " << M_unknownShape;
    }

    material_property_shape_dim_type scalarShape = std::make_pair(1,1);
    material_property_shape_dim_type matrixShape = std::make_pair(nDim,nDim);

    this->addMaterialPropertyDescription( this->convectionCoefficientName(), this->convectionCoefficientName(), { scalarShape } );
    this->addMaterialPropertyDescription( this->diffusionCoefficientName(), this->diffusionCoefficientName(), { scalarShape } );
    this->addMaterialPropertyDescription( this->reactionCoefficientName(), this->reactionCoefficientName(), { scalarShape } );
    this->addMaterialPropertyDescription( this->firstTimeDerivativeCoefficientName(), this->firstTimeDerivativeCoefficientName(), { scalarShape } );
}

template <uint16_type Dim>
ModelGenericPDEs<Dim>::ModelGenericPDEs( /*std::string const& physic*/ )
{}

template <uint16_type Dim>
void
ModelGenericPDEs<Dim>::setupGenericPDEs( std::string const& name, pt::ptree const& modelPTree )
{
    this->M_physic = name;
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
                ModelGenericPDE<nDim> mgpde( nameEqDefault, itemEq.second );
                M_pdes.push_back( std::move( mgpde ) );
            }
        }
    }


    for ( auto const& pde : M_pdes )
    {
        for (auto const& [propName,propDesc] : pde.materialPropertyDescription() )
            this->M_materialPropertyDescription[propName] = propDesc;
        this->M_physics.insert( pde.physics().begin(), pde.physics().end() );
        this->M_mapPhysicsToSubphysics[this->physic()].insert( pde.physics().begin(), pde.physics().end() );
    }
    this->M_physics.insert( this->physic() );

}

template class ModelGenericPDE<2>;
template class ModelGenericPDE<3>;
template class ModelGenericPDEs<2>;
template class ModelGenericPDEs<3>;

} // namespace FeelModels
} // namespace Feel

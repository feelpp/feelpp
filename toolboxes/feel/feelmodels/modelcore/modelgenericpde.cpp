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
ModelGenericPDE<Dim>::Infos::Infos( std::string const& name, nl::json const& jarg )
{
    M_equationName = name;

    if ( jarg.contains("name") )
    {
        auto const& j_name = jarg.at("name");
        if ( j_name.is_string() )
            M_equationName = j_name.template get<std::string>();
    }

    if ( jarg.contains("unknown") )
    {
        auto const& j_unknown = jarg.at("unknown");
        CHECK( j_unknown.contains("name") ) << "require to define unknown.name";
        auto const& j_unknown_name = j_unknown.at("name");
        if ( j_unknown_name.is_string() )
            M_unknownName = j_unknown_name.template get<std::string>();
        CHECK( !M_unknownName.empty() ) << "require to define a non empty unknown.name";

        if ( j_unknown.contains( "symbol" ) )
        {
            auto const& j_unknown_symbol = j_unknown.at( "symbol" );
            if ( j_unknown_symbol.is_string() )
                M_unknownSymbol = j_unknown_symbol.template get<std::string>();
        }
        else
             M_unknownSymbol = M_unknownName;

        if ( j_unknown.contains( "basis" ) )
        {
            auto const& j_unknown_basis = j_unknown.at( "basis" );
            if ( j_unknown_basis.is_string() )
                M_unknownBasis = j_unknown_basis.template get<std::string>();
        }
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
    auto mphysic = std::make_shared<ModelPhysic<Dim>>( this->physicType(), this->physicType(), M_infos.equationName(), *this ); // TOCHECK VINCENT

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

    this->M_physics.emplace( std::make_pair(mphysic->type(),mphysic->name()), mphysic );
}

template <uint16_type Dim>
ModelGenericPDEs<Dim>::ModelGenericPDEs( /*std::string const& physic*/ )
    :
    super_type( "GenericPDEs" ),
    ModelBase("")
{}

template <uint16_type Dim>
void
ModelGenericPDEs<Dim>::setupGenericPDEs( nl::json const& jarg )
{
    if ( jarg.is_object() )
    {
        CHECK( false ) << "TODO";
    }
    else if ( jarg.is_array() )
    {
        for ( auto const& [jargkey,jargval] : jarg.items() )
        {
            std::string nameEqDefault = (boost::format("equation%1%")%M_pdes.size()).str();
            typename ModelGenericPDE<nDim>::infos_type infos( nameEqDefault, jargval );
            M_pdes.push_back( std::make_tuple( std::move( infos ), std::shared_ptr<ModelGenericPDE<nDim>>{} ) );
        }
    }
}

template <uint16_type Dim>
void
ModelGenericPDEs<Dim>::initGenericPDEs( std::string const& name )
{
    //this->M_physicDefault = name;
    auto mphysic = std::make_shared<ModelPhysic<Dim>>( this->physicType(), this->physicType(), name, *this ); // TOCHECK VINCENT
    this->M_physics.emplace( std::make_pair(mphysic->type(),mphysic->name()), mphysic );
}

template <uint16_type Dim>
void
ModelGenericPDEs<Dim>::addGenericPDE( typename ModelGenericPDE<nDim>::infos_type const& infos )
{
    M_pdes.push_back( std::make_tuple( infos, std::shared_ptr<ModelGenericPDE<nDim>>{} ) );
}

template <uint16_type Dim>
void
ModelGenericPDEs<Dim>::updateForUseGenericPDEs( std::string const& name )
{
    //TODO VINCENT add also name in arg as initGenericPDEs or save this info??

    auto pId = std::make_pair( this->physicType(), name );
    CHECK( this->hasPhysic(pId) ) << "physic not registered";
    auto mphysic = this->physic(pId);
    //zfor ( auto & [physicId,mphysic] : this->physicsFromCurrentType() )
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


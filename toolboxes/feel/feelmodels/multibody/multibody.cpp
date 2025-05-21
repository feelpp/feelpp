/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/multibody/multibody.hpp>

namespace Feel::FeelModels {


template< typename ConvexType>
Multibody<ConvexType>::Multibody( std::string const& prefix,
                                  std::string const& keyword,
                                  worldcomm_ptr_t const& worldComm,
                                  ModelBaseRepository const& modelRep )
    :
    super_numerical_type( prefix, keyword, worldComm, "", modelRep ),
    super_physics_type( "multibody" ),
    ModelBase( prefix, keyword, worldComm, "", modelRep )
{}


template< typename ConvexType>
void
Multibody<ConvexType>::init()
{
    this->initModelProperties();

    this->initPhysics( this->shared_from_this(), this->modelProperties().models() );

    // physical properties
    if ( !M_materialsProperties )
    {
        M_materialsProperties.reset( new materialsproperties_type( this->shared_from_this() ) );
        M_materialsProperties->updateForUse( this->modelProperties().materials() );
    }

#if 0
    this->initMesh();
#endif

    this->materialsProperties()->addMesh( this->mesh() );


    for ( auto & [physicId,physicObj] : this->physicsFromCurrentType() )
    {
        auto physicMultibody = std::static_pointer_cast<ModelPhysicMultibody<nDim>>(physicObj);//->updateForUse( this->materialsProperties(), this->mesh() );
        for ( auto const& [bodyName,bodyPhysic] : physicMultibody->bodies() )
        {
            auto body = std::make_unique<body_type>( &bodyPhysic );
            body->setup( M_materialsProperties, this->mesh() );
            M_bodies.emplace( bodyName, std::move(body) );
        }
    }

    this->setIsUpdatedForUse( true );

}


template< typename ConvexType>
void
Multibody<ConvexType>::applyRemesh( mesh_ptrtype oldMesh, mesh_ptrtype newMesh, std::shared_ptr<RemeshInterpolation> remeshInterp )
{
    // material prop
    this->materialsProperties()->removeMesh( oldMesh );
    this->materialsProperties()->addMesh( newMesh );
    // bodies
    for ( auto & [bodyName,body] : M_bodies )
        body->applyRemesh( oldMesh, newMesh, remeshInterp );
}


template< typename ConvexType>
void
Multibody<ConvexType>::updateInformationObject( nl::json & p ) const
{
    if ( !this->isUpdatedForUse() )
        return;
    if ( p.contains( "Environment" ) )
        return;

    super_numerical_type::super_model_base_type::updateInformationObject( p["Environment"] );

    super_numerical_type::super_model_meshes_type::updateInformationObject( p["Meshes"] );

    super_physics_type::updateInformationObjectFromCurrentType( p["Physics"] );

    // Materials properties
    if ( this->materialsProperties() )
        this->materialsProperties()->updateInformationObject( p["Materials Properties"] );
}

template< typename ConvexType>
tabulate_informations_ptr_t
Multibody<ConvexType>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains("Environment") )
        tabInfo->add( "Environment",  super_numerical_type::super_model_base_type::tabulateInformations( jsonInfo.at("Environment"), tabInfoProp ) );

    if ( jsonInfo.contains("Physics") )
        tabInfo->add( "Physics", super_numerical_type::tabulateInformations( jsonInfo.at("Physics"), tabInfoProp ) );

    if ( this->materialsProperties() && jsonInfo.contains("Materials Properties") )
        tabInfo->add( "Materials Properties", this->materialsProperties()->tabulateInformations(jsonInfo.at("Materials Properties"), tabInfoProp ) );

    if ( jsonInfo.contains("Meshes") )
        tabInfo->add( "Meshes", super_numerical_type::super_model_meshes_type::tabulateInformations( jsonInfo.at("Meshes"), tabInfoProp ) );
    return tabInfo;
}



template class Multibody<Simplex<2,1>>;
template class Multibody<Simplex<2,2>>;
template class Multibody<Simplex<3,1>>;
template class Multibody<Simplex<3,2>>;

}


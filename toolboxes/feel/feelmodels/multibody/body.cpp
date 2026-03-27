/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/multibody/body.hpp>

namespace Feel::FeelModels {


template<typename ConvexType>
void
Body<ConvexType>::setup( materialsproperties_ptrtype materialsProperties, mesh_ptrtype mesh )
{
    M_mesh = mesh;
    M_materialsProperties = materialsProperties;


    std::set<std::string> markers;
    for ( std::string const& matName : M_physicBody->materialNames() )
    {
        auto const& matProps = M_materialsProperties->materialProperties( matName );
        for ( std::string const& m : matProps.markers() )
            markers.insert(m);
    }

    M_rangeMeshElements = markedelements(this->mesh(), markers );
    M_spaceDisplacement = space_displacement_type::New(_mesh=M_mesh,_range=M_rangeMeshElements);
    M_fieldDisplacement = M_spaceDisplacement->elementPtr();
    M_fieldDisplacementAtPreviousTime = M_spaceDisplacement->elementPtr();

    this->updateForUse();
}

#if 0
template<typename ConvexType>
void
Body<ConvexType>::setup( nl::json const& jarg, ModelMaterials const& mats, mesh_ptrtype mesh )
{
    M_mesh = mesh;
    std::set<std::string> matNames;
    if ( jarg.contains( "names" ) )
    {
        auto const& j_names = jarg.at( "names" );
        if ( j_names.is_string() )
            matNames.insert( j_names.template get<std::string>() );
        else if ( j_names.is_array() )
        {
            for ( auto const& [j_nameskey,j_namesval] : j_names.items() )
            {
                CHECK( j_namesval.is_string() ) << "should be a string";
                matNames.insert( j_namesval.template get<std::string>() );
            }
        }
    }

    ModelMarkers onlyMarkers;
    if ( jarg.contains( "markers" ) )
        onlyMarkers.setup( jarg.at( "markers" ) /*, indexes*/ );

    M_materialsProperties.reset( new materialsproperties_type( M_modelPhysics ) );
    M_materialsProperties->updateForUse( mats, matNames, onlyMarkers );
    M_materialsProperties->addMesh( M_mesh );

    // init displacement space
    auto mom = this->materialsProperties()->materialsOnMesh( this->mesh() );
    auto M_rangeMeshElements = markedelements(this->mesh(), mom->markers( M_modelPhysics->physicsAvailableFromCurrentType() ) );
    M_spaceDisplacement = space_displacement_type::New(_mesh=M_mesh,_range=M_rangeMeshElements);
    M_fieldDisplacement = M_spaceDisplacement->elementPtr();
    M_fieldDisplacementAtPreviousTime = M_spaceDisplacement->elementPtr();

    this->updateForUse();
}
#endif

template<typename ConvexType>
void
Body<ConvexType>::applyRemesh( mesh_ptrtype oldMesh, mesh_ptrtype newMesh, std::shared_ptr<RemeshInterpolation> remeshInterp )
{
    if ( M_mesh )
    {
        mesh_ptrtype oldMesh = this->mesh();

        // material prop
        this->materialsProperties()->removeMesh( oldMesh );
        this->materialsProperties()->addMesh( newMesh );

        M_mesh = newMesh;

        // function space and fields
        space_displacement_ptrtype old_spaceDisplacement = M_spaceDisplacement;
        element_displacement_ptrtype old_fieldDisplacement = M_fieldDisplacement;
        element_displacement_ptrtype old_fieldDisplacementAtPreviousTime = M_fieldDisplacementAtPreviousTime;
        element_displacement_ptrtype old_fieldElasticDisplacement = M_fieldElasticDisplacement;


        std::set<std::string> markers;
        for ( std::string const& matName : M_physicBody->materialNames() )
        {
            auto const& matProps = M_materialsProperties->materialProperties( matName );
            for ( std::string const& m : matProps.markers() )
                markers.insert(m);
        }
        M_rangeMeshElements = markedelements(this->mesh(), markers );
        M_spaceDisplacement = space_displacement_type::New(_mesh=M_mesh,_range=M_rangeMeshElements);
        M_fieldDisplacement = M_spaceDisplacement->elementPtr();
        M_fieldDisplacementAtPreviousTime = M_spaceDisplacement->elementPtr();

        // createInterpolationOp
        auto opI_displacement = opInterpolation(_domainSpace=old_spaceDisplacement,
                                                _imageSpace=M_spaceDisplacement,
                                                _range=M_rangeMeshElements
                                                );

        auto matrixInterpolation_displacement = opI_displacement->matPtr();
        matrixInterpolation_displacement->multVector( *old_fieldDisplacement, *M_fieldDisplacement );
        matrixInterpolation_displacement->multVector( *old_fieldDisplacementAtPreviousTime, *M_fieldDisplacementAtPreviousTime );

        if ( old_fieldElasticDisplacement )
        {
            M_fieldElasticDisplacement = M_spaceDisplacement->elementPtr();
            matrixInterpolation_displacement->multVector( *old_fieldElasticDisplacement, *M_fieldElasticDisplacement );
        }

        if ( M_fieldElasticVelocity )
        {
            space_velocity_ptrtype old_spaceElasticVelocity = M_spaceElasticVelocity;
            element_velocity_ptrtype old_fieldElasticVelocity = M_fieldElasticVelocity;
            M_fieldElasticVelocity.reset();
            this->initElasticVelocity();

            auto opI_elasticVelocity = opInterpolation(_domainSpace=old_spaceElasticVelocity,
                                                       _imageSpace=M_spaceElasticVelocity,
                                                       _range=M_rangeMeshElements );

            auto matrixInterpolation_elasticVelocity = opI_elasticVelocity->matPtr();
            matrixInterpolation_elasticVelocity->multVector( *old_fieldElasticVelocity, *M_fieldElasticVelocity );

        }

    }
}

template<typename ConvexType>
void
Body<ConvexType>::updateForUse()
{
    CHECK( M_materialsProperties ) << "no materialsProperties defined";

    auto mom = M_materialsProperties->materialsOnMesh(M_mesh);
    M_mass = 0;
    M_massCenter = eigen_vector_type<nRealDim>::Zero();
    double massForMassCenter = 0;
    for ( std::string const& matName : M_physicBody->materialNames() )
    {
        auto const& range = mom->rangeMeshElementsByMaterial( matName );
        auto const& density = M_materialsProperties->density( matName );
        auto const& densityExpr = density.exprScalar();
        double currentMass = integrate(_range=range,_expr=densityExpr).evaluate()(0,0);
        M_mass += currentMass;
        if ( M_physicBody->useMaterialWithMassCenterEvaluation( matName ) )
        {
            M_massCenter += integrate(_range=range,_expr=densityExpr*P()).evaluate();
            massForMassCenter += currentMass;
        }
    }
    M_massCenter /= massForMassCenter;

    if ( M_physicBody->hasMassCenterImposed() )
    {
        M_massCenter = M_physicBody->massCenterImposedExpr( /*se*/ ).evaluate();
    }

    this->computeMomentOfInertia_bodyFrame( this->massCenterExpr(), this->rigidRotationMatrix(), M_momentOfInertia_bodyFrame );
}

template<typename ConvexType>
void
Body<ConvexType>::updateDisplacementFromRigidDisplacement( eigen_vector_type<nRealDim> const& rigidTranslation, rotation_angles_type const& rigidRotationAngles )
{
    M_rigidTranslationDisplacement = rigidTranslation;
    M_rigidRotationAngles = rigidRotationAngles;
    // we compute new displacement from current mesh position + rigid body displacement update
    // TODO: check is on moving mesh
    //auto T = Feel::vf::toExpr( M_rigidTranslationDisplacement - M_rigidTranslationDisplacementAtPreviousTime ); // NOT COMPILE, should be fixed!!
    auto T = this->rigidTranslationExpr() - Feel::vf::toExpr( M_rigidTranslationDisplacementAtPreviousTime );
    auto R = Feel::vf::toExpr( Body::rigidRotationMatrix( M_rigidRotationAngles-M_rigidRotationAnglesAtPreviousTime ) );
    auto M = this->massCenterExpr();

    auto tmp = M_spaceDisplacement->element();
    tmp = this->fieldDisplacement();

    this->updateDisplacement( elements(support(M_spaceDisplacement)), idv(tmp) + R*( P() - M ) + M + T - P() );
}



template class Body<Simplex<2,1>>;
template class Body<Simplex<2,2>>;
template class Body<Simplex<3,1>>;
template class Body<Simplex<3,2>>;


}

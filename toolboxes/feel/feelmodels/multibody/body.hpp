/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_MULTIBODY_BODY_HPP
#define FEELPP_TOOLBOXES_MULTIBODY_BODY_HPP 1

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelcore/modelphysics.hpp>
#include <feel/feelmodels/modelmaterials/materialsproperties.hpp>
#include <feel/feelmodels/modelcore/remeshinterpolation.hpp>



namespace Feel::vf
{

template <int Dim>
auto toExpr( eigen_vector_type<Dim> const& ev )
{
    static_assert( Dim > 0 && Dim <=3, "toExpr only implement with Dim 1,2,3" );
    if constexpr ( Dim == 1 )
        return cst(ev(0));
    else if constexpr ( Dim == 2 )
        return vec( cst(ev(0)), cst(ev(1)) );
    else
        return vec( cst(ev(0)), cst(ev(1)), cst(ev(2)) );
}

template <int RowDim,int RowCol>
auto toExpr( eigen_matrix_type<RowDim, RowCol> const& em )
{
    static_assert( RowDim == RowCol && (RowDim == 2 || RowDim == 3), "toExpr only implement matrix 2x2 or 3x3" );
    if constexpr ( RowDim == 2 && RowCol == 2 )
    {
        return mat<2,2>( cst( em(0,0) ), cst( em(0,1) ),
                         cst( em(1,0) ), cst( em(1,1) ) );
    }
    else
    {
        return mat<3,3>( cst( em(0,0) ), cst( em(0,1) ), cst( em(0,2) ),
                         cst( em(1,0) ), cst( em(1,1) ), cst( em(1,2) ),
                         cst( em(2,0) ), cst( em(2,1) ), cst( em(2,2) ) );
    }
}

} // namespace Feel::vf


namespace Feel::FeelModels
{

// // forward declaration
// template< typename ConvexType>
// class Multibody;

/**
 * @brief Body class
 * @ingroup Multibody
 *
 */
template< typename ConvexType>
class Body
{
public :
    using self_type = Body<ConvexType>;

    // mesh
    using convex_type = ConvexType;
    static inline const uint16_type nDim = convex_type::nDim;
    static inline const uint16_type nOrderGeo = convex_type::nOrder;
    static inline const uint16_type nRealDim = convex_type::nRealDim;
    using mesh_type = Mesh<convex_type>;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;

    using mesh_range_element_type = Range<mesh_type,MESH_ELEMENTS>;

    // mesh motion
    using mesh_motion_type = MeshALE<convex_type>;
    using mesh_motion_ptrtype = std::shared_ptr<mesh_motion_type>;

    // materials properties
    using materialsproperties_type = MaterialsProperties<nRealDim>;
    using materialsproperties_ptrtype = std::shared_ptr<materialsproperties_type>;


    static constexpr int nDimRotation = (nDim==3)?3:1;
    using moment_of_inertia_type = eigen_matrix_type<nDimRotation,nDimRotation>;
    using translational_velocity_type = eigen_vector_type<nRealDim>;
    using rotation_angles_type = eigen_matrix_type<nDimRotation, 1>;
    using angular_velocity_type = rotation_angles_type;

    using space_displacement_type = typename mesh_motion_type::ale_map_functionspace_type;
    using space_displacement_ptrtype = std::shared_ptr<space_displacement_type>;
    using element_displacement_type = typename space_displacement_type::element_type;
    using element_displacement_ptrtype = std::shared_ptr<element_displacement_type>;

    using space_velocity_type = space_displacement_type;
    using space_velocity_ptrtype = space_displacement_ptrtype;
    using element_velocity_type = element_displacement_type;
    using element_velocity_ptrtype = element_displacement_ptrtype;

    using physic_body_type = typename ModelPhysicMultibody<nRealDim>::Body;

    Body( physic_body_type const* pb ) : M_physicBody( pb ) {}
    Body( Body const& ) = default;
    Body( Body && ) = default;

    void setup( materialsproperties_ptrtype materialsProperties, mesh_ptrtype mesh, mesh_motion_ptrtype meshMotion );
    void applyRemesh( mesh_ptrtype oldMesh, mesh_ptrtype newMesh, std::shared_ptr<RemeshInterpolation> remeshInterp = std::make_shared<RemeshInterpolation>() );
    void updateForUse();

    physic_body_type const* physic() const { return M_physicBody; }


    //! return the mesh containing the body mesh
    mesh_ptrtype mesh() const { return M_mesh; }

    //! return true if a MaterialsProperties has been attached to this body
    bool hasMaterialsProperties() const { return (M_materialsProperties? true : false); }
    //! return the MaterialsProperties object associated to this body
    materialsproperties_ptrtype materialsProperties() const { return M_materialsProperties; }

    //! return the current displacement (corresponding to displacement applied to the reference mesh  and including elastic displacement if enabled)
    element_displacement_type const& fieldDisplacement() const { return *M_fieldDisplacement; }
    //! return the current displacement (corresponding to displacement applied to the reference mesh and including elastic displacement if enabled)
    element_displacement_type & fieldDisplacement() { return *M_fieldDisplacement; }
    //! return the elastic displacement (corresponding to displacement applied to the current moving mesh)
    element_displacement_type const& fieldElasticDisplacement() const { return *M_fieldElasticDisplacement; }
    //! return the elastic displacement (corresponding to displacement applied to the current moving mesh)
    element_displacement_type & fieldElasticDisplacement() { return *M_fieldElasticDisplacement; }
    //! return the elastic velocity
    element_velocity_type const& fieldElasticVelocity() const { return *M_fieldElasticVelocity; }
    //! return the elastic velocity
    element_velocity_type & fieldElasticVelocity() { return *M_fieldElasticVelocity; }

    //! return the displacement field at previous time
    element_displacement_type const& fieldDisplacementAtPreviousTime() const { return *M_fieldDisplacementAtPreviousTime; }



    //! return true if an elastic displacement is defined
    bool hasElasticDisplacement() const { return M_fieldElasticDisplacement? true : false; }

    //! return true if an elastic velocity is defined
    bool hasElasticVelocity() const { return M_fieldElasticVelocity? true:false; }


    //! update the elastic displacement from an expression \e on entities \range
    template <typename RangeType, typename ExprT>
    void updateDisplacement( RangeType const& range, Expr<ExprT> const& e )
        {
            bool close = true;
            if ( !M_fieldDisplacement )
                M_fieldDisplacement = M_spaceDisplacement->elementPtr();
            M_fieldDisplacement->on(_range=range,_expr=e,_close=close);
        }

    //! return the current translation
    eigen_vector_type<nRealDim> const& rigidTranslation() const { return M_rigidTranslationDisplacement; }

    //! return the current translation as an expression
    auto rigidTranslationExpr() const { return Feel::vf::toExpr( M_rigidTranslationDisplacement ); }

    //! return rigid translation displacement at previous time
    eigen_vector_type<nRealDim> const& rigidTranslationDisplacementAtPreviousTime() const { return M_rigidTranslationDisplacementAtPreviousTime; }

    //! return the current rotation angles
    rotation_angles_type const& rigidRotationAngles() const { return M_rigidRotationAngles; }

    //! return rotation matrix from angles
    static eigen_matrix_type<nRealDim, nRealDim> rigidRotationMatrix( rotation_angles_type const& rigidRotationAngles )
        {
            eigen_matrix_type<nRealDim, nRealDim> res;
            if constexpr ( nRealDim == 2 )
            {
                double angle = rigidRotationAngles(0,0);
                res <<   std::cos(angle), -std::sin(angle),
                    /**/ std::sin(angle),  std::cos(angle);
            }
            else
            {
                double angleZ = rigidRotationAngles(2);
                double angleY = rigidRotationAngles(1);
                double angleX = rigidRotationAngles(0);
                eigen_matrix_type<3, 3> rotMatZ,rotMatY,rotMatX;
                rotMatZ << std::cos(angleZ), -std::sin(angleZ), 0,
                    /**/   std::sin(angleZ),  std::cos(angleZ), 0,
                    /**/                  0,                 0, 1;
                rotMatY << std::cos(angleY), 0, std::sin(angleY),
                    /**/                  0, 1,                0,
                    /**/  -std::sin(angleY), 0, std::cos(angleY);
                rotMatX << 1,                0,                 0,
                    /**/   0, std::cos(angleX), -std::sin(angleX),
                    /**/   0, std::sin(angleX),  std::cos(angleX);
                res = rotMatZ*rotMatY*rotMatX;
            }
            return res;
        }

    //! return the current rotation matrix
    eigen_matrix_type<nRealDim, nRealDim> rigidRotationMatrix() const { return Body::rigidRotationMatrix( M_rigidRotationAngles ); }

    //! return the current rotation matrix as an expression
    auto rigidRotationMatrixExpr() const { return toExpr( this->rigidRotationMatrix() ); }

    //! return rotation angles at previous time
    rotation_angles_type const& rigidRotationAnglesAtPreviousTime() const { return M_rigidRotationAnglesAtPreviousTime; }


    void updateDisplacementFromRigidDisplacement( eigen_vector_type<nRealDim> const& rigidTranslation, rotation_angles_type const& rigidRotationAngles );


    //! update displacement by setting disp equal to previous disp + current elastic update
    void updateDisplacementFromElasticBehavior() // DEPRECATED!!!!
        {
            M_fieldDisplacement->zero();
            M_fieldDisplacement->add( 1.0, *M_fieldDisplacementAtPreviousTime );
            M_fieldDisplacement->add( 1.0, *M_fieldElasticDisplacement );
        }

    //! update the current displacement from the rigid body displacement and eventually from the elastic behavior if enabled
    void updateDisplacementForUse();

    // void addRigidTranslationToCurrentDisplacement( eigen_vector_type<nRealDim> const& rigidTranslation )
    //     {
    //         auto tmp = M_spaceDisplacement->element();
    //         tmp = this->fieldDisplacement();
    //         this->updateDisplacement( elements(support(M_spaceDisplacement)), idv( tmp ) + Feel::vf::toExpr(rigidTranslation) );
    //     }

    template <typename ExpRotationMatrixType,typename ExprMassCenterType>
    void applyRotationToCurrentDisplacement( Expr<ExpRotationMatrixType> const& R, Expr<ExprMassCenterType> const& massCenter )
        {
            auto tmp = M_spaceDisplacement->element();
            tmp = this->fieldDisplacement();
            if ( M_meshMotionTool && !M_meshMotionTool->isMappedOntoTheInitialMesh() )
            {
                auto relativeDisp = idv(tmp) - ( P() - idv(M_meshMotionTool->fieldInitialIdentity()) );
                this->updateDisplacement( elements(support(M_spaceDisplacement)), R*(P()+relativeDisp-massCenter) + massCenter - idv(M_meshMotionTool->fieldInitialIdentity()) );
            }
            else
                this->updateDisplacement( elements(support(M_spaceDisplacement)), R*(P()+idv(tmp)-massCenter) + massCenter - P() );
        }


    //! init init elastic displacement field if not built
    void initElasticDisplacement()
        {
            if ( !M_fieldElasticDisplacement )
                M_fieldElasticDisplacement = M_spaceDisplacement->elementPtr();
        }

    //! init init elastic displacement field if not built
    void initElasticVelocity()
        {
            if ( !M_fieldElasticVelocity )
            {
                //auto mom = this->materialsProperties()->materialsOnMesh( this->mesh() );
                //auto M_rangeMeshElements = markedelements(this->mesh(), mom->markers( M_modelPhysics->physicsAvailableFromCurrentType() ) );
                M_spaceElasticVelocity = space_velocity_type::New(_mesh=M_mesh,_range=M_rangeMeshElements);
                M_fieldElasticVelocity = M_spaceElasticVelocity->elementPtr();
            }
        }

    //! update the elastic displacement from an expression \e on entities \range
    template <typename RangeType, typename ExprT>
    void updateElasticDisplacement( RangeType const& range, Expr<ExprT> const& e )
        {
            bool close = true;
            if ( !M_fieldElasticDisplacement )
                this->initElasticDisplacement();
            M_fieldElasticDisplacement->on(_range=range,_expr=e,_close=close);
        }

    //! update the elastic displacement from an expression \e on entities \range
    template <typename RangeType, typename ExprT>
    void updateElasticVelocity( RangeType const& range, Expr<ExprT> const& e )
        {
            bool close = true;
            if ( !M_fieldElasticVelocity )
                this->initElasticVelocity();
            M_fieldElasticVelocity->on(_range=range,_expr=e,_close=close);
        }



    void setMass( double m ) { M_mass = m; }
    void setMomentOfInertia_bodyFrame( moment_of_inertia_type const& m ) { M_momentOfInertia_bodyFrame = m; }
    void setMomentOfInertia_bodyFrame( double val ) { M_momentOfInertia_bodyFrame = val*moment_of_inertia_type::Identity(); }
    void setMassCenter( eigen_vector_type<nRealDim> const& massCenter ) { M_massCenter = massCenter; }

    //! return the mass of the body
    double mass() const { return M_mass; }
    //! return the mass of the body as an expression
    auto massExpr() const { return cst( M_mass ); }

    //! return the moment of inertia related to body frame
    moment_of_inertia_type const& momentOfInertia_bodyFrame() const { return M_momentOfInertia_bodyFrame; }
    //! return the moment of inertia related to body frame as an expression
    auto momentOfInertiaExpr_bodyFrame() const { return Feel::vf::toExpr(M_momentOfInertia_bodyFrame); }
    //! return the moment of inertia related to inertial frame
    moment_of_inertia_type momentOfInertia_inertialFrame() const
        {
            if constexpr ( nDim == 2 )
                return M_momentOfInertia_bodyFrame;
            else
            {
                auto R = this->rigidRotationMatrix();
                return R*M_momentOfInertia_bodyFrame*(R.transpose());
            }
        }
    //! return time derivative of moment of inertia related to body frame
    moment_of_inertia_type timeDerivativeOfMomentOfInertia_bodyFrame( double dt ) const
        {
            return (1/dt)*(M_momentOfInertia_bodyFrame - M_momentOfInertiaAtPreviousTime_bodyFrame);
        }

    //! return the center of mass of the body
    eigen_vector_type<nRealDim> const& massCenter() const { return M_massCenter; }
    //! return the center of mass of the body as an expression
    auto massCenterExpr() const { return Feel::vf::toExpr(M_massCenter); }

    //! return the mass from an expression of the density \densityExpr
    template <typename ExprType>
    double evaluateMassFromDensity( Expr<ExprType> const& densityExpr ) const;

    //! compute mass center of the body with a displacement (given as displacement of initial domain)
    //! Computation are done directly onto the current mesh (can be moving or reference or initial state)
    template <typename DispElementType>
    std::tuple<double,eigen_vector_type<nRealDim> >
    computeMassAndMassCenterFromDisplacementField( DispElementType const& d ) const;

    //!
    template <typename MassCenterExprType>
    void computeMomentOfInertia_inertialFrame( MassCenterExprType const& massCenterExpr, moment_of_inertia_type & momentOfInertia, bool addValue = false ) const;

    template <typename MassCenterExprType>
    void computeMomentOfInertia_bodyFrame( MassCenterExprType const& massCenterExpr, eigen_matrix_type<nRealDim, nRealDim> const& R, moment_of_inertia_type & momentOfInertia, bool addValue = false ) const;


    void setParameterValues( std::map<std::string,double> const& mp )
        {
            if ( M_materialsProperties )
                M_materialsProperties->setParameterValues( mp );
        }

    auto modelMeasuresQuantities( std::string const& prefix ) const
        {
            return Feel::FeelModels::modelMeasuresQuantities( modelMeasuresQuantity( prefix, "mass_center", std::bind( &self_type::massCenter, this ) ),
                                                              modelMeasuresQuantity( prefix, "rigid_rotation_angles", std::bind( &self_type::rigidRotationAngles, this ) ),
                                                              modelMeasuresQuantity( prefix, "moment_of_inertia", std::bind( &self_type::momentOfInertia_inertialFrame, this ) ),
                                                              modelMeasuresQuantity( prefix, "moment_of_inertia_body_frame", std::bind( &self_type::momentOfInertia_bodyFrame, this ) )
                                                              );
        }

    void updateTimeStep()
        {
            M_rigidTranslationDisplacementAtPreviousTime = M_rigidTranslationDisplacement;
            M_rigidRotationAnglesAtPreviousTime = M_rigidRotationAngles;
            M_momentOfInertiaAtPreviousTime_bodyFrame = M_momentOfInertia_bodyFrame;
            *M_fieldDisplacementAtPreviousTime = *M_fieldDisplacement;
        }

private:
    template <typename DispElementType>
    std::tuple<double,eigen_vector_type<nRealDim> >
    computeMassAndMassCenterFromDisplacementFieldImpl( DispElementType const& d ) const;

private:
    physic_body_type const* M_physicBody = nullptr;
    mesh_ptrtype M_mesh;
    materialsproperties_ptrtype M_materialsProperties;
    mesh_range_element_type M_rangeMeshElements;
    mesh_motion_ptrtype M_meshMotionTool;

    eigen_vector_type<nRealDim> M_massCenter;//, M_massCenterRef;
    double M_mass = 0;
    moment_of_inertia_type M_momentOfInertia_bodyFrame = moment_of_inertia_type::Zero();
    moment_of_inertia_type M_momentOfInertiaAtPreviousTime_bodyFrame = moment_of_inertia_type::Zero();

    eigen_vector_type<nRealDim> M_rigidTranslationDisplacement = eigen_vector_type<nRealDim>::Zero();
    eigen_vector_type<nRealDim> M_rigidTranslationDisplacementAtPreviousTime = eigen_vector_type<nRealDim>::Zero();
    rotation_angles_type M_rigidRotationAngles = rotation_angles_type::Zero();
    rotation_angles_type M_rigidRotationAnglesAtPreviousTime = rotation_angles_type::Zero();

    space_displacement_ptrtype M_spaceDisplacement;
    element_displacement_ptrtype M_fieldDisplacement;
    element_displacement_ptrtype M_fieldElasticDisplacement;
    element_displacement_ptrtype M_fieldDisplacementAtPreviousTime;

    space_velocity_ptrtype M_spaceElasticVelocity;
    element_velocity_ptrtype M_fieldElasticVelocity;
};



template<typename ConvexType>
template <typename ExprType>
double
Body<ConvexType>::evaluateMassFromDensity( Expr<ExprType> const& densityExpr ) const
{
    CHECK( M_materialsProperties ) << "no materialsProperties defined";
    auto mom = M_materialsProperties->materialsOnMesh(M_mesh);
    double mass = 0;
    for ( std::string const& matName : M_physicBody->materialNames() )
    {
        auto const& range = mom->rangeMeshElementsByMaterial( matName );
        mass += integrate(_range=range,_expr=densityExpr).evaluate()(0,0);
    }
    return mass;
}


template<typename ConvexType>
template <typename DispElementType>
std::tuple<double,eigen_vector_type<Body<ConvexType>::nRealDim> >
Body<ConvexType>::computeMassAndMassCenterFromDisplacementField( DispElementType const& d ) const
{
    CHECK( M_materialsProperties ) << "no materialsProperties defined";

    if ( M_meshMotionTool && !M_meshMotionTool->isMappedOntoTheInitialMesh() )
    {
        auto u = M_spaceDisplacement->element();
        auto dispApplyOnMesh = idv(d) - ( P() - idv(M_meshMotionTool->fieldInitialIdentity()) );
        u.on(_range=M_rangeMeshElements, _expr=dispApplyOnMesh);
        return computeMassAndMassCenterFromDisplacementFieldImpl( u );
    }
    else
        return computeMassAndMassCenterFromDisplacementFieldImpl( d );
}

template<typename ConvexType>
template <typename DispElementType>
std::tuple<double,eigen_vector_type<Body<ConvexType>::nRealDim> >
Body<ConvexType>::computeMassAndMassCenterFromDisplacementFieldImpl( DispElementType const& d ) const
{
    auto const Id = eye<nDim,nDim>();
    // deformation tensor
    auto F = Id+gradv(d);
    auto J = det(F);

    auto mom = M_materialsProperties->materialsOnMesh(M_mesh);
    double mass = 0;
    eigen_vector_type<nRealDim> massCenter = eigen_vector_type<nRealDim>::Zero();
    double massForMassCenter = 0;
    for ( std::string const& matName : M_physicBody->materialNames() )
    {
        auto const& range = mom->rangeMeshElementsByMaterial(matName);
        auto const& density = M_materialsProperties->density( matName );
        auto const& densityExpr = density.exprScalar();
        double currentMass = integrate(_range=range,_expr=densityExpr*J).evaluate()(0,0);
        mass += currentMass;
        if ( M_physicBody->useMaterialWithMassCenterEvaluation( matName ) )
        {
            massCenter += integrate(_range=range,_expr=densityExpr*(P()+idv(d))*J).evaluate();
            massForMassCenter += currentMass;
        }
    }
    massCenter /= massForMassCenter;
    return std::make_tuple( mass, std::move( massCenter ) );
}

template<typename ConvexType>
template <typename MassCenterExprType>
void
Body<ConvexType>::computeMomentOfInertia_inertialFrame( MassCenterExprType const& massCenterExpr, moment_of_inertia_type & momentOfInertia, bool addValue ) const
{
    auto mom = M_materialsProperties->materialsOnMesh(M_mesh);
    if ( !addValue )
        momentOfInertia = moment_of_inertia_type::Zero();
    for ( std::string const& matName : M_physicBody->materialNames() )
    {
        auto const& range = mom->rangeMeshElementsByMaterial( matName );
        auto const& density = M_materialsProperties->density( matName );
        auto const& densityExpr = density.exprScalar();

        if constexpr ( nDim == 2 )
        {
            momentOfInertia(0,0) += integrate(_range=range,_expr=densityExpr*( inner(P()-massCenterExpr) ) ).evaluate()(0,0);
        }
        else
        {
            auto rvec = P()-massCenterExpr;
            momentOfInertia += integrate(_range=range,_expr=densityExpr*( inner(rvec)*eye<nDim,nDim>() - rvec*trans(rvec) ) ).evaluate();
        }
    }
}

template<typename ConvexType>
template <typename MassCenterExprType>
void
Body<ConvexType>::computeMomentOfInertia_bodyFrame( MassCenterExprType const& massCenterExpr, eigen_matrix_type<nRealDim, nRealDim> const& R, moment_of_inertia_type & momentOfInertia, bool addValue ) const
{
    if constexpr ( nDim == 2 )
    {
        this->computeMomentOfInertia_inertialFrame( massCenterExpr, momentOfInertia, addValue );
    }
    else
    {
        moment_of_inertia_type momentOfInertia_inertialFrame;
        this->computeMomentOfInertia_inertialFrame( massCenterExpr,momentOfInertia_inertialFrame );
        if ( addValue )
            momentOfInertia += R.transpose()*momentOfInertia_inertialFrame*R;
        else
            momentOfInertia = R.transpose()*momentOfInertia_inertialFrame*R;
    }
}


} // namespace Feel::FeelModels

#endif

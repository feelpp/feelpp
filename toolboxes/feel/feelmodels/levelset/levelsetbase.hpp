/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Thibaut Metivet <thibaut.metivet@univ-grenoble-alpes.fr>
 Date: 2019-01-25

 Copyright (C) 2016 Université Joseph Fourier (Grenoble I)

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 3.0 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
/**
 \file levelsetbase.hpp
 \author Thibaut Metivet <thibaut.metivet@univ-grenoble-alpes.fr>
 \date 2019-01-25
 */
#ifndef _LEVELSETBASE_HPP
#define _LEVELSETBASE_HPP 1

#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/projector.hpp>

#include <feel/feelmodels/levelset/levelsetspacemanager.hpp>
#include <feel/feelmodels/levelset/levelsettoolmanager.hpp>
#include <feel/feelmodels/levelset/levelsetredistanciation_fm.hpp>
#include <feel/feelmodels/levelset/levelsetredistanciation_hj.hpp>

#include <feel/feelmodels/modelcore/modelbase.hpp>
#include <feel/feelmodels/modelcore/modelfields.hpp>
#include <feel/feelmodels/modelcore/modelmeasuresquantities.hpp>

#include <feel/feelmodels/levelset/parameter_map.hpp>

#if defined (MESH_ADAPTATION_LS)
 #include <levelsetmesh/meshadaptation.hpp>
// #warning MESH_ADAPTATION_LS is defined in levelset. Need to be defined identically in the application
#endif


namespace Feel {

namespace FeelModels {

/* Levelset redistanciation strategy
 * FM -> Fast-Marching
 * HJ -> Hamilton-Jacobi
 */
enum class LevelSetDistanceMethod { NONE, FASTMARCHING, HAMILTONJACOBI, RENORMALISATION };

enum class LevelSetMeasuresExported
{
    Volume, Perimeter, Position_COM, Velocity_COM
};

template<
    typename ConvexType, typename BasisType, typename PeriodicityType = NoPeriodicity, 
    typename BasisPnType = BasisType
    >
class LevelSetBase : 
    public ModelNumerical,
    public std::enable_shared_from_this< LevelSetBase<ConvexType, BasisType, PeriodicityType, BasisPnType> >
{
    typedef ModelNumerical super_type;
public:
    typedef LevelSetBase<ConvexType, BasisType, PeriodicityType, BasisPnType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    static const uint16_type Order = BasisType::nOrder;
    using size_type = typename super_type::size_type;
    using value_type = typename super_type::value_type;

    //--------------------------------------------------------------------//
    // Mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    static const uint16_type nRealDim = convex_type::nRealDim;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    typedef mesh_type mymesh_type;

    //--------------------------------------------------------------------//
    // Periodicity
    typedef PeriodicityType periodicity_type;
    //--------------------------------------------------------------------//
    // Function space manager
    typedef LevelSetSpaceManager<ConvexType, BasisType, PeriodicityType, BasisPnType> levelset_space_manager_type;
    typedef std::shared_ptr<levelset_space_manager_type> levelset_space_manager_ptrtype;
    //--------------------------------------------------------------------//
    // Function spaces and elements
    // levelset
    typedef typename levelset_space_manager_type::basis_scalar_type basis_levelset_type;
    typedef typename levelset_space_manager_type::space_scalar_type space_levelset_type;
    typedef typename levelset_space_manager_type::space_scalar_ptrtype space_levelset_ptrtype;
    typedef typename space_levelset_type::element_type element_levelset_type;
    typedef std::shared_ptr<element_levelset_type> element_levelset_ptrtype;
    // levelset PN
    typedef typename levelset_space_manager_type::basis_scalar_PN_type basis_levelset_PN_type;
    typedef typename levelset_space_manager_type::space_scalar_PN_type space_levelset_PN_type;
    typedef typename levelset_space_manager_type::space_scalar_PN_ptrtype space_levelset_PN_ptrtype;
    typedef typename space_levelset_PN_type::element_type element_levelset_PN_type;
    typedef std::shared_ptr<element_levelset_PN_type> element_levelset_PN_ptrtype;
    // scalar
    typedef typename levelset_space_manager_type::basis_scalar_type basis_scalar_type;
    typedef typename levelset_space_manager_type::space_scalar_type space_scalar_type;
    typedef typename levelset_space_manager_type::space_scalar_ptrtype space_scalar_ptrtype;
    typedef typename space_scalar_type::element_type element_scalar_type;
    typedef std::shared_ptr<element_scalar_type> element_scalar_ptrtype;
    // scalar PN
    typedef typename levelset_space_manager_type::basis_scalar_PN_type basis_scalar_PN_type;
    typedef typename levelset_space_manager_type::space_scalar_PN_type space_scalar_PN_type;
    typedef typename levelset_space_manager_type::space_scalar_PN_ptrtype space_scalar_PN_ptrtype;
    typedef typename space_scalar_PN_type::element_type element_scalar_PN_type;
    typedef std::shared_ptr<element_scalar_PN_type> element_scalar_PN_ptrtype;
    // vectorial
    typedef typename levelset_space_manager_type::basis_vectorial_type basis_vectorial_type;
    typedef typename levelset_space_manager_type::space_vectorial_type space_vectorial_type;
    typedef typename levelset_space_manager_type::space_vectorial_ptrtype space_vectorial_ptrtype;
    typedef typename space_vectorial_type::element_type element_vectorial_type;
    typedef std::shared_ptr< element_vectorial_type > element_vectorial_ptrtype;
    // markers P0
    typedef typename levelset_space_manager_type::basis_markers_type basis_markers_type;
    typedef typename levelset_space_manager_type::space_markers_type space_markers_type;
    typedef typename levelset_space_manager_type::space_markers_ptrtype space_markers_ptrtype;
    typedef typename space_markers_type::element_type element_markers_type;
    typedef std::shared_ptr<element_markers_type> element_markers_ptrtype;
    // tensor2symm
    typedef typename levelset_space_manager_type::basis_tensor2symm_type basis_tensor2symm_type;
    typedef typename levelset_space_manager_type::space_tensor2symm_type space_tensor2symm_type;
    typedef typename levelset_space_manager_type::space_tensor2symm_ptrtype space_tensor2symm_ptrtype;
    typedef typename space_tensor2symm_type::element_type element_tensor2symm_type;
    typedef std::shared_ptr< element_tensor2symm_type > element_tensor2symm_ptrtype;

    static_assert( space_levelset_type::is_scalar, "LevelSetBase function basis must be scalar" );

#if defined (MESH_ADAPTATION_LS)
    typedef MeshAdaptation<Dim, Order, 1, periodicity_type > mesh_adaptation_type;
    typedef std::shared_ptr< mesh_adaptation_type > mesh_adaptation_ptrtype;
#endif
    //--------------------------------------------------------------------//
    // Heaviside and Dirac expressions
    typedef Expr< LevelsetDeltaExpr<element_levelset_type> > levelset_delta_expr_type;

    //--------------------------------------------------------------------//
    // Range types
    typedef typename MeshTraits<mesh_type>::element_reference_wrapper_const_iterator element_reference_wrapper_const_iterator;
    typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_type elements_reference_wrapper_type;
    typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype elements_reference_wrapper_ptrtype;
    typedef elements_reference_wrapper_t<mesh_type> range_elements_type;

    typedef typename MeshTraits<mesh_type>::face_reference_wrapper_const_iterator face_reference_wrapper_const_iterator;
    typedef typename MeshTraits<mesh_type>::faces_reference_wrapper_type faces_reference_wrapper_type;
    typedef typename MeshTraits<mesh_type>::faces_reference_wrapper_ptrtype faces_reference_wrapper_ptrtype;
    typedef faces_reference_wrapper_t<mesh_type> range_faces_type;

    //--------------------------------------------------------------------//
    // Stretch and shear types
    typedef element_levelset_type element_stretch_type;
    typedef element_levelset_ptrtype element_stretch_ptrtype;

    //--------------------------------------------------------------------//
    // Tool manager
    typedef LevelSetToolManager<ConvexType, BasisType, PeriodicityType, BasisPnType> levelset_tool_manager_type;
    typedef std::shared_ptr<levelset_tool_manager_type> levelset_tool_manager_ptrtype;
    //--------------------------------------------------------------------//
    // Projectors
    typedef Projector<space_levelset_type, space_levelset_type> projector_levelset_type;
    typedef std::shared_ptr<projector_levelset_type> projector_levelset_ptrtype;
    
    typedef Projector<space_vectorial_type, space_vectorial_type> projector_levelset_vectorial_type;
    typedef std::shared_ptr<projector_levelset_vectorial_type> projector_levelset_vectorial_ptrtype;

    typedef Projector<space_tensor2symm_type, space_tensor2symm_type> projector_tensor2symm_type;
    typedef std::shared_ptr<projector_tensor2symm_type> projector_tensor2symm_ptrtype;

    //--------------------------------------------------------------------//
    // Redistanciation
    typedef LevelSetRedistanciation<space_levelset_type> redistanciation_type;
    typedef std::shared_ptr<redistanciation_type> redistanciation_ptrtype;
    typedef LevelSetRedistanciationFM<space_levelset_type> redistanciationFM_type;
    typedef std::shared_ptr<redistanciationFM_type> redistanciationFM_ptrtype;
    typedef LevelSetRedistanciationHJ<space_levelset_type> redistanciationHJ_type;
    typedef std::shared_ptr<redistanciationHJ_type> redistanciationHJ_ptrtype;

    //--------------------------------------------------------------------//
    // Initial value
    enum class ShapeType {
        SPHERE, ELLIPSE
    };
    static std::map<std::string, ShapeType> ShapeTypeMap;

    //--------------------------------------------------------------------//
    // Backend
    typedef Backend<value_type> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;

    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    //--------------------------------------------------------------------//
    // Exporter
    typedef Exporter<mymesh_type, nOrderGeo> exporter_type;
    typedef std::shared_ptr<exporter_type> exporter_ptrtype;
    typedef exporter_ptrtype exporter_manager_ptrtype;

    // Measure tools for points evaluation
    typedef MeasurePointsEvaluation<space_levelset_type> measure_points_evaluation_type;
    typedef std::shared_ptr<measure_points_evaluation_type> measure_points_evaluation_ptrtype;

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//

    //--------------------------------------------------------------------//
    // Constructor
    LevelSetBase(
            std::string const& prefix,
            std::string const& keyword = "levelset",
            worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
            std::string const& subPrefix = "",
            ModelBaseRepository const& modelRep = ModelBaseRepository() );
    LevelSetBase(
            std::string const& prefix,
            worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
            std::string const& subPrefix = "",
            ModelBaseRepository const& modelRep = ModelBaseRepository() )
        : LevelSetBase( prefix, prefix, _worldComm, subPrefix, modelRep )
    {}


    LevelSetBase( self_type const& L ) = default;

    static self_ptrtype New( 
            std::string const& prefix,
            std::string const& keyword = "levelset",
            worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
            std::string const& subPrefix = "",
            ModelBaseRepository const& modelRep = ModelBaseRepository() );

    //--------------------------------------------------------------------//
    // Initialization
    void createMesh();
    void init();
    void initLevelsetValue();
    void initUserFunctions();
    void initPostProcess() override;

    std::shared_ptr<std::ostringstream> getInfo() const override;

    //--------------------------------------------------------------------//
    // Spaces
    levelset_space_manager_ptrtype const& functionSpaceManager() const { return M_spaceManager; }
    void setFunctionSpaceManager( levelset_space_manager_ptrtype const& manager ) { M_spaceManager = manager; }

    bool buildExtendedDofSpace() const { return M_buildExtendedDofSpace; }
    void setBuildExtendedDofSpace( bool b ) { M_buildExtendedDofSpace = b; }

    space_levelset_ptrtype const& functionSpace() const { return M_spaceLevelset; }
    space_markers_ptrtype const& functionSpaceMarkers() const { return M_spaceMarkers; }
    space_scalar_ptrtype const& functionSpaceScalar() const { return M_spaceLevelset; }
    space_vectorial_ptrtype const& functionSpaceVectorial() const { return M_spaceVectorial; }

    mesh_ptrtype mesh() const { return super_type::super_model_meshes_type::mesh<mesh_type>( this->keyword() ); }
    void setMesh( mesh_ptrtype const& mesh ) { super_type::super_model_meshes_type::setMesh( this->keyword(), mesh ); }

    bool useSpaceIsoPN() const { return M_useSpaceIsoPN; }

    virtual std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"LevelsetMesh.path"); }

    mesh_ptrtype const& submeshDirac() const;
    mesh_ptrtype const& submeshOuter( double cut = 0.999 ) const;
    mesh_ptrtype const& submeshInner( double cut = 1e-3 ) const;

    periodicity_type const& periodicity() const { return M_periodicity; }

    //--------------------------------------------------------------------//
    // Mesh adaptation
#if defined (MESH_ADAPTATION_LS)
    mesh_ptrtype adaptMesh(double time, element_levelset_ptrtype elt);
    mesh_adaptation_ptrtype mesh_adapt;
    template<typename ExprT> mesh_ptrtype adaptedMeshFromExp(ExprT expr);
#endif

    //--------------------------------------------------------------------//
    // Levelset
    element_levelset_ptrtype & phiPtr() { return M_phi; }
    element_levelset_ptrtype const& phiPtr() const { return M_phi; }
    element_levelset_type & phiElt() { return *M_phi; }
    element_levelset_type const& phiElt() const { return *M_phi; }
    void setPhi( element_levelset_type const& phi, bool reinit = true ) { *M_phi = phi; M_hasRedistanciated = reinit; }
    void setPhi( element_levelset_ptrtype const& phi, bool reinit = true ) { this->setPhi( *phi, reinit ); }
    //element_levelset_ptrtype const& phinl() const { return M_phinl; }
    element_levelset_PN_ptrtype const& phiPN() const;
    element_vectorial_ptrtype const& gradPhi( bool update ) const;
    element_vectorial_ptrtype const& gradPhi() const { return this->gradPhi( M_doUpdateGradPhi ); }
    element_levelset_ptrtype const& modGradPhi( bool update ) const;
    element_levelset_ptrtype const& modGradPhi() const { return this->modGradPhi( M_doUpdateModGradPhi ); }
    element_levelset_ptrtype const& distance( bool update ) const;
    element_levelset_ptrtype const& distance() const { return this->distance( M_doUpdateDistance ); }

    element_levelset_ptrtype const& heaviside( bool update ) const;
    element_levelset_ptrtype const& heaviside() const { return this->heaviside( M_doUpdateHeaviside ); }
    element_levelset_ptrtype const& H() const { return this->heaviside(); }
    levelset_delta_expr_type diracExpr() const;
    element_levelset_ptrtype const& dirac( bool update ) const;
    element_levelset_ptrtype const& dirac() const { return this->dirac( M_doUpdateDirac ); }
    element_levelset_ptrtype const& D() const { return this->dirac(); }

    element_vectorial_ptrtype const& normal( bool update ) const;
    element_vectorial_ptrtype const& normal() const { return this->normal( M_doUpdateNormal ); }
    element_vectorial_ptrtype const& N() const { return this->normal(); }
    element_levelset_ptrtype const& curvature( bool update ) const;
    element_levelset_ptrtype const& curvature() const { return this->curvature( M_doUpdateCurvature ); }
    element_levelset_ptrtype const& K() const { return this->curvature(); }
    element_vectorial_ptrtype const& distanceNormal( bool update ) const;
    element_vectorial_ptrtype const& distanceNormal() const { return this->distanceNormal( M_doUpdateDistanceNormal ); }
    element_levelset_ptrtype const& distanceCurvature( bool update ) const;
    element_levelset_ptrtype const& distanceCurvature() const { return this->distanceCurvature( M_doUpdateDistanceCurvature ); }

    virtual void updateInterfaceQuantities();

    double thicknessInterface() const { return M_thicknessInterface; }
    void setThicknessInterface( double value ) { M_thicknessInterface = value; }

    //--------------------------------------------------------------------//
    // Interface quantities helpers
    element_vectorial_type grad( element_levelset_type const& phi, LevelSetDerivationMethod method ) const;
    element_vectorial_type grad( element_levelset_ptrtype const& phi, LevelSetDerivationMethod method ) const { return this->grad( *phi, method ); }
    element_vectorial_type grad( element_levelset_type const& phi ) const { return this->grad(phi, M_gradPhiMethod); }
    element_vectorial_type grad( element_levelset_ptrtype const& phi ) const { return this->grad(*phi); }

    element_levelset_type modGrad( element_levelset_type const& phi, LevelSetDerivationMethod method ) const;
    element_levelset_type modGrad( element_levelset_ptrtype const& phi, LevelSetDerivationMethod method ) const { return this->modGrad(*phi, method); }
    element_levelset_type modGrad( element_levelset_type const& phi ) const { return this->modGrad(phi, M_modGradPhiMethod); }
    element_levelset_type modGrad( element_levelset_ptrtype const& phi ) const { return this->modGrad(*phi); }

    //--------------------------------------------------------------------//
    // Tools
    levelset_tool_manager_ptrtype const& toolManager() const { return M_toolManager; }
    void setToolManager( levelset_tool_manager_ptrtype const& manager ) { M_toolManager = manager; }

    projector_levelset_ptrtype const& projectorL2() const { return this->projectorL2Scalar(); }
    projector_levelset_ptrtype const& projectorL2Scalar() const { return M_projectorL2Scalar; }
    projector_levelset_vectorial_ptrtype const& projectorL2Vectorial() const { return M_projectorL2Vectorial; }

    projector_levelset_ptrtype const& smoother() const { return this->projectorSMScalar(); }
    projector_levelset_ptrtype const& projectorSMScalar() const { return M_projectorSMScalar; }
    projector_levelset_vectorial_ptrtype const& smootherVectorial() const { return M_projectorSMVectorial; }
    projector_levelset_vectorial_ptrtype const& projectorSMVectorial() const { return M_projectorSMVectorial; }
    projector_levelset_ptrtype const& smootherInterface() const;
    projector_levelset_vectorial_ptrtype const& smootherInterfaceVectorial() const;

    struct FieldTag
    {
        static auto levelset_scalar( self_type const* t ) { return ModelFieldTag<self_type,0>( t ); }
        static auto levelset_vectorial( self_type const* t ) { return ModelFieldTag<self_type,1>( t ); }
    };

    //___________________________________________________________________________________//
    // toolbox fields
    //___________________________________________________________________________________//

    auto modelFields( std::string const& prefix = "" ) const
    {
        return this->modelFields( this->phiPtr(), prefix );
    }
#if 0
    auto modelFields( vector_ptrtype sol, size_type rowStartInVector = 0, std::string const& prefix = "" ) const
        {
            auto field_t = this->functionSpace()->elementPtr( *sol, rowStartInVector /*+ this->startSubBlockSpaceIndex( "field" )*/ );
            return this->modelFields( field_t, prefix );
        }
#endif
    template <typename PhiFieldType>
    auto modelFields( PhiFieldType const& field_phi, std::string const& prefix = "" ) const
    {
        typedef element_levelset_ptrtype const& (self_type::*scalar_fc_type)() const;
        typedef element_vectorial_ptrtype const& (self_type::*vec_fc_type)() const;
        return Feel::FeelModels::modelFields( 
                modelField<FieldCtx::ID|FieldCtx::GRAD|FieldCtx::GRAD_NORMAL>( FieldTag::levelset_scalar(this), prefix, "phi", field_phi, "phi", this->keyword() ),
                modelField<FieldCtx::ID>( FieldTag::levelset_scalar(this), prefix, "dirac", this->dirac(false), "dirac", this->keyword(), std::bind<scalar_fc_type>( &self_type::dirac, this ) ),
                modelField<FieldCtx::ID>( FieldTag::levelset_scalar(this), prefix, "heaviside", this->heaviside(false), "heaviside", this->keyword(), std::bind<scalar_fc_type>( &self_type::heaviside, this ) ),
                modelField<FieldCtx::ID>( FieldTag::levelset_vectorial(this), prefix, "normal", this->normal(false), "normal", this->keyword(), std::bind<vec_fc_type>( &self_type::normal, this ) ),
                modelField<FieldCtx::ID>( FieldTag::levelset_scalar(this), prefix, "curvature", this->curvature(false), "curvature", this->keyword(), std::bind<scalar_fc_type>( &self_type::curvature, this ) ),
                modelField<FieldCtx::ID>( FieldTag::levelset_vectorial(this), prefix, "gradphi", this->gradPhi(false), "gradphi", this->keyword(), std::bind<vec_fc_type>( &self_type::gradPhi, this ) ),
                modelField<FieldCtx::ID>( FieldTag::levelset_scalar(this), prefix, "modgradphi", this->modGradPhi(false), "modgradphi", this->keyword(), std::bind<scalar_fc_type>( &self_type::modGradPhi, this ) ),
                modelField<FieldCtx::ID>( FieldTag::levelset_scalar(this), prefix, "distance", this->distance(false), "distance", this->keyword(), std::bind<scalar_fc_type>( &self_type::distance, this ) ),
                modelField<FieldCtx::ID>( FieldTag::levelset_vectorial(this), prefix, "distance-normal", this->distanceNormal(false), "distance-normal", this->keyword(), std::bind<vec_fc_type>( &self_type::distanceNormal, this ) ),
                modelField<FieldCtx::ID>( FieldTag::levelset_scalar(this), prefix, "distance-curvature", this->distanceCurvature(false), "distance-curvature", this->keyword(), std::bind<scalar_fc_type>( &self_type::distanceCurvature, this ) )
                );
    }
    //___________________________________________________________________________________//
    // symbols expressions
    //___________________________________________________________________________________//

    template <typename ModelFieldsType>
    auto symbolsExpr( ModelFieldsType const& mfields ) const
    {
        //auto seToolbox = this->symbolsExprToolbox( mfields );
        auto seParam = this->symbolsExprParameter();
        //auto seMat = this->materialsProperties()->symbolsExpr();
        auto seFields = mfields.symbolsExpr();
        return Feel::vf::symbolsExpr( /*seToolbox,*/ seParam/*, seMat*/, seFields );
    }
    auto symbolsExpr( std::string const& prefix = "" ) const { return this->symbolsExpr( this->modelFields( prefix ) ); }

    // template <typename ModelFieldsType>
    // auto symbolsExprToolbox( ModelFieldsType const& mfields ) const
    //     {
    //     }
#if 0
    //--------------------------------------------------------------------//
    // Symbols expressions
    auto symbolsExpr( std::string const& prefix_symbol = "levelset_" ) const {
        return this->symbolsExpr( this->phiElt(), prefix_symbol ); 
    }
    template <typename FieldType>
    auto symbolsExpr( FieldType const& f, std::string const& prefix_symbol = "adr_" ) const
    {
        auto seField = this->symbolsExprField( f, prefix_symbol );
        return Feel::vf::symbolsExpr( seField );
    }

    auto symbolsExprField( std::string const& prefix_symbol = "levelset_" ) const { 
        return this->symbolsExprField( this->phiElt(), prefix_symbol ); 
    }
    template <typename FieldType>
    auto symbolsExprField( FieldType const& f, std::string const& prefix_symbol = "levelset_" ) const
    {
        // generate symbols levelset_phi, levelset_grad_phi(_x,_y,_z), levelset_dn_phi
        return Feel::vf::symbolsExpr( 
                symbolExpr( (boost::format("%1%phi")%prefix_symbol).str(),idv(f) ),
                symbolExpr( (boost::format("%1%grad_phi")%prefix_symbol).str(),gradv(f), SymbolExprComponentSuffix( 1,nDim ) ),
                symbolExpr( (boost::format("%1%dn_phi")%prefix_symbol).str(),dnv(f) )
                );
    }
    // Fields
    auto allFields( std::string const& prefix = "" ) const
    {
        return hana::make_tuple(
                std::make_pair( prefixvm( prefix, "phi" ), this->phiPtr() ),
                std::make_pair( prefixvm( prefix, "dirac" ), this->dirac() ),
                std::make_pair( prefixvm( prefix, "heaviside" ), this->heaviside() ),
                std::make_pair( prefixvm( prefix, "normal" ), this->normal() ),
                std::make_pair( prefixvm( prefix, "curvature" ), this->curvature() ),
                std::make_pair( prefixvm( prefix, "gradphi" ), this->gradPhi() ),
                std::make_pair( prefixvm( prefix, "modgradphi" ), this->modGradPhi() ),
                std::make_pair( prefixvm( prefix, "distance" ), this->distance() ),
                std::make_pair( prefixvm( prefix, "distance-normal" ), this->distanceNormal() ),
                std::make_pair( prefixvm( prefix, "distance-curvature" ), this->distanceCurvature() ),
                this->fieldsUserScalar(),
                this->fieldsUserVectorial()
                );
    }
#endif
    // Measures quantities
    auto allMeasuresQuantities( std::string const& prefix = "" ) const
    {
#if 0
        // TODO : we need to explicitly convert Eigen::Matrix to std::vector as
        // the begin and end iterators used in ModelNumerical::588 are only implemented from
        // Eigen v3.4 on (not released yet).
        auto eigenToVec = []( Eigen::Matrix<value_type, nDim, 1> const& m )
        {
            std::vector<value_type> v( m.size() );
            Eigen::Matrix<value_type, nDim, 1>::Map( &v[0], v.size() ) = m;
            return v;
        };
#endif
        return Feel::FeelModels::modelMeasuresQuantities( modelMeasuresQuantity( prefix, "volume", std::bind( &self_type::volume, this ) ),
                                                          modelMeasuresQuantity( prefix, "perimeter", std::bind( &self_type::perimeter, this ) ),
                                                          modelMeasuresQuantity( prefix, "position-com", std::bind( &self_type::positionCOM, this ) )
                                                          //ModelMeasuresQuantity( prefix, "position-com", std::bind( eigenToVec, std::bind( &self_type::positionCOM, this ) ) )
                                                          );
    }
    //--------------------------------------------------------------------//
    // Curvature diffusion
    bool useCurvatureDiffusion() const { return M_useCurvatureDiffusion; }
    void setUseCurvatureDiffusion( bool b ) { M_useCurvatureDiffusion = b; }

    //--------------------------------------------------------------------//
    // Markers
    element_markers_ptrtype const& markerInterface() const;
    element_markers_ptrtype const& markerDirac() const;
    element_markers_ptrtype const& markerOuter( double cut = 0.999 ) const;
    element_markers_ptrtype const& markerInner( double cut = 1e-3 ) const;
    element_markers_ptrtype const& markerHeaviside( double cut = 0.999 ) const { return this->markerOuter(cut); }
    //element_markers_ptrtype const& markerCrossedElements() const;

    //--------------------------------------------------------------------//
    // Ranges
    range_elements_type const& rangeMeshElements() const { return M_rangeMeshElements; }
    range_elements_type const& rangeInterfaceElements() const;
    range_elements_type rangeOuterElements( double cut ) const;
    range_elements_type rangeOuterElements() const { return rangeOuterElements( -this->thicknessInterface() ); }
    range_elements_type const& rangeDiracElements() const;

    range_faces_type rangeInterfaceFaces() const;

    //--------------------------------------------------------------------//
    // Utility distances
    element_levelset_ptrtype distToBoundary();
    element_levelset_ptrtype distToMarkedFaces( boost::any const& marker );
    element_levelset_ptrtype distToMarkedFaces( std::initializer_list<boost::any> marker );

    //--------------------------------------------------------------------//
    // Redistanciation
    virtual void redistanciate();
    element_levelset_type redistanciate( element_levelset_type const& phi, LevelSetDistanceMethod method ) const;

    redistanciationFM_ptrtype const& redistanciationFM( bool buildOnTheFly = true ) const;
    redistanciationHJ_ptrtype const& redistanciationHJ( bool buildOnTheFly = true ) const;

    bool hasRedistanciated() const { return M_hasRedistanciated; }

    //--------------------------------------------------------------------//
    // Initial value
    void setInitialValue(element_levelset_ptrtype const& phiv, bool doRedistanciate);
    void setInitialValue(element_levelset_ptrtype const& phiv)
    {
        this->setInitialValue(phiv, M_redistInitialValue);
    }
    template<typename ExprT>
    void setInitialValue(vf::Expr<ExprT> const& expr, bool doRedistanciate)
    {
        auto phi_init = this->functionSpace()->elementPtr();
        phi_init->on( 
                _range=this->rangeMeshElements(),
                _expr=expr
                );
        this->setInitialValue( phi_init, doRedistanciate );
    }
    template<typename ExprT>
    void setInitialValue(vf::Expr<ExprT> const& expr)
    {
        this->setInitialValue(expr, M_redistInitialValue );
    }

    element_levelset_ptrtype const& initialValue() const { return M_initialPhi; }

    double initialVolume() const { return M_initialVolume; }
    double initialPerimeter() const { return M_initialPerimeter; }

    //--------------------------------------------------------------------//
    // Export results
    exporter_ptrtype & exporter() { return M_exporter; }
    exporter_ptrtype const& exporter() const { return M_exporter; }

    virtual std::set<std::string> postProcessSaveAllFieldsAvailable() const;
    virtual std::set<std::string> postProcessExportsAllFieldsAvailable() const;

    void exportResults() { this->exportResults( this->currentTime() ); }
    virtual void exportResults( double time );
    template<typename SymbolsExpr>
    void exportResults( double time, SymbolsExpr const& symbolsExpr ) { this->exportResults( time, symbolsExpr, this->modelFields(), this->allMeasuresQuantities() ); }
    template<typename SymbolsExpr, typename TupleFieldsType, typename TupleMeasuresQuantitiesType>
    void exportResults( double time, SymbolsExpr const& symbolsExpr, TupleFieldsType const& fields, TupleMeasuresQuantitiesType const& tupleMeasuresQuantities );

    //--------------------------------------------------------------------//
    // User-defined fields
    std::map<std::string, element_scalar_ptrtype> const& fieldsUserScalar() const { return M_fieldsUserScalar; }
    std::map<std::string, element_vectorial_ptrtype> const& fieldsUserVectorial() const { return M_fieldsUserVectorial; }
    bool hasFieldUserScalar( std::string const& key ) const { return M_fieldsUserScalar.find( key ) != M_fieldsUserScalar.end(); }
    bool hasFieldUserVectorial( std::string const& key ) const { return M_fieldsUserVectorial.find( key ) != M_fieldsUserVectorial.end(); }
    element_scalar_ptrtype const& fieldUserScalarPtr( std::string const& key ) const {
        CHECK( this->hasFieldUserScalar( key ) ) << "field name " << key << " not registered"; return M_fieldsUserScalar.find( key )->second; }
    element_vectorial_ptrtype const& fieldUserVectorialPtr( std::string const& key ) const {
        CHECK( this->hasFieldUserVectorial( key ) ) << "field name " << key << " not registered"; return M_fieldsUserVectorial.find( key )->second; }
    element_scalar_type const& fieldUserScalar( std::string const& key ) const { return *this->fieldUserScalarPtr( key ); }
    element_vectorial_type const& fieldUserVectorial( std::string const& key ) const { return *this->fieldUserVectorialPtr( key ); }

    void registerCustomFieldScalar( std::string const& name )
    {
        if ( M_fieldsUserScalar.find( name ) == M_fieldsUserScalar.end() )
            M_fieldsUserScalar[name];
    }
    void registerCustomFieldVectorial( std::string const& name )
    {
        if ( M_fieldsUserVectorial.find( name ) == M_fieldsUserVectorial.end() )
            M_fieldsUserVectorial[name];
    }
    template <typename ExprT>
    void updateCustomField( std::string const& name, vf::Expr<ExprT> const& e )
    {
        this->updateCustomField( name, e, M_rangeMeshElements );
    }
    template <typename ExprT, typename OnRangeType>
    void updateCustomField( std::string const& name, vf::Expr<ExprT> const& e, OnRangeType const& range, std::enable_if_t< ExprTraits<OnRangeType, vf::Expr<ExprT>>::shape::is_scalar>* = nullptr )
    {
        if ( M_fieldsUserScalar.find( name ) == M_fieldsUserScalar.end() || !M_fieldsUserScalar[name] )
            M_fieldsUserScalar[name] = this->functionSpaceScalar()->elementPtr();
        M_fieldsUserScalar[name]->on(_range=range,_expr=e );
    }
    template <typename ExprT, typename OnRangeType>
    void updateCustomField( std::string const& name, vf::Expr<ExprT> const& e, OnRangeType const& range, std::enable_if_t< ExprTraits<OnRangeType, vf::Expr<ExprT>>::shape::is_vectorial>* = nullptr )
    {
        if ( M_fieldsUserVectorial.find( name ) == M_fieldsUserVectorial.end() || !M_fieldsUserVectorial[name] )
            M_fieldsUserVectorial[name] = this->functionSpaceVectorial()->elementPtr();
        M_fieldsUserVectorial[name]->on(_range=range,_expr=e );
    }

    void updateUserFunctions( bool onlyExprWithTimeSymbol = false );

    //--------------------------------------------------------------------//
    // Physical quantities
    double volume() const { 
        return integrate( 
                _range=this->rangeMeshElements(), 
                _expr=(1-idv(this->heaviside())) 
                ).evaluate()(0,0); 
    }
    double perimeter() const { 
        return integrate( 
                _range=this->rangeDiracElements(), 
                _expr=this->diracExpr() 
                ).evaluate()(0,0); 
    }
    Eigen::Matrix<value_type, nDim, 1> positionCOM() const { 
        return integrate( 
                _range=this->rangeMeshElements(), 
                _expr=vf::P() * (1.-idv(this->H())) 
                ).evaluate() / this->volume(); 
    }

    //--------------------------------------------------------------------//
    // Parameters
    bool redistInitialValue() const { return M_redistInitialValue; }
    void setRedistInitialValue( bool b ) { M_redistInitialValue = b; }

protected:
    //--------------------------------------------------------------------//
    void buildImpl();
    //--------------------------------------------------------------------//
    void initPostProcessExportsAndMeasures();
    //--------------------------------------------------------------------//
    // Levelset data update functions
    void updateGradPhi() const;
    void updateModGradPhi() const;
    void updateDirac() const;
    void updateHeaviside() const;

    void updateNormal() const;
    void updateCurvature() const;

    void updatePhiPN();

    void updateDistance() const;
    void updateDistanceNormal() const;
    void updateDistanceCurvature() const;

    void updateMarkerDirac();
    void markerHeavisideImpl( element_markers_ptrtype const& marker, bool invert, double cut );
    //void updateMarkerCrossedElements();
    void updateMarkerInterface();

    void updateLeftCauchyGreenTensor();

    range_elements_type rangeInterfaceElementsImpl( element_levelset_type const& phi ) const;
    range_elements_type rangeThickInterfaceElementsImpl( element_levelset_type const& phi, double thickness ) const;
    //--------------------------------------------------------------------//
    // Interface rectangular function
    element_levelset_type interfaceRectangularFunction() const { return this->interfaceRectangularFunction(this->phiElt()); }
    element_levelset_type interfaceRectangularFunction( element_levelset_ptrtype const& p ) const { return this->interfaceRectangularFunction(*p); }
    element_levelset_type interfaceRectangularFunction( element_levelset_type const& p ) const;
    //--------------------------------------------------------------------//
    // Export
    void createPostProcessExporters();
    void createPostProcessMeasures();

private:
    void loadParametersFromOptionsVm();

    void createFunctionSpaces();
    void createInterfaceQuantities();
    void createRedistanciation();
    void createRedistanciationFM();
    void createRedistanciationHJ();
    void createTools();

    //--------------------------------------------------------------------//
    void addShape( 
            pt::ptree const& shapeParameters, 
            element_levelset_type & phi );


protected:
    //--------------------------------------------------------------------//
    // Interface quantities update flags
    mutable bool M_doUpdateDirac;
    mutable bool M_doUpdateHeaviside;
    mutable bool M_doUpdateInterfaceElements;
    mutable bool M_doUpdateRangeDiracElements;
    mutable bool M_doUpdateInterfaceFaces;
    mutable bool M_doUpdateSmootherInterface;
    mutable bool M_doUpdateSmootherInterfaceVectorial;
    mutable bool M_doUpdateNormal;
    mutable bool M_doUpdateCurvature;
    mutable bool M_doUpdateGradPhi;
    mutable bool M_doUpdateModGradPhi;
    mutable bool M_doUpdatePhiPN;
    mutable bool M_doUpdateDistance;
    mutable bool M_doUpdateDistanceNormal;
    mutable bool M_doUpdateDistanceCurvature;
    mutable bool M_doUpdateMarkers;

    //--------------------------------------------------------------------//
    // Export
    exporter_ptrtype M_exporter;
    measure_points_evaluation_ptrtype M_measurePointsEvaluation;

    //--------------------------------------------------------------------//
    // User-defined fields
    std::map<std::string, element_scalar_ptrtype> M_fieldsUserScalar;
    std::map<std::string, element_vectorial_ptrtype> M_fieldsUserVectorial;

private:
    //--------------------------------------------------------------------//
    // Meshes 
    range_elements_type M_rangeMeshElements;

    mutable mesh_ptrtype M_submeshDirac;
    mutable bool M_doUpdateSubmeshDirac;
    mutable mesh_ptrtype M_submeshOuter;
    mutable bool M_doUpdateSubmeshOuter;
    mutable mesh_ptrtype M_submeshInner;
    mutable bool M_doUpdateSubmeshInner;

    //--------------------------------------------------------------------//
    // Periodicity
    periodicity_type M_periodicity;
    //--------------------------------------------------------------------//
    // Spaces
    levelset_space_manager_ptrtype M_spaceManager;
    bool M_buildExtendedDofSpace;
    bool M_useSpaceIsoPN;
    space_levelset_ptrtype M_spaceLevelset;
    space_vectorial_ptrtype M_spaceVectorial;
    space_markers_ptrtype M_spaceMarkers;

    //--------------------------------------------------------------------//
    // Fields
    element_levelset_ptrtype M_phi;
    element_levelset_ptrtype M_initialPhi;

    //--------------------------------------------------------------------//
    // Markers
    mutable element_markers_ptrtype M_markerDirac;
    mutable element_markers_ptrtype M_markerOuter;
    mutable double M_markerOuterCut;
    mutable element_markers_ptrtype M_markerInner;
    mutable double M_markerInnerCut;
    //mutable element_markers_ptrtype M_markerCrossedElements;
    mutable element_markers_ptrtype M_markerInterface;

    //--------------------------------------------------------------------//
    // Ranges
    mutable range_elements_type M_rangeDiracElements;
    //--------------------------------------------------------------------//
    // Tools (projectors)
    levelset_tool_manager_ptrtype M_toolManager;

    projector_levelset_ptrtype M_projectorL2Scalar;
    projector_levelset_vectorial_ptrtype M_projectorL2Vectorial;

    projector_levelset_ptrtype M_projectorSMScalar;
    projector_levelset_vectorial_ptrtype M_projectorSMVectorial;
    mutable projector_levelset_ptrtype M_smootherInterface;
    mutable projector_levelset_vectorial_ptrtype M_smootherInterfaceVectorial;
    //--------------------------------------------------------------------//
    // Levelset data
    mutable element_levelset_PN_ptrtype M_levelsetPhiPN;
    mutable element_vectorial_ptrtype M_levelsetGradPhi;
    mutable element_levelset_ptrtype M_levelsetModGradPhi;
    mutable element_levelset_ptrtype M_heaviside;
    mutable element_levelset_ptrtype M_dirac;
    mutable element_levelset_ptrtype M_distance;

    mutable range_elements_type M_interfaceElements;
    mutable range_faces_type M_interfaceFaces;
    //--------------------------------------------------------------------//
    // Normal, curvature
    mutable element_vectorial_ptrtype M_levelsetNormal;
    mutable element_levelset_ptrtype M_levelsetCurvature;
    mutable element_vectorial_ptrtype M_distanceNormal;
    mutable element_levelset_ptrtype M_distanceCurvature;

    //--------------------------------------------------------------------//
    // Derivation methods
    LevelSetDerivationMethod M_gradPhiMethod;
    LevelSetDerivationMethod M_modGradPhiMethod;
    LevelSetCurvatureMethod M_curvatureMethod;

    //--------------------------------------------------------------------//
    // Curvature diffusion
    bool M_useCurvatureDiffusion;

    //--------------------------------------------------------------------//
    // Redistanciation
    static const std::map<std::string, LevelSetDistanceMethod> LevelSetDistanceMethodIdMap;

    LevelSetDistanceMethod M_distanceMethod;

    redistanciationFM_ptrtype M_redistanciationFM;
    redistanciationHJ_ptrtype M_redistanciationHJ;
    bool M_redistanciationIsUpdatedForUse;

    LevelSetDistanceMethod M_redistanciationMethod;

    bool M_redistInitialValue;

    bool M_hasRedistanciated;

    //--------------------------------------------------------------------//
    // Parameters
    double M_thicknessInterface;
    double M_thicknessInterfaceRectangularFunction;
    bool M_useAdaptiveThicknessInterface;
    bool M_useRegularPhi;
    bool M_useHeavisideDiracNodalProj;

    double M_initialVolume;
    double M_initialPerimeter;

    //LevelSetTimeDiscretization M_discrMethod;

}; //class LevelSetBase

//----------------------------------------------------------------------------//
// Static member initialization
#define LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS \
    template< typename ConvexType, typename BasisType, typename PeriodicityType, typename BasisPnType > \
        /**/
#define LEVELSETBASE_CLASS_TEMPLATE_TYPE \
    LevelSetBase<ConvexType, BasisType, PeriodicityType, BasisPnType> \
        /**/
LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
std::map<std::string, typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::ShapeType>
LEVELSETBASE_CLASS_TEMPLATE_TYPE::ShapeTypeMap = {
    {"sphere", ShapeType::SPHERE},
    {"ellipse", ShapeType::ELLIPSE}
};

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
const std::map<std::string, LevelSetDistanceMethod> 
LEVELSETBASE_CLASS_TEMPLATE_TYPE::LevelSetDistanceMethodIdMap = {
    {"none", LevelSetDistanceMethod::NONE},
    {"fm", LevelSetDistanceMethod::FASTMARCHING},
    {"hj", LevelSetDistanceMethod::HAMILTONJACOBI},
    {"renormalisation", LevelSetDistanceMethod::RENORMALISATION}
};

//----------------------------------------------------------------------------//
LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
template<typename SymbolsExpr, typename ModelFieldsType, typename TupleMeasuresQuantitiesType>
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::exportResults( double time, SymbolsExpr const& symbolsExpr, ModelFieldsType const& mfields, TupleMeasuresQuantitiesType const& tupleMeasuresQuantities )
{
    this->log("LevelSetBase","exportResults", "start");
    this->timerTool("PostProcessing").start();

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().postProcess().setParameterValues( paramValues );

    this->executePostProcessExports( M_exporter, time, mfields, symbolsExpr );
    this->executePostProcessMeasures( time, this->mesh(), this->rangeMeshElements(), M_measurePointsEvaluation, symbolsExpr, mfields, tupleMeasuresQuantities );

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("LevelSetBase","exportResults", "finish");
}

} // namespace FeelModels
} // namespace Feel

#endif

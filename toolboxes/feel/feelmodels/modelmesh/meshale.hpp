/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 2011-07-17

 Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
 \file meshale.hpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2011-07-17
 */

#ifndef FEELPP_MODELS_MESHALE_H
#define FEELPP_MODELS_MESHALE_H 1

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
//#include <feel/feelvf/vf.hpp>
#include <feel/feelmesh/meshmover.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelfilters/exporter.hpp>

#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/projectors.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/geometricdata.hpp>


#include <feel/feelmodels/modelcore/feelmodelscoreconstconfig.hpp>
#include <feel/feelmodels/modelcore/modelbase.hpp>
#include <feel/feelmodels/modelmesh/ale.hpp>


#include <feel/feelmodels/modelmesh/dofrelationshipmap.hpp>

namespace Feel
{
namespace FeelModels
{

template< class Convex >
class MeshALE : public ModelBase
{

public :
    typedef ModelBase super_type;
    using self_type = MeshALE<Convex>;
    using size_type = uint32_type;
    typedef Backend<double> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;

    /*
     * Moving mesh typedefs
     */
    typedef Convex convex_type;
    typedef Mesh< convex_type > mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    using range_elements_type = elements_reference_wrapper_t<mesh_type>;

    /*
     * Reference mesh typedefs
     */
    typedef typename mesh_type::P1_mesh_type  mesh_ref_type;
    typedef std::shared_ptr<mesh_ref_type> mesh_ref_ptrtype;
    typedef typename mesh_ref_type::shape_type convex_ref_type;
    using range_elements_ref_type = elements_reference_wrapper_t<mesh_ref_type>;

    typedef ALE< convex_ref_type, mesh_type::nOrder > ale_map_type;
    typedef std::shared_ptr< ale_map_type > ale_map_ptrtype;

    typedef typename ale_map_type::ale_map_functionspace_type ale_map_functionspace_ref_type;
    typedef std::shared_ptr<ale_map_functionspace_ref_type> ale_map_functionspace_ref_ptrtype;

    typedef typename ale_map_type::ale_map_element_type ale_map_element_ref_type;
    typedef std::shared_ptr<ale_map_element_ref_type> ale_map_element_ref_ptrtype;

    typedef typename ale_map_type::ale_map_basis_type ale_map_basis_type;
    typedef FunctionSpace<mesh_type, ale_map_basis_type> ale_map_functionspace_type;
    typedef std::shared_ptr<ale_map_functionspace_type> ale_map_functionspace_ptrtype;
    typedef typename ale_map_functionspace_type::element_type ale_map_element_type;
    typedef std::shared_ptr<ale_map_element_type> ale_map_element_ptrtype;


    typedef Bdf< ale_map_functionspace_ref_type > bdf_ale_displacement_ref_type;
    typedef std::shared_ptr<bdf_ale_displacement_ref_type> bdf_ale_displacement_ref_ptrtype;

    typedef Bdf< ale_map_functionspace_type > bdf_ale_displacement_type;
    typedef std::shared_ptr<bdf_ale_displacement_type> bdf_ale_displacement_ptrtype;

    typedef Exporter<mesh_type,mesh_type::nOrder> exporter_type;
    typedef std::shared_ptr<exporter_type> exporter_ptrtype;

    typedef Exporter<mesh_ref_type,mesh_ref_type::nOrder> exporter_ref_type;
    typedef std::shared_ptr<exporter_ref_type> exporter_ref_ptrtype;

    typedef DofRelationshipMap<ale_map_functionspace_ref_type,ale_map_functionspace_type > DofRelationshipMap_type;
    typedef std::shared_ptr<DofRelationshipMap_type> DofRelationshipMap_ptrtype;


    struct ComputationalDomain
    {
    public :
        ComputationalDomain( self_type const* meshALE );
        ComputationalDomain( self_type const* meshALE, range_elements_type const& rangeElt );

        void init( bool M_isARestart );

        void generateMap( ale_map_element_type const& displacementOnMovingBoundary_HO_ref );

        void updateDisplacement( ale_map_element_type & displacement );

        //! add boundary marker with respect to type \bctype
        void addBoundaryFlags(std::string const& bctype, std::string const& marker);

        //! return the object that manage the ale computation
        ale_map_ptrtype const& aleFactory() const { return M_aleFactory; }

        //! return the object that manage the ale computation
        ale_map_ptrtype & aleFactory() { return M_aleFactory; }

        //! update meshale object to the next time step
        void updateTimeStep();
    private:
        void build();
    private :
        self_type const* M_meshALE;
        ale_map_ptrtype M_aleFactory;
        // the distance between P1_ref and HO_ref
        std::shared_ptr<ale_map_element_ref_type> M_dispP1ToHO_ref;
        std::shared_ptr<ale_map_element_ref_type> M_displacementOnMovingBoundary_P1_ref;
        std::shared_ptr<ale_map_element_ref_type> M_displacement_ref;
        bdf_ale_displacement_ref_ptrtype M_bdf_ale_displacement_ref;
        DofRelationshipMap_ptrtype M_drm;
        exporter_ref_ptrtype M_exporter_ref;
    };



    MeshALE( mesh_ptrtype mesh_moving,
             std::string const& prefix="",
             worldcomm_ptr_t const& worldcomm = Environment::worldCommPtr(),
             ModelBaseRepository const& modelRep = ModelBaseRepository() );

    void init();

    std::shared_ptr<std::ostringstream> getInfo() const;


    void setWholeMeshAsComputationalDomain( std::string const& name );
    void setComputationalDomain( std::string const& name, range_elements_type const& rangeElt );

    void addBoundaryFlags( std::string const& bctype, std::string const& marker );
    void addBoundaryFlags( std::string const& bctype, std::vector<std::string> const& markers )
        {
            for ( std::string const& marker : markers )
                this->addBoundaryFlags( bctype,marker );
        }

    /**
     * \return the reference mesh
     */
    mesh_ref_ptrtype referenceMesh() const { return M_referenceMesh; }

    /**
     * \return the moving mesh
     */
    mesh_ptrtype movingMesh() const { return M_movingMesh; }

    /**
     * \return true is reverted on reference mesh, else false
     */
    bool isOnReferenceMesh() const { return M_isOnReferenceMesh; }
    /**
     * \return true is on moving mesh, else false
     */
    bool isOnMovingMesh() const { return M_isOnMovingMesh; }
    /**
     * \return the functionspace
     */
    ale_map_functionspace_ptrtype functionSpace() const { return M_Xhmove; }

    /**
     *
     */
    ale_map_element_ptrtype identityALE() const { return M_identity_ale; }

    /**
     *
     */
    ale_map_element_ptrtype displacementOnMovingBoundary() const { return M_displacementOnMovingBoundary_HO_ref; }

    /**
     * \return the displacement
     */
    ale_map_element_ptrtype displacement() const { return M_displacement; }

    /**
     * \return the velocity mesh
     */
    ale_map_element_ptrtype velocity() const { return M_meshVelocity; }

    //! update displacement imposed from an expression \expr on entities defined by \range
    template <typename ExprType,typename RangeType >
    void updateDisplacementImposed( ExprType const& expr, RangeType const& range );

    //! update displacement imposed from a velocity expression \v on entities defined by \range (use bdf scheme)
    template <typename ExprT, typename RangeType >
    void updateDisplacementImposedFromVelocity( Expr<ExprT> const& v, RangeType const& range )
        {
            auto dispExpr = (v + idv(M_bdf_ale_identity->polyDeriv()))/M_bdf_ale_identity->polyDerivCoefficient(0) - idv(this->identityALE()) + idv( this->displacement() );
            this->updateDisplacementImposed( dispExpr, range );
        }

    //! update the field \d  from a velocity expression \v on entities defined by \range (use bdf scheme)
    template < typename FieldDispType, typename ExprT, typename RangeType >
    void updateDisplacementFieldFromVelocity( FieldDispType & d, Expr<ExprT> const& v, RangeType const& range, bool doParallelSync = false ) const
        {
            unwrap_ptr(d).on(_range=range,
                             _expr= (v + idv(M_bdf_ale_identity->polyDeriv()))/M_bdf_ale_identity->polyDerivCoefficient(0) - idv(this->identityALE()) + idv( this->displacement() ),
                             _close=doParallelSync/*,_geomap=this->geomap()*/ );
        }

    auto displacementExprAtPreviousTime() const
        {
            return idv( M_bdf_ale_identity->unknown(0) ) - idv(this->identityALE()) + idv( this->displacement() );
        }

    //! parallel sync of values given
    void updateDisplacementImposedForUse();

    /**
     * \Revert mesh in reference state, \updateMeshMeasures boolean put to false avoid to recompute some mesh measure as the aera,...
     */
    void revertReferenceMesh( bool updateMeshMeasures = true );

    /**
     * \Revert mesh in moving state, \updateMeshMeasures boolean put to false avoid to recompute some mesh measure as the aera,...
     */
    void revertMovingMesh( bool updateMeshMeasures = true );

    /**
     * update meshale object to the next time step
     */
    void updateTimeStep();

    /**
     * \Export visualisation results
     */
    void exportResults(double time=0);

    void updateForUseOnMovingMesh();

    //! update the moving mesh by using displacement imposed given
    void updateMovingMesh() { this->updateImpl(); }

#if 0
    void updateMetricMeshAdaptation( Expr<GinacExVF<2>> const& e )
        {
            M_aleFactory->updateMetricMeshAdaptation( e );
        }
#endif
    template < typename ExprT >
    void updateMetricMeshAdaptation( Expr<ExprT> const& e )
        {
            for ( auto & [name,cd] : M_computationalDomains )
                cd.aleFactory()->updateMetricMeshAdaptation( e );
        }

private :

    void updateIdentityMap();
    void initTimeStep();
    void updateImpl();

private :

    MeshMover<mesh_type> M_mesh_mover;

    mesh_ref_ptrtype M_referenceMesh;
    mesh_ptrtype M_movingMesh;

    bool M_isOnReferenceMesh, M_isOnMovingMesh;

    ale_map_functionspace_ptrtype M_Xhmove;
    std::shared_ptr<ale_map_element_type> M_identity_ale;
    std::shared_ptr<ale_map_element_type> M_displacementOnMovingBoundary_HO_ref;
    std::shared_ptr<ale_map_element_type> M_displacement;
    std::shared_ptr<ale_map_element_type> M_meshVelocity;
    std::shared_ptr<ale_map_element_type> M_fieldTmp;

    bdf_ale_displacement_ptrtype M_bdf_ale_identity;
    bdf_ale_displacement_ptrtype M_bdf_ale_velocity;

    exporter_ptrtype M_exporter;

    bool M_isARestart;
    std::string M_restartPath;

    std::set<size_type> M_dofsMultiProcessOnMovingBoundary_HO;

    std::map<std::string,ComputationalDomain> M_computationalDomains;
};

//------------------------------------------------------------------------------------------------//

template< class Convex >
template <typename ExprType,typename RangeType >
void
MeshALE<Convex>::updateDisplacementImposed( ExprType const& expr, RangeType const& range )
{
    M_displacementOnMovingBoundary_HO_ref->on(_range=range, _expr=expr );

    auto dofMultiProcessOnRange = M_displacementOnMovingBoundary_HO_ref->functionSpace()->dofs( range, ComponentType::NO_COMPONENT, true );
    M_dofsMultiProcessOnMovingBoundary_HO.insert( dofMultiProcessOnRange.begin(), dofMultiProcessOnRange.end() );
}

//------------------------------------------------------------------------------------------------//




template<typename Args>
struct compute_meshale_return
{
    typedef typename boost::remove_reference<typename parameter::binding<Args, tag::mesh>::type>::type::element_type mesh_type;
    using index_type = typename mesh_type::index_type;
    typedef typename mesh_type::shape_type convex_type;
    typedef MeshALE<convex_type> type;
    typedef std::shared_ptr<type> ptrtype;
};

BOOST_PARAMETER_FUNCTION(
    ( typename compute_meshale_return<Args>::ptrtype ), // 1. return type
    meshale,                        // 2. name of the function template
    tag,                                        // 3. namespace of tag types
    ( required
      ( mesh,    *( boost::is_convertible<mpl::_,std::shared_ptr<MeshBase<>> > ) )
      ) // required
    ( optional
      ( prefix,            (std::string), std::string("") )
      ( worldcomm,         (worldcomm_ptr_t), mesh->worldCommPtr() )
      ( directory,         (ModelBaseRepository),  ModelBaseRepository() )
      ) // optionnal
                         )
{
    typedef typename compute_meshale_return<Args>::ptrtype meshale_ptrtype;
    typedef typename compute_meshale_return<Args>::type meshale_type;
    return meshale_ptrtype( new meshale_type(mesh,prefix,worldcomm,directory) );
}


} // namespace FeelModels
} // namespace Feel
#endif // FEELPP_MODELS_MESHALE_H

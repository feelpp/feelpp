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
    using size_type = uint32_type;
    typedef Backend<double> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;

    /*
     * Moving mesh typedefs
     */
    typedef Convex convex_type;
    typedef Mesh< convex_type > mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    /*
     * Reference mesh typedefs
     */
    typedef typename mesh_type::P1_mesh_type  mesh_ref_type;
    typedef std::shared_ptr<mesh_ref_type> mesh_ref_ptrtype;
    typedef typename mesh_ref_type::shape_type convex_ref_type;

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


    MeshALE( mesh_ptrtype mesh_moving,
             std::string const& prefix="",
             worldcomm_ptr_t const& worldcomm = Environment::worldCommPtr(),
             ModelBaseRepository const& modelRep = ModelBaseRepository() );

    void init();

    std::shared_ptr<std::ostringstream> getInfo() const;


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
    bool isOnMovingMesh() const { return M_isOnMovingMesh; }
    /**
     * \return the functionspace
     */
    ale_map_functionspace_ptrtype functionSpace() const { return M_Xhmove; }

    /**
     * \return the functionspace in reference mesh
     */
    ale_map_functionspace_ref_ptrtype functionSpaceInRef() const { return M_Xhref; }

    /**
     * \return ale object define (only!) on P1_ref
     */
    ale_map_ptrtype const& aleFactory() const { return M_aleFactory; }

    /**
     * \return the ale map in P1_ref
     */
    ale_map_element_ref_type const& mapInRef() const;

    /**
     * \return the ale map in HO_ref
     */
    ale_map_element_type map();

    /**
     *
     */
    ale_map_element_ptrtype identityALE() const { return M_identity_ale; }

    /**
     * \return the distance between P1_ref and HO_ref
     */
    ale_map_element_ref_ptrtype dispP1ToHO_ref() const { return M_dispP1ToHO_ref; }

    /**
     *
     */
    ale_map_element_ptrtype displacementOnMovingBoundary() const { return M_displacementOnMovingBoundary_HO_ref; }

    /**
     *
     */
    ale_map_element_ref_ptrtype displacementOnMovingBoundaryInRef() const { return M_displacementOnMovingBoundary_P1_ref; }

    /**
     * \return the displacement
     */
    ale_map_element_ref_ptrtype displacementInRef() const { return M_displacement_ref; }

    /**
     * \return the displacement
     */
    ale_map_element_ptrtype displacement() const { return M_displacement; }

    /**
     * \return the velocity mesh
     */
    ale_map_element_ptrtype velocity() const { return M_meshVelocity; }

    /*
     * \ return dofRelationShipMap
     */
    DofRelationshipMap_ptrtype dofRelationShipMap() const { return M_drm; }


    /**
     * \Compute ale map and move the mesh
     */
    template< typename elem_type >
    void update( std::vector<elem_type> const& polyDisplacementSet );
    template< typename elem_type >
    void update( elem_type const& polyDisplacementSet );

    /**
     * \Revert mesh_ho in reference state
     */
    void revertReferenceMesh( bool updateMeshMeasures = true );

    /**
     * \Revert mesh_ho in reference state
     */
    void revertMovingMesh( bool updateMeshMeasures = true );

    /**
     * \Revert mesh_ho in moving state
     */
    void updateTimeStep();

    /**
     * \Export
     */
    void exportResults(double time=0);

    void updateMetricMeshAdaptation( Expr<GinacExVF<2>> const& e )
        {
            M_aleFactory->updateMetricMeshAdaptation( e );
        }
    template < typename ExprT >
    void updateMetricMeshAdaptation( Expr<ExprT> const& e )
        {
            M_aleFactory->updateMetricMeshAdaptation( e );
        }

    template < typename FieldDispType, typename ExprT, typename RangeType >
    void updateDisplacementFieldFromVelocity( FieldDispType & d, Expr<ExprT> const& v, RangeType const& range ) const
        {
            unwrap_ptr(d).on(_range=range,
                             _expr= (v + idv(M_bdf_ale_identity->polyDeriv()))/M_bdf_ale_identity->polyDerivCoefficient(0) - idv(this->identityALE()) + idv( this->displacement() ) );
        }

private :

    void updateIdentityMap();
    void initTimeStep();
    void updateImpl();

private :


    //backend_ptrtype M_backend;

    MeshMover<mesh_type> M_mesh_mover;

    mesh_ref_ptrtype M_referenceMesh;
    mesh_ptrtype M_movingMesh;

    bool M_isOnReferenceMesh, M_isOnMovingMesh;

    ale_map_ptrtype M_aleFactory;

    ale_map_functionspace_ref_ptrtype M_Xhref;
    ale_map_functionspace_ptrtype M_Xhmove;

    std::shared_ptr<ale_map_element_type> M_identity_ale;
    std::shared_ptr<ale_map_element_ref_type> M_dispP1ToHO_ref;

    std::shared_ptr<ale_map_element_type> M_displacementOnMovingBoundary_HO_ref;
    std::shared_ptr<ale_map_element_ref_type> M_displacementOnMovingBoundary_P1_ref;

    std::shared_ptr<ale_map_element_type> M_displacement;
    std::shared_ptr<ale_map_element_ref_type> M_displacement_ref;
    std::shared_ptr<ale_map_element_ref_type> M_map_ref;
    std::shared_ptr<ale_map_element_type> M_meshVelocity;
    std::shared_ptr<ale_map_element_type> M_fieldTmp;

    bdf_ale_displacement_ref_ptrtype M_bdf_ale_displacement_ref;
    bdf_ale_displacement_ptrtype M_bdf_ale_identity;
    bdf_ale_displacement_ptrtype M_bdf_ale_velocity;

    DofRelationshipMap_ptrtype M_drm;

    exporter_ptrtype M_exporter;
    exporter_ref_ptrtype M_exporter_ref;

    bool M_isARestart;
    std::string M_restartPath;

    std::set<size_type> M_dofsMultiProcessOnMovingBoundary_HO;
};

//------------------------------------------------------------------------------------------------//

template< class Convex >
template< typename elem_type >
void
MeshALE<Convex>::update( std::vector<elem_type> const& polyDisplacementSet )
{
    CHECK( polyDisplacementSet.size() == 1 ) << "invalid size";
    this->update( polyDisplacementSet[0] );
}
template< class Convex >
template< typename elem_type >
void
MeshALE<Convex>::update( elem_type const& polyDisplacementSet )
{
    this->log(prefixvm(this->prefix(),"MeshALE"),"update", "start");

    CHECK( this->isOnMovingMesh() ) << "meshALE must be on moving mesh\n";

    M_displacementOnMovingBoundary_HO_ref->on(_range=markedfaces(M_movingMesh,this->aleFactory()->flagSet("moving")),
        _expr=vf::idv(polyDisplacementSet) );
    sync( *M_displacementOnMovingBoundary_HO_ref, "=", M_dofsMultiProcessOnMovingBoundary_HO );

    this->updateImpl();

    this->log(prefixvm(this->prefix(),"MeshALE"),"update", "finish");
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

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
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelmesh/meshmover.hpp>
#include <feel/feelts/bdf.hpp>

#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/geometricdata.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/projectors.hpp>

#include <feel/feelmodels/modelcore/feelmodelscoreconstconfig.hpp>
#include <feel/feelmodels/modelcore/modelbase.hpp>
#include <feel/feelmodels/modelmesh/ale.hpp>

#include <feel/feelmodels/modelalg/functionSup.cpp>
#include <feel/feelmodels/modelmesh/dofrelationshipmap.hpp>

namespace Feel
{
namespace FeelModels
{

template <class Convex>
class MeshALE : public ModelBase
{

  public:
    typedef ModelBase super_type;

    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*
     * Moving mesh typedefs
     */
    typedef Convex convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*
     * Reference mesh typedefs
     */
    typedef typename mesh_type::P1_mesh_type mesh_ref_type;
    typedef boost::shared_ptr<mesh_ref_type> mesh_ref_ptrtype;
    typedef typename mesh_ref_type::shape_type convex_ref_type;

    typedef ALE<convex_ref_type, mesh_type::nOrder> ale_map_type;
    typedef boost::shared_ptr<ale_map_type> ale_map_ptrtype;

    typedef typename ale_map_type::ale_map_functionspace_type ale_map_functionspace_ref_type;
    typedef boost::shared_ptr<ale_map_functionspace_ref_type> ale_map_functionspace_ref_ptrtype;

    typedef typename ale_map_type::ale_map_element_type ale_map_element_ref_type;
    typedef boost::shared_ptr<ale_map_element_ref_type> ale_map_element_ref_ptrtype;

    typedef typename ale_map_type::ale_map_basis_type ale_map_basis_type;
    typedef FunctionSpace<mesh_type, ale_map_basis_type> ale_map_functionspace_type;
    typedef boost::shared_ptr<ale_map_functionspace_type> ale_map_functionspace_ptrtype;
    typedef typename ale_map_functionspace_type::element_type ale_map_element_type;
    typedef boost::shared_ptr<ale_map_element_type> ale_map_element_ptrtype;

    typedef FunctionSpace<mesh_type, bases<Lagrange<mesh_type::nOrder /*-1*/, Vectorial, Discontinuous>>> ale_map_functionspacedisc_type;
    typedef boost::shared_ptr<ale_map_functionspacedisc_type> ale_map_functionspacedisc_ptrtype;
    typedef typename ale_map_functionspacedisc_type::element_type ale_map_elementdisc_type;
    typedef boost::shared_ptr<ale_map_elementdisc_type> ale_map_elementdisc_ptrtype;

    typedef Bdf<ale_map_functionspace_ref_type> bdf_ale_displacement_ref_type;
    typedef boost::shared_ptr<bdf_ale_displacement_ref_type> bdf_ale_displacement_ref_ptrtype;

    typedef Bdf<ale_map_functionspace_type> bdf_ale_displacement_type;
    typedef boost::shared_ptr<bdf_ale_displacement_type> bdf_ale_displacement_ptrtype;

    typedef Exporter<mesh_type, mesh_type::nOrder> exporter_type;
    typedef boost::shared_ptr<exporter_type> exporter_ptrtype;

    typedef Exporter<mesh_ref_type, mesh_ref_type::nOrder> exporter_ref_type;
    typedef boost::shared_ptr<exporter_ref_type> exporter_ref_ptrtype;

    typedef DofRelationshipMap<ale_map_functionspace_ref_type, ale_map_functionspace_type> DofRelationshipMap_type;
    typedef boost::shared_ptr<DofRelationshipMap_type> DofRelationshipMap_ptrtype;

    MeshALE( mesh_ptrtype mesh_moving,
             //po::variables_map const& vm=Environment::vm(),
             std::string const& prefix = "",
             //std::string exportName="ExportMeshALE",
             WorldComm const& worldcomm = Environment::worldComm(),
             bool moveGhostEltFromExtendedStencil = false,
             std::string const& rootRepository = ModelBase::rootRepositoryByDefault() );

    void init();

    boost::shared_ptr<std::ostringstream> getInfo() const;

    void addBoundaryFlags( std::string __type, std::string __marker );

    /**
     * \return the reference mesh
     */
    mesh_ref_ptrtype referenceMesh();

    /**
     * \return the moving mesh
     */
    mesh_ptrtype movingMesh();

    /**
     * \return true is reverted on reference mesh, else false
     */
    bool isOnReferenceMesh() const { return M_isOnReferenceMesh; }
    bool isOnMovingMesh() const { return M_isOnMovingMesh; }
    /**
     * \return the functionspace
     */
    ale_map_functionspace_ptrtype functionSpace();

    /**
     * \return the functionspace in reference mesh
     */
    ale_map_functionspace_ref_ptrtype functionSpaceInRef();

    /**
     * \return ale object define (only!) on P1_ref
     */
    ale_map_ptrtype const& aleFactory() const;

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
    ale_map_element_ptrtype identityALE();

    /**
     * \return the distance between P1_ref and HO_ref
     */
    ale_map_element_ref_ptrtype dispP1ToHO_ref();

    /**
     *
     */
    ale_map_element_ptrtype displacementOnMovingBoundary();

    /**
     *
     */
    ale_map_element_ref_ptrtype displacementOnMovingBoundaryInRef();

    /**
     * \return the displacement
     */
    ale_map_element_ref_ptrtype displacementInRef();

    /**
     * \return the displacement
     */
    ale_map_element_ptrtype displacement();

    /**
     * \return the velocity mesh
     */
    ale_map_element_ptrtype velocity();
    ale_map_element_ptrtype const& velocity() const;

    /*
     * \ return dofRelationShipMap
     */
    DofRelationshipMap_ptrtype dofRelationShipMap();

    //bool doExport() const { return M_doExport; }

    /**
     * \Compute ale map and move the mesh
     */
    template <typename elem_type>
    void update( std::vector<elem_type> const& polyDisplacementSet );
    template <typename elem_type>
    void update( elem_type const& polyDisplacementSet );

    void updateImpl( Vector<double> const& dispInput );

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
    void updateBdf();

    /**
     * \Export
     */
    void exportResults( double time = 0 );

  private:
    void updateIdentityMap();

  private:
    //backend_ptrtype M_backend;

    MeshMover<mesh_type> M_mesh_mover;

    mesh_ref_ptrtype M_referenceMesh;
    mesh_ptrtype M_movingMesh;

    bool M_isOnReferenceMesh, M_isOnMovingMesh;

    ale_map_ptrtype M_aleFactory;

    ale_map_functionspace_ref_ptrtype M_Xhref;
    ale_map_functionspace_ptrtype M_Xhmove;

    boost::shared_ptr<ale_map_element_type> M_identity_ale;
    boost::shared_ptr<ale_map_element_ref_type> M_dispP1ToHO_ref;

    boost::shared_ptr<ale_map_element_type> M_displacementOnMovingBoundary_HO_ref;
    boost::shared_ptr<ale_map_element_ref_type> M_displacementOnMovingBoundary_P1_ref;

    boost::shared_ptr<ale_map_element_type> M_displacement;
    boost::shared_ptr<ale_map_element_ref_type> M_displacement_ref;
    boost::shared_ptr<ale_map_element_ref_type> M_map_ref;
    boost::shared_ptr<ale_map_element_type> M_meshVelocity;

    bdf_ale_displacement_ref_ptrtype M_bdf_ale_displacement_ref;
    bdf_ale_displacement_ptrtype M_bdf_ale_identity;
    bdf_ale_displacement_ptrtype M_bdf_ale_velocity;

    DofRelationshipMap_ptrtype M_drm;

    exporter_ptrtype M_exporter;
    exporter_ref_ptrtype M_exporter_ref;
    int M_cpt_export;

    bool M_isARestart;
    std::string M_restartPath;

    //bool M_doExport;
};

//------------------------------------------------------------------------------------------------//

template <class Convex>
template <typename elem_type>
void MeshALE<Convex>::update( std::vector<elem_type> const& polyDisplacementSet )
{
    CHECK( polyDisplacementSet.size() == 1 ) << "invalid size";
    this->update( polyDisplacementSet[0] );
}
template <class Convex>
template <typename elem_type>
void MeshALE<Convex>::update( elem_type /*std::vector<elem_type>*/ const& polyDisplacementSet )
{
    this->log( prefixvm( this->prefix(), "MeshALE" ), "update", "start" );

    CHECK( this->isOnMovingMesh() ) << "meshALE must be on moving mesh\n";

    //---------------------------------------------------------------------------------------------//
    // reset to 0
    //M_displacementOnMovingBoundary_HO_ref->zero();

    //---------------------------------------------------------------------------------------------//
    // interp disp to ref_ho
    auto vecDisp = backend()->newVector( M_displacementOnMovingBoundary_HO_ref->functionSpace() );
    for ( uint16_type i = 0; i < this->aleFactory()->flagSet( "moving" ).size(); ++i )
    {
#if 1
        modifVec( markedfaces( M_movingMesh, this->aleFactory()->flagSet( "moving", i ) ),
                  *M_displacementOnMovingBoundary_HO_ref,
                  vecDisp /*M_displacementOnMovingBoundary_HO_ref*/,
                  vf::idv( polyDisplacementSet /*[i]*/ ) );
#else
//M_displacementOnMovingBoundary_HO_ref->on(_range=markedfaces( M_movingMesh, this->aleFactory()->flagSet("moving",i) ),_expr=vf::idv(polyDisplacementSet) );
#endif
    }
    vecDisp->close();
    //*M_displacementOnMovingBoundary_HO_ref = *vecDisp;

    this->updateImpl( *vecDisp );

    this->log( prefixvm( this->prefix(), "MeshALE" ), "update", "finish" );
}

//------------------------------------------------------------------------------------------------//

template <typename Args>
struct compute_meshale_return
{
    typedef typename boost::remove_reference<typename parameter::binding<Args, tag::mesh>::type>::type::element_type mesh_type;
    typedef typename mesh_type::shape_type convex_type;
    typedef MeshALE<convex_type> type;
    typedef boost::shared_ptr<type> ptrtype;
};

BOOST_PARAMETER_FUNCTION(
    ( typename compute_meshale_return<Args>::ptrtype ),                                                                                                                                                            // 1. return type
    meshale,                                                                                                                                                                                                       // 2. name of the function template
    tag,                                                                                                                                                                                                           // 3. namespace of tag types
    ( required( mesh, *(boost::is_convertible<mpl::_, boost::shared_ptr<MeshBase>>)) )                                                                                                                             // required
    ( optional( prefix, ( std::string ), std::string( "" ) )( worldcomm, ( WorldComm ), Environment::worldComm() )( extended_doftable, (bool), true )( directory, ( std::string ), fs::current_path().string() ) ) // optionnal
    )
{
    typedef typename compute_meshale_return<Args>::ptrtype meshale_ptrtype;
    typedef typename compute_meshale_return<Args>::type meshale_type;
    return meshale_ptrtype( new meshale_type( mesh, prefix, worldcomm, extended_doftable, directory ) );
}

} // namespace FeelModels
} // namespace Feel
#endif // FEELPP_MODELS_MESHALE_H

/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 2011-07-17

 Copyright (C) 2011 Université Joseph Fourier (Grenoble I)
 Copyright (C) 2011-present Feel++ Consortium

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
    using size_type = uint32_type;
    typedef Backend<double> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;

    /*
     * Moving mesh typedefs
     */
    typedef Convex convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    /*
     * Reference mesh typedefs
     */
    typedef typename mesh_type::P1_mesh_type mesh_ref_type;
    typedef std::shared_ptr<mesh_ref_type> mesh_ref_ptrtype;
    typedef typename mesh_ref_type::shape_type convex_ref_type;

    typedef ALE<convex_ref_type, mesh_type::nOrder> ale_map_type;
    typedef std::shared_ptr<ale_map_type> ale_map_ptrtype;

    typedef typename ale_map_type::ale_map_functionspace_type ale_map_functionspace_ref_type;
    typedef std::shared_ptr<ale_map_functionspace_ref_type> ale_map_functionspace_ref_ptrtype;

    typedef typename ale_map_type::ale_map_element_type ale_map_element_ref_type;
    typedef std::shared_ptr<ale_map_element_ref_type> ale_map_element_ref_ptrtype;

    typedef typename ale_map_type::ale_map_basis_type ale_map_basis_type;
    typedef FunctionSpace<mesh_type, ale_map_basis_type> ale_map_functionspace_type;
    typedef std::shared_ptr<ale_map_functionspace_type> ale_map_functionspace_ptrtype;
    typedef typename ale_map_functionspace_type::element_type ale_map_element_type;
    typedef std::shared_ptr<ale_map_element_type> ale_map_element_ptrtype;

    typedef Bdf<ale_map_functionspace_ref_type> bdf_ale_displacement_ref_type;
    typedef std::shared_ptr<bdf_ale_displacement_ref_type> bdf_ale_displacement_ref_ptrtype;

    typedef Bdf<ale_map_functionspace_type> bdf_ale_displacement_type;
    typedef std::shared_ptr<bdf_ale_displacement_type> bdf_ale_displacement_ptrtype;

    typedef Exporter<mesh_type, mesh_type::nOrder> exporter_type;
    typedef std::shared_ptr<exporter_type> exporter_ptrtype;

    typedef Exporter<mesh_ref_type, mesh_ref_type::nOrder> exporter_ref_type;
    typedef std::shared_ptr<exporter_ref_type> exporter_ref_ptrtype;

    typedef DofRelationshipMap<ale_map_functionspace_ref_type, ale_map_functionspace_type> DofRelationshipMap_type;
    typedef std::shared_ptr<DofRelationshipMap_type> DofRelationshipMap_ptrtype;

    MeshALE( mesh_ptrtype mesh_moving,
             std::string const& prefix = "",
             worldcomm_ptr_t const& worldcomm = Environment::worldCommPtr(),
             ModelBaseRepository const& modelRep = ModelBaseRepository() );

    void init();

    std::shared_ptr<std::ostringstream> getInfo() const;

    void addBoundaryFlags( std::string const& bctype, std::string const& marker );
    void addBoundaryFlags( std::string const& bctype, std::vector<std::string> const& markers )
    {
        for ( std::string const& marker : markers )
            this->addBoundaryFlags( bctype, marker );
    }

    /**
     * \return the reference mesh
     */
    mesh_ref_ptrtype referenceMesh() const { return M_referenceMesh; }

    /**
     * \return the moving mesh
     */
    mesh_ptrtype movingMesh() const { return M_movingMesh; }
    mesh_ptrtype movingMeshModify() & { return M_movingMesh; }
    /**
     * \return true is reverted on reference mesh, else false
     */
    bool isOnReferenceMesh() const { return M_isOnReferenceMesh; }
    bool isOnMovingMesh() const { return M_isOnMovingMesh; }
    void switchReferenceMovingMesh()
    {
        if ( M_isOnMovingMesh )
        {
            M_isOnReferenceMesh = true;
            M_isOnMovingMesh = false;
            return;
        }
        if ( M_isOnReferenceMesh )
        {
            M_isOnReferenceMesh = false;
            M_isOnMovingMesh = true;
            return;
        }
    }
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
    template <typename elem_type>
    void update( std::vector<elem_type> const& polyDisplacementSet );
    template <typename elem_type>
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
    void exportResults( double time = 0 );

    void updateMetricMeshAdaptation( Expr<GinacExVF<2>> const& e )
    {
        M_aleFactory->updateMetricMeshAdaptation( e );
    }
    template <typename ExprT>
    void updateMetricMeshAdaptation( Expr<ExprT> const& e )
    {
        M_aleFactory->updateMetricMeshAdaptation( e );
    }

    template <typename FieldDispType, typename ExprT, typename RangeType>
    void updateDisplacementFieldFromVelocity( FieldDispType& d, Expr<ExprT> const& v, RangeType const& range ) const
    {
        unwrap_ptr( d ).on( _range = range,
                            _expr = ( v + idv( M_bdf_ale_identity->polyDeriv() ) ) / M_bdf_ale_identity->polyDerivCoefficient( 0 ) - idv( this->identityALE() ) + idv( this->displacement() ) );
    }

    // make the new mesh the new reference for the fields -> disp=0 on it
    void remeshingRevertBDF( bdf_ale_displacement_ptrtype const& old_meshALE_bdf_velocity, bdf_ale_displacement_ref_ptrtype const& old_meshALE_bdf_displacement_ref,
                             bdf_ale_displacement_ptrtype const& old_meshALE_bdf_identity )
    {
        for ( int i = old_meshALE_bdf_displacement_ref->timeOrder() + 1; i >= 0; --i )
        //for ( int  i = 0;i < old_meshALE_bdf_displacement_ref->timeOrder();++i )
        {
            old_meshALE_bdf_displacement_ref->unknown( i ).printMatlab( "bdf_disp_ref_" + std::to_string( i ) + ".m" );
            auto temp = old_meshALE_bdf_displacement_ref->functionSpace()->elementPtr();
            temp->zero();
            temp->add( 1.0, old_meshALE_bdf_displacement_ref->unknown( i ) );
            temp->add( -1.0, old_meshALE_bdf_displacement_ref->unknown( 0 ) );
            old_meshALE_bdf_displacement_ref->setUnknown( i, *temp );
            old_meshALE_bdf_displacement_ref->unknown( i ).printMatlab( "bdf_disp_ref_" + std::to_string( i ) + "_apres_soustraction.m" );
        }
#if 0
        for ( int  i = old_meshALE_bdf_velocity->timeOrder()+1;i >= 0; --i )
        //for ( int  i = 0;i < old_meshALE_bdf_velocity->timeOrder();++i )
        {
            old_meshALE_bdf_velocity->unknown(i).printMatlab("bdf_vel_"+std::to_string(i)+".m");
            auto temp = old_meshALE_bdf_velocity->functionSpace()->elementPtr();
            temp->zero();
            temp->add(1.0,old_meshALE_bdf_velocity->unknown(i));
            temp->add(-1.0,old_meshALE_bdf_velocity->unknown(0));
            
            old_meshALE_bdf_velocity->setUnknown(i,*temp);
            old_meshALE_bdf_velocity->unknown(i).printMatlab("bdf_vel_"+std::to_string(i)+"_apres_soustraction.m");
        }
        for ( int  i = old_meshALE_bdf_identity->timeOrder()+1;i >= 0; --i )
        //for ( int  i = 0;i < old_meshALE_bdf_identity->timeOrder();++i )
        {
            old_meshALE_bdf_identity->unknown(i).printMatlab("bdf_identity_"+std::to_string(i)+".m");
            auto ind_max=old_meshALE_bdf_identity->timeOrder()-1;
            auto temp = old_meshALE_bdf_identity->functionSpace()->elementPtr();
            temp->zero();
            temp->add(1.0,old_meshALE_bdf_identity->unknown(i));
            temp->add(-1.0,old_meshALE_bdf_identity->unknown(0));
            old_meshALE_bdf_identity->setUnknown(i,*temp);
            old_meshALE_bdf_velocity->unknown(i).printMatlab("bdf_identity_"+std::to_string(i)+"_apres_soustraction.m");
            
        }

#endif
#if 0        
        //M_bdf_ale_displacement_ref = unwrap_ptr(old_meshALE_bdf_displacement_ref).deepCopy();
        M_bdf_ale_displacement_ref = bdf_ale_displacement_ref_ptrtype( new bdf_ale_displacement_ref_type( *old_meshALE_bdf_displacement_ref ) );
        for ( auto it = M_bdf_ale_displacement_ref->unknowns().begin(), en = M_bdf_ale_displacement_ref->unknowns().end(); it != en; ++ it )
        {
            *it = ale_map_element_ref_ptrtype( new ale_map_element_ref_type( this->M_bdf_ale_displacement_ref->functionSpace() ) );
        }

        M_bdf_ale_velocity=unwrap_ptr(old_meshALE_bdf_velocity).deepCopy();
        M_bdf_ale_identity=unwrap_ptr(old_meshALE_bdf_identity).deepCopy();   

        for(auto itdisp = old_meshALE_bdf_displacement_ref->unknowns().begin(), endisp = old_meshALE_bdf_displacement_ref->unknowns().end(); itdisp != endisp; ++ itdisp)
           {
               //Feel::cout << itdisp << std::endl;
               /*auto element_ale = this->M_bdf_ale_displacement_ref->functionSpace()->element();
               element_ale.add(1.0,unwrap_ptr(*itdisp));
               element_ale.add(-1.0,unwrap_ptr(*endisp));*/
           }
           for(auto itdisp = old_meshALE_bdf_velocity->unknowns().begin(), endisp = old_meshALE_bdf_velocity->unknowns().end(); itdisp != endisp; ++ itdisp)
           {
               //Feel::cout << itdisp << std::endl;
               /*auto element_ale = this->M_bdf_ale_displacement_ref->functionSpace()->element();
               element_ale.add(1.0,unwrap_ptr(*itdisp));
               element_ale.add(-1.0,unwrap_ptr(*endisp));*/
           }
           for(auto itdisp = old_meshALE_bdf_identity->unknowns().begin(), endisp = old_meshALE_bdf_identity->unknowns().end(); itdisp != endisp; ++ itdisp)
           {
               //Feel::cout << itdisp << std::endl;
               /*auto element_ale = this->M_bdf_ale_displacement_ref->functionSpace()->element();
               element_ale.add(1.0,unwrap_ptr(*itdisp));
               element_ale.add(-1.0,unwrap_ptr(*endisp));*/
           }

        M_bdf_ale_displacement_ref = unwrap_ptr(old_meshALE_bdf_displacement_ref).deepCopy();
        M_bdf_ale_velocity=unwrap_ptr(old_meshALE_bdf_velocity).deepCopy();
        M_bdf_ale_identity=unwrap_ptr(old_meshALE_bdf_identity).deepCopy();

        std::copy(old_meshALE_bdf_velocity->unknowns().begin(), old_meshALE_bdf_velocity->unknowns().end(), M_bdf_ale_velocity->unknowns().begin());
        std::copy(old_meshALE_bdf_displacement_ref->unknowns().begin(), old_meshALE_bdf_displacement_ref->unknowns().end(), M_bdf_ale_displacement_ref->unknowns().begin());
        
        
        
        for (int itvel =0 ; itvel<old_meshALE_bdf_displacement_ref->timeValues().size();//, envel = this->M_bdf_ale_velocity->unknowns().end(); 
        //itvel!=envel;
        ++itvel )
        {
            Feel::cout << itvel << std::endl;
            auto end_index =old_meshALE_bdf_displacement_ref->timeValues().size();
            auto current_unknown = unwrap_ptr(old_meshALE_bdf_displacement_ref->unknown(itvel));
            auto final_unknown = unwrap_ptr(old_meshALE_bdf_displacement_ref->unknown(end_index-1));
            current_unknown.add(-1.0,final_unknown);
            this->M_bdf_ale_displacement_ref->unknown(itvel) = (current_unknown);
            //unwrap_ptr(*itvel).add(-1.0,unwrap_ptr(*envel));
        }
        for (int itvel =0 ; itvel<old_meshALE_bdf_velocity->timeValues().size();//, envel = this->M_bdf_ale_velocity->unknowns().end(); 
        //itvel!=envel;
        ++itvel )
        {
            Feel::cout << itvel << std::endl;
            auto end_index =old_meshALE_bdf_velocity->timeValues().size();Feel::cout << end_index << std::endl;
            auto current_unknown = unwrap_ptr(old_meshALE_bdf_velocity->unknown(itvel));Feel::cout << itvel << std::endl;
            auto final_unknown = unwrap_ptr(old_meshALE_bdf_velocity->unknown(itvel));Feel::cout << itvel << std::endl;
            current_unknown.add(-1.0,final_unknown);Feel::cout << itvel << std::endl;
            this->M_bdf_ale_velocity->unknown(itvel) = (current_unknown);
            //unwrap_ptr(*itvel).add(-1.0,unwrap_ptr(*envel));
        }
        /*for(auto itdisp = 0; itdisp<M_bdf_ale_displacement_ref->->timeValues().size(); ++ itdisp)
        //(auto itdisp = M_bdf_ale_displacement_ref->unknowns().begin(), endisp = M_bdf_ale_displacement_ref->unknowns().end(); itdisp != endisp; ++ itdisp)
            unwrap_ptr(*itdisp).add(-1.0,unwrap_ptr(*endisp));*/
#endif
    }

  private:
    void updateIdentityMap();
    void initTimeStep();
    void updateImpl();

  private:
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

  public:
    bdf_ale_displacement_ref_ptrtype displacementRefBDF() { return M_bdf_ale_displacement_ref; }
    bdf_ale_displacement_ref_ptrtype const& displacementRefBDF() const { return M_bdf_ale_displacement_ref; }
    bdf_ale_displacement_ptrtype aleIdentityBDF() { return M_bdf_ale_identity; }
    bdf_ale_displacement_ptrtype const& aleIdentityBDF() const { return M_bdf_ale_identity; }
    bdf_ale_displacement_ptrtype aleVelocityBDF() { return M_bdf_ale_velocity; }
    bdf_ale_displacement_ptrtype const& aleVelocityBDF() const { return M_bdf_ale_velocity; }
    exporter_ptrtype& exporterALE() { return M_exporter; };
    exporter_ref_ptrtype& exporterRef() { return M_exporter_ref; };
    mesh_ref_ptrtype& modifyReferenceMesh() { return M_referenceMesh; }

  private:
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

template <class Convex>
template <typename elem_type>
void MeshALE<Convex>::update( std::vector<elem_type> const& polyDisplacementSet )
{
    CHECK( polyDisplacementSet.size() == 1 ) << "invalid size";
    this->update( polyDisplacementSet[0] );
}
template <class Convex>
template <typename elem_type>
void MeshALE<Convex>::update( elem_type const& polyDisplacementSet )
{
    this->log( prefixvm( this->prefix(), "MeshALE" ), "update", "start" );

    CHECK( this->isOnMovingMesh() ) << "meshALE must be on moving mesh\n";

    M_displacementOnMovingBoundary_HO_ref->on( _range = markedfaces( M_movingMesh, this->aleFactory()->flagSet( "moving" ) ),
                                               _expr = vf::idv( polyDisplacementSet ) );
    sync( *M_displacementOnMovingBoundary_HO_ref, "=", M_dofsMultiProcessOnMovingBoundary_HO );

    this->updateImpl();

    this->log( prefixvm( this->prefix(), "MeshALE" ), "update", "finish" );
}

//------------------------------------------------------------------------------------------------//

template <typename Args>
struct compute_meshale_return
{
    typedef typename boost::remove_reference<typename parameter::binding<Args, tag::mesh>::type>::type::element_type mesh_type;
    using index_type = typename mesh_type::index_type;
    typedef typename mesh_type::shape_type convex_type;
    typedef MeshALE<convex_type> type;
    typedef std::shared_ptr<type> ptrtype;
};

// clang-format off
BOOST_PARAMETER_FUNCTION(
    ( typename compute_meshale_return<Args>::ptrtype ),                                                                                                                             // 1. return type
    meshale,                                                                                                                                                                        // 2. name of the function template
    tag,                                                                                                                                                                            // 3. namespace of tag types
    ( required
        ( mesh, *(boost::is_convertible<mpl::_, std::shared_ptr<MeshBase<>>>) ) 
    ) // required                                                                                             
    ( optional
        ( prefix, ( std::string ), std::string( "" ) )
        ( worldcomm, ( worldcomm_ptr_t ), mesh->worldCommPtr() )
        ( directory, ( ModelBaseRepository ), ModelBaseRepository() ) 
    ) // optionnal
)
// clang-format on
{
    typedef typename compute_meshale_return<Args>::ptrtype meshale_ptrtype;
    typedef typename compute_meshale_return<Args>::type meshale_type;
    return meshale_ptrtype( new meshale_type( mesh, prefix, worldcomm, directory ) );
}

} // namespace FeelModels
} // namespace Feel
#endif // FEELPP_MODELS_MESHALE_H

/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 This file is part of the Feel library

 Copyright (C) 2010 Université de Grenoble 1

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
 * \file doftable.hpp
 * \author Christophe Prud'homme
 */
#ifndef FEELPP_DOFTABLE_HH
#define FEELPP_DOFTABLE_HH

#include <tuple>
namespace std
{
namespace
{

// Code from boost
// Reciprocal of the golden ratio helps spread entropy
//     and handles duplicates.
// See Mike Seymour in magic-numbers-in-boosthash-combine:
//     http://stackoverflow.com/questions/4948780

template <class T>
inline void hash_combine(std::size_t& seed, T const& v)
{
    seed ^= hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

// Recursive template code derived from Matthieu M.
template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
    struct HashValueImpl
    {
        static void apply(size_t& seed, Tuple const& tuple)
            {
                HashValueImpl<Tuple, Index-1>::apply(seed, tuple);
                hash_combine(seed, get<Index>(tuple));
            }
    };

template <class Tuple>
struct HashValueImpl<Tuple,0>
{
    static void apply(size_t& seed, Tuple const& tuple)
        {
            hash_combine(seed, get<0>(tuple));
        }
};
}

template <typename ... TT>
struct hash<std::tuple<TT...>>
{
    size_t
        operator()(std::tuple<TT...> const& tt) const
    {
        size_t seed = 0;
        HashValueImpl<std::tuple<TT...> >::apply(seed, tt);
        return seed;
    }

};
}


#include <set>
#include <map>
//#include <boost/functional/hash.hpp>
//#include <boost/functional/hash/extensions.hpp>
#include <unordered_map>
#include <boost/unordered_map.hpp>
#include <vector>
#include <algorithm>
#include <cmath>


#include <feel/feelcore/boostmultiarray.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/fusion/algorithm/iteration/accumulate.hpp>
#include <boost/bimap.hpp>
#include <boost/bimap/support/lambda.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/set_of.hpp>

#include <Eigen/Core>
#include<Eigen/StdVector>

#include <feel/feelcore/feel.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelpoly/mapped.hpp>
#include <feel/feelpoly/isp0continuous.hpp>
#include <feel/feelpoly/order.hpp>
#include <feel/feelpoly/hdivpolynomialset.hpp>
#include <feel/feelpoly/hcurlpolynomialset.hpp>
#include <feel/feeldiscr/doftablebase.hpp>
#include <feel/feeldiscr/doffromelement.hpp>
#include <feel/feeldiscr/doffrommortar.hpp>
#include <feel/feeldiscr/doffromboundary.hpp>
#include <feel/feeldiscr/doffromedge.hpp>

#include <feel/feelmesh/meshsupport.hpp>

namespace Feel
{
template<class ITERATOR>
ITERATOR begin( std::pair<ITERATOR,ITERATOR> &range )
{
    return range.first;
}

template<class ITERATOR>
ITERATOR end( std::pair<ITERATOR,ITERATOR> &range )
{
    return range.second;
}

// import fusion namespace in Feel
namespace fusion = boost::fusion;
namespace bimaps = boost::bimaps;
/**
 * \class DofTable
 * \ingroup SpaceTime
 * \brief Local-to-global Degree of Freedom table
 *
 * \author Christophe Prud'homme
 * \author Goncalo Pena
 */
template<typename MeshType, typename FEType, typename MortarType>
class DofTable : public DofTableBase<typename MeshType::size_type>
{
    typedef DofTableBase<typename MeshType::size_type> super;
public:

    /**
     * mesh type
     */
    typedef MeshType mesh_type;
    typedef FEType fe_type;
    using self_type = DofTable<MeshType, FEType, MortarType>;
    using doftable_type = self_type;
    using size_type = typename mesh_type::size_type;
    typedef std::shared_ptr<FEType> fe_ptrtype;
    typedef MortarType mortar_type;
    static inline const bool is_mortar = mortar_type::is_mortar;
    typedef typename fe_type::SSpace::type mortar_fe_type;

    typedef MeshSupport<mesh_type> mesh_support_type;
    typedef std::shared_ptr<mesh_support_type> mesh_support_ptrtype;
    typedef typename super::mesh_support_base_ptrtype mesh_support_base_ptrtype;

    typedef typename mesh_type::element_const_iterator element_const_iterator;
    typedef typename mesh_type::element_type element_type;
    typedef typename mesh_type::face_type face_type;
    typedef typename mesh_type::edge_type edge_type;
    typedef typename mesh_type::gm_ptrtype gm_ptrtype;
    typedef typename mesh_type::gm_type gm_type;
    using mesh_marker_type = typename element_type::marker_type;

    typedef typename fe_type::matrix_type matrix_type;
    typedef typename fe_type::value_type value_type;
    typedef typename fe_type::reference_convex_type reference_convex_type;
    typedef typename fe_type::points_type points_type;
    //typedef ContinuityType continuity_type;
    typedef typename fe_type::continuity_type continuity_type;



    typedef typename reference_convex_type::super convex_type;

    typedef typename element_type::edge_permutation_type edge_permutation_type;
    typedef typename element_type::face_permutation_type face_permutation_type;

    using dof_from_edge_type = DofFromEdge<doftable_type,fe_type>;

    static inline const uint16_type nOrder = fe_type::nOrder;
    static inline const uint16_type nDim = mesh_type::nDim;
    static inline const uint16_type nRealDim = mesh_type::nRealDim;
    static inline const uint16_type Shape = mesh_type::Shape;
    static inline const uint16_type nComponents = fe_type::nComponents;
    static inline const uint16_type nComponents1 = fe_type::nComponents1;
    static inline const uint16_type nComponents2 = fe_type::nComponents2;


    static inline const bool is_continuous = fe_type::isContinuous;
    static inline const bool is_discontinuous_locally = fe_type::continuity_type::is_discontinuous_locally;
    static inline const bool is_discontinuous_totally = fe_type::continuity_type::is_discontinuous_totally;

    static inline const bool is_scalar = FEType::is_scalar;
    static inline const bool is_vectorial = FEType::is_vectorial;
    static inline const bool is_tensor2 = FEType::is_tensor2;
    static inline const bool is_tensor2symm = FEType::is_tensor2 && is_symm_v<FEType>;
    static inline const bool is_modal = FEType::is_modal;
    static inline const bool is_product = FEType::is_product;
    static inline const uint16_type nRealComponents = is_tensor2symm?(fe_type::nComponents1*(fe_type::nComponents1+1)/2):fe_type::nComponents;

    static inline const bool is_p0_continuous = ( ( nOrder == 0 ) && is_continuous );

    static inline const bool is_hdiv_conforming = Feel::is_hdiv_conforming<fe_type>::value;
    static inline const bool is_hcurl_conforming = Feel::is_hcurl_conforming<fe_type>::value;

    static inline const uint16_type nDofPerEdge = fe_type::nDofPerEdge;
    static inline const uint16_type nDofPerElement = mpl::if_<mpl::bool_<is_product>, mpl::int_<FEType::nLocalDof*nComponents>, mpl::int_<FEType::nLocalDof> >::type::value;

    //! @brief True if polynomial order is determined at runtime
    static constexpr bool is_order_dynamic = orderIsDynamic<fe_type>;

    static inline const bool is_periodic = false;

    static constexpr uint16_type nDofComponents() { return is_product?nComponents:1; }



    /**
     * A global dof is defined by its index in the global table
     * and by its sign.
     **/

    //typedef boost::tuple<size_type, int16_type, bool> global_dof_type;
    typedef Dof<size_type> global_dof_type;
    using globaldof_type = global_dof_type;

    //! point id to dof id relation type
    using pidtodofid_type = std::pair<std::unordered_map<size_type,size_type>,std::unordered_map<size_type,size_type> >;

    /**
     * A global dof from face is defined by
     * -its index in the global table
     * -its sign
     * -bool : periodicity
     * -local number in the element
     **/

    //typedef boost::tuple<size_type, int16_type, bool, int16_type> global_dof_fromface_type;
    typedef FaceDof<size_type> global_dof_fromface_type;
    using global_dof_from_entity_type = EntityDof<size_type>;

    /**
     * Type for the localToGlobal table.
     */
    //typedef std::unordered_map<int,std::map<int,global_dof_type> > Container;
    //typedef typename std::map<int,global_dof_type>::iterator local_map_iterator;
    typedef LocalDof<nDofComponents(),size_type> localdof_type;
    typedef boost::bimap<bimaps::set_of<localdof_type>, bimaps::multiset_of<globaldof_type> > dof_table;
    typedef typename dof_table::value_type dof_relation;
    typedef std::unordered_map<int,std::vector<global_dof_fromface_type> > Container_fromface;
    typedef typename std::vector<global_dof_fromface_type>::const_iterator face_local_dof_const_iterator;
    typedef typename dof_table::left_iterator local_dof_iterator;
    typedef typename dof_table::left_const_iterator local_dof_const_iterator;
    typedef typename dof_table::right_iterator global_dof_iterator;
    typedef typename dof_table::right_const_iterator global_dof_const_iterator;
    typedef typename std::map<int,global_dof_type> indices_per_element_type;

    typedef typename node<value_type>::type node_type;


    typedef boost::tuple<node_type, size_type, uint16_type > dof_point_type;
    typedef std::unordered_map<size_type,dof_point_type> dof_points_type;
    typedef typename std::unordered_map<size_type,dof_point_type>::iterator dof_points_iterator;
    typedef typename std::unordered_map<size_type,dof_point_type>::const_iterator dof_points_const_iterator;

    /**
     * Tuple that holds a size_type \p elt 1 uint16_type \p l and 1
     * uint16_type ent
     * \p elt shall be an element index
     * \p l shall be the local index of the dof in the element
     * \p ent shall be the entity the dof belongs to (0: vertex, 1: edge, 2: face, 3: volume)
     */
    typedef boost::tuple<size_type, uint16_type, uint16_type, uint16_type> local_dof_type;
    typedef LocalDofSet<nDofComponents(),size_type> local_dof_set_type;

    typedef std::tuple<uint16_type&,size_type&> ref_shift_type;

    /**
     * Type that hold the map between a global dof and the elements
     */
    typedef std::map<size_type, std::list<local_dof_type> >  dof_element_type;

    typedef boost::bimap<size_type,boost::bimaps::multiset_of<size_type> > dof_marker_type;
    typedef typename dof_marker_type::value_type dof2marker;

    typedef typename dof_element_type::iterator dof_iterator;
    typedef typename dof_element_type::const_iterator dof_const_iterator;

    typedef typename std::list<local_dof_type>::const_iterator ldof_const_iterator;

    // unique dof description : fist entity type (0,1,2,3: vertex, edge, face,
    // volume), then and dof id associated to the entity that is unique with
    // respect to the entity
    typedef std::tuple<uint16_type,size_type> dof_type;
    typedef std::unordered_map<dof_type, size_type> dof_map_type;
    typedef typename dof_map_type::iterator dof_map_iterator;
    typedef typename dof_map_type::const_iterator dof_map_const_iterator;
    struct global_dof_key_entry
    {
        dof_type key;
        size_type localDofIndex = invalid_v<size_type>;
        uint16_type component = 0;
    };

    typedef std::map<size_type, std::set<size_type> > dof_procset_type;
    /**
     * This type is useful to construct the sign map in the modal case
     **/

    typedef ublas::vector<bool> face_sign_info_type;


    //typedef typename mpl::if_<is_mortar,
    //mpl::identity<Eigen::Matrix<int, Eigen::Dynamic, 1> >,
    //mpl::identity<Eigen::Matrix<int, nDofPerElement, 1> > >::type::type localglobal_indices_type;
    typedef Eigen::Matrix<int, Eigen::Dynamic, 1>  localglobal_indices_type;
    using localglobal_transforms_type = std::vector<DofTransform>;

    /**
     * Type for the permutations to be done in the faces
     **/

    typedef ublas::vector<uint16_type> permutation_vector_type;

    //typedef typename std::vector<localglobal_indices_type,Eigen::aligned_allocator<localglobal_indices_type> > vector_indices_type;
    using vector_indices_type = std::unordered_map<size_type,localglobal_indices_type,
                                        std::hash<size_type>,std::equal_to<size_type>,
                                        Eigen::aligned_allocator<std::pair<const size_type,localglobal_indices_type > > >;
    using vector_transforms_type = std::unordered_map<size_type,localglobal_transforms_type>;

    DofTable( WorldComm const& _worldComm )
        :
        super( _worldComm )
        {}

    /**
     * @brief The minimal constructor
     *
     * @param _fe reference element
     *
     */
    DofTable( fe_ptrtype const& _fe, WorldComm const& _worldComm );

    /**
     * copy constructor
     *
     * @param dof2 a dof object instance
     */
    DofTable( const DofTable & dof2 );

    /**
     * @brief  Constructor accepting a mesh as parameter
     *
     *  @param mesh a RegionMesh3D
     *  @param _fe reference element
     */
    DofTable( mesh_type& mesh, fe_ptrtype const& _fe, WorldComm const& _worldComm );

    ~DofTable() override
        {
            M_el_l2g.clear();
            M_face_l2g.clear();
            M_dof_points.clear();
        }
    fe_type const& fe() const { return *M_fe; }

    size_type nRealLocalDof( bool per_component = false ) const
        {
            const auto localDofPerComponent = this->feLocalDofCount( true );
            if ( per_component )
                return localDofPerComponent;

            if constexpr ( is_tensor2symm )
                return nRealComponents * localDofPerComponent;
            else
                return this->feLocalDofCount();
        }

    DofTableInfos infos() const override
        {
            DofTableInfos infos;
            infos.nOrder = nOrder;
            infos.nDim = nDim;
            infos.nRealDim = nRealDim;
            infos.Shape = Shape;
            infos.nComponents = nComponents;
            infos.nComponents1 = nComponents1;
            infos.nComponents2 = nComponents2;
            infos.is_continuous = is_continuous;
            infos.is_discontinuous_locally = is_discontinuous_locally;
            infos.is_discontinuous_totally = is_discontinuous_totally;

            infos.is_scalar = is_scalar;
            infos.is_vectorial = is_vectorial;
            infos.is_tensor2 = is_tensor2;
            infos.is_tensor2symm = is_tensor2symm;
            infos.is_modal = is_modal;
            infos.is_product = is_product;
            infos.nRealComponents = nRealComponents;

            infos.is_p0_continuous = is_p0_continuous;

            infos.is_hdiv_conforming = is_hdiv_conforming;
            infos.is_hcurl_conforming = is_hcurl_conforming;

            infos.nDofPerEdge = nDofPerEdge;
            infos.nDofPerElement = nDofPerElement;

            infos.is_periodic = is_periodic;

            infos.nDofComponents = this->nDofComponents();

            if ( M_fe )
                infos.feFamilyName = M_fe->familyName();

            return infos;
        }

    size_type nLocalDof( bool per_component = false ) const
        {
            return this->feLocalDofCount( per_component );
        }
    size_type nLocalDofOnFace( bool per_component = false ) const
        {
            return this->feLocalDofCountOnFacet( 0, per_component );
        }
    size_type nLocalDofOnFacet( bool per_component = false ) const
        {
            return this->nLocalDofOnFace( per_component );
        }
    local_dof_set_type const&
    localDofSet( size_type eid ) const
        {
            if ( is_mortar )
                return M_local_dof_set.update( eid, getIndicesSize( eid ) );
            else
                return M_local_dof_set.update( eid );
        }
    mesh_type* mesh() { return M_mesh; }
    mesh_type* mesh() const { return M_mesh; }

    /**
     * set a mesh support where the doftable is built
     */
    void setMeshSupport( mesh_support_ptrtype const& meshSupport )
        {
            M_meshSupport = meshSupport;
        }
    /**
     * \return mesh support
     */
    mesh_support_ptrtype meshSupport() const
        {
            return M_meshSupport;
        }

    mesh_support_base_ptrtype meshSupportBase() const override
        {
            return M_meshSupport;
        }
    /**
     * \return true if a mesh support object is defined
     */
    bool hasMeshSupport() const
        {
            if ( M_meshSupport )
                return true;
            else
                return false;
        }

    /**
     * \return the number of dof for faces on the boundary
     */
    uint16_type nDofPerFaceOnBoundary() const
        {
            return M_n_dof_per_face_on_bdy;
        }

    indices_per_element_type  indices( size_type id_el ) const
        {
            indices_per_element_type ind;
            for( localdof_type const& ldof : localDofSet( id_el ) )
            {
                auto it = M_el_l2g.left.find( ldof );
                DCHECK( it != M_el_l2g.left.end() ) << "Invalid element id " << id_el;
                ind[ldof.localDof()] = *it;
            }
            return ind;
        }

    size_type getIndicesSize( int eid = 0 ) const
        {
            return getIndicesSize( eid, mpl::bool_<is_mortar>() );
        }
    size_type getIndicesSize( int eid, mpl::true_ ) const
        {
            auto itrange = localDof( eid );
            return std::distance( itrange.first, itrange.second );
        }
    size_type getIndicesSize( int eid, mpl::false_ ) const
        {
            return nLocalDof();
        }
    std::vector<size_type> getIndices( size_type id_el ) const
        {
            std::vector<size_type> ind( getIndicesSize(id_el) );
            getIndicesSet( id_el, ind );

            return ind;
        }

    std::vector<size_type> getIndices( size_type id_el, mpl::size_t<MESH_ELEMENTS> /**/ ) const
        {
            return getIndices( id_el );
        }

    void getIndicesSet( size_type id_el, std::vector<size_type>& ind ) const
        {
#if 0
            for( localdof_type const& ldof : this->localDofSet( id_el ) )
            {
                auto it = M_el_l2g.left.find( ldof );
                DCHECK(it != M_el_l2g.left.end() ) << "Invalid element id " << id_el;
                ind[ldof.localDof()] = it->second.index();
            }
#else
            for( auto const& ldof : this->localDof( id_el ) )
            {
                ind[ldof.first.localDof()] = ldof.second.index();
            }

#endif
        }

    std::vector<size_type> getIndices( size_type id_el, mpl::size_t<MESH_FACES> /**/ ) const
        {
            std::vector<size_type> ind;
            ind.reserve( nLocalDofOnFace() );

            auto eit = M_face_l2g.find( id_el );
            DCHECK( eit != M_face_l2g.end() ) << "Invalid face id " << id_el;
            std::for_each( eit->second.begin(), eit->second.end(),
                           [&ind]( FaceDof<size_type> const& f ) { ind.push_back( f.index() ); } );
            return ind;
        }


    bool getIndicesSetOnGlobalCluster( size_type id_el, std::vector<size_type>& ind ) const
        {
            bool is_empty = false;
            for( localdof_type const& ldof: this->localDofSet( id_el ) )
            {
                auto it = M_el_l2g.left.find( ldof );
                DLOG_IF( WARNING, it == M_el_l2g.left.end() ) << "Invalid element id " << id_el;
                is_empty = is_empty || it == M_el_l2g.left.end();
                if ( it != M_el_l2g.left.end() )
                    ind[ldof.localDof()] =this->mapGlobalProcessToGlobalCluster()[ it->second.index() ];

            }
            return is_empty;
        }

    std::vector<size_type> getIndicesOnGlobalCluster( size_type id_el ) const
        {
            const size_type s = getIndicesSize( id_el );
            std::vector<size_type> ind(s);
            bool is_empty = getIndicesSetOnGlobalCluster( id_el, ind );
            if ( is_empty ) ind.clear();
            return ind;
        }

    /**
     * @return the coordinates of the nodal dofs associated with the
     * element \p el
     */
    const dof_point_type& dofPoint( size_type i ) const
        {
            if (!hasDofPoints()) this->generateDofPoints(*M_mesh);
            auto itFindDp = M_dof_points.find( i );
            CHECK( itFindDp != M_dof_points.end() )
                << "dof index " << i << " has no representative point in FE "
                << ( M_fe ? M_fe->familyName() : std::string( "<null>" ) )
                << ". This is expected for moment/modal finite elements; use FE functional metadata instead of dofPoint().";
            return itFindDp->second;
        }

    /**
     * \return true if dof index \p i has representative point coordinates.
     */
    bool hasDofPoint( size_type i ) const
        {
            if (!hasDofPoints()) this->generateDofPoints(*M_mesh);
            return M_dof_points.find( i ) != M_dof_points.end();
        }

    /**
     * @return the dof points data structure
     * it allows for example to do:
     * \code
     * for( auto const& pt: dofPoints())
     * {
     *   // do something on pt
     * }
     * \endcode
     */
    dof_points_type const& dofPoints() const
        {
            if (!hasDofPoints()) this->generateDofPoints(*M_mesh);
            return M_dof_points;
        }
    /**
     * @return an iterator at the beginning of dof points
     */
    dof_points_const_iterator dofPointBegin() const
        {
            if (!hasDofPoints()) this->generateDofPoints(*M_mesh);
            return M_dof_points.begin();
        }

    /**
     * @return an iterator at the beginning of dof points
     */
    dof_points_iterator dofPointBegin()
        {
            if (!hasDofPoints()) this->generateDofPoints(*M_mesh);
            return M_dof_points.begin();
        }

    /**
     * @return an iterator at the end of dof points
     */
    dof_points_const_iterator dofPointEnd() const
        {
            if (!hasDofPoints()) this->generateDofPoints(*M_mesh);
            return M_dof_points.end();
        }

    /**
     * @return an iterator at the end of dof points
     */
    dof_points_iterator dofPointEnd()
        {
            if (!hasDofPoints()) this->generateDofPoints(*M_mesh);
            return M_dof_points.end();
        }

    /**
     * insted of creating the dof indices on the fly, get them from a
     * vector. The situation typically arises when we want to have dof
     * correspondence between two spaces
     *
     * \see OperatorLagrangeP1
     */
#if 0
    void setDofIndices( std::vector<globaldof_type> const& dof )
        {
            M_dof_indices.resize( dof.size() );
            std::copy( dof.begin(), dof.end(), M_dof_indices.begin() );

            if ( dof.empty() )
                return ;

#if 1
            std::set<size_type> eltid;
            std::set<size_type> dofs;

            for( globaldof_type const& thedof:  dof )
            {
                eltid.insert( std::get<0>( thedof ) );
                dofs.insert( std::get<0>( thedof ) );
            }
#endif

            for( globaldof_type const& thedof:  dof )
            {
                M_el_l2g.insert( dof_relation( localdof_type( std::get<0>( thedof ), std::get<0>( thedof ) ),
                                               Dof( std::get<0>( thedof ), 0, false ) ) );

            }
            int processor = this->worldComm().localRank();

            this->M_first_df[processor] = 0;
            this->M_last_df[processor] = dofs.size()-1;
            this->M_n_dofs = dofs.size();

            this->M_n_localWithGhost_df[processor] = this->M_last_df[processor] - this->M_first_df[processor] + 1;
            this->M_n_localWithoutGhost_df[processor]=this->M_n_localWithGhost_df[processor];
            this->M_first_df_globalcluster[processor]=this->M_first_df[processor];
            this->M_last_df_globalcluster[processor]=this->M_last_df[processor];
        }
#endif

    /**
     * \return the dof index
     */
    size_type dofIndex( size_type dof ) const
        {
            return dof;
#if 0

            if ( M_dof_indices.empty() )
                return dof;

            FEELPP_ASSERT( dof < M_dof_indices.size() )( dof )( M_dof_indices.size() ).warn( "invalid dof index" );
            return M_dof_indices[dof];
#endif
        }

    /**
     * \return the local to global indices
     */
    vector_indices_type const& localToGlobalIndices() const
        {
            return M_locglob_indices;
        }

    /**
     * \return the local to global indices
     */
    localglobal_indices_type const& localToGlobalIndices( size_type ElId ) const
        {
            auto itFindElt = M_locglob_indices.find( ElId );
            DCHECK( itFindElt != M_locglob_indices.end() ) << "no locglob_indices in elt : " << ElId;
            return itFindElt->second;
        }

    /**
     * \return the local to global indices
     */
    localglobal_indices_type localToGlobalIndices( size_type ElId, std::vector<size_type> const& basisToContainerGlobalProcess ) const
        {
            auto const& basisIndices = this->localToGlobalIndices( ElId );
            int nLocalDof = basisIndices.size();
            localglobal_indices_type res = localglobal_indices_type::Zero( basisIndices.size() );
            for ( int j=0 ; j<nLocalDof ; ++j )
                res( j ) = basisToContainerGlobalProcess[ basisIndices(j) ];
            return res;
        }

    /**
     * \return the signs of the global dof (=1 in nodal case, +-1 in modal case)
     */
    localglobal_indices_type const& localToGlobalSigns( size_type ElId ) const
        {
            if ( is_hdiv_conforming || is_hcurl_conforming )
                return M_locglob_signs.find( ElId )->second;
            else
                return M_locglob_nosigns;
        }

    /**
     * \return transform descriptors of the local-to-global map.
     *
     * H(div)/H(curl) spaces store FE-owned orientation transforms here. The
     * legacy sign vector is the sign projection of this view.
     */
    localglobal_transforms_type const& localToGlobalTransforms( size_type ElId ) const
        {
            if constexpr ( is_hdiv_conforming || is_hcurl_conforming )
            {
                auto itFindElt = M_locglob_transforms.find( ElId );
                DCHECK( itFindElt != M_locglob_transforms.end() ) << "no locglob_transforms in elt : " << ElId;
                return itFindElt->second;
            }
            else
                return M_locglob_notransforms;
        }

    /**
     * \return transform descriptor for one local dof.
     */
    DofTransform const& localToGlobalTransform( size_type ElId, uint16_type localNode ) const
        {
            auto const& transforms = this->localToGlobalTransforms( ElId );
            DCHECK_LT( static_cast<std::size_t>( localNode ), transforms.size() )
                << "invalid transform local dof " << localNode << " in element " << ElId;
            return transforms[localNode];
        }

    [[nodiscard]] static int16_type dofTransformSignProjection( DofTransform const& transform ) noexcept
        {
            return transform.kind == DofTransformKind::Sign ? transform.sign : 1;
        }

    void setLocalToGlobalTransform( size_type ElId, uint16_type localDof, DofTransform transform )
        {
            if constexpr ( is_hdiv_conforming || is_hcurl_conforming )
            {
                auto itTransform = M_locglob_transforms.try_emplace( ElId, localglobal_transforms_type( runtimeNDofPerElement() ) ).first;
                if ( itTransform->second.size() < runtimeNDofPerElement() )
                    itTransform->second.resize( runtimeNDofPerElement() );
                DCHECK_LT( static_cast<std::size_t>( localDof ), itTransform->second.size() )
                    << "invalid transform local dof " << localDof << " in element " << ElId;
                itTransform->second[localDof] = std::move( transform );

                auto const sign = dofTransformSignProjection( itTransform->second[localDof] );
                auto signIt = M_locglob_signs.find( ElId );
                if ( signIt != M_locglob_signs.end() )
                {
                    DCHECK_LT( localDof, signIt->second.size() )
                        << "invalid sign local dof " << localDof << " in element " << ElId;
                    signIt->second[localDof] = sign;
                }
            }
        }

    void setLocalToGlobalSign( size_type ElId, uint16_type localDof, int16_type sign )
        {
            this->setLocalToGlobalTransform( ElId, localDof, DofTransform::signedOrientation( sign ) );
        }

    /**
     * \return the specified entries of the localToGlobal table
     *
     * \param ElId the element ID
     * \param localNode the local DOF numbering (starting from 1)
     * \param c the component index, default is 0-th component
     *
     * \return the global numbering of a DOF, given an element and the local numbering
     */
    size_type localToGlobalId( const size_type ElId,
                               const uint16_type id ) const
        {
            auto it = M_el_l2g.left.find( localdof_type(ElId, id ) );
            DCHECK( it != M_el_l2g.left.end() ) << "Invalid dof entry ( " << ElId << ", " << id << ")";
            DCHECK( it->second.index() < this->nDof() ) << "Invalid Dof Entry: " << it->second.index() << " > " << this->nDof();
            return it->second.index();
        }

    std::pair<global_dof_const_iterator,global_dof_const_iterator>  globalDof()  const
        {
            return std::make_pair( M_el_l2g.right.begin(), M_el_l2g.right.end() );
        }
    std::pair<global_dof_const_iterator,global_dof_const_iterator> globalDof( size_type GlobalDofId ) const
        {
            auto lower = M_el_l2g.right.lower_bound( globaldof_type(GlobalDofId,-1) );
            auto upper = M_el_l2g.right.upper_bound( globaldof_type(GlobalDofId,2) );
            return std::make_pair( lower, upper );
        }
    //!
    //! get the neighbor dofs
    //!
    std::set<globaldof_type> globalNeighbors( size_type GlobalDofId, bool add_self = true ) const
        {
            std::set<globaldof_type> neigh;
            auto const& [beg,end] = globalDof( GlobalDofId );
            //neigh.reserve( std::distance( beg, end ) * nLocalDof() );
            for( auto it = beg; it != end; ++ it )
            {
                for( auto const& [lid,gid]  : localDof( it->second.elementId() ) )
                    if ( gid.index() != GlobalDofId || add_self )
                        neigh.emplace( gid.index() );
            }

            return neigh;
        }
    /**
     * \return the specified entries of the globalToLocal table
     *
     * \param DofId the Dof ID
     *
     * \return the element id and local dof id
     */
    localdof_type const& globalToLocal( size_type dof )  const
        {
            auto it = M_el_l2g.right.find( Dof( dof ) );
            DCHECK( it != M_el_l2g.right.end() ) << "Invalid global dof entry ( " << dof << ")";
            return it->second;
        }

    uint16_type localDofId( uint16_type const lid, uint16_type const c = 0 ) const
        {
            return runtimeLocalDofId( lid, c );
        }

    std::pair<local_dof_const_iterator,local_dof_const_iterator> localDof() const
        {
            return std::make_pair( M_el_l2g.left.begin(), M_el_l2g.left.end() );
        }
    std::pair<local_dof_const_iterator,local_dof_const_iterator> localDof( size_type ElId ) const
        {
            auto lower = M_el_l2g.left.lower_bound( localdof_type(ElId) );
            auto upper = M_el_l2g.left.upper_bound( localdof_type(ElId,invalid_uint16_type_value) );
            //DCHECK( it.first != M_el_l2g.left.end() ) << "Invalid element dof entry " << ElId;
            return std::make_pair( lower, upper );
        }

    std::pair<face_local_dof_const_iterator,face_local_dof_const_iterator> faceLocalDof( size_type ElId ) const
        {
            auto it = M_face_l2g.find( ElId );
            if (  it == M_face_l2g.end() )
                return std::make_pair( face_local_dof_const_iterator(), face_local_dof_const_iterator() );
            auto be = it->second.begin();
            auto en = it->second.end();
            return std::make_pair( be, en );
        }
    std::pair<face_local_dof_const_iterator,face_local_dof_const_iterator> facetLocalDof( size_type ElId ) const
        {
            return this->faceLocalDof( ElId );
        }

    std::vector<global_dof_from_entity_type> edgeLocalDof( size_type elid, uint16_type edge_id ) const
        {
            return M_dfe( elid, edge_id );
        }

    template <typename MeshEntityType>
        requires std::is_same_v<MeshEntityType,typename mesh_type::element_type>
    auto localDof( MeshEntityType const& elt ) const { return this->localDof( elt.id() ); }
    template <typename MeshEntityType>
        requires std::is_same_v<MeshEntityType,typename mesh_type::face_type> && (mesh_type::nDim > 0)
    auto localDof( MeshEntityType const& face ) const { return this->faceLocalDof( face.id() ); }
#if 0
    template <typename MeshEntityType>
        requires std::is_same_v<MeshEntityType,typename mesh_type::edge_type> && (mesh_type::nDim == 3)
    auto localDof( MeshEntityType const& edge ) const { this->edgeLocalDof( edge.id() ); }
#endif


    template<typename ElemTest,typename ElemTrial>
    std::vector<uint16_type> const& localIndices( ElemTest const& eltTest, ElemTrial const& eltTrial  ) const
        {
            return localIndices( eltTest,eltTrial,mpl::bool_<ElemTest::nDim==1 && ElemTrial::nDim==1>() );
        }
    template<typename ElemTest,typename ElemTrial>
    std::vector<uint16_type> const& localIndices( ElemTest const& eltTest, ElemTrial const& eltTrial, mpl::false_ ) const
        {
            return M_localIndicesIdentity;
        }
    template<typename ElemTest,typename ElemTrial>
    std::vector<uint16_type> const& localIndices( ElemTest const& eltTest, ElemTrial const& eltTrial, mpl::true_ ) const
        {
            auto const& gTest = eltTest.G();
            auto const& gTrial = eltTrial.G();
            Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> gTestMap( gTest.data().begin(), gTest.size1(), gTest.size2() );
            Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> gTrialMap( gTrial.data().begin(), gTrial.size1(), gTrial.size2() );
            double dotVec = ( gTestMap.col( 1 ) - gTestMap.col( 0 ) ).dot( gTrialMap.col( 1 ) - gTrialMap.col( 0 ) );
            CHECK( std::abs( dotVec ) > 1e-9 ) << " inner_prod is null " << dotVec << "\n";

            if ( dotVec > 0 ) // identity permutation
                return M_localIndicesIdentity;
            else // reverse permutation
                return  M_localIndicesPerm;
        }


    global_dof_type const& localToGlobal( const size_type ElId,
                                          const uint16_type localNode,
                                          const uint16_type c = 0 ) const override
        {
            const uint16_type localDof = runtimeLocalDofId( localNode, c );
            auto it = M_el_l2g.left.find( localdof_type(ElId, localDof ) );
            DCHECK( it != M_el_l2g.left.end() ) << "Invalid dof entry ( " << ElId << ", " << localDof << ")";
            //DCHECK( it->second.index() < nDof() && nDof() > 0 ) << "Invalid Dof Entry: " << it->second.index() << " > " << this->nDof();
            return it->second;
        }

    global_dof_type localToGlobalOnCluster( const size_type ElId,
                                            const uint16_type localNode,
                                            const uint16_type c = 0 ) const
        {
            Dof resloc = M_el_l2g.left.find( localdof_type( ElId, runtimeLocalDofId( localNode, c ) ) )->second;
            resloc.setIndex( this->mapGlobalProcessToGlobalCluster()[resloc.index()] );
            return resloc;
        }

  private:
    [[nodiscard]] uint16_type componentFromGlobalProcessDof( size_type localDofIndex ) const noexcept
        {
            if constexpr ( FiniteElementDofLayoutProvider<fe_type> )
            {
                if ( M_fe )
                {
                    auto it = M_el_l2g.right.find( globaldof_type( localDofIndex ) );
                    if ( it != M_el_l2g.right.end() )
                        return M_fe->component( it->second.localDof() );
                }
            }
            return 0;
        }

    [[nodiscard]] bool localDofHasRepresentativePoint( uint16_type localDofId ) const
        {
            return M_fe && M_fe->dofHasRepresentativePoint( localDofId );
        }

    [[nodiscard]] uint16_type localDofRepresentativePointIndex( uint16_type localDofId ) const
        {
            if ( M_fe )
                return M_fe->dofRepresentativePointIndex( localDofId );
            return localDofId;
        }

    [[nodiscard]] size_type legacyLocalDofPerComponent() const noexcept
        {
            if constexpr ( is_order_dynamic )
            {
                const uint16_type nDofPerVertex = runtimeDofPerVertex();
                const uint16_type nDofPerEdge = runtimeDofPerEdge();
                const uint16_type nDofPerFace = runtimeDofPerFace();
                const uint16_type nDofPerVolume = runtimeDofPerVolume();
                return nDofPerVolume * element_type::numVolumes +
                       nDofPerFace * element_type::numGeometricFaces +
                       nDofPerEdge * element_type::numEdges +
                       nDofPerVertex * element_type::numVertices;
            }
            else
            {
                return fe_type::nDofPerVolume * element_type::numVolumes +
                       fe_type::nDofPerFace * element_type::numGeometricFaces +
                       fe_type::nDofPerEdge * element_type::numEdges +
                       fe_type::nDofPerVertex * element_type::numVertices;
            }
        }

    [[nodiscard]] size_type legacyLocalDofCount( bool perComponent = false ) const noexcept
        {
            const size_type localDofPerComponent = this->legacyLocalDofPerComponent();
            if ( perComponent )
                return localDofPerComponent;
            if constexpr ( is_product )
                return nComponents * localDofPerComponent;
            else
                return localDofPerComponent;
        }

    [[nodiscard]] size_type legacyLocalDofCountOnFacet( bool perComponent = false ) const noexcept
        {
            size_type localDofOnFacet;
            if constexpr ( is_order_dynamic )
            {
                const uint16_type nDofPerVertex = runtimeDofPerVertex();
                const uint16_type nDofPerEdge = runtimeDofPerEdge();
                const uint16_type nDofPerFace = runtimeDofPerFace();
                localDofOnFacet = face_type::numVertices * nDofPerVertex +
                                  face_type::numEdges * nDofPerEdge +
                                  face_type::numFaces * nDofPerFace;
            }
            else
            {
                localDofOnFacet = face_type::numVertices * fe_type::nDofPerVertex +
                                  face_type::numEdges * fe_type::nDofPerEdge +
                                  face_type::numFaces * fe_type::nDofPerFace;
            }

            if ( perComponent )
                return localDofOnFacet;
            if constexpr ( is_product )
                return nComponents * localDofOnFacet;
            else
                return localDofOnFacet;
        }

    [[nodiscard]] size_type feLocalDofCount( bool perComponent = false ) const noexcept
        {
            if ( M_fe )
            {
                if constexpr ( requires( fe_ptrtype const& fe, bool value ) { fe->localDofCount( value ); } )
                    return M_fe->localDofCount( perComponent );
            }
            return this->legacyLocalDofCount( perComponent );
        }

    [[nodiscard]] size_type feLocalDofCountOnFacet( uint16_type localFacet = 0,
                                                    bool perComponent = false ) const noexcept
        {
            if ( M_fe )
            {
                if constexpr ( requires( fe_ptrtype const& fe, uint16_type facet, bool value ) { fe->localDofCountOnFacet( facet, value ); } )
                    return M_fe->localDofCountOnFacet( localFacet, perComponent );
            }
            return this->legacyLocalDofCountOnFacet( perComponent );
        }

    [[nodiscard]] bool initializeDescriptorLocalIndexPermutation( uint16_type nFlatLocalDof ) noexcept
        {
            if constexpr ( !FiniteElementDofLayoutProvider<fe_type> )
                return false;
            else
            {
                if ( !M_fe )
                    return false;

                using key_type = std::tuple<int,uint16_type,uint16_type,uint16_type>;
                std::map<key_type,uint16_type> localDofByAttachment;
                for ( uint16_type localDof = 0; localDof < nFlatLocalDof; ++localDof )
                {
                    auto const layout = M_fe->localDofLayout( localDof );
                    if ( !layout.attachment.isValid() )
                        return false;
                    localDofByAttachment.emplace( key_type{ layout.attachment.entityDim,
                                                            layout.attachment.entityId,
                                                            layout.attachment.ordinal,
                                                            layout.component },
                                                  localDof );
                }

                for ( uint16_type localDof = 0; localDof < nFlatLocalDof; ++localDof )
                {
                    auto const layout = M_fe->localDofLayout( localDof );
                    auto attachment = layout.attachment;

                    if ( attachment.entityDim == 0 )
                    {
                        if ( attachment.entityId >= element_type::numVertices )
                            return false;
                        attachment.entityId = static_cast<uint16_type>( element_type::numVertices - 1 - attachment.entityId );
                    }
                    else if ( attachment.entityDim == 1 )
                    {
                        const uint16_type nEntityDof = static_cast<uint16_type>(
                            M_fe->localDofCountOnEntity( 1, attachment.entityId, true ) );
                        if ( nEntityDof == 0 || attachment.ordinal >= nEntityDof )
                            return false;
                        attachment.ordinal = static_cast<uint16_type>( nEntityDof - 1 - attachment.ordinal );
                    }
                    else
                        continue;

                    auto it = localDofByAttachment.find( key_type{ attachment.entityDim,
                                                                   attachment.entityId,
                                                                   attachment.ordinal,
                                                                   layout.component } );
                    if ( it == localDofByAttachment.end() )
                        return false;
                    M_localIndicesPerm[localDof] = it->second;
                }
                return true;
            }
        }

    void initializeLegacyLocalIndexPermutation( uint16_type nFlatLocalDof ) noexcept
        {
            const uint16_type nLocalDofRt = runtimeNLocalDof();
            const uint16_type dofPerVertex = runtimeDofPerVertex();
            const uint16_type dofPerEdge = runtimeDofPerEdge();
            const uint16_type ncdof = nLocalDofRt > 0
                ? static_cast<uint16_type>( nFlatLocalDof / nLocalDofRt )
                : 0;

            for ( uint16_type i = 0; i < nLocalDofRt; ++i )
                for ( uint16_type c = 0; c < ncdof; ++c )
                {
                    const uint16_type flatLocalDof = static_cast<uint16_type>( nLocalDofRt * c + i );
                    M_localIndicesIdentity[flatLocalDof] = flatLocalDof;

                    if ( i < dofPerVertex * element_type::numVertices )
                        M_localIndicesPerm[flatLocalDof] = static_cast<uint16_type>(
                            nLocalDofRt * c + dofPerVertex * element_type::numVertices - 1 - i );
                    else if ( i < dofPerVertex * element_type::numVertices + dofPerEdge * element_type::numEdges )
                        M_localIndicesPerm[flatLocalDof] = static_cast<uint16_type>(
                            nLocalDofRt * c + 2 * dofPerVertex * element_type::numVertices +
                            dofPerEdge * element_type::numEdges - 1 - i );
                }
        }

    void initializeLocalIndexPermutations()
        {
            const uint16_type nFlatLocalDof = static_cast<uint16_type>( runtimeNDofPerElement() );
            if ( M_localIndicesIdentity.size() < nFlatLocalDof )
                M_localIndicesIdentity.resize( nFlatLocalDof );
            if ( M_localIndicesPerm.size() < nFlatLocalDof )
                M_localIndicesPerm.resize( nFlatLocalDof );

            for ( uint16_type localDof = 0; localDof < nFlatLocalDof; ++localDof )
            {
                M_localIndicesIdentity[localDof] = localDof;
                M_localIndicesPerm[localDof] = localDof;
            }

            if ( nFlatLocalDof == 0 )
                return;

            if ( this->initializeDescriptorLocalIndexPermutation( nFlatLocalDof ) )
                return;

            this->initializeLegacyLocalIndexPermutation( nFlatLocalDof );
        }

    [[nodiscard]] uint16_type runtimeLocalDofId( uint16_type localNode, uint16_type c = 0 ) const noexcept
        {
            if ( !M_fe )
                return localNode;

            if constexpr ( requires( fe_type const& fe, uint16_type lid, uint16_type comp ) { fe.localDofId( lid, comp ); } )
                return M_fe->localDofId( localNode, c );
            else
                return static_cast<uint16_type>( runtimeNLocalDof() * c + localNode );
        }

    //! @brief Get nLocalDof value for runtime use (avoids code duplication)
    [[nodiscard]] uint16_type runtimeNLocalDof() const noexcept
        {
            return static_cast<uint16_type>( this->feLocalDofCount( true ) );
        }

    [[nodiscard]] uint16_type runtimeDofPerVertex() const noexcept
        {
            if constexpr ( is_order_dynamic )
            {
                if ( !M_fe )
                    return 0;
                if constexpr ( requires( fe_ptrtype const& fe ) { fe->dofPerVertex(); } )
                    return M_fe->dofPerVertex();
                else if constexpr ( requires( fe_ptrtype const& fe ) { fe->runtimeDofPerVertex(); } )
                    return M_fe->runtimeDofPerVertex();
                else
                    return 0;
            }
            else
                return fe_type::nDofPerVertex;
        }

    [[nodiscard]] uint16_type runtimeDofPerEdge() const noexcept
        {
            if constexpr ( is_order_dynamic )
            {
                if ( !M_fe )
                    return 0;
                if constexpr ( requires( fe_ptrtype const& fe ) { fe->dofPerEdge(); } )
                    return M_fe->dofPerEdge();
                else if constexpr ( requires( fe_ptrtype const& fe ) { fe->runtimeDofPerEdge(); } )
                    return M_fe->runtimeDofPerEdge();
                else
                    return 0;
            }
            else
                return fe_type::nDofPerEdge;
        }

    [[nodiscard]] uint16_type runtimeDofPerFace() const noexcept
        {
            if constexpr ( is_order_dynamic )
            {
                if ( !M_fe )
                    return 0;
                if constexpr ( requires( fe_ptrtype const& fe ) { fe->dofPerFace(); } )
                    return M_fe->dofPerFace();
                else if constexpr ( requires( fe_ptrtype const& fe ) { fe->runtimeDofPerFace(); } )
                    return M_fe->runtimeDofPerFace();
                else
                    return 0;
            }
            else
                return fe_type::nDofPerFace;
        }

    [[nodiscard]] uint16_type runtimeDofPerVolume() const noexcept
        {
            if constexpr ( is_order_dynamic )
            {
                if ( !M_fe )
                    return 0;
                if constexpr ( requires( fe_ptrtype const& fe ) { fe->dofPerVolume(); } )
                    return M_fe->dofPerVolume();
                else if constexpr ( requires( fe_ptrtype const& fe ) { fe->runtimeDofPerVolume(); } )
                    return M_fe->runtimeDofPerVolume();
                else
                    return 0;
            }
            else
                return fe_type::nDofPerVolume;
        }

    //! @brief Get nDofPerElement value for runtime use
    [[nodiscard]] size_type runtimeNDofPerElement() const noexcept
        {
            return this->feLocalDofCount();
        }
  public:

    global_dof_fromface_type const& faceLocalToGlobal( const size_type ElId,
                                                       const uint16_type localNode,
                                                       const uint16_type c = 0 ) const override
        {
            const size_type nDofF = nLocalDofOnFace( true );
            return M_face_l2g.find( ElId )->second[ nDofF*c+localNode ];
        }
    global_dof_fromface_type const& facetLocalToGlobal( const size_type ElId,
                                                        const uint16_type localNode,
                                                        const uint16_type c = 0 ) const
        {
            return this->faceLocalToGlobal( ElId, localNode, c );
        }

    struct element_access
    {
        element_access( DofTable const& __d )
            :
            M_d( __d )
            {}
        global_dof_type const& operator()( size_type __id, uint16_type __loc, uint16_type c = 0 ) const
            {
                return M_d.M_el_l2g.left.find( localdof_type( __id, M_d.localDofId( __loc, c ) ) )->second;
            }
        uint16_type localDofInElement( size_type __id, uint16_type __loc, uint16_type c = 0 ) const
            {
                return M_d.localDofId( __loc, c );
            }
        DofTable const& M_d;
    };
    friend struct element_access;

    struct face_access
    {

        face_access( DofTable const& __d )
            :
            M_d( __d )
            {}
        global_dof_fromface_type operator()( size_type __id, uint16_type __loc, uint16_type c = 0 ) const
            {
                return M_d.M_face_l2g.find( __id)->second[M_d.nLocalDofOnFace( true )*c+__loc];
            }

        uint16_type localDofInElement( size_type __id, uint16_type __loc, uint16_type c = 0 ) const
            {
                auto const& ldofFace = M_d.M_face_l2g.find( __id )->second[M_d.nLocalDofOnFace( true )*c+__loc];
                const uint16_type ldof = ldofFace.localDof();
                if ( M_d.fe().component( ldof ) == c )
                    return ldof;
                return M_d.localDofId( M_d.fe().dofParent( ldof ), c );
            }

        DofTable const& M_d;
    };
    friend struct face_access;

    /**
     * @brief local to global mapping
     */
    template<typename Elem>
    typename mpl::if_<mpl::equal_to<mpl::int_<Elem::nDim>,mpl::int_<nDim> >,
                      global_dof_type,
                      global_dof_fromface_type >::type
    localToGlobal( Elem const& El, const uint16_type localNode, const uint16_type c = 0 ) const
        {
            typedef typename mpl::if_<mpl::equal_to<mpl::int_<Elem::nDim>,mpl::int_<nDim> >,mpl::identity<element_access>,mpl::identity<face_access> >::type::type access_type;
            //DVLOG(2) << "dof:(" << El.id() << ", " << localNode << ")= "
            //<< access_type(*this)( El.id(), localNode, c ) << "\n";
            return access_type( *this )( El.id(), localNode, c );
        }

    template<typename Elem>
    uint16_type
    localDofInElement( Elem const& El, const uint16_type localNode, const uint16_type c = 0 ) const
        {
            typedef typename mpl::if_<mpl::equal_to<mpl::int_<Elem::nDim>,mpl::int_<nDim> >,mpl::identity<element_access>,mpl::identity<face_access> >::type::type access_type;

            return ( access_type( *this ) ).localDofInElement( El.id(), localNode, c );
        }

    /**
     * Number of elements in mesh
     */
    size_type numElements() const
        {
            return M_n_el;
        }

    /**
     * Number of local vertices (in an element)
     */
    uint16_type numLocalVertices() const
        {
            return fe_type::numVertices;
        }

    /**
     * Number of local edges (in an element)
     */
    uint16_type numLocalEdges() const
        {
            return fe_type::numEdges;
        }

    /**
     * Number of local faces (in an element)
     */
    uint16_type numLocalFaces() const
        {
            return fe_type::numFaces;
        }

    /**
     * show some information about the dof instance
     */
    void showMe() const;

    void dump() const
        {
#if 0

            for ( size_type __i = 0; __i < M_face_l2g.nrows(); ++__i )
            {
                for ( size_type __l = 0; __l < M_face_l2g.ncols(); ++__l )
                {
                    std::cout << "face " << __i << " local " << __l
                              << " to global " << M_face_l2g[ __i][ __l ] << "\n";
                }
            }

#endif // 0
        }

    /**
     * @param elt id the of element
     * @param c component of the dof
     *
     * @return true if element dof have all been computed, false otherwise
     */
    bool isElementDone( size_type elt, int c = 0 ) const
        {
            Feel::detail::ignore_unused_variable_warning( c );
            bool done = true;
            for( auto const& local_dof : this->localDofSet( elt ) )
            {
                auto it = M_el_l2g.left.find( local_dof );
                if ( it == M_el_l2g.left.end() )
                    return false;
            }
            return done;
        }

    /**
     * Initialize the dof map table
     */
    void initDofMap( mesh_type& M );

    /**
     * build the dof map
     */
    void build( mesh_type* M )
        {
            this->build( *M );
        }

    /**
     * build the dof map
     */
    void build( std::shared_ptr<mesh_type>& M )
        {
            this->build( *M );
        }

    /**
     * build the dof map
     */
    void build( mesh_type& M );

private :

    /**
     * build dof associated to local discontinuities
     */
    size_type buildLocallyDiscontinuousDofMap( mesh_type& M, size_type start_next_free_dof );

    /**
     * @brief Build the localToGlobal table
     *
     *  \param mesh mesh
     */
    void buildDofMap( mesh_type& mesh, size_type start_next_free_dof = 0 );

    /**
     * @brief Build the localToGlobal table for the boundary
     *
     * \param mesh A mesh
     */
    void buildBoundaryDofMap( mesh_type& mesh );

    /**
     * build the GlobalProcessToGlobalClusterDof table
     */
    void buildGhostDofMap( mesh_type& mesh );

    /**
     * subroutines
     */
    bool buildGlobalProcessToGlobalClusterDofMapDescriptorKeys( mesh_type& mesh );
    void buildGlobalProcessToGlobalClusterDofMapOthersMesh( mesh_type& mesh );
    void buildGlobalProcessToGlobalClusterInterprocessDofs( mesh_type& mesh,
                                                            std::map<rank_type, std::map<size_type,std::vector<uint16_type> > > & dataToSend,
                                                            std::map<rank_type, std::map<size_type,std::vector<size_type> > > & dataMemory );
    void buildGhostDofMapExtended( mesh_type& mesh, Range<mesh_type,MESH_ELEMENTS> const& ghostEltRange );

    void updateMultiprocessDofForUse();

public:
    DofTableExtendedType dofTableExtended() const noexcept { return M_buildDofTableMPIExtended; }
    bool hasDofTableExtended() const { return M_buildDofTableMPIExtended == DofTableExtendedType::VERTICES; }
    void setDofTableExtended( DofTableExtendedType b )
        {
            if ( b == DofTableExtendedType::DEFAULT )
                b = DofTableExtendedType::VERTICES;
            M_buildDofTableMPIExtended = b;
        }


    size_type nGhostDofAddedInExtendedDofTable() const { return M_nGhostDofAddedInExtendedDofTable; }




    /**
     * \return the dictionary for the global dof
     */
    dof_map_type const& mapGDof() const
        {
            return map_gdof;
        }

    /**
     * \return the dictionary for the global dof
     */
    dof_map_type& mapGDof()
        {
            return map_gdof;
        }

    [[nodiscard]] bool mapGDofUsesFlatLocalDofEntries() const noexcept
        {
            return M_mapGDofHasFlatEntries && !M_mapGDofHasLegacyEntries;
        }

    [[nodiscard]] uint16_type legacyMapGDofComponentCount() const noexcept
        {
            if constexpr ( is_tensor2symm )
                return nRealComponents;
            else if constexpr ( is_product )
                return nComponents;
            else
                return uint16_type( 1 );
        }

    template<typename Visitor>
    void forEachGlobalDofKeyEntry( Visitor&& visitor ) const
        {
            if ( this->mapGDofUsesFlatLocalDofEntries() )
            {
                for ( auto const& [dofKey,localDofIndex] : map_gdof )
                {
                    visitor( global_dof_key_entry{ .key = dofKey,
                                                   .localDofIndex = localDofIndex,
                                                   .component = this->componentFromGlobalProcessDof( localDofIndex ) } );
                }
            }
            else
            {
                const uint16_type nCompPerDof = this->legacyMapGDofComponentCount();
                for ( auto const& [dofKey,baseLocalDof] : map_gdof )
                {
                    for ( uint16_type c = 0; c < nCompPerDof; ++c )
                    {
                        visitor( global_dof_key_entry{ .key = dofKey,
                                                       .localDofIndex = baseLocalDof + c,
                                                       .component = c } );
                    }
                }
            }
        }

    /**
     * clear the dictionary
     */
    void clearMapGDof()
        {
            map_gdof.clear();
            M_mapGDofHasLegacyEntries = false;
            M_mapGDofHasFlatEntries = false;
        }

    /**
     * set the dictionary for the dictionary of the global dof
     */
    void setMapGDof( dof_map_type const& mapdof )
        {
            map_gdof = mapdof;
            M_mapGDofHasFlatEntries = false;
            M_mapGDofHasLegacyEntries = !map_gdof.empty();
        }

    typename dof_marker_type::right_range_type
    markerToDof( boost::any const& marker )
        {
            using namespace boost::bimaps;
            int id = M_mesh->markerId( marker );
            return M_dof_marker.right.range( id <= _key, _key<id+1 );
        }

    typename dof_marker_type::right_range_type
    markerToDofLessThan( boost::any const& marker )
        {
            using namespace boost::bimaps;
            int id = M_mesh->markerId( marker );
            return M_dof_marker.right.range( unbounded, _key<id );
        }
    typename dof_marker_type::right_range_type
    markerToDofGreaterThan( boost::any const& marker )
        {
            using namespace boost::bimaps;
            int id = M_mesh->markerId( marker );
            return M_dof_marker.right.range( id<_key, unbounded );
        }

    void printDofMarker(std::string const& filename )
        {
            // std::ofstream ofs( filename.c_str() );
            // for( auto dof : _M_dof_marker )
            // {
            //     //ofs << dof.first << " " << dof.second << "\n";
            // }
            std::ofstream ofs( filename.c_str() );
            for( auto dofleft : M_dof_marker.left )
            {
                ofs << dofleft.first << " " << dofleft.second << "\n";
            }
        }
    /**
     * The dof are ordered such that they are contiguous per element
     * and components. This way an extraction of the dof indices in
     * one element allows to extract a view of the corresponding
     * coefficient in a given basis which is then very helpful for
     * interpolation for example.
     *
     * \param ie index of the element
     * \param lc_dof index of the dof in the element
     * \param lc local index of the entity associated with the dof in the element
     * \param gDof global dof index
     * \param pDof dof index in the processor
     *
     * \return the index of the next free dof in the processor
     */
    bool insertDof( size_type ie,
                    uint16_type l_dof,
                    uint16_type lc,
                    dof_type && gDof,
                    rank_type processor,
                    size_type& pDof,
                    int32_type sign = 1,
                    bool is_dof_periodic = false,
                    size_type shift = 0,
                    mesh_marker_type const& marker = mesh_marker_type{} )
        {
            bool res = true;
            M_mapGDofHasLegacyEntries = true;
            const int ncdof = is_product?nComponents:1;
            const uint16_type nldof = runtimeNLocalDof();

            //for ( int c = 0; c < ncdof; ++c )
            {
                uint16_type lc_dof = nldof*0+l_dof;
                Feel::detail::ignore_unused_variable_warning( lc );
                dof_map_iterator endof = map_gdof.end();

                auto [ itdof, __inserted ] = map_gdof.try_emplace( gDof, dofIndex( pDof ) );

                if ( __inserted )
                {
                    if constexpr ( is_tensor2symm )
                        pDof += nRealComponents;
                    else
                        pDof += ncdof;
                }
                M_ldof.set(ie,lc_dof);
                auto eit = M_el_l2g.left.find( M_ldof  );
                // make sure that no already created dof is overwritten here (may be done elsewhere)
                if ( eit == M_el_l2g.left.end() )
                {
                    M_gdof.set( itdof->second+shift, sign, is_dof_periodic );
                    DCHECK( itdof->first == gDof ) << "very bad logical error in insertDof";
#if 0
                    DCHECK( lc_dof >= fe_type::nLocalDof*(std::get<1>(itdof->first)) &&
                            lc_dof < fe_type::nLocalDof*( std::get<1>(itdof->first)+1 ) )
                        << "invalid local dof index"
                        <<  lc_dof << ", " << fe_type::nLocalDof*std::get<1>(itdof->first);
#endif
                    if constexpr ( is_tensor2symm )
                    {
                        for( int c1 = 0; c1 < nComponents1; ++c1 )
                        {
                            for( int c2 = 0; c2 < c1; ++c2 )
                            {
                                const int k = Feel::detail::symmetricIndex(c1,c2,nComponents1);
                                M_ldof.setLocalDof( nldof*(nComponents1*c1+c2)+l_dof );
                                M_gdof.setIndex( itdof->second+shift+k );
                                auto res = M_el_l2g.insert( dof_relation( M_ldof, M_gdof ) );
                                DCHECK( res.second ) << "global dof " << itdof->second+shift+k << " not inserted in local dof (" <<
                                    ie << "," << lc_dof << ")";

                                M_ldof.setLocalDof( nldof*(nComponents1*c2+c1)+l_dof );
                                res = M_el_l2g.insert( dof_relation( M_ldof, M_gdof ) );
                                DCHECK( res.second ) << "global dof " << itdof->second+shift+k << " not inserted in local dof (" <<
                                    ie << "," << lc_dof << ")";
                                //(Dof  itdof->second+shift, sign, is_dof_periodic, 0, 0, marker.value() ) ) );

                                if ( !marker.empty() )
                                    M_dof_marker.insert( dof2marker( itdof->second+shift+k,  marker.value() ) );
                            }
                            M_ldof.setLocalDof( nldof*(nComponents1*c1+c1)+l_dof );
                            const int k = Feel::detail::symmetricIndex(c1,c1,nComponents1);
                            M_gdof.setIndex( itdof->second+shift+k);
                            auto res = M_el_l2g.insert( dof_relation( M_ldof, M_gdof ) );
                            DCHECK( res.second ) << "global dof " << itdof->second+shift+k << " not inserted in local dof (" <<
                                ie << "," << lc_dof << ")";
                            if ( !marker.empty() )
                                M_dof_marker.insert( dof2marker( itdof->second+shift+k,  marker.value() ) );
                        }
                    }
                    else
                    {
                        for( int c = 0; c < ncdof; ++c )
                        {
                            M_ldof.setLocalDof( nldof*c+l_dof );
                            M_gdof.setIndex( itdof->second+shift+c );
                            auto res = M_el_l2g.insert( dof_relation( M_ldof, M_gdof ) );
                            //(Dof  itdof->second+shift, sign, is_dof_periodic, 0, 0, marker.value() ) ) );
                            DCHECK( res.second ) << "global dof " << itdof->second+shift << " not inserted in local dof (" <<
                                ie << "," << lc_dof << ")";
                            if ( !marker.empty() )
                                M_dof_marker.insert( dof2marker( itdof->second+shift+c,  marker.value() ) );
                        }
                    }



#if 0// !defined(NDEBUG)
                    M_dof2elt[itdof->second+shift].push_back( boost::make_tuple( ie, lc_dof, lc, std::get<0>(itdof->first) ) );
#endif
#if 0
                    M_dof_view.insert( Dof( itdof->second+shift,      // global index
                                            sign,                     // sign
                                            std::get<0>(itdof->first),    // entity type
                                            false,                    // is on boundary ?
                                            0                         // marker
                                            ) );
#endif
#if 0
                    for ( index i2 = 0; i2 < nLocalDof(); ++i2 )
                        VLOG(1) << "dof table( " << ie << ", " << lc  << ")=" << M_el_l2g.left.find(localdof_type(ie,i2))->second.index() << "\n";

#endif
                }
#if 0
                else
                {
                    size_type _dof = M_el_l2g.left.find(localdof_type(ie,lc_dof))->second.index();

                    CHECK(  M_dof_marker[_dof] == marker.value() ) << "Invalid dof marker, element id: " <<  ie
                                                                   << ", local dof id: " << lc_dof
                                                                   << ", global dof id: "<< _dof
                                                                   << ", dof marker: " <<  M_dof_marker[_dof]
                                                                   << ", marker: " << marker.value() << "\n";
                }
#endif

                res = res && ( __inserted || ( ( M_el_l2g.left.find(localdof_type(ie,lc_dof)) != M_el_l2g.left.end() ) && shift ) );
            }

            return res;
        }

    /**
     * Insert exactly one flattened local DOF described by the FE layout.
     *
     * This is the descriptor-backed insertion path. It intentionally does not
     * expand product components; callers pass the already flattened local slot
     * and the FE-owned component metadata needed to build a unique global key.
     */
    bool insertFlatDof( size_type ie,
                        uint16_type localDof,
                        uint16_type lc,
                        dof_type gDof,
                        uint16_type component,
                        uint16_type componentStride,
                        rank_type processor,
                        size_type& pDof,
                        int32_type sign = 1,
                        bool is_dof_periodic = false,
                        size_type shift = 0,
                        mesh_marker_type const& marker = mesh_marker_type{} )
        {
            Feel::detail::ignore_unused_variable_warning( lc );
            Feel::detail::ignore_unused_variable_warning( processor );
            M_mapGDofHasFlatEntries = true;

            const uint16_type stride = componentStride > 0 ? componentStride : 1;
            FEELPP_ASSERT( component < stride )
                ( component )( stride )( localDof ).error( "invalid flattened dof component" );

            if ( stride > 1 )
                std::get<1>( gDof ) = std::get<1>( gDof ) * stride + component;

            auto [itdof, inserted] = map_gdof.try_emplace( gDof, dofIndex( pDof ) );
            if ( inserted )
                ++pDof;

            M_ldof.set( ie, localDof );
            auto eit = M_el_l2g.left.find( M_ldof );
            if ( eit == M_el_l2g.left.end() )
            {
                M_gdof.set( itdof->second + shift, sign, is_dof_periodic );
                DCHECK( itdof->first == gDof ) << "invalid descriptor-backed dof insertion";
                auto res = M_el_l2g.insert( dof_relation( M_ldof, M_gdof ) );
                DCHECK( res.second ) << "global dof " << itdof->second + shift
                                     << " not inserted in flat local dof ("
                                     << ie << "," << localDof << ")";
                if ( !marker.empty() )
                    M_dof_marker.insert( dof2marker( itdof->second + shift, marker.value() ) );
            }

            return inserted || M_el_l2g.left.find( localdof_type( ie, localDof ) ) != M_el_l2g.left.end();
        }

    /**
     * rebuild dof points
     */
    void rebuildDofPoints( mesh_type& M )
        {
            M_dof_points.clear();
            M_hasBuiltDofPoints = false;
            this->generateDofPoints(M);

            if ( this->worldComm().localSize()>1 && this->hasDofTableExtended() )
            {
                Range<mesh_type,MESH_ELEMENTS> rangeExtendedElements;
                if (this->hasMeshSupport())
                {
                    rangeExtendedElements = elements(this->meshSupport(), entity_process_t::GHOST_ONLY );
                }
                else
                {
                    rangeExtendedElements =  elements( M, entity_process_t::GHOST_ONLY );
                }
                this->generateDofPoints( rangeExtendedElements );
                //this->generateDofPointsExtendedGhostMap(M);
            }
        }

    /**
     * build point id to dof id relationship
     * if \p dof2pid is true then generate dof to point id relation,
     * if \p pid2dof is true then generate point id to dof relation,
     */
    std::pair<std::unordered_map<size_type,size_type>,std::unordered_map<size_type,size_type> >
    pointIdToDofRelation( std::string fname="", bool dof2pid = true, bool pid2dof = true ) const;
private:
    template<typename, typename > friend class DofFromElement;
    template<typename, typename, typename > friend class DofFromMortar;
    template<typename, typename > friend class DofFromBoundary;
    template<typename, typename > friend class DofFromEdge;

    void addSubstructuringDofMap( mesh_type const& M, size_type next_free_dof );
    void addSubstructuringDofVertex( mesh_type const& M, size_type next_free_dof );
    void addSubstructuringDofEdge( mesh_type const& M, size_type next_free_dof, mpl::int_<1> );
    void addSubstructuringDofEdge( mesh_type const& M, size_type next_free_dof, mpl::int_<2> );
    void addSubstructuringDofEdge( mesh_type const& M, size_type next_free_dof, mpl::int_<3> );
    void addSubstructuringDofFace( mesh_type const& M, size_type next_free_dof, mpl::int_<1> );
    void addSubstructuringDofFace( mesh_type const& M, size_type next_free_dof, mpl::int_<2> );
    void addSubstructuringDofFace( mesh_type const& M, size_type next_free_dof, mpl::int_<3> );

    /**
     * @brief Checks if the dofs associated with entity_id N are continuous
     *
     * \param M A mesh
     */
    template<int N>
    void checkDofEntity ( mesh_type& M )
        {
            Feel::detail::ignore_unused_variable_warning( M );

            if ( !is_scalar )
                return;

#if !defined(NDEBUG)

            using namespace Feel;

            gm_ptrtype M_gm_ptr = M.gm();

            fe_type M_basis;

            //value_type tol = value_type(100.0)*type_traits<value_type>::epsilon();
            value_type tol = value_type( 10.0 )*type_traits<double>::epsilon();

            bool global_signs_good = 1;

            std::vector<size_type> bad_dof;

            for ( uint16_type gDof = 0; gDof < this->nDof(); ++gDof )
            {
                uint16_type _numEntities = M_dof2elt[gDof].size();
                uint16_type _ent = M_dof2elt[gDof].begin()->template get<3>();

                if ( _numEntities > 1 && _ent == mpl::int_<N>() )
                {
                    bool signs_good = 1;

                    std::vector< ublas::vector<value_type> > basis_eval;

                    std::vector< points_type > real_coordinates;

                    ldof_const_iterator __ldofit = M_dof2elt[gDof].begin();
                    ldof_const_iterator __ldofen = M_dof2elt[gDof].end();

                    while ( __ldofit != __ldofen )
                    {
                        size_type entity_element_id = __ldofit->template get<0>();
                        uint16_type entity_local_dof_id = __ldofit->template get<1>();
                        uint16_type entity_local_id = __ldofit->template get<2>();

                        PointSetMapped<element_type, convex_type, nOrder> test( M.element( entity_element_id ) );

                        points_type Gt = test.pointsBySubEntity( N, entity_local_id );

                        real_coordinates.push_back( test.pointsBySubEntity( N, entity_local_id, 0, 1 ) );

                        int sign = boost::get<1>( localToGlobal( entity_element_id, entity_local_dof_id ) );

                        basis_eval.push_back( value_type( sign )*ublas::row( M_basis.evaluate( Gt ), entity_local_dof_id ) );

                        ++__ldofit;
                    }

                    for ( uint16_type i=1; i < _numEntities; i++ )
                    {

                        FEELPP_ASSERT( ublas::norm_inf( real_coordinates[i] - real_coordinates[0] ) < tol  )
                            ( gDof )
                            ( real_coordinates[0] )
                            ( real_coordinates[i] ).error( "Reference points aren't being mapped to the same real one's" );

                        if ( ublas::norm_inf( basis_eval[i] - basis_eval[0] ) > tol )
                        {
                            signs_good = 0;
                            global_signs_good = 0;
                        }
                    }

                    basis_eval.empty();
                    real_coordinates.empty();

                    if ( signs_good == 0 )
                        bad_dof.push_back( gDof );
                }
            }

            if ( !bad_dof.empty() )
            {
                for ( uint16_type i = 0; i < bad_dof.size(); ++i )
                    LOG(WARNING) << bad_dof[i] << "\n";

                if ( mpl::int_<N>() == 1 )
                    LOG(WARNING) << "Edges: ";

                else
                    LOG(WARNING) << "Faces: ";

                LOG(WARNING) << "Bad dof signs. \n";
            }

#endif
        }

    void checkDofContinuity( mesh_type& /*mesh*/, mpl::int_<1> ) {}

    void checkDofContinuity( mesh_type& mesh, mpl::int_<2> )
        {
            checkDofEntity<1>( mesh );
        }

    void checkDofContinuity( mesh_type& mesh, mpl::int_<3> )
        {
            checkDofContinuity( mesh, mpl::int_<2>() );
            checkDofEntity<2>( mesh );
        }

    [[nodiscard]] static bool
    hasTriangularSidePointCount( uint16_type nFaceDof, uint16_type& nSidePoints )
        {
            if ( nFaceDof == 0 )
                return false;

            const int disc = 1 + 8 * int( nFaceDof );
            const int root = static_cast<int>( std::lround( std::sqrt( static_cast<double>( disc ) ) ) );
            if ( root * root != disc )
                return false;

            const int sidePoints = ( root - 1 ) / 2;
            if ( sidePoints <= 0 || sidePoints * ( sidePoints + 1 ) / 2 != int( nFaceDof ) )
                return false;

            nSidePoints = static_cast<uint16_type>( sidePoints );
            return true;
        }

    [[nodiscard]] static bool
    hasQuadrangularSidePointCount( uint16_type nFaceDof, uint16_type& nSidePoints )
        {
            if ( nFaceDof == 0 )
                return false;

            const int root = static_cast<int>( std::lround( std::sqrt( static_cast<double>( nFaceDof ) ) ) );
            if ( root <= 0 || root * root != int( nFaceDof ) )
                return false;

            nSidePoints = static_cast<uint16_type>( root );
            return true;
        }

    [[nodiscard]] static permutation_vector_type
    composeFacePermutation( permutation_vector_type const& first,
                            permutation_vector_type const& second )
        {
            FEELPP_ASSERT( first.size() == second.size() )
                ( first.size() )( second.size() ).error( "invalid permutation composition sizes" );
            permutation_vector_type out( first.size() );
            for ( uint16_type i = 0; i < first.size(); ++i )
            {
                FEELPP_ASSERT( first( i ) < second.size() )
                    ( i )( first( i ) )( second.size() ).error( "invalid permutation composition index" );
                out( i ) = second( first( i ) );
            }
            return out;
        }

    void
    clearFacePermutationTables()
        {
            vector_permutation.clear();
            vector_permutation_dense.clear();
            vector_permutation_dense.resize( face_permutation_type::N_PERMUTATIONS );
        }

    void
    setFacePermutationVector( face_permutation_type permutation,
                              permutation_vector_type perm )
        {
            const auto idx = static_cast<size_type>( permutation.value() );
            FEELPP_ASSERT( idx < vector_permutation_dense.size() )
                ( int( permutation.value() ) )
                ( vector_permutation_dense.size() )
                .error( "invalid face permutation index" );

            vector_permutation[permutation] = perm;
            vector_permutation_dense[idx] = std::move( perm );
        }

    bool
    tryGenerateSimplexTriangularFacePermutations( uint16_type nDofPerFace )
        {
            if constexpr ( !( nDim == 3 && convex_type::is_simplex ) )
                return false;

            uint16_type nSidePoints = 0;
            if ( !hasTriangularSidePointCount( nDofPerFace, nSidePoints ) )
                return false;

            permutation_vector_type reverseHypotenuse( nDofPerFace );
            permutation_vector_type reverseBase( nDofPerFace );

            // Reverse along the hypotenuse direction.
            for ( uint16_type i = 0; i < nSidePoints; ++i )
            {
                for ( uint16_type j = 0; j <= i; ++j )
                {
                    const uint16_type first = i + j * ( j - 1 ) / 2 + j * ( nSidePoints - j );
                    const uint16_type last = i +
                                             ( i - j ) * ( i - j - 1 ) / 2 +
                                             ( i - j ) * ( nSidePoints - ( i - j ) );
                    FEELPP_ASSERT( first < nDofPerFace && last < nDofPerFace )
                        ( first )( last )( nDofPerFace ).error( "invalid triangular permutation index" );
                    reverseHypotenuse( first ) = last;
                }
            }

            // Reverse along the base direction.
            for ( uint16_type i = 0; i < nSidePoints; ++i )
            {
                const uint16_type begin = i * nSidePoints - i * ( i - 1 ) / 2;
                const uint16_type end = ( i + 1 ) * nSidePoints - ( i + 1 ) * i / 2 - 1;
                for ( uint16_type j = 0; j <= nSidePoints - i - 1; ++j )
                {
                    FEELPP_ASSERT( begin + j < nDofPerFace && end >= j && end - j < nDofPerFace )
                        ( begin + j )( end - j )( nDofPerFace ).error( "invalid triangular permutation index" );
                    reverseBase( begin + j ) = end - j;
                }
            }

            auto const rotationAntiClock = composeFacePermutation( reverseBase, reverseHypotenuse );
            auto const rotationClockwise = composeFacePermutation( reverseHypotenuse, reverseBase );
            auto const reverseHeight = composeFacePermutation( reverseBase, rotationClockwise );

            this->setFacePermutationVector( face_permutation_type( face_permutation_type::REVERSE_HYPOTENUSE ),
                                            std::move( reverseHypotenuse ) );
            this->setFacePermutationVector( face_permutation_type( face_permutation_type::REVERSE_BASE ),
                                            std::move( reverseBase ) );
            this->setFacePermutationVector( face_permutation_type( face_permutation_type::ROTATION_ANTICLOCK ),
                                            std::move( rotationAntiClock ) );
            this->setFacePermutationVector( face_permutation_type( face_permutation_type::ROTATION_CLOCKWISE ),
                                            std::move( rotationClockwise ) );
            this->setFacePermutationVector( face_permutation_type( face_permutation_type::REVERSE_HEIGHT ),
                                            std::move( reverseHeight ) );
            return true;
        }

    bool
    tryGenerateHypercubeQuadrangularFacePermutations( uint16_type nDofPerFace )
        {
            if constexpr ( !( nDim == 3 && !convex_type::is_simplex ) )
                return false;

            uint16_type nSidePoints = 0;
            if ( !hasQuadrangularSidePointCount( nDofPerFace, nSidePoints ) )
                return false;

            permutation_vector_type reverseBase( nDofPerFace );
            permutation_vector_type rotationAntiClock( nDofPerFace );

            uint16_type p = 0;
            for ( int16_type i = static_cast<int16_type>( nSidePoints ) - 1; i >= 0; --i )
            {
                const uint16_type first = static_cast<uint16_type>( i * nSidePoints );
                for ( uint16_type j = 0; j < nSidePoints; ++j, ++p )
                    reverseBase( p ) = static_cast<uint16_type>( first + j );
            }

            p = 0;
            for ( int16_type i = static_cast<int16_type>( nSidePoints ) - 1; i >= 0; --i )
            {
                for ( uint16_type j = 0; j < nSidePoints; ++j, ++p )
                    rotationAntiClock( p ) = static_cast<uint16_type>( i + nSidePoints * j );
            }

            auto const secondDiagonal = composeFacePermutation( reverseBase, rotationAntiClock );
            auto const reverseHeight = composeFacePermutation( secondDiagonal, rotationAntiClock );
            auto const rotationTwiceClockwise = composeFacePermutation( reverseHeight, reverseBase );
            auto const principalDiagonal = composeFacePermutation( reverseHeight, rotationAntiClock );
            auto const rotationClockwise = composeFacePermutation( rotationAntiClock, rotationTwiceClockwise );

            this->setFacePermutationVector( face_permutation_type( face_permutation_type::REVERSE_BASE ),
                                            std::move( reverseBase ) );
            this->setFacePermutationVector( face_permutation_type( face_permutation_type::ROTATION_ANTICLOCK ),
                                            std::move( rotationAntiClock ) );
            this->setFacePermutationVector( face_permutation_type( face_permutation_type::SECOND_DIAGONAL ),
                                            std::move( secondDiagonal ) );
            this->setFacePermutationVector( face_permutation_type( face_permutation_type::REVERSE_HEIGHT ),
                                            std::move( reverseHeight ) );
            this->setFacePermutationVector( face_permutation_type( face_permutation_type::ROTATION_TWICE_CLOCKWISE ),
                                            std::move( rotationTwiceClockwise ) );
            this->setFacePermutationVector( face_permutation_type( face_permutation_type::PRINCIPAL_DIAGONAL ),
                                            std::move( principalDiagonal ) );
            this->setFacePermutationVector( face_permutation_type( face_permutation_type::ROTATION_CLOCKWISE ),
                                            std::move( rotationClockwise ) );
            return true;
        }

    [[nodiscard]] bool
    hasValidFacePermutation( face_permutation_type permutation,
                             uint16_type nDofPerFace ) const
        {
            const auto idx = static_cast<size_type>( permutation.value() );
            if ( idx < vector_permutation_dense.size() )
                return vector_permutation_dense[idx].size() == nDofPerFace;
            auto const it = vector_permutation.find( permutation );
            return it != vector_permutation.end() && it->second.size() == nDofPerFace;
        }

    [[nodiscard]] permutation_vector_type const&
    facePermutationVector( face_permutation_type permutation,
                           uint16_type nDofPerFace ) const
        {
            FEELPP_ASSERT( hasValidFacePermutation( permutation, nDofPerFace ) )
                ( int( permutation.value() ) )( int( nDofPerFace ) )
                .error( "missing/invalid face permutation vector" );

            const auto idx = static_cast<size_type>( permutation.value() );
            if ( idx < vector_permutation_dense.size() )
                return vector_permutation_dense[idx];

            auto const it = vector_permutation.find( permutation );
            FEELPP_ASSERT( it != vector_permutation.end() )
                ( int( permutation.value() ) ).error( "missing face permutation vector" );
            return it->second;
        }

    void generateFacePermutations ( mesh_type& /*mesh*/, mpl::bool_<false> ) {}

    void generateFacePermutations ( mesh_type& mesh, mpl::bool_<true> )
        {
            if (! mesh.numElements() )
            {
                this->clearFacePermutationTables();
                return;
            }

            const uint16_type nDofPerFace = runtimeDofPerFace();
            if ( nDofPerFace <= 1 )
            {
                this->clearFacePermutationTables();
                return;
            }

            this->clearFacePermutationTables();

            if constexpr ( nDim < 3 )
            {
                return;
            }
            else if constexpr ( is_order_dynamic )
            {
                bool generated = false;
                if constexpr ( convex_type::is_simplex )
                    generated = this->tryGenerateSimplexTriangularFacePermutations( nDofPerFace );
                else
                    generated = this->tryGenerateHypercubeQuadrangularFacePermutations( nDofPerFace );

                FEELPP_ASSERT( generated )
                    ( int( nDofPerFace ) )
                    .error( "dynamic-order face permutations require triangular or quadrangular face cardinality; provide FE-owned facet ordinal permutations for this element" );
            }
            else
            {
                element_type const& _elt = mesh.beginElement()->second;
                PointSetMapped<element_type, convex_type, nOrder> pts( _elt );

                for ( uint16_type i = 2; i < face_permutation_type::N_PERMUTATIONS; i++ )
                {
                    auto perm = pts.getVectorPermutation( face_permutation_type( i ) );
                    if ( perm.size() == nDofPerFace )
                        this->setFacePermutationVector( face_permutation_type( i ), std::move( perm ) );
                }

                if constexpr ( convex_type::is_simplex )
                {
                    bool missing = false;
                    for ( uint16_type i = 2; i < face_permutation_type::N_PERMUTATIONS; ++i )
                    {
                        if ( !hasValidFacePermutation( face_permutation_type( i ), nDofPerFace ) )
                        {
                            missing = true;
                            break;
                        }
                    }
                    if ( missing )
                        this->tryGenerateSimplexTriangularFacePermutations( nDofPerFace );
                }
            }

            for ( uint16_type i = 2; i < face_permutation_type::N_PERMUTATIONS; ++i )
            {
                FEELPP_ASSERT( hasValidFacePermutation( face_permutation_type( i ), nDofPerFace ) )
                    ( int( i ) )( int( nDofPerFace ) )
                    ( int( face_permutation_type::N_PERMUTATIONS ) )
                    .error( "missing/invalid face permutation table" );
            }
        }
    /**
     * @return true if dof points are computed, false otherwise
     * @attention dof points are n
     */
    bool hasDofPoints() const { return M_hasBuiltDofPoints;/*!M_dof_points.empty();*/ }
    void generateDofPoints( mesh_type& M, bool buildMinimalParallel = false ) const;
    void generateDofPointsExtendedGhostMap( mesh_type& M ) const;
    void generateDofPoints( Range<mesh_type,MESH_ELEMENTS> const& range ) const;

private:
    //void generateDofPoints( mesh_type& M, bool buildMinimalParallel, mpl::bool_<true> ) const;
    //void generateDofPoints( mesh_type& M, bool buildMinimalParallel, mpl::bool_<false> ) const;
private:

    mesh_type* M_mesh;
    mesh_support_ptrtype M_meshSupport;

    fe_ptrtype M_fe;

    reference_convex_type M_convex_ref;

    size_type M_n_el;
    uint16_type M_n_dof_per_face_on_bdy;
    uint16_type M_n_dof_per_face;

    dof_table M_el_l2g;
    Container_fromface M_face_l2g;

    mutable local_dof_set_type M_local_dof_set;
    dof_element_type M_dof2elt;
    dof_marker_type M_dof_marker;

    dof_map_type map_gdof;
    bool M_mapGDofHasLegacyEntries = false;
    bool M_mapGDofHasFlatEntries = false;
    localdof_type M_ldof;
    global_dof_type M_gdof;

    std::map<face_permutation_type, permutation_vector_type> vector_permutation;
    // Dense permutation cache indexed by face_permutation_type::value() for hot-path access.
    std::vector<permutation_vector_type> vector_permutation_dense;

    /**
     * coordinates of the nodal dofs
     */
    mutable dof_points_type M_dof_points;
    mutable bool M_hasBuiltDofPoints;

    std::vector<globaldof_type> M_dof_indices;

    /// a view of the dof container
    //dof_container_type M_dof_view;

    vector_indices_type M_locglob_indices;
    vector_indices_type M_locglob_signs;
    localglobal_indices_type M_locglob_nosigns;
    vector_transforms_type M_locglob_transforms;
    localglobal_transforms_type M_locglob_notransforms;

    DofTableExtendedType M_buildDofTableMPIExtended = DofTableExtendedType::VERTICES;
    size_type M_nGhostDofAddedInExtendedDofTable;
    bool M_hasDescriptorKeyClusterDofMap = false;

    std::vector<uint16_type> M_localIndicesPerm, M_localIndicesIdentity;

    dof_from_edge_type M_dfe;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<typename MeshType, typename FEType, typename MortarType>
DofTable<MeshType, FEType, MortarType>::DofTable( mesh_type& mesh,
                                                                   fe_ptrtype const& _fe,
                                                                   WorldComm const& _worldComm )
    :
    super( _worldComm ),
    M_fe( _fe ),
    M_n_el( invalid_v<size_type> ),
    M_n_dof_per_face_on_bdy( invalid_uint16_type_value ),
    M_n_dof_per_face( invalid_uint16_type_value ),
    M_el_l2g(),
    M_face_l2g(),
    M_local_dof_set( 0, is_order_dynamic ? 0 : nLocalDof() ),
    map_gdof(),
    M_mapGDofHasLegacyEntries( false ),
    M_mapGDofHasFlatEntries( false ),
    M_hasBuiltDofPoints( false ),
    M_dof_indices(),
    M_buildDofTableMPIExtended( DofTableExtendedType::VERTICES ),
    M_nGhostDofAddedInExtendedDofTable( 0 ),
    M_localIndicesPerm( is_order_dynamic ? 0 : nDofPerElement ),
    M_localIndicesIdentity( is_order_dynamic ? 0 : nDofPerElement ),
    M_dfe( this )
{
    // For dynamic order, resize vectors after M_fe is available
    if constexpr ( is_order_dynamic )
    {
        const size_type ndpe = runtimeNDofPerElement();
        M_local_dof_set = local_dof_set_type( 0, nLocalDof() );
        M_localIndicesPerm.resize( ndpe );
        M_localIndicesIdentity.resize( ndpe );
    }

    size_type start_next_free_dof = 0;

    buildDofMap( mesh, start_next_free_dof );
    if ( !is_mortar )
        buildBoundaryDofMap( mesh );
    map_gdof.clear();
}

template<typename MeshType, typename FEType, typename MortarType>
DofTable<MeshType, FEType, MortarType>::DofTable( fe_ptrtype const& _fe,
                                                                   WorldComm const& _worldComm )
    :
    super( _worldComm ),
    M_fe( _fe ),
    M_n_el( 0 ),
    M_n_dof_per_face_on_bdy( invalid_uint16_type_value ),
    M_n_dof_per_face( invalid_uint16_type_value ),
    M_el_l2g(),
    M_face_l2g(),
    M_local_dof_set( 0, is_order_dynamic ? 0 : nLocalDof() ),
    map_gdof(),
    M_mapGDofHasLegacyEntries( false ),
    M_mapGDofHasFlatEntries( false ),
    M_hasBuiltDofPoints( false ),
    M_dof_indices(),
    M_buildDofTableMPIExtended( DofTableExtendedType::VERTICES ),
    M_nGhostDofAddedInExtendedDofTable( 0 ),
    M_localIndicesPerm( is_order_dynamic ? 0 : nDofPerElement ),
    M_localIndicesIdentity( is_order_dynamic ? 0 : nDofPerElement ),
    M_dfe( this )
{
    // For dynamic order, resize vectors after M_fe is available
    if constexpr ( is_order_dynamic )
    {
        const size_type ndpe = runtimeNDofPerElement();
        M_local_dof_set = local_dof_set_type( 0, nLocalDof() );
        M_localIndicesPerm.resize( ndpe );
        M_localIndicesIdentity.resize( ndpe );
    }
}

template<typename MeshType, typename FEType, typename MortarType>
DofTable<MeshType, FEType, MortarType>::DofTable( const self_type & dof2 )
    :
    super( dof2 ),
    M_fe( dof2.M_fe ),
    M_n_el( dof2.M_n_el ),
    M_n_dof_per_face_on_bdy( dof2.M_n_dof_per_face_on_bdy ),
    M_n_dof_per_face( dof2.M_n_dof_per_face ),
    M_el_l2g( dof2.M_el_l2g ),
    M_face_l2g( dof2.M_face_l2g ),
    M_local_dof_set( dof2.M_local_dof_set),
    map_gdof( dof2.map_gdof ),
    M_mapGDofHasLegacyEntries( dof2.M_mapGDofHasLegacyEntries ),
    M_mapGDofHasFlatEntries( dof2.M_mapGDofHasFlatEntries ),
    M_hasBuiltDofPoints( false ),
    M_dof_indices( dof2.M_dof_indices ),
    M_locglob_indices( dof2.M_locglob_indices ),
    M_locglob_signs( dof2.M_locglob_signs ),
    M_locglob_nosigns( dof2.M_locglob_nosigns ),
    M_locglob_transforms( dof2.M_locglob_transforms ),
    M_locglob_notransforms( dof2.M_locglob_notransforms ),
    M_buildDofTableMPIExtended( dof2.M_buildDofTableMPIExtended ),
    M_nGhostDofAddedInExtendedDofTable( dof2.M_nGhostDofAddedInExtendedDofTable ),
    M_localIndicesPerm( dof2.M_localIndicesPerm ),
    M_localIndicesIdentity( dof2.M_localIndicesIdentity ),
    M_dfe( dof2.M_dfe )
{
}

template<typename MeshType, typename FEType, typename MortarType>
void
DofTable<MeshType, FEType, MortarType>::showMe() const
{
    LOG(INFO)  << " Degree of Freedom (DofTable) Object" << "\n";
    //if ( verbose )
    {
        LOG(INFO) <<  "* nDof = " << this->nLocalDof() << "\n";
        LOG(INFO)  << "************************************************************" << "\n";
        LOG(INFO)  << "           Local to Global DOF table" << "\n";
        LOG(INFO)  << "************************************************************" << "\n";
        LOG(INFO)  << "Element Id    Loc. N.    Global N.   Sign#    Element Id   Loc. N.  Global N.  Sign" << "\n";

        for ( size_type i = 0; i < M_n_el; ++i )
        {

            for ( size_type j = 0; j < nDofPerElement; ++j )
            {

                LOG(INFO)<< "elt id " << i << " : "
                         << "(local/global : " << j << " : "
                         << std::get<0>( localToGlobal( i  , j ) ) << "  ";
            }

        }

        LOG(INFO)  << "\n";

        LOG(INFO)  << "************************************************************" << "\n";
        LOG(INFO)  << " Boundary  Local to Global DOF table" << "\n";
        LOG(INFO)  << "************************************************************" << "\n";

        auto it = M_face_l2g.begin();
        auto en = M_face_l2g.end();

        for ( size_type f = 0; it!=en; ++it,++f )
        {
            std::ostringstream ostr;
            ostr  << "face id " << it->first << " : ";

            auto it2 = it->second.begin();
            auto en2 = it->second.end();

            for ( auto const& facedof : it->second )
            {
                ostr << "(local/global/sign dof : " << facedof.localDof() << " : "
                     << facedof.index()  << "\n";
            }

            LOG(INFO) << ostr.str() << "\n";
        }
    }

}

template<typename MeshType, typename FEType, typename MortarType>
void
DofTable<MeshType, FEType, MortarType>::initDofMap( mesh_type& M )
{
    size_type numMeshElements = (this->hasMeshSupport())? this->meshSupport()->numElements() : M.numElements();
    M_n_el = numMeshElements;

    const size_type nldof = this->feLocalDofCount( true );
    const size_type nFlatLocalDof = this->feLocalDofCount();
    if constexpr ( is_order_dynamic )
    {
        const uint16_type nDofPerVertex = runtimeDofPerVertex();
        const uint16_type nDofPerEdge = runtimeDofPerEdge();
        const uint16_type nDofPerFace = runtimeDofPerFace();
        const uint16_type nDofPerVolume = runtimeDofPerVolume();

        VLOG(2) << "==============================\n";
        VLOG(2) << "[initDofMap] (dynamic order)\n";
        VLOG(2) << "is_hdiv_conforming     = "  << is_hdiv_conforming << "\n";
        VLOG(2) << "is_hcurl_conforming    = "  << is_hcurl_conforming << "\n";
        VLOG(2) << "nldof                  = "  << int( nldof ) << "\n";
        VLOG(2) << "nFlatLocalDof          = "  << int( nFlatLocalDof ) << "\n";
        VLOG(2) << "runtimeLocalDof        = "  << int( runtimeNLocalDof() ) << "\n";
        VLOG(2) << "nDofPerVolume          = "  << int( nDofPerVolume ) << "\n";
        VLOG(2) << "nDofPerFace            = "  << int( nDofPerFace ) << "\n";
        VLOG(2) << "nDofPerEdge            = "  << int( nDofPerEdge ) << "\n";
        VLOG(2) << "nDofPerVertex          = "  << int( nDofPerVertex ) << "\n";
        VLOG(2) << "element_type::numVolumes= "  << int( element_type::numVolumes ) << "\n";
        VLOG(2) << "element_type::numFaces= "    << int( element_type::numFaces ) << "\n";
        VLOG(2) << "element_type::numEdges= "    << int( element_type::numEdges ) << "\n";
        VLOG(2) << "element_type::numVertices= " << int( element_type::numVertices ) << "\n";
        VLOG(2) << "==============================\n";
    }
    else
    {
        VLOG(2) << "==============================\n";
        VLOG(2) << "[initDofMap]\n";
        VLOG(2) << "is_hdiv_conforming     = "  << is_hdiv_conforming << "\n";
        VLOG(2) << "is_hcurl_conforming    = "  << is_hcurl_conforming << "\n";
        VLOG(2) << "nldof                  = "  << int( nldof ) << "\n";
        VLOG(2) << "nFlatLocalDof          = "  << int( nFlatLocalDof ) << "\n";
        VLOG(2) << "fe_type::nLocalDof     = "  << int( fe_type::nLocalDof ) << "\n";
        VLOG(2) << "fe_type::nDofPerVolume = "  << int( fe_type::nDofPerVolume ) << "\n";
        VLOG(2) << "fe_type::nDofPerFace   = "  << int( fe_type::nDofPerFace ) << "\n";
        VLOG(2) << "fe_type::nDofPerEdge   = "  << int( fe_type::nDofPerEdge ) << "\n";
        VLOG(2) << "fe_type::nDofPerVertex = "  << int( fe_type::nDofPerVertex ) << "\n";
        VLOG(2) << "element_type::numVolumes= "  << int( element_type::numVolumes ) << "\n";
        VLOG(2) << "element_type::numFaces= "    << int( element_type::numFaces ) << "\n";
        VLOG(2) << "element_type::numEdges= "    << int( element_type::numEdges ) << "\n";
        VLOG(2) << "element_type::numVertices= " << int( element_type::numVertices ) << "\n";
        VLOG(2) << "==============================\n";
    }

    // initialize the local to global map and fill it with invalid
    // values that will allow to check whether we have a new dof or
    // not when building the table
    const size_type nV = numMeshElements;
    //M_locglob_indices.resize( nV, localglobal_indices_type::Zero( nDofPerElement ) );
    M_locglob_indices.reserve( nV );
    this->initNumberOfDofIdToContainerId( 1 );

    if ( is_hdiv_conforming || is_hcurl_conforming )
    {
        //M_locglob_signs.resize( nV, localglobal_indices_type::Ones( nDofPerElement ) );
        M_locglob_signs.reserve( nV );
        M_locglob_transforms.reserve( nV );
    }
    else
    {
        M_locglob_nosigns = localglobal_indices_type::Ones( runtimeNDofPerElement() );
        M_locglob_notransforms = localglobal_transforms_type( runtimeNDofPerElement() );
    }

    if ( this->hasMeshSupport() && this->meshSupport()->isPartialSupport() )
    {
        //for ( size_type eltId : this->meshSupport()->rangeMeshElementsIdsPartialSupport() )
        for ( auto const& eltWrap : elements(this->meshSupport(), entity_process_t::ALL ) )
        {
            size_type eltId = unwrap_ref( eltWrap ).id();
            M_locglob_indices[eltId] = localglobal_indices_type::Zero( runtimeNDofPerElement() );
            if ( is_hdiv_conforming || is_hcurl_conforming )
            {
                M_locglob_signs[eltId] = localglobal_indices_type::Ones( runtimeNDofPerElement() );
                M_locglob_transforms[eltId] = localglobal_transforms_type( runtimeNDofPerElement() );
            }
        }
    }
    else
    {
        for ( auto const& elt : allelements( M ) )
        {
            size_type eltId = unwrap_ref( elt ).id();
            M_locglob_indices[eltId] = localglobal_indices_type::Zero( runtimeNDofPerElement() );
            if ( is_hdiv_conforming || is_hcurl_conforming )
            {
                M_locglob_signs[eltId] = localglobal_indices_type::Ones( runtimeNDofPerElement() );
                M_locglob_transforms[eltId] = localglobal_transforms_type( runtimeNDofPerElement() );
            }
        }
    }

    constexpr bool doperm = is_hdiv_conforming || is_hcurl_conforming ||
                            ( ( ( Shape == SHAPE_TETRA ) && ( nOrder > 2 ) ) ||
                              ( ( Shape == SHAPE_HEXA ) && ( nOrder > 1 ) ) );
    DVLOG(2) << "generateFacePermutations: " << doperm << "\n";
    generateFacePermutations( M, mpl::bool_<doperm>() );
}
template<typename MeshType, typename FEType, typename MortarType>
void
DofTable<MeshType, FEType, MortarType>::build( mesh_type& M )
{
    tic();
    M_mesh = boost::addressof( M );
    wc( this )->print( fmt::format( "[DofTable::build] starts, has mesh support: {}", this->hasMeshSupport() ),
                       Environment::logVerbosityLevel() > 1, Environment::logVerbosityLevel() > 0, Environment::logVerbosityLevel() > 1 );

#if 0
    if ( this->hasMeshSupport() )
    {
        tic();
        this->meshSupport()->updateParallelData();
#if 0
        this->meshSupport()->updateBoundaryInternalFaces();
#endif
        toc("DofTable::meshSupport", Environment::logVerbosityLevel()>1);
    }
#endif

    tic();
    VLOG(2) << "[Dof::build] initDofMap\n";
    this->initDofMap( M );

    VLOG(2) << "[Dof::build] start building dof map\n";
    size_type start_next_free_dof = 0;
    VLOG(2) << "[Dof::build] start_next_free_dof = " << start_next_free_dof << "\n";
    toc("DofTable::init", Environment::logVerbosityLevel()>1);
    tic();
    if ( is_discontinuous_locally )
    {
        VLOG(2) << "[build] call buildLocallyDiscontinuousDofMap()\n";
        start_next_free_dof = this->buildLocallyDiscontinuousDofMap( M, start_next_free_dof );
        VLOG(2) << "[Dof::build] start_next_free_dof(after local discontinuities) = " << start_next_free_dof << "\n";
    }
    toc("DofTable::buildLocalDiscon", Environment::logVerbosityLevel()>1);
    tic();
    VLOG(2) << "[build] call buildDofMap()\n";
    this->buildDofMap( M, start_next_free_dof );
    //std::cout << "[build] callFINISH buildDofMap() with god rank " << this->worldComm().godRank() <<"\n";
    toc("DofTable::call buildDofMap", Environment::logVerbosityLevel()>1);
    tic();

#if !defined(NDEBUG)
    VLOG(2) << "[build] check that all elements dof were assigned()\n";
    element_const_iterator fit, fen;
    boost::tie( fit, fen ) = M.elementsRange();
    std::vector<boost::tuple<size_type,uint16_type,size_type> > em;

    for ( ; fit != fen; ++fit )
    {
        auto const& elt = fit->second;
        if ( !this->isElementDone( elt.id() ) )
            em.push_back( boost::make_tuple( elt.id(), uint16_type( 0 ), (elt.hasMarker())? elt.marker().value() : 0 ) );
    }
    if ( !em.empty() )
    {
        VLOG(2) << "[build] some element dof were not assigned\n";

        for ( size_type i = 0; i < em.size(); ++i )
        {
            VLOG(3) << " - element " << boost::get<0>( em[i] ) << " c=" << boost::get<1>( em[i] )
                    << " m=" << boost::get<2>( em[i] ) << "\n";
        }
    }

    else
    {
        VLOG(2) << "[build] check that all elements dof were assigned: OK\n";
    }

#endif // NDEBUG
    VLOG(2) << "[Dof::build] n_dof = " << this->nLocalDofWithGhost() << "\n";

    toc("DofTable::checki dof element assignement",Environment::logVerbosityLevel()>1);
    // if ( !is_mortar )
    // {
    //     VLOG(2) << "[build] call buildBoundaryDofMap()\n";
    //     this->buildBoundaryDofMap( M );
    // }

    tic( );
    // multi process
    if ( this->worldComm().localSize()>1 )
    {
        bool isP0continuous = isP0Continuous<fe_type>::result;
        if ( !isP0continuous )
        {
            // add neighbor partition
            this->setNeighborSubdomains(M.neighborSubdomains());

            VLOG(2) << "[build] call buildGhostDofMap () with god rank " << this->worldComm().godRank()  << "\n";
            this->buildGhostDofMap( M );
            VLOG(2) << "[build] callFINISH buildGhostDofMap () with god rank " << this->worldComm().godRank()  << "\n";
        }
        else
        {
            // add all partition as neighbor (if has localdof)
            if ( this->nLocalDofWithGhost() > 0 )
                for ( rank_type proc=0; proc<this->worldComm().localSize(); ++proc )
                    if ( proc!=this->worldComm().rank() && this->nLocalDofWithGhost(proc) > 0 )
                        this->addNeighborSubdomain( proc );

            rank_type themasterRank = 0;
            bool findMasterProc=false;
            uint16_type nDofP0 = (fe_type::is_product)? fe_type::nComponents : 1;
            for ( rank_type proc=0; proc<this->worldComm().localSize(); ++proc )
            {
                if (!findMasterProc && this->nLocalDofWithGhost(proc) > 0)
                {
                    CHECK( nDofP0 == this->nLocalDofWithGhost(proc) ) << "invalid number of dofs" << nDofP0 << " vs " << this->nLocalDofWithGhost(proc);
                    this->M_n_localWithoutGhost_df[proc] = nDofP0;
                    this->M_first_df_globalcluster[proc] = 0;
                    this->M_last_df_globalcluster[proc] = nDofP0-1;
                    themasterRank=proc;
                    findMasterProc=true;
                }
                else
                {
                    this->M_n_localWithoutGhost_df[proc] = 0;
                    this->M_first_df_globalcluster[proc] = 25;// 0;
                    this->M_last_df_globalcluster[proc] = 25; //0;
                }
            }

            if (this->nLocalDofWithGhost() >0 )
            {
                this->M_mapGlobalProcessToGlobalCluster.resize( nDofP0 );
                std::iota( this->M_mapGlobalProcessToGlobalCluster.begin(),
                           this->M_mapGlobalProcessToGlobalCluster.end(),
                           0 );
            }
            this->M_n_dofs = nDofP0;

            if ( themasterRank == this->worldComm().localRank() )
            {
                for ( size_type k = 0; k < this->nLocalDofWithGhost() ; ++k )
                {
                    for ( rank_type proc=0; proc<this->worldComm().localSize(); ++proc )
                    {
                        if ( proc == themasterRank ) continue;
                        if ( this->nLocalDofWithGhost(proc) == 0 ) continue;
                        this->M_activeDofSharedOnCluster[k].insert(proc);
                    }
                }
            }
        }
    }
    else
    {
    toc("DofTable::multi process", Environment::logVerbosityLevel()>1);
    tic();
        // in sequential : identity map
        const size_type s = this->M_n_localWithGhost_df[this->comm().rank()];
        this->M_mapGlobalProcessToGlobalCluster.resize( s );

        std::iota( this->M_mapGlobalProcessToGlobalCluster.begin(),
                   this->M_mapGlobalProcessToGlobalCluster.end(),
                   0 );
    }

    if ( !is_mortar )
    {
        VLOG(2) << "[build] call buildBoundaryDofMap()\n";
        this->buildBoundaryDofMap( M );
    }


    toc("DofTable::sequential map", Environment::logVerbosityLevel()>1);
    tic();
    // reordoring of global process id in doftable (active dofs before and ghost dofs after)
    if ( this->worldComm().localSize()>1 )
    {
        this->updateMultiprocessDofForUse();
    }

    this->initDofIdToContainerIdIdentity( 0,this->nLocalDofWithGhost() );
    toc("DofTable::reordering global id in doftable", Environment::logVerbosityLevel()>1);
    tic();
    EntityProcessType entityProcess = this->hasDofTableExtended()? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
    Range<mesh_type,MESH_ELEMENTS> rangeMeshElt;
    if ( this->hasMeshSupport() )
        rangeMeshElt = elements(this->meshSupport(), entityProcess );
    else
        rangeMeshElt = elements( M, entityProcess );

    for ( auto const& eltWrap : rangeMeshElt )
    {
        auto const& elt = boost::unwrap_ref(eltWrap);
        size_type elid= elt.id();
        if ( is_mortar && elt.isOnBoundary() )
        {
            VLOG(1) << "resizing indices and signs for mortar...";
            auto const& ldof = this->localDof( elid );
            size_type ne = std::distance( ldof.first, ldof.second );
            VLOG(1) << "resizing indices and signs for mortar:  " << ne;
            M_locglob_indices[elid].resize( ne );
            //M_locglob_signs[elid].resize( ne );
        }
        for( auto const& dof: this->localDof( elid ) )
        {
            M_locglob_indices[elid][dof.first.localDof()] = dof.second.index();
            //M_locglob_signs[elid][dof.first.localDof()] = dof.second.sign();
        }
    }
    toc("DofTable::build - locglob indices", Environment::logVerbosityLevel()>1);

    this->buildIndexSplit();

    // build splits with components
    if ( is_product && nRealComponents > 1 )
        this->buildIndexSplitWithComponents( nRealComponents );

    toc("DofTable::build", Environment::logVerbosityLevel()>1);
}

template<typename MeshType, typename FEType, typename MortarType>
typename DofTable<MeshType, FEType, MortarType>::size_type
DofTable<MeshType, FEType, MortarType>::buildLocallyDiscontinuousDofMap( mesh_type& M, size_type start_next_free_dof )
{
    typedef typename continuity_type::template apply<MeshType, self_type> builder;
    return fusion::accumulate( typename continuity_type::discontinuity_markers_type(), start_next_free_dof,  builder( M, *this ) );
}
template<typename MeshType, typename FEType, typename MortarType>
void
DofTable<MeshType, FEType, MortarType>::buildDofMap( mesh_type& M, size_type start_next_free_dof )
{
    wc( this )->print( fmt::format( "[DofTable::buildDofMap] starts, dof_indices empty: {}", M_dof_indices.empty() ), Environment::logVerbosityLevel() > 1, Environment::logVerbosityLevel() > 0, Environment::logVerbosityLevel() > 1 );

    if ( !M_dof_indices.empty() )
    {
        return;
    }

    tic();
    tic();
    const uint16_type dofPerVertex = runtimeDofPerVertex();
    const uint16_type dofPerEdge = runtimeDofPerEdge();
    const uint16_type dofPerFace = runtimeDofPerFace();
    const uint16_type dofPerVolume = runtimeDofPerVolume();
    const uint16_type nLocalDofRt = runtimeNLocalDof();

    size_type legacyEntityLocalDof =
        dofPerVolume * element_type::numVolumes +
        dofPerFace * element_type::numGeometricFaces +
        dofPerEdge * element_type::numEdges +
        dofPerVertex * element_type::numVertices;

    if ( legacyEntityLocalDof != nLocalDofRt )
        VLOG(1) << "[DofTable::buildDofMap] FE-owned local dof count differs from legacy entity arithmetic: "
                << legacyEntityLocalDof << " != " << nLocalDofRt << "\n";

    this->initializeLocalIndexPermutations();
    toc( "DofTable buildDofMap allocation", Environment::logVerbosityLevel() > 1 );
    wc( this )->print( fmt::format( "[DofTable::buildDofMap] allocation done" ), Environment::logVerbosityLevel() > 1, Environment::logVerbosityLevel() > 0, Environment::logVerbosityLevel() > 1 );

    tic();

    // compute the number of dof on current processor
    //Range<MeshType,MESH_ELEMENTS> rangeElements = (this->hasMeshSupport())? elements(this->meshSupport()) : elements(M);
    entity_process_t ept = isP0Continuous<fe_type>::result? entity_process_t::ALL : entity_process_t::LOCAL_ONLY; // WARNING, special case with P0 continuous
    Range<MeshType,MESH_ELEMENTS> rangeElements = (this->hasMeshSupport())? elements(this->meshSupport(),ept) : elements(M,ept);



    auto it_elt = rangeElements.begin();
    auto en_elt = rangeElements.end();
    bool hasNoElt = ( it_elt == en_elt );

    //size_type n_elts = std::distance( it_elt, en_elt);
    wc( this )->print( fmt::format( "[DofTable::buildDofMap]  n_elts =  {} on processor {}", std::distance( it_elt, en_elt ), this->worldComm().localRank() ), Environment::logVerbosityLevel() > 1, Environment::logVerbosityLevel() > 0, Environment::logVerbosityLevel() > 1 );

    size_type theFirstDf = start_next_free_dof;

    if ( is_discontinuous_locally )
        theFirstDf = 0;

    //if ( is_periodic || is_discontinuous_locally )
    //    this->M_first_df[processor] =  0;
    size_type next_free_dof = start_next_free_dof;
    DofFromElement<self_type,fe_type> dfe( this, *M_fe );
    mortar_fe_type mfe;
    if ( nDim == 1 && is_mortar )
        CHECK( mfe.nLocalDof == M_fe->nLocalDof-1 ) << "Invalid number of dof : "
                                                    << " mortar : " << mfe.nLocalDof
                                                    << " fe : " << M_fe->nLocalDof;

    DofFromMortar<self_type,mortar_fe_type,fe_type> dfe_mortar( this, mfe, *M_fe );
    tic();
    for ( ; it_elt!=en_elt; ++it_elt )
    {
        auto const& elt = boost::unwrap_ref( *it_elt );
        //if ( !this->isElementDone( it_elt->id() ) )
        {
            if ( is_mortar )
            {

                if ( !elt.isOnBoundary() )
                {
                    VLOG(1) << "add standard element " << elt.id() << " ndof : " << M_fe->nLocalDof;
                    dfe.add( elt, next_free_dof, this->worldComm().localRank() );
                }
                else
                {
                    VLOG(1) << "add mortar element " << elt.id() << " ndof : " << mfe.nLocalDof;
                    dfe_mortar.add( elt, next_free_dof, this->worldComm().localRank() );
                }
            }
            else
            {
                dfe.add( elt, next_free_dof, this->worldComm().localRank() );
            }
        }
    } // elements loop
    toc("DofTable buildDofMap element loop", Environment::logVerbosityLevel()>1);
    // update extended doftable for P0 continuous
#if 0
    if ( isP0Continuous<fe_type>::result && this->hasDofTableExtended() )
    {
        for (auto const& ghostEltWrap : elements(M,EntityProcessType::GHOST_ONLY ) )
        {
            auto const& ghostElt = boost::unwrap_ref( ghostEltWrap );
            dfe.add( ghostElt, next_free_dof, this->worldComm().localRank() );
        }
    }
#endif

    toc( "DofTable buildDofMap dof generation", Environment::logVerbosityLevel() > 1 );
    size_type mynDofWithGhost = next_free_dof;//next_free_dof - start_next_free_dof;

    //const size_type thelastDof = ( !hasNoElt )?next_free_dof-1:0;
    const rank_type myrank = this->worldComm().localRank();
    wc( this )->print( fmt::format( "[builddofmap - {}] dof generation nLocalDof : {}",
                                    rank(this), mynDofWithGhost ), Environment::logVerbosityLevel() > 1, Environment::logVerbosityLevel() > 0, Environment::logVerbosityLevel() > 1 );
    tic();
#if 0
    for ( auto mit = M_dof_marker.right.begin(), men = M_dof_marker.right.end() ; mit != men ; ++mit )
    {
        LOG(INFO) << "marker " << mit->first << " dof id " << mit->second;
    }
#endif

#if 0
    LOG(INFO) << "local to global view";
    for( auto it = M_el_l2g.left.begin(), en = M_el_l2g.left.end();
         it != en; ++it )
    {
        LOG(INFO) << "local dof (" << it->first.elementId()<< ","<< it->first.localDof() << ") --> global dof " << it->second.index();
    }
    LOG(INFO) << "global to local  view";
    for( auto it = M_el_l2g.right.begin(), en = M_el_l2g.right.end();
         it != en; ++it )
    {
        LOG(INFO) << "global dof " << it->first.index() << " --> local dof (" << it->second.elementId()<< ","<< it->second.localDof() << ")";
    }
#endif


    if ( isP0Continuous<fe_type>::result /*|| !is_continuous*/ )
    {
        tic();
        mpi::all_gather( this->worldComm().localComm(),
                         mynDofWithGhost,
                         this->M_n_localWithGhost_df );
        toc("DofTable buildDofMap all_gather", Environment::logVerbosityLevel()>1);
    }
    else
    {
        // up only with myrank (completed in buildGhostDofMap)
        this->M_n_localWithGhost_df[myrank] = mynDofWithGhost;
    }

#if 0
    std::cout << "\n build Dof Map --2---with god rank " << this->worldComm().godRank()
              << " local rank DofT " << this->worldComm().localRank()
              << " local rank mesh " << M.worldComm().localRank()
              << std::endl;
    LOG( INFO ) << fmt::format( "[builddofmap] localrank {}", this->worldComm().localRank() ) << std::endl;
    #if !defined(FEELPP_HAS_SPDLOG)

    google::FlushLogFiles(google::INFO);

    #else

    Logger::flush();

    #endif
#endif

    // only true in sequential, redefine in buildDofGhostMap
    this->M_n_localWithoutGhost_df[myrank] = this->M_n_localWithGhost_df[myrank];
    this->M_first_df_globalcluster[myrank] = this->firstDof(); //this->M_first_df[myrank];
    this->M_last_df_globalcluster[myrank] = this->lastDof(); //this->M_last_df[myrank];
    this->M_n_dofs = next_free_dof;

#if 0
    it_elt = M.beginElementWithProcessId();

    for ( ; it_elt != en_elt; ++it_elt )
    {
        size_type elid= it_elt->id();
        if ( is_mortar && it_elt->isOnBoundary() )
        {
            VLOG(1) << "resizing indices and signs for mortar...";
            auto const& ldof = this->localDof( elid );
            size_type ne = std::distance( ldof.first, ldof.second );
            VLOG(1) << "resizing indices and signs for mortar:  " << ne;
            M_locglob_indices[elid].resize( ne );
            //M_locglob_signs[elid].resize( ne );
            for( auto const& dof: this->localDof( elid ) )
            {
                M_locglob_indices[elid][dof.first.localDof()] = dof.second.index();
                //M_locglob_signs[elid][dof.first.localDof()] = dof.second.sign();
            }
        }
        else
            for ( int i = 0; i < FEType::nLocalDof; ++i )
            {
                int nc1 = ( is_product?nComponents:1 );

                for ( int c1 =0; c1 < nc1; ++c1 )
                {
                    int ind = FEType::nLocalDof*c1+i;
                    auto const& dof = localToGlobal( elid, i, c1 );
                    M_locglob_indices[elid][ind] = dof.index();
                    //M_locglob_signs[elid][ind] = dof.sign();
                }
            }
    }
#endif
    tic();
    // the dof points are necessary to build the parallel dof table
    // if ( this->worldComm().localSize() > 1 )
    //     this->generateDofPoints( M, true );
    toc("DofTable generateDofPoints", Environment::logVerbosityLevel()>1);

    toc( "DofTable buildDofMap done", Environment::logVerbosityLevel()>1);
}

template<typename MeshType, typename FEType, typename MortarType>
void
DofTable<MeshType, FEType, MortarType>::buildBoundaryDofMap( mesh_type& mesh )
{
    tic();
    size_type nDofF = nLocalDofOnFace(true);
    size_type nDofFFlat = nLocalDofOnFace();
    M_n_dof_per_face_on_bdy = nDofF;
    DVLOG(2) << "vertex dof : " <<  face_type::numVertices * runtimeDofPerVertex() << "\n";
    DVLOG(2) << "edge dof : " <<  face_type::numEdges * runtimeDofPerEdge() << "\n";
    DVLOG(2) << "face dof : " << face_type::numFaces * runtimeDofPerFace()  << "\n";
    DVLOG(2) << "number of Dof on an Element Face : " << nDofF << "\n";

    if ( nDofFFlat == 0 ) return;

    // Face dof
    DofFromBoundary<self_type, fe_type> dfb( this, *M_fe );

    auto rangeFaces = this->hasMeshSupport() && this->meshSupport()->isPartialSupport()?
        faces( this->meshSupport(), entity_process_t::ALL ) : faces( mesh, entity_process_t::ALL );
    std::map<size_type, std::map<rank_type,size_type>> isolatedFaces;
    for ( auto const& faceWrap : rangeFaces )
    {
        auto const& face = unwrap_ref( faceWrap );
        LOG_IF(WARNING, !face.isConnectedTo0() )
            << "face " << face.id() << " not connected"
            << " hasMarker : " << face.hasMarker()
            << " connectedTo0 : " << face.isConnectedTo0()
            << " connectedTo1 : " << face.isConnectedTo1();

        if ( !face.isConnectedTo0() ) continue;

#if !defined(NDEBUG)
        if (  face.isOnBoundary() )
        {
            DVLOG(4) << "[buildBoundaryDofMap] boundary global face id : " << face.id()
                     << " hasMarker: " << face.hasMarker()<< "\n";
        }
        else
        {
            DVLOG(4) << "[buildBoundaryDofMap] global face id : " << face.id() << "\n";
        }
#endif
        M_face_l2g[ face.id()].resize( nDofFFlat );
        if ( !dfb.add( face ) )
            isolatedFaces.emplace( face.id(), face.idInOthersPartitions() );
    }

    //DCHECK( isolatedFaces.empty() ) << "TODO: finish implementation of this case below";
    LOG_IF(WARNING, isolatedFaces.empty() ) << "TODO: finish implementation of this case below";
    if ( false && this->worldComm().localSize()>1 && this->hasMeshSupport() && this->meshSupport()->isPartialSupport() )
    {
        int nbMaxRequest = 2*mesh.neighborSubdomains().size();
        std::vector<mpi::request> reqs( nbMaxRequest );
        int countRequest = 0;

        std::map<rank_type, std::vector<size_type> > dataToSend;
        std::map<rank_type, std::vector<size_type> > dataToRecv;

        for ( auto const& [faceId,mapProcessToFaceId] : isolatedFaces )
        {
            std::cout << "mapProcessToFaceId.size:" << mapProcessToFaceId.size() << std::endl;
            for ( auto const& [rank,faceIdOtherProcess] : mapProcessToFaceId )
                dataToSend[rank].push_back(faceIdOtherProcess);
        }

        // get size of data to transfer
        std::map<rank_type,std::size_t> sizeRecv;
        std::map<rank_type,std::size_t> sizeSend;
        for ( rank_type neighborRank : mesh.neighborSubdomains() )
        {
            sizeSend[neighborRank] = dataToSend[neighborRank].size();
            reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank, 0, sizeSend[neighborRank] );
            reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank, 0, sizeRecv[neighborRank] );
        }
        // wait all requests
        mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
        countRequest = 0;


        // step 1 :send/recv of data
        for ( rank_type neighborRank : mesh.neighborSubdomains() )
        {
            std::size_t nSendData = dataToSend[neighborRank].size();
            if ( nSendData > 0 )
                reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank , 0, dataToSend[neighborRank].data(), nSendData );
            std::size_t nRecvData = sizeRecv[neighborRank];
            dataToRecv[neighborRank].resize( nRecvData );
            if ( nRecvData > 0 )
                reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank , 0, dataToRecv[neighborRank].data(), nRecvData );
        }
        // step 1 :wait all requests
        mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
        countRequest = 0;


        // step2 : send globl face from active elts, prepare mpi data of subentities required (from ghost elts)
        std::map< rank_type, std::vector<std::vector<size_type>> > dataToSendStep2, dataToRecvStep2;
        for ( auto const& [rank,faceIds] : dataToRecv )
        {
            dataToSendStep2[rank].reserve( faceIds.size() );
            for ( auto const& faceId : faceIds )
            {
                std::vector<size_type> ind;
                auto eit = M_face_l2g.find( faceId );
                if ( eit != M_face_l2g.end() )
                {
                    ind.reserve( nLocalDofOnFace() );
                    std::for_each( eit->second.begin(), eit->second.end(),
                                   [this,&ind]( FaceDof<size_type> const& f ) { ind.push_back( this->mapGlobalProcessToGlobalCluster().at( f.index() ) ); } );
                }
                dataToSendStep2[rank].push_back( std::move(ind) );
            }
            std::cout << "dataToSendStep2[rank].size:"<<dataToSendStep2[rank].size()<<std::endl;
        }
        // step 2 :send/recv of data
        for ( rank_type neighborRank : mesh.neighborSubdomains() )
        {
            std::size_t nSendData = sizeRecv[neighborRank]; // dataToSendStep2[neighborRank].size();
            if ( nSendData > 0 )
                reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank , 0, dataToSendStep2[neighborRank].data(), nSendData );
            std::size_t nRecvData = sizeSend[neighborRank];
            dataToRecvStep2[neighborRank].resize( nRecvData );
            if ( nRecvData > 0 )
                reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank , 0, dataToRecvStep2[neighborRank].data(), nRecvData );
        }
        // step 2 :wait all requests
        mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
        countRequest = 0;


        // step 3 : fetch info about ghost dofs required
        const rank_type myRank = this->worldComm().localRank();
        const rank_type nProc = this->worldComm().localSize();
        size_type start_next_free_dof = this->M_n_localWithGhost_df[myRank];
        size_type next_free_dof = start_next_free_dof;
        std::map<size_type,size_type> dofGlobalClusterToGlobalProcess;
        for ( auto const& [rank,faceData] : dataToRecvStep2 )
        {
            for ( uint16_type k=0;k<faceData.size();++k )
            {
                std::cout << "use faceData.size:" << faceData.size()<<std::endl;
                if ( faceData[k].empty() )
                    continue;

                for ( uint16_type q=0;q<faceData[k].size();++q )
                {
                    size_type gc = faceData[k][q];
                    auto itFind = dofGlobalClusterToGlobalProcess.find( gc );
                    if ( itFind == dofGlobalClusterToGlobalProcess.end() )
                        dofGlobalClusterToGlobalProcess.emplace( gc, next_free_dof++ );
                }
            }
        }

        //------------------------------------------------------------------------------//
        // step 3 : update local datamap
#if 0 // TODO
        this->M_nGhostDofAddedInExtendedDofTable = next_free_dof-start_next_free_dof;
        std::vector<size_type> dataRecvFromGather;
        mpi::all_gather( this->worldComm().localComm(),
                         this->M_nGhostDofAddedInExtendedDofTable,
                         dataRecvFromGather );
        for (rank_type p=0;p<nProc;++p)
        {
            this->M_last_df[p] += dataRecvFromGather[p];
            this->M_n_localWithGhost_df[p] += dataRecvFromGather[p];
        }
        this->M_mapGlobalProcessToGlobalCluster.resize( this->M_n_localWithGhost_df[myRank],invalid_v<size_type> );


        //
        for ( auto const& [dofGlobalClusterIndex,dofIndex] : dofGlobalClusterToGlobalProcess )
        {
            std::cout << "add dofIndex:" << dofIndex << " dofGlobalClusterIndex:"<<dofGlobalClusterIndex<<std::endl;
            this->M_mapGlobalProcessToGlobalCluster[dofIndex] = dofGlobalClusterIndex;
        }
#endif

        // TODO face mapping
    }


    toc( "DofTable::buildBoundaryDofMap", Environment::logVerbosityLevel()>1 );

}    // updateBoundaryDof


template<typename MeshType, typename FEType, typename MortarType>
void
DofTable<MeshType, FEType, MortarType>::updateMultiprocessDofForUse()
{
    size_type _nLocalDofWithGhost = this->nLocalDofWithGhost();
    size_type _nLocalDofWithoutGhost = this->nLocalDofWithoutGhost();

    //! reordering of global dof : actives first, then ghosts
    std::vector<size_type> previousGlobalIdToNewGlobalId( _nLocalDofWithGhost );
    size_type currentActiveDof=0,currentGhostDof=_nLocalDofWithoutGhost;
    std::vector<size_type> newMapGlobalProcessToGlobalCluster( _nLocalDofWithGhost );
    size_type firstGlobIndex = this->firstDofGlobalCluster();

    for ( size_type k=0;k<_nLocalDofWithGhost;++k )
    {
        size_type gcdof = this->M_mapGlobalProcessToGlobalCluster[k];
        if ( this->dofGlobalProcessIsGhost(k) )
            previousGlobalIdToNewGlobalId[k]=currentGhostDof++;
        else
            previousGlobalIdToNewGlobalId[k]=currentActiveDof++;
        DCHECK( previousGlobalIdToNewGlobalId[k] < newMapGlobalProcessToGlobalCluster.size() )
            << fmt::format("index out of range : id: {} vs  size:{}  isghost:{}",previousGlobalIdToNewGlobalId[k], newMapGlobalProcessToGlobalCluster.size(), this->dofGlobalProcessIsGhost(k) );
        newMapGlobalProcessToGlobalCluster[previousGlobalIdToNewGlobalId[k]] = gcdof;
    }
    this->M_mapGlobalProcessToGlobalCluster = std::move( newMapGlobalProcessToGlobalCluster );
    this->updateWorldIndexForUse();

    for( auto it = M_el_l2g.left.begin(), en = M_el_l2g.left.end(); it != en; ++it )
    {
        auto const& previousGDof=it->second;
        Dof newGDof( previousGDof );
        CHECK( previousGDof.index() < previousGlobalIdToNewGlobalId.size() ) << fmt::format("index out of range index: {}  size: {}",
                                                                                            previousGDof.index(), previousGlobalIdToNewGlobalId.size() );
        newGDof.setIndex( previousGlobalIdToNewGlobalId[previousGDof.index()] );
        bool successfulModify = M_el_l2g.left.modify_data( it, boost::bimaps::_data = newGDof );
        CHECK( successfulModify ) << "modify global dof id fails";
    }
    for ( auto & faceDataElt : M_face_l2g )
        for ( FaceDof<size_type> & faceDataDof : faceDataElt.second )
            faceDataDof.setIndex( previousGlobalIdToNewGlobalId[faceDataDof.index()] );

    dof_points_type newDofPoints;
    for ( auto const& dofPt : M_dof_points )
    {
        size_type newDofId = previousGlobalIdToNewGlobalId[ dofPt.first ];
        auto const& dofPtData = dofPt.second;
        newDofPoints[newDofId] = boost::make_tuple( boost::get<0>( dofPtData ),newDofId,boost::get<2>( dofPtData ) );
    }
    M_dof_points.clear();
    M_dof_points.swap( newDofPoints );

    dof_marker_type newDofMarker;
    for ( auto it = M_dof_marker.left.begin(), en = M_dof_marker.left.end(); it != en; ++it )
        newDofMarker.insert( dof2marker(previousGlobalIdToNewGlobalId[it->first],it->second) );
    M_dof_marker.clear();
    M_dof_marker.swap( newDofMarker );

    auto reorderActiveDofSharedOnCluster = [&previousGlobalIdToNewGlobalId,this]()
    {
        std::map<size_type, std::set<rank_type> > newActiveDofSharedOnCluster;
        for ( auto const& activeDof : this->M_activeDofSharedOnCluster )
        {
            DCHECK( activeDof.first < previousGlobalIdToNewGlobalId.size() ) << fmt::format("activeDof.first {} vs size{}",activeDof.first,previousGlobalIdToNewGlobalId.size());
            newActiveDofSharedOnCluster.emplace( std::make_pair( previousGlobalIdToNewGlobalId[activeDof.first], activeDof.second ) );
        }
        this->M_activeDofSharedOnCluster = std::move( newActiveDofSharedOnCluster );
    };

    // ---------------------------------
    // update activeDofSharedOnCluster
    if constexpr ( isP0Continuous<fe_type>::result )
    {
        // in that case, activeDofSharedOnCluster is already built, just apply reordering
        reorderActiveDofSharedOnCluster();
    }
    else if ( this->M_hasDescriptorKeyClusterDofMap )
    {
        reorderActiveDofSharedOnCluster();
    }
    else
    {
        // clear
        this->M_activeDofSharedOnCluster.clear();

        const rank_type myRank = this->worldComm().localRank();
        const rank_type nProc = this->worldComm().localSize();

        int nbMaxRequest = 2*this->neighborSubdomains().size();
        std::vector<mpi::request> reqs( nbMaxRequest );
        int countRequest = 0;

        // send global process cluster
        std::map<rank_type, std::vector<size_type> > dataToSend;
        std::map<rank_type, std::vector<size_type> > dataToRecv;

        for ( size_type k=_nLocalDofWithoutGhost;k<_nLocalDofWithGhost;++k )
        {
            size_type gcdof = this->M_mapGlobalProcessToGlobalCluster[k];
            rank_type activeProcId = this->procOnGlobalCluster( gcdof );
            DCHECK( activeProcId != myRank ) << "should be a ghost dof";
            dataToSend[activeProcId].push_back( gcdof );
        }

        std::map<size_type,size_type> mapActiveGcToGp;
        for ( size_type k=0;k<_nLocalDofWithoutGhost;++k )
            mapActiveGcToGp.emplace( this->M_mapGlobalProcessToGlobalCluster[k], k );

        // get size of data to transfer
        std::map<rank_type,std::size_t> sizeRecv;
        std::map<rank_type,std::size_t> sizeSend;
        for ( rank_type neighborRank : this->neighborSubdomains() )
        {
            sizeSend[neighborRank] = dataToSend[neighborRank].size();
            reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank, 0, sizeSend[neighborRank] );
            reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank, 0, sizeRecv[neighborRank] );
        }
        // wait all requests
        mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
        countRequest = 0;

        // step 1 :send/recv of data
        for ( rank_type neighborRank : this->neighborSubdomains() )
        {
            std::size_t nSendData = dataToSend[neighborRank].size();
            if ( nSendData > 0 )
                reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank , 0, dataToSend[neighborRank].data(), nSendData );
            std::size_t nRecvData = sizeRecv[neighborRank];
            dataToRecv[neighborRank].resize( nRecvData );
            if ( nRecvData > 0 )
                reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank , 0, dataToRecv[neighborRank].data(), nRecvData );
        }
        // step 1 :wait all requests
        mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
        countRequest = 0;

        for ( auto const& [rank,gcDofs] : dataToRecv )
        {
            for ( size_type gcDofIndex : gcDofs )
            {
                //this->addNeighborSubdomain( rank );
                auto itDofIndex = mapActiveGcToGp.find( gcDofIndex );
                CHECK( itDofIndex != mapActiveGcToGp.end() )
                    << fmt::format( "missing active gc dof {} received from rank {} (active map size={}, nLocalWithoutGhost={})",
                                    gcDofIndex, rank, mapActiveGcToGp.size(), _nLocalDofWithoutGhost );
                size_type dofIndex = itDofIndex->second;
                this->M_activeDofSharedOnCluster[dofIndex].insert(rank);
            }
        }
    } // !isP0continuous
}


template<typename MeshType, typename FEType, typename MortarType>
void
DofTable<MeshType, FEType, MortarType>::generateDofPoints(  mesh_type& M, bool __buildMinimalParallel/*, mpl::bool_<false>*/ ) const
{
    if ( M_hasBuiltDofPoints )// !M_dof_points.empty() )
        return;

    if ( fe_type::is_modal )
    {
        M_hasBuiltDofPoints = true;
        return;
    }

    DVLOG(2) << "[Dof::generateDofPoints] generating dof coordinates\n";

    auto rangeElements = (this->hasMeshSupport())? elements( this->meshSupport(), entity_process_t::ALL ) : elements( M, entity_process_t::ALL );
    auto it_elt = rangeElements.begin();
    auto en_elt = rangeElements.end();

    if ( it_elt == en_elt )
    {
        M_hasBuiltDofPoints = true;
        return;
    }

    auto const& fe = this->fe();
    auto gm = M.gm();
    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    typename gm_type::precompute_ptrtype __geopc( new typename gm_type::precompute_type( gm, fe.points() ) );

    using gm_context_type = typename gm_type::template Context<element_type>;
    using gm_context_ptrtype = std::shared_ptr<gm_context_type>;
    gm_context_ptrtype ctx = gm->template context<vm::POINT>( unwrap_ref( *it_elt ), __geopc );
    gm_context_ptrtype mctx;
    if constexpr( is_mortar )
    {
        mortar_fe_type mfe;
        typename gm_type::precompute_ptrtype __mgeopc( new typename gm_type::precompute_type( gm, mfe.points() ) );
        mctx = gm->template context<vm::POINT>( unwrap_ref( *it_elt ), __mgeopc );
    }

    gm_context_ptrtype ctxCurrent;

    for ( size_type dof_id = 0; it_elt!=en_elt ; ++it_elt )
    {
        auto const& elt = boost::unwrap_ref( *it_elt );

        if constexpr( is_mortar )
        {
            if ( elt.isOnBoundary() )
                ctxCurrent = mctx;
            else
                ctxCurrent = ctx;
        }
        else
            ctxCurrent = ctx;

        ctxCurrent->template update<vm::POINT>( elt );

        for ( auto const& ldof : this->localDof( elt.id() ) )
        {
            size_type thedof = ldof.second.index();
            uint16_type ldofId = ldof.first.localDof();
            if ( !this->localDofHasRepresentativePoint( ldofId ) )
                continue;

            uint16_type pointId = this->localDofRepresentativePointIndex( ldofId );
            FEELPP_ASSERT( pointId < static_cast<uint16_type>( fe.points().size2() ) )
                ( int( pointId ) )( int( fe.points().size2() ) )( int( ldofId ) )
                .error( "invalid FE representative point index" );
            if ( ( thedof >= this->firstDof() ) && ( thedof <= this->lastDof() ) )
            {
                DCHECK( thedof < this->nLocalDofWithGhost() )
                    << "invalid local dof index "
                    <<  thedof << ", " << this->nLocalDofWithGhost() << "," << this->firstDof()  << ","
                    <<  this->lastDof() << "," << elt.id() << "," << ldofId << "," << pointId;

                if ( M_dof_points.find( thedof ) == M_dof_points.end() )
                {
                    uint16_type c1 = fe.component( ldofId );
                    M_dof_points[thedof] = boost::make_tuple( ctxCurrent->xReal( pointId ), thedof, c1 );
                }
#if !defined( NDEBUG )
                else if ( !isP0Continuous<fe_type>::result && !is_mortar )
                {
                    auto dofpointFromGmc = ctxCurrent->xReal( pointId );
                    auto dofpointStored = M_dof_points[thedof].template get<0>();
                    bool find2=true;
                    for (uint16_type d=0;d< nRealDim;++d)
                    {
                        find2 = find2 && (std::abs( dofpointFromGmc[d]-dofpointStored[d] )<1e-9);
                    }
                    CHECK(find2) << " error localToGlobal for "<< pointId <<" with " << dofpointFromGmc << " and " << dofpointStored <<"\n" ;
                }
#endif
            }

        }
    }

    M_hasBuiltDofPoints = true;
#if !defined( NDEBUG )
        if ( !hasDofTableExtended() )
            for ( auto const& dofPointEntry : M_dof_points )
            {
                auto const dof_id = dofPointEntry.first;
                auto const& dof_point = dofPointEntry.second;
                CHECK( boost::get<1>( dof_point ) >= this->firstDof() &&
                       boost::get<1>( dof_point ) <= this->lastDof() )
                    <<  "invalid dof point "
                    <<  dof_id << ", " <<  this->firstDof() << ", " <<  this->lastDof() << ", " <<  this->nLocalDofWithGhost()
                    << ", " << boost::get<1>( dof_point )
                    << ", " <<  boost::get<0>( dof_point ) ;
            }
#endif
    DVLOG(2) << "[Dof::generateDofPoints] generating dof coordinates done\n";
}
template<typename MeshType, typename FEType, typename MortarType>
void
DofTable<MeshType, FEType, MortarType>::addSubstructuringDofMap( mesh_type const& M, size_type next_free_dof )
{
    addSubstructuringDofVertex( M, next_free_dof );
    addSubstructuringDofEdge( M, next_free_dof, mpl::int_<nDim>() );
    addSubstructuringDofFace( M, next_free_dof, mpl::int_<nDim>() );
}

template<typename MeshType, typename FEType, typename MortarType>
void
DofTable<MeshType, FEType, MortarType>::addSubstructuringDofVertex(mesh_type const& M,
                                                                                    size_type next_free_dof )
{
    std::cout << "found CrossPoints and WireBasket\n";
    std::cout << "n cp: " << std::distance( M.beginPointWithMarker( M.markerName("CrossPoints") ), M.endPointWithMarker( M.markerName("CrossPoints") ) ) << "\n";
#if 0
    std::cout << "n wb: " << std::distance( M.beginEdgeWithMarker( M.markerName("WireBasket") ), M.endEdgeWithMarker( M.markerName("WireBasket") ) ) << "\n";
#endif
    // go through all the crosspoints and add them to the dof table

    for( auto pit = M.beginPointWithMarker( M.markerName("CrossPoints") ),
             pen = M.endPointWithMarker( M.markerName("CrossPoints") );
         pit!=pen; ++pit )
    {
        const uint16_type nDofPerVertexRt = runtimeDofPerVertex();
        if ( nDofPerVertexRt == 0 )
            continue;

        // get one element
        auto __elt = M.element( *pit->elements().begin() );
        size_type ie = __elt.id();
        int lc = 0;
        for ( uint16_type i = 0; i < element_type::numVertices; ++i )
        {
            for ( uint16_type l = 0; l < nDofPerVertexRt; ++l, ++lc )
            {
                if (__elt.point( i ).id()==pit->id() )
                {
                    const size_type gDof = ( __elt.point( i ).id() ) * nDofPerVertexRt + l;
                    this->insertDof( ie, lc, i, boost::make_tuple( 0, 0, gDof ),
                                     M.worldComm().localRank(), next_free_dof, 1, false, 0 );
                    std::cout << "Adding crosspoint " << pit->id() << " with dof " << next_free_dof << "\n";
                }
            }
        }
    }
}

template<typename MeshType, typename FEType, typename MortarType>
void
DofTable<MeshType, FEType, MortarType>::addSubstructuringDofEdge( mesh_type const& M,
                                                                                   size_type next_free_dof,
                                                                                   mpl::int_<1> )
{}

template<typename MeshType, typename FEType, typename MortarType>
void
DofTable<MeshType, FEType, MortarType>::addSubstructuringDofEdge( mesh_type const& M,
                                                                                   size_type next_free_dof,
                                                                                   mpl::int_<2> )
{}

template<typename MeshType, typename FEType, typename MortarType>
void
DofTable<MeshType, FEType, MortarType>::addSubstructuringDofEdge( mesh_type const& M,
                                                                                   size_type next_free_dof,
                                                                                   mpl::int_<3> )
{
    const uint16_type nDofPerEdgeRt = runtimeDofPerEdge();
    if ( nDofPerEdgeRt == 0 )
        return;

    // go through all Wirebasket edges
    for( auto pit = M.beginEdgeWithMarker( M.markerName("WireBasket") ),
             pen = M.endEdgeWithMarker( M.markerName("WireBasket") );
         pit!=pen; ++pit )
    {
        auto __elt = M.element( *pit->elements().begin() );
        std::cout << "Adding wirebasket edge " << pit->id() << " using element "  << __elt.id() << "\n";
        size_type ie = __elt.id();
        uint16_type lc = 0;

        for ( uint16_type i = 0; i < element_type::numEdges; ++i )
        {
            for ( uint16_type l = 0; l < nDofPerEdgeRt; ++l, ++lc )
            {
                if (__elt.edge( i ).id()==pit->id() )
                {
                    size_type gDof = __elt.edge( i ).id() * nDofPerEdgeRt;
                    int32_type sign = 1;

                    if ( __elt.edgePermutation( i ).value()  == edge_permutation_type::IDENTITY )
                    {
                        gDof += l ; // both nodal and modal case
                    }
                    else if ( __elt.edgePermutation( i ).value()  == edge_permutation_type::REVERSE_PERMUTATION )
                    {

                        if ( fe_type::is_modal )
                        {
                            //only half of the modes (odd polynomial order) are negative.
                            sign = ( l%2 )?( -1 ):( 1 );
                            gDof += l;
                        }

                        else
                            gDof += nDofPerEdgeRt - 1 - l ;
                    }
                    else
                        FEELPP_ASSERT( 0 ).error ( "invalid edge permutation" );

                    this->insertDof( ie, lc, i, boost::make_tuple( 1, 0, gDof ), M.worldComm().localRank(), next_free_dof, sign, false, 0 );
                    std::cout << "Adding wirebasket edge " << pit->id() << " with dof " << next_free_dof << "\n";
                }
            }
        }

    }
}
template<typename MeshType, typename FEType, typename MortarType>
void
DofTable<MeshType, FEType, MortarType>::addSubstructuringDofFace( mesh_type const& M,
                                                                                   size_type next_free_dof,
                                                                                   mpl::int_<1> )
{}

template<typename MeshType, typename FEType, typename MortarType>
void
DofTable<MeshType, FEType, MortarType>::addSubstructuringDofFace( mesh_type const& M,
                                                                                   size_type next_free_dof,
                                                                                   mpl::int_<2> )
{}

template<typename MeshType, typename FEType, typename MortarType>
void
DofTable<MeshType, FEType, MortarType>::addSubstructuringDofFace( mesh_type const& M,
                                                                                   size_type next_free_dof,
                                                                                   mpl::int_<3> )
{
    const uint16_type nDofPerFaceRt = runtimeDofPerFace();
    if ( nDofPerFaceRt == 0 )
        return;

    std::vector<std::string> faces = assign::list_of("TOP")("BOTTOM")("NORTH")("EAST")("WEST")("SOUTH");
    for( auto face : faces )
    {
        auto faces = markedfaces( &M, face );

        for( auto pit = faces.template get<1>(), pen = faces.template get<2>(); pit!=pen; ++pit )
        {
            auto __elt = M.element( *pit->elements().begin() );
            std::cout << "Adding face " << pit->id() << " with marker " << face << " using element "  << __elt.id() << "\n";
            size_type ie = __elt.id();
            uint16_type lc = 0;

            for ( uint16_type i = 0; i < element_type::numFaces; ++i )
            {
                face_permutation_type permutation = __elt.facePermutation( i );
                FEELPP_ASSERT( permutation != face_permutation_type( 0 ) ).error ( "invalid face permutation" );

                // Polynomial order in each direction
                uint16_type p=1;
                uint16_type q=0;

                // MaxOrder = Order - 2
                int MaxOrder = int( ( 3 + std::sqrt( 1 + 8 * nDofPerFaceRt ) ) / 2 ) - 2;

                for ( uint16_type l = 0; l < nDofPerFaceRt; ++l, ++lc )
                {
                    if (__elt.face( i ).id()==pit->id() )
                    {
                        // TODO: orient the dof indices such
                        // that they match properly the faces
                        // dof of the connected faces. There
                        // are a priori many permutations of
                        // the dof face indices
                        size_type gDof = __elt.face( i ).id() * nDofPerFaceRt;
                        int32_type sign = 1;

                        q=q+1;

                        if ( q > MaxOrder )
                        {
                            q = 1;
                            p = p+1;
                            MaxOrder = MaxOrder-1;
                        }

                        if ( !fe_type::is_modal )
                        {
                            // no need of permutation if identity or only one dof on face
                            if ( permutation == face_permutation_type( face_permutation_type::IDENTITY ) || nDofPerFaceRt == 1 )
                                gDof += l;
                            else
                            {
                                auto const& perm = facePermutationVector( permutation, nDofPerFaceRt );
                                gDof += perm( l );
                            }
                        }

                        else
                        {
                            gDof += l;

                            if ( permutation == face_permutation_type( 2 ) )
                            {
                                // Reverse sign if polynomial order in
                                // eta_1 direction is odd

                                if ( p%2 == 0 )
                                    sign = -1;

                            }
                        }

                        this->insertDof( ie, lc, i, boost::make_tuple( 2, 0, gDof ), M.worldComm().localRank(), next_free_dof, sign, false, 0 );
                        std::cout << "Adding face " << pit->id() << " with dof " << next_free_dof << "\n";
                    }
                }
            }

        }
    }

}

template<typename MeshType, typename FEType, typename MortarType>
typename DofTable<MeshType, FEType, MortarType>::pidtodofid_type
DofTable<MeshType, FEType, MortarType>::pointIdToDofRelation(std::string fname, bool dof2pid, bool pid2dof ) const
{
    std::unordered_map<size_type,size_type> pidtodof,doftopid;
    auto rangeElements = M_mesh->elementsWithProcessId( M_mesh->worldComm().localRank() );
    auto it_elt = std::get<0>( rangeElements );
    auto en_elt = std::get<1>( rangeElements );

    if ( it_elt == en_elt )
        return std::make_pair(doftopid,pidtodof);
    int ncdof  = is_product?nComponents:1;

    if ( dof2pid )
        doftopid.reserve( this->nLocalDof() );
    if ( pid2dof )
        pidtodof.reserve( ncdof*std::distance( it_elt, en_elt )*M_mesh->numLocalVertices() );
    for ( size_type dof_id = 0; it_elt!=en_elt ; ++it_elt )
    {
        auto const& elt = boost::unwrap_ref( *it_elt );
        for ( uint16_type i = 0; i < M_mesh->numLocalVertices(); ++i )
        {
            for ( uint16_type c1 = 0; c1 < ncdof; ++c1 )
            {
                const size_type gDof = ( elt.point( i ).id() );
                size_type thedof = localToGlobal( elt.id(), i, c1 ).index();
                //pidtodof[ncdof*it_elt->point(l).id()+c1] = thedof;
                if ( pid2dof )
                    pidtodof[ncdof*gDof+c1] = thedof;
                if ( dof2pid )
                    doftopid[thedof] = ncdof*gDof+c1;

            }
        }
    }
    if ( !fname.empty() )
    {
        std::ostringstream os1,os2;
        os1 << fs::path( fname ).stem().string() << "_pidtodof" << fs::path( fname ).extension().string();
        os2 << fs::path( fname ).stem().string() << "_doftopid" << fs::path( fname ).extension().string();
        if ( pid2dof )
        {
            std::ofstream ofs( os1.str().c_str() );
            auto it = pidtodof.begin();
            auto en = pidtodof.end();
            std::for_each( it, en,
                           [&ofs]( std::pair<size_type, size_type> const& p )
                               {
                                   ofs << p.first << " " << p.second << "\n";
                               });
        }
        if ( dof2pid )
        {
            std::ofstream ofs2( os2.str().c_str() );
            auto it = doftopid.begin();
            auto en = doftopid.end();
            std::for_each( it, en,
                           [&ofs2]( std::pair<size_type, size_type> const& p ) {
                               ofs2 << p.first << " " << p.second << "\n";
                           } );
        }
    }
    return std::make_pair(doftopid,pidtodof);
}
} // namespace Feel



#include <feel/feeldiscr/doftablempi.hpp>

#endif //FEELPP_DOFTABLE_HH

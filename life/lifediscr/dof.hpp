/*
  This file is part of the Life library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
 * \file dof.hpp
 * \author Christophe Prud'homme
 */
#ifndef _DOF_HH
#define _DOF_HH

#include <set>
#include <map>
#include <vector>
#include <algorithm>

#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifealg/glas.hpp>
#include <life/lifepoly/mapped.hpp>
#include <life/lifealg/datamap.hpp>

namespace Life
{
/**
 * \class Dof
 * \ingroup SpaceTime
 * \brief Local-to-global Degree of Freedom table
 *
 * \author Christophe Prud'homme
 * \author Goncalo Pena
 */
template<typename MeshType,  typename FEType, typename PeriodicityType>
class Dof : public DataMap
{
    typedef DataMap super;
public:

    /**
     * mesh type
     */
    typedef MeshType mesh_type;
    typedef FEType fe_type;
    typedef boost::shared_ptr<FEType> fe_ptrtype;


    typedef typename mesh_type::pid_element_const_iterator pid_element_const_iterator;
    typedef typename mesh_type::element_const_iterator element_const_iterator;
    typedef typename mesh_type::element_type element_type;
    typedef typename mesh_type::face_type face_type;
    typedef typename mesh_type::gm_ptrtype gm_ptrtype;
    typedef typename mesh_type::gm_type gm_type;

    typedef typename fe_type::matrix_type matrix_type;
    typedef typename fe_type::value_type value_type;
    typedef typename fe_type::reference_convex_type reference_convex_type;
    typedef typename fe_type::points_type points_type;

    typedef typename reference_convex_type::super convex_type;

    typedef typename element_type::edge_permutation_type edge_permutation_type;
    typedef typename element_type::face_permutation_type face_permutation_type;

    static const uint16_type nOrder = fe_type::nOrder;
    static const uint16_type nDim = mesh_type::nDim;
    static const uint16_type Shape = mesh_type::Shape;
    static const uint16_type nComponents = fe_type::nComponents;
    static const uint16_type nComponents1 = fe_type::nComponents1;
    static const uint16_type nComponents2 = fe_type::nComponents2;

    static const bool is_continuous = FEType::is_continuous;
    static const bool is_scalar = FEType::is_scalar;
    static const bool is_vectorial = FEType::is_vectorial;
    static const bool is_tensor2 = FEType::is_tensor2;
    static const bool is_modal = FEType::is_modal;

    static const bool is_p0_continuous = ( (nOrder == 0) && is_continuous );

    typedef PeriodicityType periodicity_type;
    static const bool is_periodic = periodicity_type::is_periodic;


    /**
     * A global dof is defined by its index in the global table
     * and by its sign.
     **/

    typedef boost::tuple<size_type, int16_type, bool> global_dof_type;

    /**
     * Type for the localToGlobal table.
     */
    typedef boost::multi_array<global_dof_type,2> Container;
    typedef typename Container::array_view<1>::type indices_per_element_type;

    typedef typename node<value_type>::type node_type;

    typedef boost::array<Container::index,2> dims_type;
    typedef typename Container::index_range range_type;
    typedef boost::tuple<node_type, size_type, uint16_type > dof_point_type;
    typedef std::vector<dof_point_type> dof_points_type;
    typedef typename std::vector<dof_point_type>::iterator dof_points_iterator;
    typedef typename std::vector<dof_point_type>::const_iterator dof_points_const_iterator;

    /**
     * Tuple that holds a size_type \p elt 1 uint16_type \p l and 1
     * uint16_type ent
     * \p elt shall be an element index
     * \p l shall be the local index of the dof in the element
     * \p ent shall be the entity the dof belongs to (0: vertex, 1: edge, 2: face, 3: volume)
     */
    typedef boost::tuple<size_type, uint16_type, uint16_type, uint16_type> local_dof_type;

    typedef boost::tuple<uint16_type&,size_type&> ref_shift_type;

    /**
     * Type that hold the map between a global dof and the elements
     */
    typedef std::map<size_type, std::list<local_dof_type> >  dof_element_type;
    typedef typename dof_element_type::iterator dof_iterator;
    typedef typename dof_element_type::const_iterator dof_const_iterator;

    typedef std::list<local_dof_type>::const_iterator ldof_const_iterator;

    typedef boost::tuple<uint16_type,uint16_type,size_type> dof_type;
    typedef std::map<dof_type, size_type> dof_map_type;
    typedef std::map<dof_type, size_type>::iterator dof_map_iterator;
    typedef std::map<dof_type, size_type, size_type>::const_iterator dof_map_const_iterator;

    typedef std::map<size_type, std::set<size_type> > dof_procset_type;
    /**
     * This type is useful to construct the sign map in the modal case
     **/

    typedef ublas::vector<bool> face_sign_info_type;

    /**
     * Type for the permutations to be done in the faces
     **/

    typedef ublas::vector<uint16_type> permutation_vector_type;

    typedef boost::tuple<element_type const*, face_type const*> element_face_pair_type;
    typedef std::list<element_face_pair_type> periodic_element_list_type;
    typedef typename periodic_element_list_type::iterator periodic_element_list_iterator;
    typedef typename periodic_element_list_type::const_iterator periodic_element_list_const_iterator;
    typedef boost::tuple<size_type /*element id*/, uint16_type /*lid*/, size_type /*gDof*/> periodic_dof_type;
    typedef std::multimap<size_type /*gid*/, periodic_dof_type> periodic_dof_map_type;

    /**
     * @brief The minimal constructor
     *
     * @param _fe reference element
     *
     */
    Dof( fe_ptrtype const& _fe, periodicity_type const& periodicity );

    /**
     * copy constructor
     *
     * @param dof2 a dof object instance
     */
    Dof( const Dof & dof2 );

    /**
     * @brief  Constructor accepting a mesh as parameter
     *
     *  @param mesh a RegionMesh3D
     *  @param _fe reference element
     */
    Dof( mesh_type& mesh, fe_ptrtype const& _fe, periodicity_type const& periodicity );

    /**
     * The number of local Dof (nodes) in the finite element
     */
    uint16_type nDofPerElement() const
    {
        return _M_fe->nLocalDof;
    }

    /**
     * \return the number of dof for faces on the boundary
     */
    uint16_type nDofPerFaceOnBoundary() const
    {
        return _M_n_dof_per_face_on_bdy;
    }

    indices_per_element_type indices( size_type id_el ) const
    {
        return _M_el_l2g[ boost::indices[id_el][range_type()] ];
    }

    std::vector<size_type> getIndices( size_type id_el ) const
    {
        size_type nldof =
            fe_type::nDofPerVolume * element_type::numVolumes +
            fe_type::nDofPerFace * element_type::numGeometricFaces +
            fe_type::nDofPerEdge * element_type::numEdges +
            fe_type::nDofPerVertex * element_type::numVertices;
        std::vector<size_type> ind( nComponents*nldof );
        for( size_type i = 0; i < nComponents*nldof; ++i )
            ind[i] = boost::get<0>( _M_el_l2g[ id_el][ i ] );
        return ind;
    }

    /**
     * \return the number of processors \c dof belongs to
     */
    size_type dofNProc( size_type dof ) const { return _M_dof_procset.find( dof )->second.size(); }

    /**
     * \return the set of processors \c dof belongs to
     */
    std::set<size_type> dofProcSet( size_type dof ) const { return _M_dof_procset.find( dof )->second; }

    /**
     * @return the coordinates of the nodal dofs associated with the
     * element \p el
     */
    const dof_point_type& dofPoint( size_type i ) const { return M_dof_points[i]; }

    /**
     * @return an iterator at the beginning of dof points
     */
    dof_points_const_iterator dofPointBegin() const { return M_dof_points.begin(); }

    /**
     * @return an iterator at the beginning of dof points
     */
    dof_points_iterator dofPointBegin() { return M_dof_points.begin(); }

    /**
     * @return an iterator at the end of dof points
     */
    dof_points_const_iterator dofPointEnd() const { return M_dof_points.end(); }

    /**
     * @return an iterator at the end of dof points
     */
    dof_points_iterator dofPointEnd() { return M_dof_points.end(); }


    /**
     * insted of creating the dof indices on the fly, get them from a
     * vector. The situation typically arises when we want to have dof
     * correspondance between two spaces
     *
     * \see OperatorLagrangeP1
     */
    void setDofIndices( std::vector<boost::tuple<size_type, uint16_type, size_type> > const& dof )
    {
        M_dof_indices.resize( dof.size() );
        std::copy( dof.begin(), dof.end(), M_dof_indices.begin() );

        if ( dof.empty() )
            return ;
        size_type nldof =
            fe_type::nDofPerVolume * element_type::numVolumes +
            fe_type::nDofPerFace * element_type::numGeometricFaces +
            fe_type::nDofPerEdge * element_type::numEdges +
            fe_type::nDofPerVertex * element_type::numVertices;
        typedef boost::tuple<size_type, uint16_type,size_type> thedof_type;
#if 1
        std::set<size_type> eltid;
        std::set<size_type> dofs;

        BOOST_FOREACH( thedof_type thedof,  dof )
            {
                eltid.insert( boost::get<0>( thedof ) );
                dofs.insert( boost::get<2>( thedof ) );
            }
#endif
        _M_el_l2g.resize( boost::extents[eltid.size()][nComponents*nldof] );


        BOOST_FOREACH( thedof_type thedof,  dof )
            {
                _M_el_l2g[ boost::get<0>( thedof ) ][ boost::get<1>( thedof ) ] = boost::get<2>( thedof );

            }
        this->_M_first_df[0] = 0;
        this->_M_last_df[0] = dofs.size()-1;
        this->_M_n_dofs = dofs.size();
    }

    /**
     * \return the dof index
     */
    size_type dofIndex( size_type dof ) const
    {
        return dof;
#if 0
        if ( M_dof_indices.empty() )
            return dof;
        LIFE_ASSERT( dof < M_dof_indices.size() )( dof )( M_dof_indices.size() ).warn( "invalid dof index" );
        return M_dof_indices[dof];
#endif
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

    global_dof_type localToGlobal( const size_type ElId,
                                   const uint16_type localNode,
                                   const uint16_type c = 0 ) const
    {
        return _M_el_l2g[ ElId][ fe_type::nLocalDof * c  + localNode ];
    }

    global_dof_type faceLocalToGlobal( const size_type ElId,
                                       const uint16_type localNode,
                                       const uint16_type c = 0 ) const
    {
        const size_type nDofF = ( face_type::numVertices * fe_type::nDofPerVertex +
                                  face_type::numEdges * fe_type::nDofPerEdge +
                                  face_type::numFaces * fe_type::nDofPerFace );
        return _M_face_l2g[ ElId][ nDofF*c+localNode ];
    }

    struct element_access
    {
        element_access( Dof const& __d )
            :
            _M_d( __d )
        {}
        global_dof_type operator()( size_type __id, uint16_type __loc, uint16_type c = 0 ) const
        {
            return _M_d._M_el_l2g[ __id][ fe_type::nLocalDof * c+ __loc ];
        }
        Dof const& _M_d;
    };
    friend class element_access;

    struct face_access
    {
        face_access( Dof const& __d )
            :
            _M_d( __d )
        {}
        global_dof_type operator()( size_type __id, size_type __loc, uint16_type c = 0 ) const
        {
            return _M_d._M_face_l2g[ __id][ __loc ];
        }
        Dof const& _M_d;
    };
    friend class face_access;

    /**
     * @brief local to global mapping
     */
    template<typename Elem>
    global_dof_type
    localToGlobal( Elem const& El, const uint16_type localNode, const uint16_type c = 0 ) const
    {
        typedef typename mpl::if_<mpl::equal_to<mpl::int_<Elem::nDim>,mpl::int_<nDim> >,mpl::identity<element_access>,mpl::identity<face_access> >::type::type access_type;
        //Debug( 5005 ) << "dof:(" << El.id() << ", " << localNode << ")= "
        //<< access_type(*this)( El.id(), localNode, c ) << "\n";
        return access_type(*this)( El.id(), localNode, c );
    }



    /**
     * Number of elements in mesh
     */
    size_type numElements() const
    {
        return _M_n_el;
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
        for ( size_type __i = 0;__i < _M_face_l2g.nrows();++__i )
            {
                for ( size_type __l = 0;__l < _M_face_l2g.ncols();++__l )
                    {
                        std::cout << "face " << __i << " local " << __l
                                  << " to global " << _M_face_l2g[ __i][ __l ] << "\n";
                    }
            }
#endif // 0
    }

    /**
     * \return the local to global map
     */
    std::map<size_type, std::pair<size_type, size_type> >
    localToGlobalMap() const
    { return _M_local2global;}

    /**
     * local(on processor) dof to global dof
     */
    size_type localToGlobal( size_type l ) const
    {
        return _M_local2global.find( l )->second.first;
    }

    /**
     * Add a new periodic dof to the dof map and to the list of periodic dof \p
     * periodic_dof for a given tag \p tag
     *
     * \arg __elt geometric element
     * \arg __face  face of the element that holds the periodic dof
     * \arg next_free_dof integer that will get possibilty incremented and holds the current number of dof
     * \arg periodic_dof table of periodic dof
     * \arg tag the boundary tag associated with the periodic dof
     */
    void addPeriodicDof( element_type const& __elt,
                         face_type const& __face,
                         size_type& next_free_dof,
                         std::map<size_type,periodic_dof_map_type>& periodic_dof,
                         size_type tag );
    /**
     * build dof map associated to the periodic dof, must be called
     * before buildDofMap
     */
    size_type buildPeriodicDofMap( mesh_type& M );

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

private:


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
                    uint16_type lc_dof,
                    uint16_type lc,
                    dof_type const& gDof,
                    uint16_type processor,
                    size_type& pDof,
                    int32_type sign = 1,
                    bool is_dof_periodic = false )
    {
        Life::detail::ignore_unused_variable_warning(lc);
        dof_map_iterator itdof = map_gdof.find( gDof );
        dof_map_iterator endof = map_gdof.end();
        bool __inserted = false;
        if ( itdof == endof )
            {
                boost::tie( itdof, __inserted ) = map_gdof.insert( std::make_pair( gDof, dofIndex( pDof ) ) );
                pDof += 1;

                LIFE_ASSERT( __inserted == true )( ie )( lc_dof )
                    (gDof.get<0>())(gDof.get<1>())( gDof.get<2>() )
                    ( processor )( itdof->second ).error( "dof should have been inserted");
            }
#if !defined( NDEBUG )
        Debug( 5005 ) << "global dof = " << itdof->second
                      << " local dof = " << fe_type::nLocalDof*itdof->first.get<1>() + lc_dof
                      << " element = " << ie
                      << " entity = " << itdof->first.get<0>()
                      << " component = " << itdof->first.get<1>()
                      << " index = " << itdof->first.get<2>() << "\n";
#endif

        LIFE_ASSERT( itdof->first == gDof ).error( "very bad logical error in insertDof" );

        LIFE_ASSERT( lc_dof >= fe_type::nLocalDof*itdof->first.get<1>() &&
                      lc_dof < fe_type::nLocalDof*(itdof->first.get<1>()+1) )
            ( lc_dof )
            ( fe_type::nLocalDof*itdof->first.get<1>() ).error( "invalid local dof index" );

        // add processor to the set of processors this dof belongs to
        _M_dof_procset[ itdof->second ].insert( processor );

        _M_el_l2g[ ie][ lc_dof ] = boost::make_tuple(itdof->second, sign, is_dof_periodic );

#if !defined(NDEBUG)
        _M_dof2elt[itdof->second].push_back( boost::make_tuple( ie, lc_dof, lc, itdof->first.get<0>() ) );
#endif

#if 0
        typedef Container::index index;

        for ( index i2 = 0; i2 < nComponents*fe_type::nLocalDof; ++i2 )
            Debug() << "dof table( " << ie  << ")=" << boost::get<0>(_M_el_l2g[ ie][ i2 ]) << "\n";
#endif
        return __inserted;
    }

    void addVertexDof( element_type const& __elt, uint16_type processor, uint16_type c, size_type& next_free_dof,
                       ref_shift_type& shifts  )
    {
        addVertexDof( __elt, processor, c, next_free_dof, shifts, mpl::bool_<fe_type::nDofPerVertex>() );
    }
    void addVertexDof( element_type const& /*M*/, uint16_type /*processor*/,  uint16_type /*c*/, size_type& /*next_free_dof*/,
                       ref_shift_type& /*shifts*/, mpl::bool_<false> )
    {}
    void addVertexDof( element_type const& __elt, uint16_type processor, uint16_type c, size_type& next_free_dof,
                       ref_shift_type& shifts, mpl::bool_<true> )
    {
        uint16_type local_shift;
        size_type global_shift;
        boost::tie( local_shift, global_shift ) = shifts;


        size_type ie = __elt.id();

        uint16_type lc = local_shift;
        for ( uint16_type i = 0; i < element_type::numVertices; ++i )
            {
                for ( uint16_type l = 0; l < fe_type::nDofPerVertex; ++l, ++lc )
                    {
                        //const size_type gDof = global_shift + ( __elt.point( i ).id() ) * fe_type::nDofPerVertex + l;
                        const size_type gDof = ( __elt.point( i ).id() ) * fe_type::nDofPerVertex + l;

                        this->insertDof( ie, lc, i, boost::make_tuple(0, c, gDof), processor, next_free_dof );
                    }
            }
        // update shifts
        shifts.get<0>() = lc;

#if !defined(NDEBUG)
        Debug( 5005 ) << "[Dof::updateVolumeDof(addVertexDof] vertex proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
    }
    void addEdgeDof( element_type const& __elt, uint16_type processor, uint16_type c, size_type& next_free_dof,
                     ref_shift_type& shifts )
    {
        return addEdgeDof( __elt,
                           processor,
                           c,
                           next_free_dof,
                           shifts,
                           mpl::int_<fe_type::nDim>(),
                           mpl::bool_<fe_type::nDofPerEdge>() );
    }
    void addEdgeDof( element_type const& /*M*/, uint16_type /*processor*/, uint16_type /*c*/, size_type& /*next_free_dof*/,
                     ref_shift_type& /*shifts*/, mpl::int_<1>, mpl::bool_<false> )
    {}

    void addEdgeDof( element_type const& __elt, uint16_type processor, uint16_type c, size_type& next_free_dof,
                     ref_shift_type& shifts, mpl::int_<1>, mpl::bool_<true> )
    {
        uint16_type local_shift;
        size_type global_shift;
        boost::tie( local_shift, global_shift ) = shifts;

        size_type ie = __elt.id();
        uint16_type lc = local_shift;
        for ( uint16_type l = 0; l < fe_type::nDofPerEdge; ++l, ++lc )
            {
                const size_type gDof = is_p0_continuous? l:ie * fe_type::nDofPerEdge + l;
                this->insertDof( ie, lc, l, boost::make_tuple(1, c, gDof), processor, next_free_dof );
            }
        // update shifts
        shifts.get<0>() = lc;
#if !defined(NDEBUG)
        Debug( 5005 ) << "[Dof::addEdgeDof(1)] element proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
    }
    void addEdgeDof( element_type const& /*__elt*/, uint16_type /*processor*/, uint16_type /*c*/, size_type& /*next_free_dof*/,
                     ref_shift_type& /*shifts*/, mpl::int_<2>, mpl::bool_<false> )
    {}
    void addEdgeDof( element_type const& __elt, uint16_type processor, uint16_type c, size_type& next_free_dof,
                     ref_shift_type& shifts, mpl::int_<2>, mpl::bool_<true> )
    {
        uint16_type local_shift;
        size_type global_shift;
        boost::tie( local_shift, global_shift ) = shifts;

        size_type ie = __elt.id();
        uint16_type lc = local_shift;

        /** The boundary dofs are constructed in the same way if the basis is modal **/

        for ( uint16_type i = 0; i < element_type::numEdges; ++i )
            {
                for ( uint16_type l = 0; l < fe_type::nDofPerEdge; ++l, ++lc )
                    {
                        size_type gDof = __elt.edge( i ).id() * fe_type::nDofPerEdge;
                        int32_type sign = 1;


                        if ( __elt.edgePermutation(i).value()  == edge_permutation_type::IDENTITY )
                            {
                                gDof += l ; // both nodal and modal case
                            }
                        else if ( __elt.edgePermutation(i).value()  == edge_permutation_type::REVERSE_PERMUTATION)
                            {

                                if(fe_type::is_modal)
                                    {
                                        //only half of the modes (odd polynomial order) are negative.
                                        sign = (l%2)?(-1):(1);
                                        gDof += l;
                                    }
                                else
                                    gDof += fe_type::nDofPerEdge - 1 - l ;
                            }
                        else
                            LIFE_ASSERT( 0 ).error ( "invalid edge permutation" );
                        this->insertDof( ie, lc, i, boost::make_tuple(1, c, gDof), processor, next_free_dof, sign );
                    }
            }

        // update shifts
        shifts.get<0>() = lc;
#if !defined(NDEBUG)
        Debug( 5005 ) << "[Dof::addEdgeDof] edge proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
    }

    void addEdgeDof( element_type const& /*__elt*/, uint16_type /*processor*/, uint16_type /*c*/, size_type& /*next_free_dof*/,
                     ref_shift_type& /*shifts*/, mpl::int_<3>, mpl::bool_<false> )
    {}

    void addEdgeDof( element_type const& __elt, uint16_type processor, uint16_type c, size_type& next_free_dof,
                     ref_shift_type& shifts, mpl::int_<3>, mpl::bool_<true> )
    {
        uint16_type local_shift;
        size_type global_shift;
        boost::tie( local_shift, global_shift ) = shifts;

        size_type ie = __elt.id();
        uint16_type lc = local_shift;
        for ( uint16_type i = 0; i < element_type::numEdges; ++i )
            {
                for ( uint16_type l = 0; l < fe_type::nDofPerEdge; ++l, ++lc )
                    {
                        size_type gDof = __elt.edge( i ).id() * fe_type::nDofPerEdge;

                        int32_type sign = 1;

                        if ( __elt.edgePermutation(i).value()  == edge_permutation_type::IDENTITY )
                            {
                                gDof += l ; // both nodal and modal case
                            }
                        else if ( __elt.edgePermutation(i).value()  == edge_permutation_type::REVERSE_PERMUTATION)
                            {

                                if(fe_type::is_modal)
                                    {
                                        //only half of the modes (odd polynomial order) are negative.
                                        sign = (l%2)?(-1):(1);
                                        gDof += l;
                                    }
                                else
                                    gDof += fe_type::nDofPerEdge - 1 - l ;
                            }
                        else
                            LIFE_ASSERT( 0 ).error ( "invalid edge permutation" );

                        this->insertDof( ie, lc, i, boost::make_tuple(1, c, gDof), processor, next_free_dof, sign );
                    }
            }
        // update shifts
        shifts.get<0>() = lc;
#if !defined(NDEBUG)
        Debug( 5005 ) << "[Dof::addEdgeDof] edge proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
    }


    void addFaceDof( element_type const& __elt, uint16_type processor, uint16_type c, size_type& next_free_dof,
                     ref_shift_type& shifts )
    {
        return addFaceDof( __elt, processor, c, next_free_dof, shifts, mpl::int_<fe_type::nDim>(), mpl::bool_<fe_type::nDofPerFace>() );
    }
    void addFaceDof( element_type const& /*M*/, uint16_type /*processor*/, uint16_type /*c*/, size_type& /*next_free_dof*/,
                     ref_shift_type& /*shifts*/, mpl::int_<1>, mpl::bool_<false> )
    {}
    void addFaceDof( element_type const& /*M*/, uint16_type /*processor*/, uint16_type /*c*/, size_type& /*next_free_dof*/,
                     ref_shift_type& /*shifts*/, mpl::int_<2>, mpl::bool_<false> )
    {}
    void addFaceDof( element_type const& __elt, uint16_type processor, uint16_type c, size_type& next_free_dof,
                     ref_shift_type& shifts, mpl::int_<2>, mpl::bool_<true> )
    {
        uint16_type local_shift;
        size_type global_shift;
        boost::tie( local_shift, global_shift ) = shifts;

        size_type ie = __elt.id();
        uint16_type lc = local_shift;
        for ( uint16_type l = 0; l < fe_type::nDofPerFace; ++l, ++lc )
            {
                const size_type gDof = is_p0_continuous? l:ie * fe_type::nDofPerFace + l;
                this->insertDof( ie, lc, l, boost::make_tuple(2, c, gDof), processor, next_free_dof );
            }
        // update shifts
        shifts.get<0>() = lc;
#if !defined(NDEBUG)
        Debug( 5005 ) << "[Dof::addFaceDof(2,true)] face proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
    }
    void addFaceDof( element_type const& /*M*/, uint16_type /*processor*/, uint16_type /*c*/, size_type& /*next_free_dof*/,
                     ref_shift_type& /*shifts*/, mpl::int_<3>, mpl::bool_<false> )
    {}
    void addFaceDof( element_type const& __elt, uint16_type processor, uint16_type c, size_type& next_free_dof,
                     ref_shift_type& shifts, mpl::int_<3>, mpl::bool_<true> )
    {
        uint16_type local_shift;
        size_type global_shift;
        boost::tie( local_shift, global_shift ) = shifts;

        size_type ie = __elt.id();

        uint16_type lc = local_shift;
        for ( uint16_type i = 0; i < element_type::numFaces; ++i )
            {
                face_permutation_type permutation = __elt.facePermutation(i);
                LIFE_ASSERT( permutation != face_permutation_type(0) ).error ( "invalid face permutation" );

                // Polynomial order in each direction
                uint16_type p=1;
                uint16_type q=0;

                // MaxOrder = Order - 2
                int MaxOrder = int((3 + std::sqrt(1+8*fe_type::nDofPerFace))/2) - 2;

                for ( uint16_type l = 0; l < fe_type::nDofPerFace; ++l, ++lc )
                    {

                        // TODO: orient the dof indices such
                        // that they match properly the faces
                        // dof of the connected faces. There
                        // are a priori many permutations of
                        // the dof face indices
                        size_type gDof = __elt.face( i ).id() * fe_type::nDofPerFace;
                        int32_type sign = 1;

                        q=q+1;
                        if(q > MaxOrder)
                            {
                                q = 1;
                                p = p+1;
                                MaxOrder = MaxOrder-1;
                            }

                        if ( !fe_type::is_modal )
                            {
                                if ( permutation  == face_permutation_type(1) )
                                    gDof += l;
                                else
                                    gDof += vector_permutation[permutation][l];
                            }
                        else
                            {
                                gDof += l;
                                if ( permutation == face_permutation_type(2) )
                                    {
                                        // Reverse sign if polynomial order in
                                        // eta_1 direction is odd

                                        if(p%2 == 0)
                                            sign = -1;

                                    }
                            }
                        this->insertDof( ie, lc, i, boost::make_tuple(2, c, gDof), processor, next_free_dof, sign );

                    }
            }
        // update shifts
        shifts.get<0>() = lc;
#if !defined(NDEBUG)
        Debug( 5005 ) << "[Dof::addFaceDof<3>] face proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
    }
    void addVolumeDof( element_type const& __elt, uint16_type processor, uint16_type c, size_type& next_free_dof,
                       ref_shift_type& shifts )
    {
        return addVolumeDof( __elt, processor, c, next_free_dof, shifts, mpl::bool_<fe_type::nDofPerVolume>() );
    }
    void addVolumeDof( element_type const& /*M*/, uint16_type /*processor*/, uint16_type /*c*/, size_type& /*next_free_dof*/,
                       ref_shift_type& /*shifts*/, mpl::bool_<false> )
    {}
    void addVolumeDof( element_type const& __elt, uint16_type processor, uint16_type c, size_type& next_free_dof,
                       ref_shift_type& shifts, mpl::bool_<true> )
    {
        BOOST_STATIC_ASSERT( element_type::numVolumes );
        uint16_type local_shift;
        size_type global_shift;
        boost::tie( local_shift, global_shift ) = shifts;

        size_type ie = __elt.id();
        uint16_type lc = local_shift;
        for ( uint16_type l = 0; l < fe_type::nDofPerVolume; ++l, ++lc )
            {
                const size_type gDof = is_p0_continuous? l:ie * fe_type::nDofPerVolume + l;
                this->insertDof( ie, lc, l, boost::make_tuple(3, c, gDof), processor, next_free_dof );
            }
        // update shifts
        shifts.get<0>() = lc;
#if !defined(NDEBUG)
        Debug( 5005 ) << "[Dof::updateVolumeDof(<2>)] element proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
    }

    template<typename FaceIterator>
    void addVertexBoundaryDof( FaceIterator __face_it, uint16_type c, uint16_type& lc )
    {
        addVertexBoundaryDof( __face_it, c, lc, mpl::bool_<fe_type::nDofPerVertex>() );
    }
    template<typename FaceIterator>
    void addVertexBoundaryDof( FaceIterator /*__face_it*/, uint16_type /*c*/, uint16_type& /*lc*/, mpl::bool_<false> )
    {
    }
    template<typename FaceIterator>
    void addVertexBoundaryDof( FaceIterator __face_it, uint16_type c, uint16_type& lc, mpl::bool_<true>  )
    {
        BOOST_STATIC_ASSERT( face_type::numVertices );

        // id of the element adjacent to the face
        // \warning NEED TO INVESTIGATE THIS
        size_type iElAd = __face_it->ad_first();
        LIFE_ASSERT( iElAd != invalid_size_type_value )( __face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

        // local id of the face in its adjacent element
        uint16_type iFaEl = __face_it->pos_first();
        LIFE_ASSERT( iFaEl != invalid_uint16_type_value ).error ("invalid element index in face");

        //_M_dof2elt[gDof].push_back( boost::make_tuple( iElAd, lc-1, 48, 0 ) );
        // loop on face vertices
        for ( uint16_type iVeFa = 0; iVeFa < face_type::numVertices; ++iVeFa )
            {
                // local vertex number (in element)
                uint16_type iVeEl = element_type::fToP( iFaEl, iVeFa );

                LIFE_ASSERT( iVeEl != invalid_uint16_type_value ).error( "invalid local dof" );

                // Loop number of Dof per vertex
                for ( uint16_type l = 0; l < fe_type::nDofPerVertex; ++l )
                    {
                        _M_face_l2g[ __face_it->id()][ lc++ ] = this->localToGlobal( iElAd,
                                                                                     iVeEl * fe_type::nDofPerVertex + l,
                                                                                     c );
                    }
            }
    }

    template<typename FaceIterator>
    void addEdgeBoundaryDof( FaceIterator __face_it, uint16_type c, uint16_type& lc )
    {
        addEdgeBoundaryDof( __face_it, c, lc, mpl::bool_<fe_type::nDofPerEdge*face_type::numEdges>(), mpl::int_<nDim>() );
    }
    template<typename FaceIterator>
    void addEdgeBoundaryDof( FaceIterator /*__face_it*/, uint16_type, uint16_type& /*lc*/, mpl::bool_<false>, mpl::int_<1> ){}
    template<typename FaceIterator>
    void addEdgeBoundaryDof( FaceIterator /*__face_it*/, uint16_type, uint16_type& /*lc*/, mpl::bool_<false>, mpl::int_<2> ){}
    template<typename FaceIterator>
    void addEdgeBoundaryDof( FaceIterator /*__face_it*/, uint16_type, uint16_type& /*lc*/, mpl::bool_<false>, mpl::int_<3> ){}
    template<typename FaceIterator>
    void addEdgeBoundaryDof( FaceIterator __face_it, uint16_type c, uint16_type& lc, mpl::bool_<true>, mpl::int_<2> )
    {
        // id of the element adjacent to the face
        // \warning NEED TO INVESTIGATE THIS
        size_type iElAd = __face_it->ad_first();
        LIFE_ASSERT( iElAd != invalid_size_type_value )
            ( __face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

        // local id of the face in its adjacent element
        uint16_type iFaEl = __face_it->pos_first();
        LIFE_ASSERT( iFaEl != invalid_uint16_type_value ).error ("invalid element index in face");
#if !defined(NDEBUG)
        Debug( 5005 ) << " local face id : " << iFaEl << "\n";
#endif
        // Loop number of Dof per edge
        for ( uint16_type l = 0; l < fe_type::nDofPerEdge; ++l )
            {
                _M_face_l2g[ __face_it->id()][ lc++ ] = this->localToGlobal( iElAd,
                                                                             element_type::numVertices*fe_type::nDofPerVertex +
                                                                             iFaEl * fe_type::nDofPerEdge + l,
                                                                             c );
            }
    }
    template<typename FaceIterator>
    void addEdgeBoundaryDof( FaceIterator __face_it, uint16_type c, uint16_type& lc, mpl::bool_<true>, mpl::int_<3> )
    {
        //BOOST_STATIC_ASSERT( face_type::numEdges );

        // id of the element adjacent to the face
        // \warning NEED TO INVESTIGATE THIS
        size_type iElAd = __face_it->ad_first();
        LIFE_ASSERT( iElAd != invalid_size_type_value )( __face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

        // local id of the face in its adjacent element
        uint16_type iFaEl = __face_it->pos_first();
        LIFE_ASSERT( iFaEl != invalid_uint16_type_value ).error ("invalid element index in face");
#if !defined(NDEBUG)
        Debug( 5005 ) << " local face id : " << iFaEl << "\n";
#endif

        // loop on face vertices
        for ( uint16_type iEdFa = 0; iEdFa < face_type::numEdges; ++iEdFa )
            {
                // local edge number (in element)
                uint16_type iEdEl = element_type::fToE( iFaEl, iEdFa );

                LIFE_ASSERT( iEdEl != invalid_uint16_type_value ).error( "invalid local dof" );

                // Loop number of Dof per edge
                for ( uint16_type l = 0; l < fe_type::nDofPerEdge; ++l )
                    {
                        _M_face_l2g[ __face_it->id()][ lc++ ] = this->localToGlobal( iElAd,
                                                                                     element_type::numVertices*fe_type::nDofPerVertex +
                                                                                     iEdEl * fe_type::nDofPerEdge + l,
                                                                                     c );
                    }
            }
    }


    template<typename FaceIterator>
    void addFaceBoundaryDof( FaceIterator __face_it, uint16_type c, uint16_type lc )
    {
        addFaceBoundaryDof( __face_it, c, lc, mpl::bool_<face_type::numFaces*fe_type::nDofPerFace>() );
    }
    template<typename FaceIterator>
    void addFaceBoundaryDof( FaceIterator /*__face_it*/, uint16_type, uint16_type /*lc*/, mpl::bool_<false> )
    {
    }
    template<typename FaceIterator>
    void addFaceBoundaryDof( FaceIterator __face_it, uint16_type c, uint16_type lc, mpl::bool_<true> )
    {
        // id of the element adjacent to the face
        // \warning NEED TO INVESTIGATE THIS
        size_type iElAd = __face_it->ad_first();
        LIFE_ASSERT( iElAd != invalid_size_type_value )( __face_it->id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

        // local id of the face in its adjacent element
        uint16_type iFaEl = __face_it->pos_first();
        LIFE_ASSERT( iFaEl != invalid_uint16_type_value ).error ("invalid element index in face");
#if !defined(NDEBUG)
        Debug( 5005 ) << " local face id : " << iFaEl << "\n";
#endif
        // Loop on number of Dof per face
        for ( uint16_type l = 0; l < fe_type::nDofPerFace; ++l )
            {
                _M_face_l2g[ __face_it->id()][ lc++ ] = this->localToGlobal( iElAd,
                                                                             element_type::numVertices*fe_type::nDofPerVertex +
                                                                             element_type::numEdges*fe_type::nDofPerEdge +
                                                                             iFaEl * fe_type::nDofPerFace + l,
                                                                             c );
            }
    }

    /**
     * @brief Checks if the dofs associated with entity_id N are continuous
     *
     * \param M A mesh
     */
    template<int N>
    void checkDofEntity (mesh_type& M)
    {
        Life::detail::ignore_unused_variable_warning(M);
        if ( !is_scalar )
            return;

#if !defined(NDEBUG)

using namespace Life;

        gm_ptrtype _M_gm_ptr = M.gm();

        fe_type _M_basis;

        //value_type tol = value_type(100.0)*type_traits<value_type>::epsilon();
        value_type tol = value_type(10.0)*type_traits<double>::epsilon();

        bool global_signs_good = 1;

        std::vector<size_type> bad_dof;

        for (uint16_type gDof = 0; gDof < nDof(); ++gDof)
            {
                uint16_type _numEntities = _M_dof2elt[gDof].size();
                uint16_type _ent = _M_dof2elt[gDof].begin()->get<3>();

                if ( _numEntities > 1 && _ent == mpl::int_<N>() )
                    {
                        bool signs_good = 1;

                        std::vector< ublas::vector<value_type> > basis_eval;

                        std::vector< points_type > real_coordinates;

                        ldof_const_iterator __ldofit = _M_dof2elt[gDof].begin();
                        ldof_const_iterator __ldofen = _M_dof2elt[gDof].end();

                        while ( __ldofit != __ldofen )
                            {
                                size_type entity_element_id = __ldofit->get<0>();
                                uint16_type entity_local_dof_id = __ldofit->get<1>();
                                uint16_type entity_local_id = __ldofit->get<2>();

                                PointSetMapped<element_type, convex_type, nOrder> test( M.element(entity_element_id) );

                                points_type Gt = test.pointsBySubEntity(N, entity_local_id);

                                real_coordinates.push_back( test.pointsBySubEntity(N, entity_local_id, 0, 1) );

                                int sign = boost::get<1>(localToGlobal( entity_element_id, entity_local_dof_id ));

                                basis_eval.push_back( value_type(sign)*ublas::row( _M_basis.evaluate(Gt), entity_local_dof_id ) );

                                ++__ldofit;
                            }

                        for (uint16_type i=1; i < _numEntities; i++)
                            {

                                LIFE_ASSERT( ublas::norm_inf( real_coordinates[i] - real_coordinates[0] ) < tol  )
                                    ( gDof )
                                    ( real_coordinates[0] )
                                    ( real_coordinates[i] ).error( "Reference points aren't being mapped to the same real one's" );

                                if ( ublas::norm_inf(basis_eval[i] - basis_eval[0]) > tol )
                                    {
                                        signs_good = 0;
                                        global_signs_good = 0;
                                    }
                            }

                        basis_eval.empty();
                        real_coordinates.empty();

                        if (signs_good == 0)
                            bad_dof.push_back(gDof);
                    }
            }

        if ( !bad_dof.empty() )
            {
                for (uint16_type i = 0; i < bad_dof.size(); ++i)
                    Warning() << bad_dof[i] << "\n";
                if (mpl::int_<N>() == 1)
                    Warning() << "Edges: ";
                else
                    Warning() << "Faces: ";
                Warning() << "Bad dof signs. \n";
            }
#endif
    }

    void checkDofContinuity( mesh_type& /*mesh*/, mpl::int_<1> ) {}

    void checkDofContinuity( mesh_type& mesh, mpl::int_<2> ) { checkDofEntity<1>(mesh); }

    void checkDofContinuity( mesh_type& mesh, mpl::int_<3> )
    {
        checkDofContinuity(mesh, mpl::int_<2>() );
        checkDofEntity<2>(mesh);
    }

    void generateFacePermutations ( mesh_type& /*mesh*/, mpl::bool_<false> ) {}

    void generateFacePermutations ( mesh_type& mesh, mpl::bool_<true> )
    {
        PointSetMapped<element_type, convex_type, nOrder> pts( mesh.element(0) );

        for (uint16_type i = 2; i < face_permutation_type::N_PERMUTATIONS; i++)
            vector_permutation[face_permutation_type(i)] = pts.getVectorPermutation(face_permutation_type(i));
    }
    void generateDofPoints( mesh_type& M );
    void generatePeriodicDofPoints( mesh_type& M, periodic_element_list_type const& periodic_elements, dof_points_type& periodic_dof_points );
private:

    fe_ptrtype _M_fe;

    reference_convex_type _M_convex_ref;

    size_type _M_n_el;

    uint16_type _M_n_dof_per_face_on_bdy;
    uint16_type _M_n_dof_per_face;

    Container _M_el_l2g;
    Container _M_face_l2g;

    dof_element_type _M_dof2elt;

    dof_map_type map_gdof;

    std::map<face_permutation_type, permutation_vector_type> vector_permutation;

    /**
     * map between local(on processor) dof and global dofs
     */
    std::map<size_type, std::pair<size_type, size_type> > _M_local2global;


    face_sign_info_type _M_face_sign;

    /**
     * for each dof, store the set of processor it belongs to
     */
    dof_procset_type _M_dof_procset;

    /**
     * coordinates of the nodal dofs
     */
    dof_points_type M_dof_points;

    std::vector<boost::tuple<size_type, uint16_type, size_type> > M_dof_indices;

    periodicity_type M_periodicity;
};

template<typename MeshType, typename FEType, typename PeriodicityType>
const uint16_type Dof<MeshType, FEType, PeriodicityType>::nComponents;

template<typename MeshType, typename FEType, typename PeriodicityType>
Dof<MeshType, FEType, PeriodicityType>::Dof( mesh_type& mesh, fe_ptrtype const& _fe, periodicity_type const& periodicity )
    :
    super(),
    _M_fe( _fe ),
    _M_n_el( invalid_size_type_value ),
    _M_n_dof_per_face_on_bdy( invalid_uint16_type_value ),
    _M_n_dof_per_face( invalid_uint16_type_value ),
    _M_el_l2g(),
    _M_face_l2g(),
    map_gdof(),
    _M_face_sign(),
    M_dof_indices(),
    M_periodicity( periodicity )
{
    Debug( 5015 ) << "[dof] is_periodic = " << is_periodic << "\n";
    size_type start_next_free_dof = 0;
    if ( is_periodic )
        start_next_free_dof = buildPeriodicDofMap( mesh );
    buildDofMap( mesh, start_next_free_dof );
    buildBoundaryDofMap( mesh );
}

template<typename MeshType, typename FEType, typename PeriodicityType>
Dof<MeshType, FEType, PeriodicityType>::Dof( fe_ptrtype const& _fe, periodicity_type const& periodicity )
    :
    super(),
    _M_fe( _fe ),
    _M_n_el( 0 ),
    _M_n_dof_per_face_on_bdy( invalid_uint16_type_value ),
    _M_n_dof_per_face( invalid_uint16_type_value ),
    _M_el_l2g(),
    _M_face_l2g(),
    map_gdof(),
    _M_face_sign(),
    M_dof_indices(),
    M_periodicity( periodicity )

{
}

template<typename MeshType, typename FEType, typename PeriodicityType>
Dof<MeshType, FEType, PeriodicityType>::Dof( const Dof<MeshType, FEType, PeriodicityType> & dof2 )
    :
    super( dof2 ),
    _M_fe( dof2.fe ),
    _M_n_el( dof2._M_n_el ),
    _M_n_dof_per_face_on_bdy( dof2._M_n_dof_per_face_on_bdy ),
    _M_n_dof_per_face( dof2._M_n_dof_per_face ),
    _M_el_l2g( dof2._M_el_l2g ),
    _M_face_l2g( dof2._M_face_l2g ),
    map_gdof( dof2.map_gdof ),
    _M_face_sign( dof2._M_face_sign),
    M_dof_indices( dof2.M_dof_indices ),
    M_periodicity( dof2.M_periodicity )
{
}

template<typename MeshType, typename FEType, typename PeriodicityType>
void
Dof<MeshType, FEType, PeriodicityType>::showMe() const
{
    Debug(5005)  << " Degree of Freedom (Dof) Object" << "\n";
    //if ( verbose )
    {
        Debug(5005)  << "************************************************************" << "\n";
        Debug(5005)  << "           Local to Global DOF table" << "\n";
        Debug(5005)  << "************************************************************" << "\n";
        Debug(5005)  << "Element Id    Loc. N.    Global N.   Sign#    Element Id   Loc. N.  Global N.  Sign" << "\n";

        for ( size_type i = 0; i < _M_n_el;++i )
            {

                for ( size_type j = 0; j < nDofPerElement();++j )
                    {

                        Debug(5005)<< "elt id " << i << " : "
                                   << "(local/global/sign dof : " << j << " : "
                                   << boost::get<0>(localToGlobal( i  , j )) << " : "
                                   << boost::get<1>(localToGlobal( i  , j )) << "\n";
                    }

            }
        Debug(5005)  << "\n";

        Debug(5005)  << "************************************************************" << "\n";
        Debug(5005)  << " Boundary  Local to Global DOF table" << "\n";
        Debug(5005)  << "************************************************************" << "\n";
        typedef typename Container::const_iterator const_iterator;
        const_iterator it = _M_face_l2g.begin();
        const_iterator en = _M_face_l2g.end();
        for(size_type f = 0;it!=en;++it,++f)
            {
                std::ostringstream ostr;
                ostr  << "face id " << f << " : ";
                typedef typename Container::template subarray<1>::type::const_iterator const_iterator2;
                const_iterator2 it2 = it->begin();
                const_iterator2 en2 = it->end();
                for(size_type l = 0;it2!=en2;++it2,++l)
                    {
                        ostr << "(local/global/sign dof : " << l << " : "
                             << boost::get<0>(*it2)<< " : "
                             << boost::get<1>(*it2) << "\n";
                    }
                Debug(5005) << ostr.str() << "\n";
            }
    }

}

template<typename MeshType, typename FEType, typename PeriodicityType>
void
Dof<MeshType, FEType, PeriodicityType>::addPeriodicDof( element_type const& __elt,
                                                        face_type const& __face,
                                                        size_type& next_free_dof,
                                                        std::map<size_type,periodic_dof_map_type>& periodic_dof,
                                                        size_type tag )
{
    // store the element and local dof id for further
    // reference when inserting the associated global dof
    // id of the element adjacent to the face

    size_type iElAd = __face.ad_first();
    LIFE_ASSERT( iElAd != invalid_size_type_value )( __face.id() ).error( "[periodic]invalid face/element in face" );
    Life::detail::ignore_unused_variable_warning(iElAd);

    // local id of the face in its adjacent element
    uint16_type iFaEl = __face.pos_first();
    LIFE_ASSERT( iFaEl != invalid_uint16_type_value ).error ("invalid element index in face");

    // loop on face vertices
    for ( uint16_type iVeFa = 0; iVeFa < face_type::numVertices; ++iVeFa )
        {
            // local vertex number (in element)
            uint16_type iVeEl = element_type::fToP( iFaEl, iVeFa );
            Life::detail::ignore_unused_variable_warning(iVeEl);

            LIFE_ASSERT( iVeEl != invalid_uint16_type_value ).error( "invalid local dof" );

            // Loop number of Dof per vertex
            for ( uint16_type l = 0; l < fe_type::nDofPerVertex; ++l )
                {
                    uint16_type lid = iVeEl * fe_type::nDofPerVertex + l;
                    //const size_type gDof = global_shift + ( __elt.point( i ).id() ) * fe_type::nDofPerVertex + l;
                    const size_type gDof = ( __elt.point( iVeEl ).id() ) * fe_type::nDofPerVertex + l;

                    Debug( 5015 ) << "add periodic doc " << next_free_dof << " in element " << __elt.id() << " lid = " << lid << "\n";
                    size_type dof_id = next_free_dof;
                    // next_free_dof might be incremented if a new dof is created
                    bool inserted = this->insertDof( __elt.id(), lid, iVeEl, boost::make_tuple(0, 0, gDof), 0, next_free_dof, 1, true );
                    Debug( 5015 ) << "periodic dof inserted : " << inserted << "\n";

                    // add the pair (elt, lid) to the map associated
                    // with dof_id, one dof can be shared by several
                    // elements
                    dof_id = boost::get<0>(localToGlobal( __elt.id(), lid, 0 ));
                    periodic_dof[tag].insert( std::make_pair( dof_id, boost::make_tuple( __elt.id(), lid, gDof ) ) );

                    Debug( 5015 ) << "added " <<  __elt.id() << ", " <<  lid << ", " << boost::get<0>(localToGlobal( __elt.id(), lid, 0 )) << "\n";
                }

        }

}
template<typename MeshType, typename FEType, typename PeriodicityType>
size_type
Dof<MeshType, FEType, PeriodicityType>::buildPeriodicDofMap( mesh_type& M )
{
    _M_n_el = M.numElements();

    size_type nldof =
        fe_type::nDofPerVolume * element_type::numVolumes +
        fe_type::nDofPerFace * element_type::numGeometricFaces +
        fe_type::nDofPerEdge * element_type::numEdges +
        fe_type::nDofPerVertex * element_type::numVertices;

    Debug(5015) << "==============================\n";
    Debug(5015) << "[buildPeriodicDofMap]\n";
    Debug(5015) << "nldof                   = "  << int(nldof) << "\n";
    Debug(5015) << "fe_type::nLocalDof     = "  << int(fe_type::nLocalDof) << "\n";
    Debug(5015) << "fe_type::nDofPerVolume = "  << int(fe_type::nDofPerVolume) << "\n";
    Debug(5015) << "fe_type::nDofPerFace   = "  << int(fe_type::nDofPerFace) << "\n";
    Debug(5015) << "fe_type::nDofPerEdge   = "  << int(fe_type::nDofPerEdge) << "\n";
    Debug(5015) << "fe_type::nDofPerVertex = "  << int(fe_type::nDofPerVertex) << "\n";
    Debug(5015) << "element_type::numVolumes= "  << int(element_type::numVolumes) << "\n";
    Debug(5015) << "element_type::numFaces= "    << int(element_type::numFaces) << "\n";
    Debug(5015) << "element_type::numEdges= "    << int(element_type::numEdges) << "\n";
    Debug(5015) << "element_type::numVertices= " << int(element_type::numVertices) << "\n";
    Debug(5015) << "==============================\n";

    LIFE_ASSERT( nldof == fe_type::nLocalDof )
        ( nldof )
        ( fe_type::nLocalDof ).error( "Something wrong in FE specification" ) ;

    // initialize the local to global map and fill it with invalid
    // values that will allow to check whether we have a new dof or
    // not when building the table
    const size_type nV = M.numElements();
    _M_el_l2g.resize( boost::extents[nV][nComponents*nldof] );
    typedef Container::index index;
    for ( index i1 = 0; i1 < index(nV); ++i1 )
        for ( index i2 = 0; i2 < index(nComponents*nldof); ++i2 )
            _M_el_l2g[i1][i2] = boost::make_tuple(invalid_size_type_value,0,false); // 0 is the invalid value for the sign !

    const size_type n_proc  = Application::nProcess();

    _M_face_sign = ublas::scalar_vector<bool>(M.numFaces(), false);

    generateFacePermutations( M, mpl::bool_< ((Shape == SHAPE_TETRA && nOrder > 2 ) || (Shape == SHAPE_HEXA && nOrder > 1 ))>() );

    //! list of elements which have a periodic face Tag2
    periodic_element_list_type periodic_elements;

    for (size_type processor=0; processor<n_proc; processor++)
        {
            // compute the number of dof on current processor
            element_const_iterator it_elt = M.beginElementWithProcessId( processor );
            element_const_iterator en_elt = M.endElementWithProcessId( processor );
            size_type n_elts = std::distance( it_elt, en_elt);
            Debug( 5015 ) << "[buildDofMap] n_elts =  " << n_elts << " on processor " << processor << "\n";
            //this->_M_first_df[processor] = next_free_dof;

            it_elt = M.beginElementWithProcessId( processor );

            it_elt = M.beginElementWithProcessId( processor );
            Debug( 5015 ) << "[buildDofMap] starting with elt " << it_elt->id() << "\n";
            for (;it_elt!=en_elt; ++it_elt )
                {
                    element_type const& __elt = *it_elt;
                    //Debug( 5015 ) << "next_free_dof " << next_free_dof  << "\n";
                    //Debug( 5015 ) << "current dof " << dofIndex( next_free_dof ) << "\n";

                    typename element_type::face_const_iterator it, en;
                    boost::tie( it, en ) = it_elt->faces();

                    bool found_periodic_face_in_element = false;
                    for( ;it != en; ++it )
                        {
                            if ( (*it)->marker().value() == periodicity_type::tag2 ||
                                 (*it)->marker().value() == periodicity_type::tag1 )
                                {
                                    // store the element reference for the end, the associated
                                    // dof on the periodic face is in fact already taken care of.
                                    // the "internal" dof or on not periodic face will be added
                                    periodic_elements.push_back( boost::make_tuple( boost::addressof(__elt), *it ) );
                                    found_periodic_face_in_element = true;
                                    break;
                                }
                        }
                }
        }
    Debug(5015) << "[buildPeriodicDofMap] built periodic_elements " << periodic_elements.size() << "\n";
    std::map<size_type,periodic_dof_map_type> periodic_dof;
    /*
     * Generate the periodic dof, assign a gid to the tag1 dof and set
     * the tag2 dof to invalid_size_type_value for now.
     */
    periodic_element_list_iterator it_periodic = periodic_elements.begin();
    periodic_element_list_iterator en_periodic = periodic_elements.end();
    size_type next_free_dof = 0;
    while( it_periodic != en_periodic )
        {
            element_type const& __elt = *it_periodic->template get<0>();
            face_type const& __face = *it_periodic->template get<1>();
            if ( __face.marker().value() == periodicity_type::tag1 )
                addPeriodicDof( __elt, __face, next_free_dof, periodic_dof, __face.marker().value() );

            ++it_periodic;
        }
    it_periodic = periodic_elements.begin();
    while( it_periodic != en_periodic )
        {
            element_type const& __elt = *it_periodic->template get<0>();
            face_type const& __face = *it_periodic->template get<1>();
            if ( __face.marker().value() == periodicity_type::tag2 )
                addPeriodicDof( __elt, __face, next_free_dof, periodic_dof, __face.marker().value() );

            ++it_periodic;
        }

    Debug( 5015 ) << "[periodic dof table] next_free_dof : " << next_free_dof << "\n";
    Debug( 5015 ) << "[periodic dof table] number of periodic dof : " << periodic_dof[periodicity_type::tag1].size() << "\n";

    dof_points_type periodic_dof_points( next_free_dof );
    generatePeriodicDofPoints( M, periodic_elements, periodic_dof_points );

    Debug( 5015 ) << "[periodic dof table] generated dof points\n";
    Debug( 5015 ) << "[periodic dof table] start matching the dof points\n";

    size_type max_gid = 0;
    std::pair<size_type,periodic_dof_type> dof;
    BOOST_FOREACH( dof, periodic_dof[periodicity_type::tag1] )
        {
            size_type gid = dof.first;
            max_gid = (max_gid > gid)?max_gid:gid;
        }

    size_type max_gid2 = 0;
    std::pair<size_type,periodic_dof_type> dof2;
    BOOST_FOREACH( dof2, periodic_dof[periodicity_type::tag2] )
        {
            size_type gid2 = dof2.first;
            max_gid2 = (max_gid2 > gid2)?max_gid2:gid2;
        }
    LIFE_ASSERT( (max_gid+1) == (max_gid2+1-(max_gid+1)) )( max_gid )( max_gid2 ).error( "[periodic] invalid periodic setup" );

    std::vector<bool> periodic_dof_done( max_gid+1 );
    std::fill( periodic_dof_done.begin(), periodic_dof_done.end(), false );

    BOOST_FOREACH( dof, periodic_dof[periodicity_type::tag1] )
        {
            size_type gid = dof.first;
            if ( periodic_dof_done[gid] )
                continue;
            node_type x1 = periodic_dof_points[gid].template get<0>();
            bool match = false;
            typename periodic_dof_map_type::iterator it_dof2 = periodic_dof[periodicity_type::tag2].begin();
            typename periodic_dof_map_type::iterator en_dof2 = periodic_dof[periodicity_type::tag2].end();

            for( ; it_dof2 != en_dof2; ++ it_dof2 )
                {
                    size_type gid2 = it_dof2->first;
                    LIFE_ASSERT( gid2 < next_free_dof )( gid )( gid2 )( next_free_dof ).error( "[periodic] invalid dof id" );
                    node_type x2 = periodic_dof_points[gid2].template get<0>();
                    //LIFE_ASSERT( math::abs( x2[0]-M_periodicity.translation()[0]) < 1e-10 )
                    //( x1 )( x2 )( M_periodicity.translation() ).error( "[periodic] invalid periodic setup");
                }
            it_dof2 = periodic_dof[periodicity_type::tag2].begin();
            size_type corresponding_gid = invalid_size_type_value;
            for( ; it_dof2 != en_dof2; ++ it_dof2 )
                {
                    size_type gid2 = it_dof2->first;
                    LIFE_ASSERT( gid2 < next_free_dof )( gid )( gid2 )( next_free_dof ).error( "[periodic] invalid dof id" );
                    node_type x2 = periodic_dof_points[gid2].template get<0>();
                    LIFE_ASSERT( ( x1.size() == x2.size() ) &&
                                 ( x1.size() == M_periodicity.translation().size() ) )
                        ( gid )( dof.second.template get<0>() )( dof.second.template get<1>())
                        ( gid2 )( it_dof2->second.template get<0>() )( it_dof2->second.template get<1>())
                        ( x1 )( x2 )( M_periodicity.translation() ).error( "invalid point size" );
                    //LIFE_ASSERT( math::abs( x1[0]-(x2[0]-M_periodicity.translation()[0])) < 1e-10 )
                    //( x1 )( x2 )( M_periodicity.translation() ).error( "[periodic] invalid periodic setup");
                    if ( ublas::norm_2( x1-(x2-M_periodicity.translation()) ) < 1e-10 )
                        {
                            // loop on each pair (element, lid) which
                            // has a global id gid2 and set it to gid
                            corresponding_gid = gid2;
                            match = true;
                            break;
                        }
                }
            // if we have --- actually we must have one --- a match, remove the
            // iterator from dof2 to quicken the search for the next dof1 match
            if ( match )
                {

                    it_dof2 = periodic_dof[periodicity_type::tag2].lower_bound( corresponding_gid );
                    en_dof2 = periodic_dof[periodicity_type::tag2].upper_bound( corresponding_gid );
                    while( it_dof2 != en_dof2 )
                        {

                            size_type ie = it_dof2->second.template get<0>();
                            size_type lid = it_dof2->second.template get<1>();
                            size_type gDof = it_dof2->second.template get<2>();

                            Debug( 5015 ) << "link " <<  boost::get<0>( _M_el_l2g[ ie][ lid ] )  << " -> " << gid << "\n";

                            // gid is given by dof1
                            _M_el_l2g[ ie][ lid ] = boost::make_tuple(gid, 1, true );


                            // warning: must modify the data structure that allows to
                            // generate unique global dof ids
                            map_gdof[ boost::make_tuple(0, 0, gDof) ] = gid;

                            ++it_dof2;
                        }
                    periodic_dof[periodicity_type::tag2].erase( it_dof2, en_dof2 );
                    periodic_dof_done[gid] =  true;
                }
            else
                {
                    // we have a problem, no match was found, this should not happen
                    Debug( 5015 ) << "[periodic] invalid point/dof matching\n";
                    Debug( 5015 ) << "[periodic] n = " << x1 << "\n";
                }

        }
    Debug( 5015 ) << "[periodic dof table] done matching the dof points\n";
    Debug( 5015 ) << "[periodic dof table] is empty : " << periodic_dof[periodicity_type::tag2].empty() << "\n";

    // ensure that periodic_dof[periodicity_type::tag2] is empty
    if ( !periodic_dof[periodicity_type::tag2].empty() )
        {
            Debug( 5015 ) << "[periodic] periodic conditions not set properly, some periodic dof were not assigned\n";
            typename periodic_dof_map_type::iterator it_dof2 = periodic_dof[periodicity_type::tag2].begin();
            typename periodic_dof_map_type::iterator en_dof2 = periodic_dof[periodicity_type::tag2].end();
            while( it_dof2 != en_dof2 )
                {

                    size_type ie = it_dof2->second.template get<0>();
                    size_type lid = it_dof2->second.template get<1>();

                    Debug( 5015 ) << "[periodic] dof " << it_dof2->first << " not assigned, "
                          << "x = " << periodic_dof_points[it_dof2->first].template get<0>() << " "
                          << "elt = " << ie << ", lid= " << lid << "\n";



                    ++it_dof2;
                }
        }
    else
        {
            Debug( 5015 ) << "[periodic] periodic condition done\n";
        }

    return max_gid+1;
}
template<typename MeshType, typename FEType, typename PeriodicityType>
void
Dof<MeshType, FEType, PeriodicityType>::buildDofMap( mesh_type& M, size_type start_next_free_dof )
{
    if ( !M_dof_indices.empty() )
        {
            generateDofPoints( M );
            return;
        }

    _M_n_el = M.numElements();

    size_type nldof =
        fe_type::nDofPerVolume * element_type::numVolumes +
        fe_type::nDofPerFace * element_type::numGeometricFaces +
        fe_type::nDofPerEdge * element_type::numEdges +
        fe_type::nDofPerVertex * element_type::numVertices;

    Debug(5005) << "==============================\n";
    Debug(5005) << "[buildDofMap]\n";
    Debug(5005) << "nldof                   = "  << int(nldof) << "\n";
    Debug(5005) << "fe_type::nLocalDof     = "  << int(fe_type::nLocalDof) << "\n";
    Debug(5005) << "fe_type::nDofPerVolume = "  << int(fe_type::nDofPerVolume) << "\n";
    Debug(5005) << "fe_type::nDofPerFace   = "  << int(fe_type::nDofPerFace) << "\n";
    Debug(5005) << "fe_type::nDofPerEdge   = "  << int(fe_type::nDofPerEdge) << "\n";
    Debug(5005) << "fe_type::nDofPerVertex = "  << int(fe_type::nDofPerVertex) << "\n";
    Debug(5005) << "element_type::numVolumes= "  << int(element_type::numVolumes) << "\n";
    Debug(5005) << "element_type::numFaces= "    << int(element_type::numFaces) << "\n";
    Debug(5005) << "element_type::numEdges= "    << int(element_type::numEdges) << "\n";
    Debug(5005) << "element_type::numVertices= " << int(element_type::numVertices) << "\n";
    Debug(5005) << "==============================\n";

    LIFE_ASSERT( nldof == fe_type::nLocalDof )
        ( nldof )
        ( fe_type::nLocalDof ).error( "Something wrong in FE specification" ) ;

    // initialize the local to global map and fill it with invalid
    // values that will allow to check whether we have a new dof or
    // not when building the table
    const size_type nV = M.numElements();
    _M_el_l2g.resize( boost::extents[nV][nComponents*nldof] );
    typedef Container::index index;
    for ( index i1 = 0; i1 < index(nV); ++i1 )
        for ( index i2 = 0; i2 < index(nComponents*nldof); ++i2 )
            _M_el_l2g[i1][i2] = boost::make_tuple(invalid_size_type_value,0,false); // 0 is the invalid value for the sign !

    const size_type n_proc  = Application::nProcess();

    _M_face_sign = ublas::scalar_vector<bool>(M.numFaces(), false);

    generateFacePermutations( M, mpl::bool_< ((Shape == SHAPE_TETRA && nOrder > 2 ) || (Shape == SHAPE_HEXA && nOrder > 1 ))>() );
    /* counter that stores the number dof per entities
       already registered this counter is used to shift
       the global dof index for each topological
       entities. We start from 0 and each topological
       entity (point,edge,face,element) will add there
       contribution to this counter to be then used the
       next topological dofs*/
    size_type gdofcount = 0;

    //! list of elements which have a periodic face Tag2
    std::list<std::pair<element_type const*, face_type const*> > periodic_elements;

    size_type next_free_dof = start_next_free_dof;
    for (size_type processor=0; processor<n_proc; processor++)
        {
            // compute the number of dof on current processor
            element_const_iterator it_elt = M.beginElementWithProcessId( processor );
            element_const_iterator en_elt = M.endElementWithProcessId( processor );
            size_type n_elts = std::distance( it_elt, en_elt);
            Debug( 5005 ) << "[buildDofMap] n_elts =  " << n_elts << " on processor " << processor << "\n";
            this->_M_first_df[processor] = next_free_dof;
            if ( is_periodic )
                this->_M_first_df[processor] =  0;
            it_elt = M.beginElementWithProcessId( processor );

            for( int c = 0; c < nComponents; ++c )
                {
                    Debug( 5005 ) << "[buildDofMap] component " << c << "\n";
                    it_elt = M.beginElementWithProcessId( processor );
                    Debug( 5005 ) << "[buildDofMap] starting with elt " << it_elt->id() << "\n";
                    for (;it_elt!=en_elt; ++it_elt )
                        {
                            element_type const& __elt = *it_elt;
                            Debug( 5005 ) << "next_free_dof " << next_free_dof  << "\n";
                            Debug( 5005 ) << "current dof " << dofIndex( next_free_dof ) << "\n";

                            /*
                             * Only in the continuous , we need to have the ordering [vertex,edge,face,volume]
                             */
                            if ( is_continuous )
                                {

                                    /* idem as above but for local element
                                       numbering except that it is
                                       reset to 0 after each element */
                                    uint16_type ldofcount = c*nldof;

                                    /* pack the shifts into a tuple */
                                    boost::tuple<uint16_type&,size_type&> shifts = boost::make_tuple( boost::ref(ldofcount),
                                                                                                      boost::ref(gdofcount) );

                                    /* \warning: the order of function calls is
                                       crucial here we order the degrees of freedom
                                       wrt the topological entities of the mesh
                                       elements from lowest dimension (vertex) to
                                       highest dimension (element)
                                    */
                                    addVertexDof( __elt, processor, c, next_free_dof, shifts  );
                                    addEdgeDof( __elt, processor, c, next_free_dof, shifts );
                                    addFaceDof( __elt, processor, c, next_free_dof, shifts );
                                    addVolumeDof( __elt, processor, c, next_free_dof, shifts );
                                }
                            else
                                {

                                    size_type ie = __elt.id();

                                    for ( uint16_type l = 0; l < nldof; ++l, ++next_free_dof )
                                        {
                                            _M_el_l2g[ ie][ fe_type::nLocalDof*c + l ] = boost::make_tuple(( dofIndex(next_free_dof)) , 1, false );
                                        }
                                }

                            // construct on the fly the associated dof points
#if 0
                            __c->update( __elt, __geopc );
                            for( uint16_type l =0; l < fe_type::nLocalDof; ++l )
                                {
                                    size_type thedof = boost::get<0>(localToGlobal( __elt->id(), l, c ));
                                    bool is_dof_periodic = boost::get<2>(localToGlobal( __elt->id(), l, c ));
                                    if ( ( thedof >= firstDof() ) && ( thedof <= lastDof() ) )
                                        {
                                            // get only the local dof
                                            //size_type thedofonproc = thedof - firstDof();
                                            thedof -= firstDof();
                                            if ( dof_done[ thedof ] == false )
                                                {
                                                    M_dof_points[thedof] = boost::make_tuple( __c->xReal( l ), firstDof()+thedof, c );

                                                    // store the periodic dofs
                                                    if ( is_dof_periodic )
                                                        M_dof_points_periodic[thedof] = boost::make_tuple( __c->xReal( l ), firstDof()+thedof, c );

                                                    dof_done[thedof] = true;
                                                }

                                        }
                                }
#endif // 0

                        } // elements loop
                } // nComponents loop

            // printing Dof table only in debug mode
#if !defined( NDEBUG )
            for( int c = 0; c < nComponents; ++c )
                {
                    it_elt = M.beginElementWithProcessId( processor );
                    for (;it_elt!=en_elt; ++it_elt )
                        {
                            element_type const& __elt = *it_elt;
                            std::ostringstream ostr;
                            ostr << "element " << __elt.id() << " : ";
                            for ( uint16_type l = 0; l < nldof; ++l )
                                {
                                    ostr << _M_el_l2g[ __elt.id()][ fe_type::nLocalDof*c + l ] << " ";
                                }
                            Debug( 5005 ) << ostr.str() << "\n";
                        }
                }
#endif
            this->_M_last_df[processor] = next_free_dof-1;
        }
    this->_M_n_dofs = next_free_dof;

    Debug( 5005 ) << " n global dof " << nDof() << "\n";
    Debug( 5005 ) << " n local dof " << nLocalDof() << "\n";
    for (size_type processor=0; processor<Application::nProcess(); processor++)
        {
            Debug( 5005 ) << "o processor " << processor << "\n";
            Debug( 5005 ) << "  - n dof on proc " << nDofOnProcessor(processor) << "\n";
            Debug( 5005 ) << "  - first dof " << firstDof(processor) << "\n";
            Debug( 5005 ) << "  - last dof " << lastDof(processor) << "\n";
        }

    //if ( is_continuous )
    //checkDofContinuity( M, mpl::int_<fe_type::nDim>() );

    generateDofPoints( M );
}

template<typename MeshType, typename FEType, typename PeriodicityType>
void
Dof<MeshType, FEType, PeriodicityType>::buildBoundaryDofMap( mesh_type& M )
{
    size_type nDofF = ( face_type::numVertices * fe_type::nDofPerVertex +
                        face_type::numEdges * fe_type::nDofPerEdge +
                        face_type::numFaces * fe_type::nDofPerFace );
    _M_n_dof_per_face_on_bdy = nDofF;
    Debug( 5005 ) << "number of Dof on an Element Face : " << nDofF << "\n";

    //
    // Face dof
    //
    typename mesh_type::face_const_iterator __face_it = M.beginFace();
    typename mesh_type::face_const_iterator __face_en = M.endFace();

    const size_type nF = M.faces().size();
    _M_face_l2g.resize( boost::extents[nF][nComponents*nDofF] );
    typedef Container::index index;
    global_dof_type default_dof = boost::make_tuple(invalid_size_type_value,0,false);
    for ( index i1 = 0; i1 < index(nF); ++i1 )
        for ( index i2 = 0; i2 < index(nComponents*nDofF); ++i2 )
            // 0 is the invalid value for the sign !
            _M_face_l2g[i1][i2] = default_dof;

    Debug( 5005 ) << "[buildBoundaryDofMap] nb faces : " << nF << "\n";
    Debug( 5005 ) << "[buildBoundaryDofMap] nb dof faces : " << nDofF*nComponents << "\n";

    for( int c = 0; c < nComponents; ++c )
        {
            __face_it = M.beginFace();
            for ( size_type nf = 0; __face_it != __face_en; ++__face_it, ++nf )
                {
                    LIFE_ASSERT( __face_it->isConnectedTo0() )
                        ( __face_it->id() )
                        ( __face_it->marker() )
                        ( __face_it->isConnectedTo0() )
                        ( __face_it->isConnectedTo1() ).error( "[Dof::buildFaceDofMap] face not connected" );
#if !defined(NDEBUG)
                    if (  __face_it->isOnBoundary() )
                        Debug( 5005 ) << "[buildBoundaryDofMap] boundary global face id : " << __face_it->id()
                                      << " marker: " << __face_it->marker()<< "\n";
                    else
                        Debug( 5005 ) << "[buildBoundaryDofMap] global face id : " << __face_it->id() << "\n";
#endif
                    uint16_type lc = c*nDofF;


                    addVertexBoundaryDof( __face_it, c, lc );
                    addEdgeBoundaryDof( __face_it, c, lc );
                    addFaceBoundaryDof( __face_it, c, lc );

                    LIFE_ASSERT( lc == (c+1)*nDofF )( lc )( c )( nDofF )( (c+1)*nDofF ).error( "invalid face local dof construction");
                }
        }
#if !defined(NDEBUG)
    for ( index i1 = 0; i1 < index(nF); ++i1 )
        for ( index i2 = 0; i2 < index(nComponents*nDofF); ++i2 )
            LIFE_ASSERT( boost::get<0>(_M_face_l2g[i1][i2]) != invalid_size_type_value )( i1 )( i2 ).error( "invalid dof table: initialized dof entries" );
#endif
}    // updateBoundaryDof

template<typename MeshType, typename FEType, typename PeriodicityType>
void
Dof<MeshType, FEType, PeriodicityType>::generateDofPoints(  mesh_type& M )
{
    if ( !M_dof_points.empty() )
        return;
    if ( fe_type::is_modal )
        return;
    Debug( 5005 ) << "[Dof::generateDofPoints] generating dof coordinates\n";
    typedef typename gm_type::template Context<vm::POINT, element_type> gm_context_type;
    typedef boost::shared_ptr<gm_context_type> gm_context_ptrtype;

    typedef typename fe_type::template Context<vm::POINT, fe_type, gm_type, element_type> fecontext_type;

    gm_ptrtype gm( new gm_type );
    fe_type fe;
    //
    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    //
    typename gm_type::precompute_ptrtype __geopc( new typename gm_type::precompute_type( gm, fe.points() ) );

    //const uint16_type ndofv = fe_type::nDof;

    element_const_iterator it_elt = M.beginElementWithProcessId( Application::processId() );
    element_const_iterator en_elt = M.endElementWithProcessId( Application::processId() );

    if ( it_elt == en_elt )
        return;

    gm_context_ptrtype __c( new gm_context_type( gm, *it_elt, __geopc ) );

    std::vector<bool> dof_done( nLocalDof() );
    M_dof_points.resize( nLocalDof() );
    std::fill( dof_done.begin(), dof_done.end(), false );
    for ( size_type dof_id = 0; it_elt!=en_elt ; ++it_elt )
    {
        __c->update( *it_elt, __geopc );
        for( uint16_type l =0; l < fe_type::nLocalDof; ++l )
            {
                for( uint16_type c1 = 0; c1 < nComponents; ++c1 )
                    {
                        size_type thedof = boost::get<0>(localToGlobal( it_elt->id(), l, c1 ));
                        if ( ( thedof >= firstDof() ) && ( thedof <= lastDof() ) )
                            {
                                // get only the local dof
                                //size_type thedofonproc = thedof - firstDof();
                                thedof -= firstDof();
                                LIFE_ASSERT( thedof < nLocalDof() )( thedof )( nLocalDof() )( firstDof() )( lastDof() )( it_elt->id() )( l )( c1 ).error ( "invalid local dof index");
                                if ( dof_done[ thedof ] == false )
                                    {
                                        std::set<uint16_type> lid;
                                        lid.insert( l );
                                        //M_dof_points[dof_id] = boost::make_tuple( thedof, __c->xReal( l ) );
                                        M_dof_points[thedof] = boost::make_tuple( __c->xReal( l ), firstDof()+thedof, c1 );
                                        dof_done[thedof] = true;
                                        ++dof_id;
                                    }
                            }
                    }
            }
    }
    for ( size_type dof_id = 0; dof_id < nLocalDof() ; ++dof_id )
        {
            LIFE_ASSERT( boost::get<1>(M_dof_points[dof_id]) >= firstDof() &&
                         boost::get<1>(M_dof_points[dof_id]) <= lastDof() )
                ( dof_id )( firstDof() )( lastDof() )( nLocalDof() )
                ( boost::get<1>(M_dof_points[dof_id]) )
                ( boost::get<0>(M_dof_points[dof_id]) ).error( "invalid dof point" );
            LIFE_ASSERT( dof_done[dof_id] == true )( dof_id ).error( "invalid dof point" );
        }
    Debug( 5005 ) << "[Dof::generateDofPoints] generating dof coordinates done\n";
}
template<typename MeshType, typename FEType, typename PeriodicityType>
void
Dof<MeshType, FEType, PeriodicityType>::generatePeriodicDofPoints(  mesh_type& M,
                                                                    periodic_element_list_type const& periodic_elements,
                                                                    dof_points_type& periodic_dof_points )
{
    if ( fe_type::is_modal )
        return;
    Debug( 5005 ) << "[Dof::generateDofPoints] generating dof coordinates\n";
    typedef typename gm_type::template Context<vm::POINT, element_type> gm_context_type;
    typedef boost::shared_ptr<gm_context_type> gm_context_ptrtype;

    typedef typename fe_type::template Context<vm::POINT, fe_type, gm_type, element_type> fecontext_type;

    gm_ptrtype gm( new gm_type );
    fe_type fe;
    //
    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    //
    typename gm_type::precompute_ptrtype __geopc( new typename gm_type::precompute_type( gm, fe.points() ) );

    //const uint16_type ndofv = fe_type::nDof;

    periodic_element_list_const_iterator it_elt = periodic_elements.begin();
    periodic_element_list_const_iterator en_elt = periodic_elements.end();

    if ( it_elt == en_elt )
        return;

    gm_context_ptrtype __c( new gm_context_type( gm, *it_elt->template get<0>(), __geopc ) );

    std::vector<bool> dof_done( periodic_dof_points.size() );
    std::fill( dof_done.begin(), dof_done.end(), false );
    for ( size_type dof_id = 0; it_elt!=en_elt ; ++it_elt )
    {
        __c->update( *it_elt->template get<0>(), __geopc );

        face_type const& __face = *it_elt->template get<1>();

        size_type iElAd = __face.ad_first();
        LIFE_ASSERT( iElAd != invalid_size_type_value )( __face.id() ).error( "[periodic]invalid face/element in face" );
        Life::detail::ignore_unused_variable_warning(iElAd);

        // local id of the face in its adjacent element
        uint16_type iFaEl = __face.pos_first();
        LIFE_ASSERT( iFaEl != invalid_uint16_type_value ).error ("invalid element index in face");

        // loop on face vertices
        for ( uint16_type iVeFa = 0; iVeFa < face_type::numVertices; ++iVeFa )
            {
                // local vertex number (in element)
                uint16_type iVeEl = element_type::fToP( iFaEl, iVeFa );
                Life::detail::ignore_unused_variable_warning(iVeEl);

                LIFE_ASSERT( iVeEl != invalid_uint16_type_value ).error( "invalid local dof" );

                // Loop number of Dof per vertex
                for ( uint16_type l = 0; l < fe_type::nDofPerVertex; ++l )
                    {
                        uint16_type lid = iVeEl * fe_type::nDofPerVertex + l;

                        size_type thedof = boost::get<0>(localToGlobal( it_elt->template get<0>()->id(), lid, 0 ));
                        LIFE_ASSERT( thedof < dof_done.size() )
                            ( thedof )
                            ( dof_done.size() )
                            ( it_elt->template get<0>()->id() )
                            ( lid ).error ("[generatePeriodicDofPoints] invalid dof id" );
                        if ( dof_done[ thedof ] == false )
                            {
                                periodic_dof_points[thedof] = boost::make_tuple( __c->xReal( lid ), thedof, 0 );
                                // these tests are problem specific x=0 and x=translation
#if 0
                                if ( __face.marker().value() == periodicity_type::tag1 )
                                    LIFE_ASSERT( math::abs( __c->xReal( lid )[0] ) < 1e-10 )( __c->xReal( lid ) ).warn( "[periodic] invalid p[eriodic point tag1");
                                if ( __face.marker().value() == periodicity_type::tag2 )
                                    LIFE_ASSERT( math::abs( __c->xReal( lid )[0] - M_periodicity.translation()[0] ) < 1e-10 )
                                        ( __c->xReal( lid ) )( M_periodicity.translation()).warn( "[periodic] invalid p[eriodic point tag1");
#endif
                                dof_done[thedof] = true;
                                ++dof_id;
                            }
                    }
            }
    }
    for ( size_type dof_id = 0; dof_id < periodic_dof_points.size() ; ++dof_id )
        {
            LIFE_ASSERT( boost::get<1>(periodic_dof_points[dof_id]) >= 0 &&
                         boost::get<1>(periodic_dof_points[dof_id]) < periodic_dof_points.size() )
                ( dof_id )( periodic_dof_points.size() )
                ( boost::get<1>(periodic_dof_points[dof_id]) )
                ( boost::get<0>(periodic_dof_points[dof_id]) ).error( "invalid dof point" );
            LIFE_ASSERT( dof_done[dof_id] == true )( dof_id ).error( "invalid dof point" );
        }
}

} // namespace Life
#endif

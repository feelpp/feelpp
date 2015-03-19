/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2015-03-06

  Copyright (C) 2015 Feel++ Consortium

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
#ifndef FEELPP_DOFFROMEDGE_HPP
#define FEELPP_DOFFROMEDGE_HPP 1

namespace Feel
{
/**
 * \brief Local Dof contribution from edge dof
 *
 * @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * @see DofTable, Dof, DofFromElement
 */
template <typename DofTableType, typename FEType>
class DofFromBoundary
{
public:

    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{

    typedef DofTableType doftable_type;
    using mesh_type = typename doftable_type::mesh_type;
    typedef typename doftable_type::element_type element_type;
    typedef typename doftable_type::face_type face_type;
    typedef typename doftable_type::ref_shift_type ref_shift_type;
    typedef typename doftable_type::localdof_type localdof_type;
    typedef FEType fe_type;

    typedef typename element_type::edge_permutation_type edge_permutation_type;

    static const uint16_type nOrder = fe_type::nOrder;
    static const uint16_type nDim = mesh_type::nDim;
    static const uint16_type nRealDim = mesh_type::nRealDim;
    static const uint16_type Shape = mesh_type::Shape;
    static const uint16_type nComponents = fe_type::nComponents;
    static const uint16_type nComponents1 = fe_type::nComponents1;
    static const uint16_type nComponents2 = fe_type::nComponents2;


    static const bool is_continuous = fe_type::isContinuous;
    static const bool is_discontinuous_locally = fe_type::continuity_type::is_discontinuous_locally;
    static const bool is_discontinuous_totally = fe_type::continuity_type::is_discontinuous_totally;

    static const bool is_scalar = fe_type::is_scalar;
    static const bool is_vectorial = fe_type::is_vectorial;
    static const bool is_tensor2 = fe_type::is_tensor2;
    static const bool is_modal = fe_type::is_modal;
    static const bool is_product = fe_type::is_product;

    static const bool is_p0_continuous = ( ( nOrder == 0 ) && is_continuous );

    static const uint16_type nDofPerElement = mpl::if_<mpl::bool_<is_product>, mpl::int_<fe_type::nLocalDof*nComponents1>, mpl::int_<fe_type::nLocalDof> >::type::value;

    //@}

    /** @name Constructors, destructor
     */
    //@{
    DofFromBoundary() = delete;
    //! copy constructor
    DofFromBoundary( DofFromBoundary const & ) = delete;
    //! copy operator
    DofFromBoundary& operator=( DofFromBoundary const & o) = delete;

    
    DofFromBoundary( doftable_type* doftable, fe_type const& fe )
        :
        M_doftable( doftable ),
        M_fe( fe )
        {}

    //! destructor
    ~DofFromBoundary() {}

    //@}

    /** @name Operator overloads
     */
    //@{

    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    template<typename Iterator>
    std::vector<size_type> edgeDof( Iterator it )
        {
            uint16_type lcVertex = 0;
            uint16_type lcEdge = 0;
            uint16_type lcFace = 0;

            std::vector<size_type> edge_dof;
            addVertexBoundaryDof( it, lcVertex, edge_dof);
            addEdgeBoundaryDof( it, lcEdge, edge_dof );
        }

    //@}

protected:

private:
    doftable_type * M_doftable;
    fe_type const& M_fe;

private:


    template<typename Iterator>
    void addVertexBoundaryDof( Iterator face_it, uint16_type& lc, std::vector<size_type>& d )
    {
        addVertexBoundaryDof( face_it, lc, d, mpl::bool_<(fe_type::nDofPerVertex>0)>(), mpl::int_<nDim>() );
    }
    template<typename Iterator> void addVertexBoundaryDof( Iterator /*face_it*/, bool /*useConnection0*/, uint16_type& /*lc*/, mpl::bool_<false>, mpl::int_<1> ) {}
    template<typename Iterator> void addVertexBoundaryDof( Iterator /*face_it*/, bool /*useConnection0*/, uint16_type& /*lc*/, mpl::bool_<false>, mpl::int_<2> ) {}
    template<typename Iterator> void addVertexBoundaryDof( Iterator /*face_it*/, bool /*useConnection0*/, uint16_type& /*lc*/, mpl::bool_<false>, mpl::int_<3> ) {}
    template<typename Iterator>
    void addVertexBoundaryDof( Iterator face_it, uint16_type& lc, std::vector<EdgeDof>& edof, mpl::bool_<true>, mpl::int_<1>  )
    {
        BOOST_STATIC_ASSERT( face_type::numVertices );

        uint16_type iFaEl = edge_it->element().id();
        size_type iElAd = edge_it->element().edgeId();

        // Loop number of Dof per vertex
        const int ncdof = is_product?nComponents:1;

        for ( int c = 0; c < ncdof; ++c )
        {
            for ( uint16_type l = 0; l < fe_type::nDofPerVertex; ++l, ++lc )
            {
                uint16_type ldinelt = iFaEl * fe_type::nDofPerVertex + l;
                auto const& temp= M_doftable->localToGlobal( iElAd, ldinelt, c );
                edof.push_back( EdgeDof( temp, lc, ldinelt ) );
            }
        }
    }
    template<typename Iterator>
    void addVertexBoundaryDof( Iterator face_it, uint16_type& lc, std::vector<EdgeDof>& edof, mpl::bool_<true>, mpl::int_<2>  )
    {
        BOOST_STATIC_ASSERT( face_type::numVertices );

        uint16_type iFaEl = edge_it->element().id();
        size_type iElAd = edge_it->element().edgeId();

        size_type ndofF = ( edge_type::numVertices * fe_type::nDofPerVertex +
                            edge_type::numEdges * fe_type::nDofPerEdge );
                            
        // loop on face vertices
        const int ncdof = is_product?nComponents:1;

        for ( int c = 0; c < ncdof; ++c )
        {
            for ( uint16_type iVeFa = 0; iVeFa < edge_type::numVertices; ++iVeFa )
            {
                // local vertex number (in element)
                uint16_type iVeEl = element_type::fToE( iFaEl, iVeFa );

                DCHECK( iVeEl != invalid_uint16_type_value ) <<  "invalid local dof";

                // Loop number of Dof per vertex
                for ( uint16_type l = 0; l < fe_type::nDofPerVertex; ++l )
                {
                    uint16_type ldinelt = iVeEl * fe_type::nDofPerVertex + l;
                    auto const& temp = M_doftable->localToGlobal( iElAd, ldinelt, c );
                    
                }
            }
        }
    }
    template<typename Iterator>
    void addVertexBoundaryDof( Iterator face_it, uint16_type& lc, std::vector<EdgeDof>& edof, mpl::bool_<true>, mpl::int_<3>  )
    {
        addVertexBoundaryDof( face_it, useConnection0, lc, mpl::bool_<true>(), mpl::int_<2>() );
    }
    template<typename Iterator>
    void addEdgeBoundaryDof( Iterator face_it, bool useConnection0, uint16_type& lc )
    {
        static const bool cond = fe_type::nDofPerEdge*face_type::numEdges > 0;
        addEdgeBoundaryDof( face_it, lc, std::vector<EdgeDof>& edof, mpl::bool_<cond>() );
    }
    template<typename Iterator>
    void addEdgeBoundaryDof( Iterator face_it, bool useConnection0, uint16_type& lc, mpl::bool_<false> )
        {}
    template<typename Iterator>
    void addEdgeBoundaryDof( Iterator face_it, bool useConnection0, uint16_type& lc, mpl::bool_<true> )
    {
        uint16_type iFaEl;
        size_type iElAd;
        DCHECK( iFaEl != invalid_uint16_type_value ) << "invalid element index in face";

        size_type nVerticesF = edge_type::numVertices * fe_type::nDofPerVertex;
        size_type ndofF = ( edge_type::numVertices * fe_type::nDofPerVertex +
                            edge_type::numEdges * fe_type::nDofPerEdge );

        const int ncdof = is_product?nComponents:1;

        for ( int c = 0; c < ncdof; ++c )
        {
            uint16_type lcc=nVerticesF+c*ndofF;

            // Loop number of Dof per edge
            for ( uint16_type l = 0; l < fe_type::nDofPerEdge; ++l, ++lcc )
            {
                uint16_type ldinelt = element_type::numVertices*fe_type::nDofPerVertex +
                    iFaEl * fe_type::nDofPerEdge + l ;
                auto const& temp = M_doftable->localToGlobal( iElAd,ldinelt, c );
                edof.push_back( EdgeDof( temp, lcc, ldinelt ) );
            }
        }
    }

}; 
}
#endif /* FEELPP_DofFromBoundary_H */

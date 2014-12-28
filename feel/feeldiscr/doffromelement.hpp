/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-10-23

  Copyright (C) 2013 Université de Strasbourg

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
   \file doflocal.hpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-10-23
 */
#ifndef FEELPP_DofFromElement_H
#define FEELPP_DofFromElement_H 1

#include <feel/feelpoly/hdivpolynomialset.hpp>
#include <feel/feelpoly/hcurlpolynomialset.hpp>
#include <feel/feeldiscr/dof.hpp>

namespace Feel
{
/**
 * \brief local dof contribution from an element
 *
 * @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * @see DofTable, Dof
 */
template <typename DofTableType, typename FEType>
class DofFromElement
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
    typedef typename doftable_type::mesh_type mesh_type;
    typedef typename doftable_type::element_type element_type;
    typedef typename doftable_type::face_type face_type;
    typedef typename doftable_type::ref_shift_type ref_shift_type;
    typedef typename doftable_type::localdof_type localdof_type;
    typedef FEType fe_type;
    typedef typename doftable_type::dof_relation dof_relation;
    typedef typename element_type::edge_permutation_type edge_permutation_type;
    typedef typename element_type::face_permutation_type face_permutation_type;


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

    DofFromElement( doftable_type* doftable, fe_type const& fe ): M_doftable(doftable), M_fe( fe ) {}

    //! destructor
    ~DofFromElement() {}

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

    void add( element_type const& __elt,
              size_type& next_free_dof,
              rank_type processor = 0,
              size_type shift = 0 );

    //@}



protected:

private:
    doftable_type* M_doftable;
    fe_type const& M_fe;
private:
    //! default constructor
    DofFromElement();
    //! copy constructor
    DofFromElement( DofFromElement const & );
    //! copy operator
    DofFromElement& operator=( DofFromElement const & o)
        {
            if (this != &o )
            {
            }
            return *this;
        }

    void addVertexDof( element_type const& __elt, rank_type processor, size_type& next_free_dof,
                       ref_shift_type& shifts  )
    {
        addVertexDof( __elt, processor, next_free_dof, shifts, mpl::bool_<(fe_type::nDofPerVertex>0)>() );
    }
    void addVertexDof( element_type const& /*M*/, rank_type /*processor*/,  size_type& /*next_free_dof*/,
                       ref_shift_type& /*shifts*/, mpl::bool_<false> )
    {}
    void addVertexDof( element_type const& __elt, rank_type processor, size_type& next_free_dof,
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
                M_doftable->insertDof( ie, lc, i, std::make_tuple(  0, gDof ),
                                 processor, next_free_dof, 1, false, global_shift, __elt.point( i ).marker() );
            }
        }

        // update shifts
        shifts.template get<0>() = lc;

#if !defined(NDEBUG)
        DVLOG(4) << "[Dof::updateVolumeDof(addVertexDof] vertex proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
    }
    void addEdgeDof( element_type const& __elt, rank_type processor, size_type& next_free_dof,
                     ref_shift_type& shifts )
    {
        static const bool cond = fe_type::nDofPerEdge > 0;
        return addEdgeDof( __elt,
                           processor,
                           next_free_dof,
                           shifts,
                           mpl::int_<fe_type::nDim>(),
                           mpl::bool_<cond>() );
    }
    void addEdgeDof( element_type const& /*M*/, rank_type /*processor*/, size_type& /*next_free_dof*/,
                     ref_shift_type& /*shifts*/, mpl::int_<1>, mpl::bool_<false> )
    {}

    void addEdgeDof( element_type const& __elt, rank_type processor, size_type& next_free_dof,
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
            M_doftable->insertDof( ie, lc, l, std::make_tuple(  1, gDof ), processor, next_free_dof, 1, false, global_shift, __elt.marker() );
        }

        // update shifts
        shifts.template get<0>() = lc;
#if !defined(NDEBUG)
        DVLOG(4) << "[Dof::addEdgeDof(1)] element proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
    }
    void addEdgeDof( element_type const& /*__elt*/, rank_type /*processor*/, size_type& /*next_free_dof*/,
                     ref_shift_type& /*shifts*/, mpl::int_<2>, mpl::bool_<false> )
    {}
    void addEdgeDof( element_type const& __elt, rank_type processor, size_type& next_free_dof,
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


                if ( __elt.edgePermutation( i ).value()  == edge_permutation_type::IDENTITY )
                {
                    gDof += l ; // both nodal and modal case
                    if ( is_hdiv_conforming<fe_type>::value || is_hcurl_conforming<fe_type>::value )
                    {
                        
                        M_doftable->M_locglob_signs[ie][lc] = 1;
                    }
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
                        gDof += fe_type::nDofPerEdge - 1 - l ;
                    if ( is_hdiv_conforming<fe_type>::value || is_hcurl_conforming<fe_type>::value )
                    {
                        sign = -1;
                        M_doftable->M_locglob_signs[ie][lc] = -1;
                    }

                }

                else
                    FEELPP_ASSERT( 0 ).error ( "invalid edge permutation" );

                M_doftable->insertDof( ie, lc, i, std::make_tuple(  1, gDof ), processor, next_free_dof, sign, false, global_shift, __elt.edge( i ).marker() );
            }
        }

        // update shifts
        shifts.template get<0>() = lc;
#if !defined(NDEBUG)
        DVLOG(4) << "[Dof::addEdgeDof] edge proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
    }

    void addEdgeDof( element_type const& /*__elt*/, rank_type /*processor*/, size_type& /*next_free_dof*/,
                     ref_shift_type& /*shifts*/, mpl::int_<3>, mpl::bool_<false> )
    {}

    void addEdgeDof( element_type const& __elt, rank_type processor, size_type& next_free_dof,
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

                if ( __elt.edgePermutation( i ).value()  == edge_permutation_type::IDENTITY )
                {
                    gDof += l ; // both nodal and modal case
                    if ( is_hcurl_conforming<fe_type>::value )
                    {
                        M_doftable->M_locglob_signs[ie][lc] = 1;
                    }
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
                        gDof += fe_type::nDofPerEdge - 1 - l ;
                    if( is_hcurl_conforming<fe_type>::value )
                    {
                        sign = -1;
                        M_doftable->M_locglob_signs[ie][lc] = -1;
                    }
                }

                else
                    FEELPP_ASSERT( 0 ).error ( "invalid edge permutation" );

                M_doftable->insertDof( ie, lc, i, std::make_tuple(  1, gDof ), processor, next_free_dof, sign, false, global_shift, __elt.edge( i ).marker() );
            }
        }

        // update shifts
        shifts.template get<0>() = lc;
#if !defined(NDEBUG)
        DVLOG(4) << "[Dof::addEdgeDof] edge proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
    }


    void addFaceDof( element_type const& __elt, rank_type processor, size_type& next_free_dof,
                     ref_shift_type& shifts )
    {
        return addFaceDof( __elt, processor, next_free_dof, shifts, mpl::int_<fe_type::nDim>(), mpl::bool_<(fe_type::nDofPerFace > 0)>() );
    }
    void addFaceDof( element_type const& /*M*/, rank_type /*processor*/, size_type& /*next_free_dof*/,
                     ref_shift_type& /*shifts*/, mpl::int_<1>, mpl::bool_<false> )
    {}
    void addFaceDof( element_type const& /*M*/, rank_type /*processor*/, size_type& /*next_free_dof*/,
                     ref_shift_type& /*shifts*/, mpl::int_<2>, mpl::bool_<false> )
    {}
    void addFaceDof( element_type const& __elt, rank_type processor, size_type& next_free_dof,
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
            M_doftable->insertDof( ie, lc, l, std::make_tuple(  2, gDof ), processor, next_free_dof, 1, false, global_shift, __elt.marker() );
        }

        // update shifts
        shifts.template get<0>() = lc;
#if !defined(NDEBUG)
        DVLOG(4) << "[Dof::addFaceDof(2,true)] face proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
    }
    void addFaceDof( element_type const& /*M*/, rank_type /*processor*/, size_type& /*next_free_dof*/,
                     ref_shift_type& /*shifts*/, mpl::int_<3>, mpl::bool_<false> )
    {}
    void addFaceDof( element_type const& __elt, rank_type processor, size_type& next_free_dof,
                     ref_shift_type& shifts, mpl::int_<3>, mpl::bool_<true> )
    {
        uint16_type local_shift;
        size_type global_shift;
        boost::tie( local_shift, global_shift ) = shifts;

        size_type ie = __elt.id();

        uint16_type lc = local_shift;

        for ( uint16_type i = 0; i < element_type::numFaces; ++i )
        {
            face_permutation_type permutation = __elt.facePermutation( i );
            FEELPP_ASSERT( permutation != face_permutation_type( 0 ) ).error ( "invalid face permutation" );

            // Polynomial order in each direction
            uint16_type p=1;
            uint16_type q=0;

            // MaxOrder = Order - 2
            int MaxOrder = int( ( 3 + std::sqrt( 1+8*fe_type::nDofPerFace ) )/2 ) - 2;

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

                if ( q > MaxOrder )
                {
                    q = 1;
                    p = p+1;
                    MaxOrder = MaxOrder-1;
                }

                if ( !fe_type::is_modal )
                {
                    if( is_hdiv_conforming<fe_type>::value || is_hcurl_conforming<fe_type>::value )
                    {
                        if( fe_type::nDofPerFace == 1 )
                            gDof += l;
                        else
                            gDof += M_doftable->vector_permutation[permutation][l];

                        if (permutation  == face_permutation_type( 1 ))
                            M_doftable->M_locglob_signs[ie][l] = 1;
                        else
                        {
                            sign=-1;
                            M_doftable->M_locglob_signs[ie][l] = -1;
                        }
                    }
                    else
                    {
                        // no need of permutation is identity or only one dof on face
                        if ( permutation  == face_permutation_type( 1 ) || fe_type::nDofPerFace == 1 )
                            gDof += l;
                        else
                            gDof += M_doftable->vector_permutation[permutation][l];
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

                M_doftable->insertDof( ie, lc, i, std::make_tuple(  2, gDof ), processor, next_free_dof, sign, false, global_shift,__elt.face( i ).marker() );

            }
        }

        // update shifts
        shifts.template get<0>() = lc;
#if !defined(NDEBUG)
        DVLOG(4) << "[Dof::addFaceDof<3>] face proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
    }
    void addVolumeDof( element_type const& __elt, rank_type processor, size_type& next_free_dof,
                       ref_shift_type& shifts )
    {
        return addVolumeDof( __elt, processor, next_free_dof, shifts, mpl::bool_<(fe_type::nDofPerVolume>0)>() );
    }
    void addVolumeDof( element_type const& /*M*/, rank_type /*processor*/, size_type& /*next_free_dof*/,
                       ref_shift_type& /*shifts*/, mpl::bool_<false> )
    {}
    void addVolumeDof( element_type const& __elt, rank_type processor, size_type& next_free_dof,
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
            M_doftable->insertDof( ie, lc, l, std::make_tuple(  3, gDof ), processor, next_free_dof, 1, false, global_shift, __elt.marker() );
        }

        // update shifts
        shifts.template get<0>() = lc;
#if !defined(NDEBUG)
        DVLOG(4) << "[Dof::updateVolumeDof(<2>)] element proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
    }

};

template <typename DofTableType,typename FEType>
void
DofFromElement<DofTableType,FEType>::add( element_type const& __elt,
                                          size_type& next_free_dof,
                                          rank_type processor,
                                          size_type shift )
{

    size_type nldof = M_doftable->nLocalDof( true );

    DVLOG(3) << "adding dof from element " << __elt.id() << "\n";
    size_type gdofcount = shift;
    DVLOG(3) << "next_free_dof " << next_free_dof  << "\n";
    DVLOG(3) << "current dof " << M_doftable->dofIndex( next_free_dof ) << "\n";


    /*
     * Only in the continuous , we need to have the ordering [vertex,edge,face,volume]
     */
    if ( is_continuous || is_discontinuous_locally )
    {

        /* idem as above but for local element
           numbering except that it is
           reset to 0 after each element */
        uint16_type ldofcount = 0;

        /* pack the shifts into a tuple */
        boost::tuple<uint16_type&,size_type&> shifts = boost::make_tuple( boost::ref( ldofcount ),
                                                                          boost::ref( gdofcount ) );

        /* \warning: the order of function calls is
           crucial here we order the degrees of freedom
           wrt the topological entities of the mesh
           elements from lowest dimension (vertex) to
           highest dimension (element)
        */
        addVertexDof( __elt, processor, next_free_dof, shifts  );
        addEdgeDof( __elt, processor, next_free_dof, shifts );
        addFaceDof( __elt, processor, next_free_dof, shifts );
        addVolumeDof( __elt, processor, next_free_dof, shifts );
    }

    else
    {

        size_type ie = __elt.id();

        const int ncdof = is_product?nComponents:1;

        for ( uint16_type l = 0; l < nldof; ++l )
        {
            for ( int c = 0; c < ncdof; ++c, ++next_free_dof )
            {
                M_doftable->M_el_l2g.insert( dof_relation( localdof_type( ie, fe_type::nLocalDof*c + l ),
                                                           Dof( ( M_doftable->dofIndex( next_free_dof ) ) , 1, false ) ) );
            }
        }
    }
}
}
#endif /* FEELPP_DofFromElement_H */

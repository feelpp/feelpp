/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-10-23

  Copyright (C) 2013 Feel++ Consortium

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
#ifndef FEELPP_DofFromMortar_H
#define FEELPP_DofFromMortar_H 1

namespace Feel
{
/**
 * \brief local dof contribution from a mortar element
 *
 * @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * @see DofTable, Dof
 */
template <typename DofTableType, typename MortarFEType, typename FEType>
class DofFromMortar
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
    typedef MortarFEType mortar_fe_type;
    typedef FEType fe_type;
    typedef typename doftable_type::dof_relation dof_relation;
    typedef typename element_type::edge_permutation_type edge_permutation_type;
    typedef typename element_type::face_permutation_type face_permutation_type;

    static const uint16_type nOrder = mortar_fe_type::nOrder;
    static const uint16_type nDim = mesh_type::nDim;
    static const uint16_type nRealDim = mesh_type::nRealDim;
    static const uint16_type Shape = mesh_type::Shape;
    static const uint16_type nComponents = mortar_fe_type::nComponents;
    static const uint16_type nComponents1 = mortar_fe_type::nComponents1;
    static const uint16_type nComponents2 = mortar_fe_type::nComponents2;


    static const bool is_continuous = mortar_fe_type::isContinuous;
    static const bool is_discontinuous_locally = mortar_fe_type::continuity_type::is_discontinuous_locally;
    static const bool is_discontinuous_totally = mortar_fe_type::continuity_type::is_discontinuous_totally;

    static const bool is_scalar = mortar_fe_type::is_scalar;
    static const bool is_vectorial = mortar_fe_type::is_vectorial;
    static const bool is_tensor2 = mortar_fe_type::is_tensor2;
    static const bool is_modal = mortar_fe_type::is_modal;
    static const bool is_product = mortar_fe_type::is_product;

    static const bool is_p0_continuous = ( ( nOrder == 0 ) && is_continuous );

    static const uint16_type nDofPerElement = mpl::if_<mpl::bool_<is_product>, mpl::int_<mortar_fe_type::nLocalDof*nComponents1>, mpl::int_<mortar_fe_type::nLocalDof> >::type::value;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    DofFromMortar( doftable_type* doftable, mortar_fe_type const& mfe, fe_type const& fe )
        :
        M_doftable(doftable),
        M_mortar_fe( mfe ),
        M_fe( fe)
        {}

    //! destructor
    ~DofFromMortar() {}

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

    void add( element_type const& elt,
              size_type& next_free_dof,
              rank_type processor = 0,
              size_type shift = 0 );

    //@}



protected:

private:
    doftable_type* M_doftable;
    mortar_fe_type const& M_mortar_fe;
    fe_type const& M_fe;

private:
    //! default constructor
    DofFromMortar();
    //! copy constructor
    DofFromMortar( DofFromMortar const & );
    //! copy operator
    DofFromMortar& operator=( DofFromMortar const & o)
        {
            if (this != &o )
            {
            }
            return *this;
        }

    void addVertexDof( element_type const& elt, rank_type processor, size_type& next_free_dof,
                       ref_shift_type& shifts  )
        {
            uint16_type local_shift;
            size_type global_shift;
            boost::tie( local_shift, global_shift ) = shifts;

            size_type ie = elt.id();

            uint16_type lc = local_shift;
            auto n = elt.neighbor( 0 );
            if ( n.first  != invalid_size_type_value )
            {
                size_type gDof = M_doftable->mesh()->element( n.first ).point(1).id()*fe_type::nDofPerVertex;
                DVLOG(2) << "inserting vertex dof " << gDof << "," << next_free_dof << "," << ie;
                M_doftable->insertDof( ie, lc++, 0, std::make_tuple(  0, gDof ), processor, next_free_dof, 1, false, global_shift, elt.point(0).marker() );
                gDof = ( elt.point( 1 ).id() ) * mortar_fe_type::nDofPerVertex;
                M_doftable->insertDof( ie, lc++, 1, std::make_tuple(  0, gDof ),
                                       processor, next_free_dof, 1, false, global_shift, elt.point( 1 ).marker() );
                DVLOG(2) << "inserting vertex dof " << gDof << "," << next_free_dof << "," << ie;
            }
            else
            {
                n = elt.neighbor( 1 );
                CHECK( n.first != invalid_size_type_value ) << "the element should be connected to at least one other element, it is not the case";

                size_type gDof = ( elt.point( 0 ).id() ) * mortar_fe_type::nDofPerVertex;
                DVLOG(2) << "inserting vertex dof " << gDof << "," << next_free_dof << "," << ie;
                M_doftable->insertDof( ie, lc++, 0, std::make_tuple(  0, gDof ),
                                       processor, next_free_dof, 1, false, global_shift, elt.point( 0 ).marker() );

                gDof = M_doftable->mesh()->element( n.first ).point(0).id()*fe_type::nDofPerVertex;
                DVLOG(2) << "inserting vertex dof " << gDof << "," << next_free_dof << "," << ie;
                M_doftable->insertDof( ie, lc++, 1, std::make_tuple(  0, gDof ), processor, next_free_dof, 1, false, global_shift, elt.point(1).marker() );
            }





            // update shifts
            shifts.template get<0>() = lc;

#if !defined(NDEBUG)
            DVLOG(4) << "[Dof::updateVolumeDof(addVertexDof] vertex proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
        }
    void addEdgeDof( element_type const& elt, rank_type processor, size_type& next_free_dof,
                     ref_shift_type& shifts )
        {
            //static const bool cond = mortar_fe_type::nDofPerEdge > 0 && mortar_fe_type::nOrder > 1;
            return addEdgeDof( elt,
                               processor,
                               next_free_dof,
                               shifts,
                               mpl::int_<mortar_fe_type::nDim>(),
                               mpl::bool_<true>() );
        }
    void addEdgeDof( element_type const& /*M*/, rank_type /*processor*/, size_type& /*next_free_dof*/,
                     ref_shift_type& /*shifts*/, mpl::int_<1>, mpl::bool_<false> )
        {}

    void addEdgeDof( element_type const& elt, rank_type processor, size_type& next_free_dof,
                     ref_shift_type& shifts, mpl::int_<1>, mpl::bool_<true> )
        {
            uint16_type local_shift;
            size_type global_shift;
            boost::tie( local_shift, global_shift ) = shifts;

            size_type ie = elt.id();
            uint16_type lc = local_shift;
            DVLOG(2) << "adding mortar dof on edge " << mortar_fe_type::nDofPerEdge << " nOrder = " << nOrder;
            for ( uint16_type l = 0; l < mortar_fe_type::nDofPerEdge; ++l, ++lc )
            {
                //size_type gDof = ie * mortar_fe_type::nDofPerEdge + l;
                size_type gDof = ie * fe_type::nDofPerEdge + l;
                if ( nOrder == 0 )//is_p0_continuous )
                {
                    auto n = elt.neighbor( 0 );
                    if ( n.first  != invalid_size_type_value )
                    {
                        gDof = M_doftable->mesh()->element( n.first ).point(1).id()*fe_type::nDofPerVertex;
                    }
                    else
                    {
                        n = elt.neighbor( 1 );
                        CHECK( n.first != invalid_size_type_value ) << "the element should be connected to at least one other element, it is not the case";
                        gDof = M_doftable->mesh()->element( n.first ).point(0).id()*fe_type::nDofPerVertex;
                    }
                    DVLOG(2) << "inserting dof " << gDof << "," << next_free_dof << "," << ie;
                    M_doftable->insertDof( ie, lc, l, std::make_tuple(  0, gDof ), processor, next_free_dof, 1, false, global_shift, elt.marker() );
                }
                else
                {
                    M_doftable->insertDof( ie, lc, l, std::make_tuple(  1, gDof ), processor, next_free_dof, 1, false, global_shift, elt.marker() );
                }
            }

            // update shifts
            shifts.template get<0>() = lc;
#if !defined(NDEBUG)
            DVLOG(4) << "[Dof::addEdgeDof(1)] element proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
        }
    void addEdgeDof( element_type const& /*elt*/, rank_type /*processor*/, size_type& /*next_free_dof*/,
                     ref_shift_type& /*shifts*/, mpl::int_<2>, mpl::bool_<false> )
        {}
    void addEdgeDof( element_type const& elt, rank_type processor, size_type& next_free_dof,
                     ref_shift_type& shifts, mpl::int_<2>, mpl::bool_<true> )
        {
            uint16_type local_shift;
            size_type global_shift;
            boost::tie( local_shift, global_shift ) = shifts;

            size_type ie = elt.id();
            uint16_type lc = local_shift;

            /** The boundary dofs are constructed in the same way if the basis is modal **/

            for ( uint16_type i = 0; i < element_type::numEdges; ++i )
            {
                for ( uint16_type l = 0; l < mortar_fe_type::nDofPerEdge; ++l, ++lc )
                {
                    size_type gDof = elt.edge( i ).id() * mortar_fe_type::nDofPerEdge;
                    int32_type sign = 1;


                    if ( elt.edgePermutation( i ).value()  == edge_permutation_type::IDENTITY )
                    {
                        gDof += l ; // both nodal and modal case
                    }

                    else if ( elt.edgePermutation( i ).value()  == edge_permutation_type::REVERSE_PERMUTATION )
                    {

                        if ( mortar_fe_type::is_modal )
                        {
                            //only half of the modes (odd polynomial order) are negative.
                            sign = ( l%2 )?( -1 ):( 1 );
                            gDof += l;
                        }

                        else
                            gDof += mortar_fe_type::nDofPerEdge - 1 - l ;
                    }

                    else
                        FEELPP_ASSERT( 0 ).error ( "invalid edge permutation" );

                    M_doftable->insertDof( ie, lc, i, std::make_tuple(  1, gDof ), processor, next_free_dof, sign, false, global_shift, elt.edge( i ).marker() );
                }
            }

            // update shifts
            shifts.template get<0>() = lc;
#if !defined(NDEBUG)
            DVLOG(4) << "[Dof::addEdgeDof] edge proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
        }

    void addEdgeDof( element_type const& /*elt*/, rank_type /*processor*/, size_type& /*next_free_dof*/,
                     ref_shift_type& /*shifts*/, mpl::int_<3>, mpl::bool_<false> )
        {}

    void addEdgeDof( element_type const& elt, rank_type processor, size_type& next_free_dof,
                     ref_shift_type& shifts, mpl::int_<3>, mpl::bool_<true> )
        {
            uint16_type local_shift;
            size_type global_shift;
            boost::tie( local_shift, global_shift ) = shifts;

            size_type ie = elt.id();
            uint16_type lc = local_shift;

            for ( uint16_type i = 0; i < element_type::numEdges; ++i )
            {
                for ( uint16_type l = 0; l < mortar_fe_type::nDofPerEdge; ++l, ++lc )
                {
                    size_type gDof = elt.edge( i ).id() * mortar_fe_type::nDofPerEdge;

                    int32_type sign = 1;

                    if ( elt.edgePermutation( i ).value()  == edge_permutation_type::IDENTITY )
                    {
                        gDof += l ; // both nodal and modal case
                    }

                    else if ( elt.edgePermutation( i ).value()  == edge_permutation_type::REVERSE_PERMUTATION )
                    {

                        if ( mortar_fe_type::is_modal )
                        {
                            //only half of the modes (odd polynomial order) are negative.
                            sign = ( l%2 )?( -1 ):( 1 );
                            gDof += l;
                        }

                        else
                            gDof += mortar_fe_type::nDofPerEdge - 1 - l ;
                    }

                    else
                        FEELPP_ASSERT( 0 ).error ( "invalid edge permutation" );

                    M_doftable->insertDof( ie, lc, i, std::make_tuple(  1, gDof ), processor, next_free_dof, sign, false, global_shift, elt.edge( i ).marker() );
                }
            }

            // update shifts
            shifts.template get<0>() = lc;
#if !defined(NDEBUG)
            DVLOG(4) << "[Dof::addEdgeDof] edge proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
        }


    void addFaceDof( element_type const& elt, rank_type processor, size_type& next_free_dof,
                     ref_shift_type& shifts )
        {
            return addFaceDof( elt, processor, next_free_dof, shifts, mpl::int_<mortar_fe_type::nDim>(), mpl::bool_<(mortar_fe_type::nDofPerFace > 0)>() );
        }
    void addFaceDof( element_type const& /*M*/, rank_type /*processor*/, size_type& /*next_free_dof*/,
                     ref_shift_type& /*shifts*/, mpl::int_<1>, mpl::bool_<false> )
        {}
    void addFaceDof( element_type const& /*M*/, rank_type /*processor*/, size_type& /*next_free_dof*/,
                     ref_shift_type& /*shifts*/, mpl::int_<2>, mpl::bool_<false> )
        {}
    void addFaceDof( element_type const& elt, rank_type processor, size_type& next_free_dof,
                     ref_shift_type& shifts, mpl::int_<2>, mpl::bool_<true> )
        {
            uint16_type local_shift;
            size_type global_shift;
            boost::tie( local_shift, global_shift ) = shifts;

            size_type ie = elt.id();
            uint16_type lc = local_shift;

            for ( uint16_type l = 0; l < mortar_fe_type::nDofPerFace; ++l, ++lc )
            {
                const size_type gDof = is_p0_continuous? l:ie * mortar_fe_type::nDofPerFace + l;
                M_doftable->insertDof( ie, lc, l, std::make_tuple(  2, gDof ), processor, next_free_dof, 1, false, global_shift, elt.marker() );
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
    void addFaceDof( element_type const& elt, rank_type processor, size_type& next_free_dof,
                     ref_shift_type& shifts, mpl::int_<3>, mpl::bool_<true> )
        {
            uint16_type local_shift;
            size_type global_shift;
            boost::tie( local_shift, global_shift ) = shifts;

            size_type ie = elt.id();

            uint16_type lc = local_shift;

            for ( uint16_type i = 0; i < element_type::numFaces; ++i )
            {
                face_permutation_type permutation = elt.facePermutation( i );
                FEELPP_ASSERT( permutation != face_permutation_type( 0 ) ).error ( "invalid face permutation" );

                // Polynomial order in each direction
                uint16_type p=1;
                uint16_type q=0;

                // MaxOrder = Order - 2
                int MaxOrder = int( ( 3 + std::sqrt( 1+8*mortar_fe_type::nDofPerFace ) )/2 ) - 2;

                for ( uint16_type l = 0; l < mortar_fe_type::nDofPerFace; ++l, ++lc )
                {

                    // TODO: orient the dof indices such
                    // that they match properly the faces
                    // dof of the connected faces. There
                    // are a priori many permutations of
                    // the dof face indices
                    size_type gDof = elt.face( i ).id() * mortar_fe_type::nDofPerFace;
                    int32_type sign = 1;

                    q=q+1;

                    if ( q > MaxOrder )
                    {
                        q = 1;
                        p = p+1;
                        MaxOrder = MaxOrder-1;
                    }

                    if ( !mortar_fe_type::is_modal )
                    {
                        // no need of permutation is identity or only one dof on face
                        if ( permutation  == face_permutation_type( 1 ) || mortar_fe_type::nDofPerFace == 1 )
                            gDof += l;

                        else
                            gDof += M_doftable->vector_permutation[permutation][l];
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

                    M_doftable->insertDof( ie, lc, i, std::make_tuple(  2, gDof ), processor, next_free_dof, sign, false, global_shift,elt.face( i ).marker() );

                }
            }

            // update shifts
            shifts.template get<0>() = lc;
#if !defined(NDEBUG)
            DVLOG(4) << "[Dof::addFaceDof<3>] face proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
        }
    void addVolumeDof( element_type const& elt, rank_type processor, size_type& next_free_dof,
                       ref_shift_type& shifts )
        {
            return addVolumeDof( elt, processor, next_free_dof, shifts, mpl::bool_<(mortar_fe_type::nDofPerVolume>0)>() );
        }
    void addVolumeDof( element_type const& /*M*/, rank_type /*processor*/, size_type& /*next_free_dof*/,
                       ref_shift_type& /*shifts*/, mpl::bool_<false> )
        {}
    void addVolumeDof( element_type const& elt, rank_type processor, size_type& next_free_dof,
                       ref_shift_type& shifts, mpl::bool_<true> )
        {
            BOOST_STATIC_ASSERT( element_type::numVolumes );
            uint16_type local_shift;
            size_type global_shift;
            boost::tie( local_shift, global_shift ) = shifts;

            size_type ie = elt.id();
            uint16_type lc = local_shift;

            for ( uint16_type l = 0; l < mortar_fe_type::nDofPerVolume; ++l, ++lc )
            {
                const size_type gDof = is_p0_continuous? l:ie * mortar_fe_type::nDofPerVolume + l;
                M_doftable->insertDof( ie, lc, l, std::make_tuple(  3, gDof ), processor, next_free_dof, 1, false, global_shift, elt.marker() );
            }

            // update shifts
            shifts.template get<0>() = lc;
#if !defined(NDEBUG)
            DVLOG(4) << "[Dof::updateVolumeDof(<2>)] element proc" << processor << " next_free_dof = " << next_free_dof << "\n";
#endif
        }

};

template <typename DofTableType, typename MortarFEType, typename FEType>
void
DofFromMortar<DofTableType,MortarFEType,FEType>::add( element_type const& elt,
                                                      size_type& next_free_dof,
                                                      rank_type processor,
                                                      size_type shift )
{
#if 0
    CHECK( elt.nNeighbors() == 1 ) << "Invalid mortar element, the number of neighbors should be exactly 1 in 1D"
                                   << " number of neighbors : " << elt.nNeighbors()
                                   << " element id : " << elt.id();
#endif
    // now get the face id

    size_type nldof = M_doftable->nLocalDof( true );

    DVLOG(3) << "adding dof from element " << elt.id() << "\n";
    size_type gdofcount = shift;
    DVLOG(3) << "next_free_dof " << next_free_dof  << "\n";
    DVLOG(3) << "current dof " << M_doftable->dofIndex( next_free_dof ) << "\n";


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
    if ( mortar_fe_type::nDofPerVertex > 0 )
        addVertexDof( elt, processor, next_free_dof, shifts  );
    if ( mortar_fe_type::nDofPerEdge > 0 )
        addEdgeDof( elt, processor, next_free_dof, shifts );
    if ( mortar_fe_type::nDofPerFace > 0 )
        addFaceDof( elt, processor, next_free_dof, shifts );
}
}
#endif /* FEELPP_DofFromMortar_H */

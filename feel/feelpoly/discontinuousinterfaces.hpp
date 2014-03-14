/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-01-23

  Copyright (C) 2009 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2013 Feel++ Consortium

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
   \file discontinuousinterfaces.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-01-23
 */
#ifndef __DiscontinuousInterfaces_H
#define __DiscontinuousInterfaces_H 1

#include <feel/feelpoly/continuity.hpp>
#include <feel/feeldiscr/doffromelement.hpp>

namespace Feel
{
/**
 * \class DiscontinuousInterfaces
 * \brief describes discontinuous interfaces and provide dof creation
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename A0>
class DiscontinuousInterfaces : public Feel::detail::continuity_base
{
public:


    /** @name Constants
     */
    //@{

    static const bool is_continuous = false;
    static const bool is_discontinuous_locally = true;
    static const bool is_discontinuous_totally = false;

    static const uint16_type n_discontinuities = fusion::result_of::size<A0>::type::value;




    //@}

    /** @name Typedefs
     */
    //@{
    typedef A0 discontinuity_markers_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    DiscontinuousInterfaces()
        :
        M_d_faces()
    {}

    //! copy constructor
    DiscontinuousInterfaces( DiscontinuousInterfaces const & d )
        :
        M_d_faces( d.M_d_faces )
    {}

    //! destructor
    ~DiscontinuousInterfaces()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    DiscontinuousInterfaces & operator=( DiscontinuousInterfaces const &  d )
    {
        if ( this != &d )
        {
            M_d_faces = d.M_d_faces;
        }

        return *this;
    }

    //@}

    /** @name Accessors
     */
    //@{

    discontinuity_markers_type const& discontinuityMarkers() const
    {
        return M_d_faces;
    }


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{


    //@}



    template<typename MeshType, typename DofType>
    class apply
    {
    public:
        typedef size_type result_type;
        typedef MeshType mesh_type;
        typedef DofType dof_type;
        typedef typename dof_type::fe_type fe_type;

        apply( MeshType& M, DofType& D )
            :
            M_mesh( M ),
            M_dof( D ),
            M_fe( D.fe() )
        {}
        template<typename T>
#if BOOST_VERSION < 104200
        result_type operator()( const T& t, const size_type& start )
#else
        result_type operator()( const size_type& start, const T& t )
#endif
        {
            boost::tuple<size_type,size_type,size_type> disc = boost::make_tuple( mpl::at<T,mpl::int_<0> >::type::value,
                    mpl::at<T,mpl::int_<1> >::type::value,
                    mpl::at<T,mpl::int_<2> >::type::value );

            return build( disc, start );
        }
    private:

        // -loop on faces marked as discontinuous
        //    - verify that the face is connected to two elements
        //    - for each element associated to the face
        //        - extract the id of the element
        //        - loop on dof associated to the face in the current element
        //        - insert dof by "just" increment next_free_dof, be careful that
        //          map_gdof must be updated to ensure that these will _never_ be revisited (don't use insertDof)
        //    - return the next_free_dof
        size_type
        build( boost::tuple<size_type,size_type,size_type> const& marker, size_type start )
        {
            size_type next_free_dof = start;
            size_type n_dof = 0;

            typedef typename mesh_type::element_type element_type;
            //typedef typename mesh_type::marker_element_iterator element_iterator;
            typedef typename mesh_type::element_const_iterator element_const_iterator;

            element_const_iterator fit, fen;
            //boost::tie( fit, fen ) = M_mesh.elementsWithMarker( boost::get<1>( marker ), M_mesh.rank() );
            boost::tie( fit, fen ) = M_mesh.elementsRange();
            DVLOG(2) << "[DiscontinuousInterfaces::build] n_elements = " << std::distance( fit, fen )
                          << " with marker " << boost::get<1>( marker ) << "\n";
#if 0

            while ( fit != fen )
            {
                if ( fit->marker().value() != boost::get<1>( marker ) )
                {
                    ++fit;
                    continue;
                }

                DVLOG(2) << "found element with marker " << fit->marker().value() << "\n";
                typename element_type::face_const_iterator it, en;
                boost::tie( it, en ) = fit->faces();


                for ( ; it != en; ++it )
                {
                    DVLOG(2) << "face with marker " << ( *it )->marker().value() << "\n";
                    //if ( (*it)->marker().value() == boost::get<0>( marker ) )
                    {
                        DVLOG(2) << "------------------------------------------------------------\n";
                        DVLOG(2) << "face " << ( *it )->id()
                                      << " marker = " << boost::get<0>( marker )
                                      << " elt marker 0 = " << boost::get<1>( marker )
                                      << " elt marker 1 = " << boost::get<2>( marker ) << "\n";
                        DVLOG(2) << "element marker " << fit->marker() << "\n";

                        addVertexDof( *fit, *( *it ), next_free_dof, 0, mpl::bool_<fe_type::nDofPerVertex>() );
                        addEdgeDof( *fit, *( *it ), next_free_dof, 0, mpl::bool_<fe_type::nDofPerEdge>(), mpl::int_<mesh_type::nDim>() );


                    }
                }

                ++fit;
            }

#else
            boost::tie( fit, fen ) = M_mesh.elementsRange();

            DofFromElement<dof_type,fe_type> dfe( &M_dof, M_fe );
            while ( fit != fen )
            {
                if ( fit->marker().value() == ( int )boost::get<1>( marker ) )
                {

                    //M_dof.addDofFromElement( *fit, next_free_dof, 0 );
                    dfe.add( *fit, next_free_dof, M_dof.worldComm().localRank()  );
                }

                ++fit;
            }

#endif
            //n_dof = next_free_dof-start;
            n_dof = next_free_dof;

            DVLOG(2) << "[DiscontinuousInterfaces::build] n_dof = " << n_dof << "\n";

            //boost::tie( fit, fen ) = M_mesh.elementsWithMarker( boost::get<2>( marker ), M_mesh.rank() );
            boost::tie( fit, fen ) = M_mesh.elementsRange();
            DVLOG(2) << "[DiscontinuousInterfaces::build] n_elements = " << std::distance( fit, fen )
                          << " with marker " << boost::get<2>( marker ) << "\n";
#if 0

            while ( fit != fen )
            {
                if ( fit->marker().value() != boost::get<2>( marker ) )
                {
                    ++fit;
                    continue;
                }

                DVLOG(2) << "found element with marker " << fit->marker().value() << "\n";
                typename element_type::face_const_iterator it, en;
                boost::tie( it, en ) = fit->faces();


                for ( ; it != en; ++it )
                {
                    DVLOG(2) << "face with marker " << ( *it )->marker().value() << "\n";
                    //if ( (*it)->marker().value() == boost::get<0>( marker ) )
                    {
                        DVLOG(2) << "------------------------------------------------------------\n";
                        DVLOG(2) << "face " << ( *it )->id()
                                      << " marker = " << boost::get<0>( marker )
                                      << " elt marker 0 = " << boost::get<1>( marker )
                                      << " elt marker 1 = " << boost::get<2>( marker ) << "\n";
                        DVLOG(2) << "element marker " << fit->marker() << "\n";

                        addVertexDof( *fit, *( *it ), next_free_dof, n_dof, mpl::bool_<fe_type::nDofPerVertex>() );
                        addEdgeDof( *fit, *( *it ), next_free_dof, n_dof, mpl::bool_<fe_type::nDofPerEdge>(), mpl::int_<mesh_type::nDim>() );


                    }

                }

                ++fit;
            }

#else

            typedef typename DofType::dof_map_type dof_map_type;
            dof_map_type m1( M_dof.mapGDof() );
            M_dof.clearMapGDof();
            boost::tie( fit, fen ) = M_mesh.elementsRange();

            while ( fit != fen )
            {
                if ( fit->marker().value() == boost::get<2>( marker ) )
                {
                    //M_dof.addDofFromElement( *fit, next_free_dof );
                    dfe.add( *fit, next_free_dof, M_dof.worldComm().localRank()  );
                }

                ++fit;
            }

#if 0
            dof_map_type m2( M_dof.mapGDof() );
            M_dof.clearMapGDof();
            std::merge( m1.begin(), m1.end(), m2.begin(), m2.end(), M_dof.mapGDof().begin() );
#else
            dof_map_type m2( M_dof.mapGDof() );
            M_dof.clearMapGDof();
            typename dof_map_type::iterator it = m1.begin();
            typename dof_map_type::iterator en = m1.end();

            for ( ; it != en; ++it )
            {
                M_dof.mapGDof().insert( *it );
            }

            it = m2.begin();
            en = m2.end();

            for ( ; it != en; ++it )
            {
                M_dof.mapGDof().insert( *it );
            }

            DVLOG(2) << "size dictionnary = " << M_dof.mapGDof().size() << " next_free_dof = " << next_free_dof+n_dof << "\n";
#endif
#endif

            return next_free_dof;
        }
        template<typename element_type, typename face_type>
        void
        addVertexDof( element_type const& elt, face_type const& face, size_type& next_free_dof, size_type shift, mpl::bool_<false> )
        {

        }
        template<typename element_type, typename face_type>
        void
        addVertexDof( element_type const& elt, face_type const& face, size_type& next_free_dof, size_type shift, mpl::bool_<true> )
        {
            uint16_type iFaEl = invalid_uint16_type_value;

            if ( face.ad_first() == elt.id() )
                iFaEl = face.pos_first();

            else if ( face.ad_second() == elt.id() )
                iFaEl = face.pos_second();

            for ( uint16_type iVeFa = 0; iVeFa < face_type::numVertices; ++iVeFa )
            {
                // local vertex number (in element)
                uint16_type iVeEl = element_type::fToP( iFaEl, iVeFa );
                Feel::detail::ignore_unused_variable_warning( iVeEl );

                FEELPP_ASSERT( iVeEl != invalid_uint16_type_value ).error( "invalid local dof" );

                // Loop number of Dof per vertex
                for ( uint16_type l = 0; l < fe_type::nDofPerVertex; ++l )
                {
                    uint16_type lid = iVeEl * fe_type::nDofPerVertex + l;
                    const size_type gDof = ( elt.point( iVeEl ).id() ) * fe_type::nDofPerVertex + l;

                    DVLOG(2) << "add vertex discontinuous dof " << next_free_dof << " in element " << elt.id() << " lid = " << lid << "\n";
                    bool inserted = M_dof.insertDof( elt.id(), lid, iVeEl, boost::make_tuple( 0, 0, gDof ), 0, next_free_dof, 1, false, shift );

                    if ( shift )
                    {
                        FEELPP_ASSERT( inserted == false )( elt.id() )
                        ( lid )( gDof )( next_free_dof ).error( "should have inserted unique dof" );
                    }

                    DVLOG(2) << "vertex discontinuous dof inserted : " << inserted << "\n";

                    DVLOG(2) << "added vertex discontinuous dof " <<  elt.id() << ", "
                                  <<  lid << ", "
                                  << boost::get<0>( M_dof.localToGlobal( elt.id(), lid, 0 ) ) << "\n";
                }

            }


        }
        template<typename element_type, typename face_type>
        void
        addEdgeDof( element_type const& elt, face_type const& face, size_type& next_free_dof, size_type shift, mpl::bool_<false>, mpl::int_<1> )
        {}
        template<typename element_type, typename face_type>
        void
        addEdgeDof( element_type const& elt, face_type const& face, size_type& next_free_dof, size_type shift, mpl::bool_<true>, mpl::int_<1> )
        {}

        template<typename element_type, typename face_type>
        void
        addEdgeDof( element_type const& elt, face_type const& face, size_type& next_free_dof, size_type shift, mpl::bool_<false>, mpl::int_<2> )
        {}
        template<typename element_type, typename face_type>
        void
        addEdgeDof( element_type const& elt, face_type const& face, size_type& next_free_dof, size_type shift, mpl::bool_<true>, mpl::int_<2> )
        {
            uint16_type iFaEl = invalid_uint16_type_value;

            if ( face.ad_first() == elt.id() )
                iFaEl = face.pos_first();

            else if ( face.ad_second() == elt.id() )
                iFaEl = face.pos_second();

            // Loop number of Dof per edge
            for ( uint16_type l = 0; l < fe_type::nDofPerEdge; ++l )
            {
                uint16_type lid = element_type::numVertices*fe_type::nDofPerVertex + iFaEl * fe_type::nDofPerEdge + l;
                const size_type gDof = ( elt.edge( iFaEl ).id() ) * fe_type::nDofPerEdge + l;

                DVLOG(2) << "add edge discontinuous dof " << next_free_dof << " in element " << elt.id() << " lid = " << lid << "\n";
                bool inserted = M_dof.insertDof( elt.id(), lid, iFaEl, boost::make_tuple( 1, 0, gDof ), 0, next_free_dof, 1, false, shift );
                DVLOG(2) << "edge discontinuous dof inserted (1 or 0) : " << inserted << "\n";

                DVLOG(2) << "added edge discontinuous dof "
                              <<  elt.id() << ", "
                              <<  lid << ", "
                              << boost::get<0>( M_dof.localToGlobal( elt.id(), lid, 0 ) ) << "\n";

            }

        }
    private:
        MeshType& M_mesh;
        DofType& M_dof;
        fe_type const& M_fe;
    };

private:
    discontinuity_markers_type M_d_faces;
};
} // Feel
#endif /* __DiscontinuousInterfaces_H */

/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-03-31

  Copyright (C) 2007 Universit�� Joseph Fourier (Grenoble I)

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
   \file MeshMover.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-03-31
 */
#ifndef __MeshMover_H
#define __MeshMover_H 1


namespace Feel
{
/**
 * \class MeshMover
 *  \brief Move mesh according to a given map
 *
 *  @author Christophe Prud'homme
 */
template<typename MeshType>
class MeshMover
{
public:


    /** @name Typedefs
     */
    //@{

    typedef MeshType mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef typename mesh_type::value_type value_type;

    typedef typename mesh_type::gm_type gm_type;
    typedef boost::shared_ptr<gm_type> gm_ptrtype;

    typedef typename mesh_type::element_iterator element_iterator;
    typedef typename mesh_type::element_const_iterator element_const_iterator;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    MeshMover()
        :
        M_updateMeshMeasures( true )
    {}

    MeshMover( mesh_ptrtype const& mesh )
        :
        M_mesh( mesh ),
        M_updateMeshMeasures( true )
    {}

    MeshMover( MeshMover const & )
    {}

    virtual ~MeshMover()
    {}

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
    void setUpdateMeshMeasures( bool b ) { M_updateMeshMeasures = b; }


    //@}

    /** @name  Methods
     */
    //@{

    template<typename DisplType>
    value_type tryit( DisplType const& u );

    template<typename DisplType>
    //boost::tuple<mesh_ptrtype,value_type> apply( mesh_ptrtype const& imesh, DisplType const& u );
    void apply( mesh_ptrtype& imesh, DisplType const& u );

    //@}
private:
    mesh_ptrtype M_mesh;
    bool M_updateMeshMeasures;
};



template<typename MeshType>
template<typename DisplType>
typename MeshMover<MeshType>::value_type
MeshMover<MeshType>::tryit( DisplType const& u )
{
    // check that the displacement function is vectorial
    //BOOST_MPL_ASSERT_MSG( DisplType::is_vectorial, INVALID_DISPLACEMENT_TYPE, DisplType );


}
template<typename MeshType>
template<typename DisplType>
//boost::tuple<typename MeshMover<MeshType>::mesh_ptrtype,
//             typename MeshMover<MeshType>::value_type>
void
MeshMover<MeshType>::apply( mesh_ptrtype& imesh, DisplType const& u )
{
    // check that the displacement function is vectorial
    //BOOST_MPL_ASSERT_MSG( DisplType::is_vectorial, INVALID_DISPLACEMENT_TYPE, DisplType );

    DVLOG(2) << "[Dof::generateDofPoints] generating dof coordinates\n";
    typedef typename mesh_type::element_type element_type;
    typedef typename mesh_type::face_type face_type;
    typedef typename gm_type::template Context<vm::POINT, element_type> gm_context_type;
    typedef boost::shared_ptr<gm_context_type> gm_context_ptrtype;

    typedef typename DisplType::functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::POINT, fe_type, gm_type, element_type> fecontext_type;


    gm_ptrtype gm( new gm_type );
    fe_type fe;

    //
    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    //
    typename gm_type::precompute_ptrtype __geopc( new typename gm_type::precompute_type( gm, gm->points() ) );

    //const uint16_type ndofv = fe_type::nDof;
    //mesh_ptrtype omesh( new mesh_type );
    //*omesh = *imesh;

    bool addExtendedMPIElt =  (imesh->worldComm().localSize() > 1) && u.functionSpace()->dof()->buildDofTableMPIExtended();
    EntityProcessType entityProcess = (addExtendedMPIElt)? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
    auto rangeElt = elements( imesh, entityProcess );
    auto it_elt = rangeElt.template get<1>();
    auto en_elt = rangeElt.template get<2>();
    if ( std::distance(it_elt,en_elt)==0 )
    {
        // call updateMeasures in parallel here because this function is call ( at the end of this function)
        // by others proc which have elements and need collective comm
        if ( imesh->worldComm().localSize() > 1 ) {
            //imesh->updateForUse();
            imesh->updateMeasures();
        }
        return;
    }

    typedef typename DisplType::pc_type pc_type;
    typedef boost::shared_ptr<pc_type> pc_ptrtype;
    pc_ptrtype __pc( new pc_type( u.functionSpace()->fe(), gm->points() ) );
    gm_context_ptrtype __c( new gm_context_type( gm, *it_elt, __geopc ) );

    typedef typename mesh_type::element_type geoelement_type;
    typedef typename fe_type::template Context<0, fe_type, gm_type, geoelement_type,gm_context_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( u.functionSpace()->fe(),
                                         __c,
                                         __pc ) );
    typedef typename fectx_type::id_type m_type;
    typedef boost::multi_array<m_type,1> array_type;
    array_type uvalues( u.idExtents( *__ctx ) );
    std::fill( uvalues.data(), uvalues.data()+uvalues.num_elements(), m_type::Zero() );
    u.id( *__ctx, uvalues );

    //const uint16_type ndofv = fe_type::nDof;


    std::map<int,bool> points_done;

    //if ( !u.areGlobalValuesUpdated() )
    u.updateGlobalValues();


    uint16_type nptsperelem = gm->points().size2();
    ublas::vector<value_type> val( fe_type::nComponents );

    for ( ; it_elt != en_elt; ++it_elt )
    {
        element_type const& curElt = *it_elt;

        __c->update( *it_elt );
        __ctx->update( __c );
        std::fill( uvalues.data(), uvalues.data()+uvalues.num_elements(), m_type::Zero() );
        u.id( *__ctx, uvalues );

        for ( uint16_type l =0; l < nptsperelem; ++l )
        {
            for ( uint16_type comp = 0; comp < fe_type::nComponents; ++comp )
            {
                val[ comp ] = uvalues[l]( comp,0 );
            }

            if ( points_done.find( curElt.point( l ).id() ) == points_done.end() )
            {
                //std::cout << "Pt: " << thedof << "Elem " << curElt.id() << " G=" << curElt.G() << "\n";
                imesh->elements().modify( imesh->elementIterator( curElt ), // it_elt,
                                          lambda::bind( &element_type::applyDisplacement,
                                                        lambda::_1,
                                                        l,
                                                        val ) );
                points_done[curElt.point( l ).id()] = true;
                //std::cout << "Pt: " << thedof << " Moved Elem " << curElt.id() << " G=" << curElt.G() << "\n";
            }

            else
            {
                imesh->elements().modify( imesh->elementIterator( curElt ), //it_elt,
                                          lambda::bind( &element_type::applyDisplacementG,
                                                        lambda::_1,
                                                        l,
                                                        val ) );
            }
        }

        // update internal data point of faces attached on this elt
        for ( size_type j = 0; j < imesh->numLocalFaces(); j++ )
        {
            if ( !curElt.facePtr( j ) ) continue;
            face_type const& curFace = curElt.face( j );

            for ( int f = 0; f < face_type::numPoints; ++f )
            {
                uint16_type ptLocalId = ( MeshType::nDim==1 )?j:curElt.fToP( j, f );
                auto const& curPoint = curElt.point( ptLocalId );
                for ( uint16_type comp = 0; comp < fe_type::nComponents; ++comp )
                {
                    val[ comp ] = curPoint( comp );
                }
                imesh->faces().modify( imesh->faceIterator( curFace ),
                                       lambda::bind( &face_type::setPointCoordG,
                                                     lambda::_1,
                                                     f,
                                                     val ) );
            }
        }

        // Todo : edges
    }

    //imesh->updateForUse();

    if ( M_updateMeshMeasures )
        imesh->updateMeasures();

    // reset geomap cache
    imesh->gm()->initCache( imesh.get() );
    imesh->gm1()->initCache( imesh.get() );

    // reset localisation tool
    imesh->tool_localization()->reset();

#if !defined( __INTEL_COMPILER )
    // notify observers that the mesh has changed
    LOG(INFO) << "Notify observers that the mesh has changed\n"; 
    imesh->meshChanged(MESH_CHANGES_POINTS_COORDINATES );
#endif
}

template<typename MeshType, typename DisplType>
boost::shared_ptr<MeshType>
meshMove( boost::shared_ptr<MeshType>& m, DisplType const& u )
{
    MeshMover<MeshType> mover( m );
    mover.apply( m, u );
    return m;
}
} // Feel
#endif /* __MeshMover_H */

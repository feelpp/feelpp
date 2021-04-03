/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    typedef typename mesh_type::value_type value_type;

    typedef typename mesh_type::gm_type gm_type;
    typedef std::shared_ptr<gm_type> gm_ptrtype;

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

private:
    void applyGhostElements( mesh_ptrtype& imesh,
                             std::map<rank_type,std::vector< boost::tuple<size_type,std::vector< ublas::vector<value_type> > > > > const& dataToSend,
                             std::unordered_set<size_type> & points_done );

    void updateForUse( mesh_ptrtype& imesh );
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
    typedef typename gm_type::template Context<element_type> gm_context_type;
    typedef std::shared_ptr<gm_context_type> gm_context_ptrtype;

    typedef typename DisplType::functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::POINT, fe_type, gm_type, element_type> fecontext_type;

    std::map<rank_type,std::vector< boost::tuple<size_type,std::vector< ublas::vector<value_type> > > > > dataToSend;
    std::unordered_set<size_type> points_done;

    auto rangeElt = elements( imesh );
    auto it_elt = rangeElt.template get<1>();
    auto en_elt = rangeElt.template get<2>();
    if ( std::distance(it_elt,en_elt)==0 )
    {
        if ( imesh->worldComm().localSize() > 1 )
            this->applyGhostElements( imesh, dataToSend, points_done );

        this->updateForUse( imesh );
        return;
    }

    //gm_ptrtype gm( new gm_type );
    gm_ptrtype gm = imesh->gm();

    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    typename gm_type::precompute_ptrtype __geopc( new typename gm_type::precompute_type( gm, gm->points() ) );

    typedef typename DisplType::pc_type pc_type;
    typedef std::shared_ptr<pc_type> pc_ptrtype;
    pc_ptrtype __pc( new pc_type( u.functionSpace()->fe(), gm->points() ) );
    gm_context_ptrtype __c = gm->template context<vm::POINT>( unwrap_ref( *it_elt ), __geopc );

    typedef typename mesh_type::element_type geoelement_type;
    typedef typename fe_type::template Context<0, fe_type, gm_type, geoelement_type,vm::POINT> fectx_type;
    typedef std::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( u.functionSpace()->fe(),
                                         __c,
                                         __pc ) );
    typedef typename fectx_type::id_type m_type;
    typedef boost::multi_array<m_type,1> array_type;
    array_type uvalues( u.idExtents( *__ctx ) );
    m_type m;
    std::fill( uvalues.data(), uvalues.data()+uvalues.num_elements(), m.constant(0.) );
    u.id( *__ctx, uvalues );


    uint16_type nptsperelem = gm->points().size2();
    ublas::vector<value_type> val( fe_type::nComponents );
    std::vector< ublas::vector<value_type> > dataEltToSend( nptsperelem, ublas::vector<value_type>( fe_type::nComponents ) );
    for ( ; it_elt != en_elt; ++it_elt )
    {
        element_type const& curElt = unwrap_ref( *it_elt );
        auto & eltModified = imesh->elementIterator( curElt )->second;
        __c->template update<vm::POINT>( *it_elt );
        __ctx->update( __c );
        std::fill( uvalues.data(), uvalues.data()+uvalues.num_elements(), m.constant(0.) );
        u.id( *__ctx, uvalues );

        bool isGhostInOtherPart = !curElt.idInOthersPartitions().empty();
        for ( uint16_type l =0; l < nptsperelem; ++l )
        {
            for ( uint16_type comp = 0; comp < fe_type::nComponents; ++comp )
            {
                val[ comp ] = uvalues[l]( comp,0 );
            }

            size_type ptId = curElt.point( l ).id();
            if ( points_done.find( ptId ) == points_done.end() )
            {
                //std::cout << "Pt: " << thedof << "Elem " << curElt.id() << " G=" << curElt.G() << "\n";
                eltModified.applyDisplacement( l, val );
                points_done.insert( ptId );
                //std::cout << "Pt: " << thedof << " Moved Elem " << curElt.id() << " G=" << curElt.G() << "\n";
            }
            if ( isGhostInOtherPart )
                dataEltToSend[l] = val;
        }
        if ( isGhostInOtherPart )
            for ( auto const& idOtherPart : curElt.idInOthersPartitions() )
                dataToSend[idOtherPart.first].push_back( boost::make_tuple( idOtherPart.second, dataEltToSend ) );
    }

    if ( imesh->worldComm().localSize() > 1 )
        this->applyGhostElements( imesh, dataToSend, points_done );

    this->updateForUse( imesh );
}

template<typename MeshType>
void
MeshMover<MeshType>::applyGhostElements( mesh_ptrtype& imesh,
                                         std::map<rank_type,std::vector< boost::tuple<size_type,std::vector< ublas::vector<value_type> > > > > const& dataToSend,
                                         std::unordered_set<size_type> & points_done )
{
    // mpi comm
    int neighborSubdomains = imesh->neighborSubdomains().size();
    int nbRequest = 2*neighborSubdomains;
    mpi::request * reqs = new mpi::request[nbRequest];
    int cptRequest=0;
    std::map<rank_type,std::vector< boost::tuple<size_type,std::vector< ublas::vector<value_type> > > > > dataToRecv;
    for ( rank_type neighborRank : imesh->neighborSubdomains() )
    {
        CHECK( dataToSend.find( neighborRank ) != dataToSend.end() ) << "something wrong in parallel datastructure of mesh";
        reqs[cptRequest++] = imesh->worldComm().localComm().isend( neighborRank , 0, dataToSend.find( neighborRank )->second );
        reqs[cptRequest++] = imesh->worldComm().localComm().irecv( neighborRank , 0, dataToRecv[neighborRank] );
    }
    mpi::wait_all(reqs, reqs + nbRequest);
    delete [] reqs;

    // move points in ghost elements
    for ( auto const& dataRecvByProc : dataToRecv )
    {
        for ( auto const& dataRecvByElt : dataRecvByProc.second )
        {
            size_type eltId = boost::get<0>( dataRecvByElt );
            auto const& pointsData =  boost::get<1>( dataRecvByElt );
            auto & eltModified = imesh->elementIterator( eltId )->second;
            for ( uint16_type p=0;p<mesh_type::element_type::numPoints;++p )
            {
                size_type ptId = eltModified.point( p ).id();
                if ( points_done.find( ptId ) == points_done.end() )
                {
                    eltModified.applyDisplacement( p, pointsData[p] );
                    points_done.insert( ptId );
                }
            }
        }
    }
}

template<typename MeshType>
void
MeshMover<MeshType>::updateForUse( mesh_ptrtype& imesh )
{
    for ( auto m : imesh->meshesWithNodesShared() )
        m->updateForUseAfterMovingNodes( M_updateMeshMeasures );
}


template<typename MeshType, typename DisplType>
std::shared_ptr<MeshType>
meshMove( std::shared_ptr<MeshType>& m, DisplType const& u )
{
    MeshMover<MeshType> mover( m );
    mover.apply( m, u );
    return m;
}
} // Feel
#endif /* __MeshMover_H */

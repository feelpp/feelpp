/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-01-31

  Copyright (C) 2008-2010 Université Joseph Fourier (Grenoble I)

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
   \file lagrangep1.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-01-31
 */
#ifndef __OperatorLagrangeP1_H
#define __OperatorLagrangeP1_H 1

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/fekete.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelfilters/pointsettomesh.hpp>
#include <feel/feelfilters/exporterquick.hpp>

namespace Feel
{
namespace detail
{
/**
 * \internal
 * class that allows to compute the basis function type
 */
template<typename SpaceType>
struct SpaceToLagrangeP1Space
{
    typedef SpaceType domain_space_type;
    typedef typename domain_space_type::fe_type::convex_type convex_type;

    template< template<uint16_type Dim> class PsetType>
    struct SelectConvex
    {
        //typedef typename PsetType<convex_type::nDim>::value_type value_type;
        typedef typename SpaceType::value_type value_type;

        typedef Lagrange<1,PsetType> type;
    };
    typedef typename convex_type::template shape<convex_type::nDim, 1, convex_type::nDim>::type image_convex_type;
    typedef Mesh<image_convex_type> image_mesh_type;

    typedef typename mpl::if_<mpl::bool_<domain_space_type::is_scalar>,
                              mpl::identity<typename SelectConvex<Scalar>::type>,
                              typename mpl::if_<mpl::bool_<domain_space_type::is_vectorial>,
                                                mpl::identity<typename SelectConvex<Vectorial>::type>,
                                                mpl::identity<typename SelectConvex<Tensor2>::type> >::type>::type::type basis_type;


    typedef FunctionSpace<image_mesh_type, fusion::vector<basis_type> > image_space_type;
};

}
/**
 * \class OperatorLagrangeP1
 * \ingroup SpaceTime
 * @author build a P1 Lagrange dof-equivalent space from a \c FunctionSpace
 *
 * \see FunctionSpace
 */
template<typename SpaceType>
class OperatorLagrangeP1
    :
    public OperatorInterpolation<SpaceType, typename detail::SpaceToLagrangeP1Space<SpaceType>::image_space_type>
{
    typedef OperatorInterpolation<SpaceType, typename detail::SpaceToLagrangeP1Space<SpaceType>::image_space_type> super;
public:


    /** @name Typedefs
     */
    //@{
    /**
     * domain space
     */
    typedef typename super::domain_space_type domain_space_type;
    typedef typename super::domain_space_ptrtype domain_space_ptrtype;
    typedef typename domain_space_type::mesh_type domain_mesh_type;
    typedef typename domain_space_type::mesh_ptrtype domain_mesh_ptrtype;
    typedef typename domain_space_type::value_type value_type;
    typedef typename domain_space_type::fe_type domain_fe_type;


    typedef typename super::backend_ptrtype backend_ptrtype;
    /**
     * image space
     */
    typedef typename super::dual_image_space_type dual_image_space_type;
    typedef typename super::dual_image_space_ptrtype dual_image_space_ptrtype;
    typedef typename dual_image_space_type::mesh_type image_mesh_type;
    typedef typename dual_image_space_type::mesh_ptrtype image_mesh_ptrtype;
    typedef typename dual_image_space_type::fe_type image_fe_type;

    /**
     * type which instances will hold the correspondance between the
     * mesh elements associated to the domain space and the ones
     * associated to the image space
     */
    typedef std::vector<std::list<size_type> > el2el_type;

    typedef typename detail::SpaceToLagrangeP1Space<SpaceType>::convex_type domain_convex_type;
    typedef typename detail::SpaceToLagrangeP1Space<SpaceType>::image_convex_type image_convex_type;
    //typedef PointSetWarpBlend<domain_convex_type, domain_space_type::basis_type::nOrder, value_type> pset_type;
    typedef PointSetFekete<domain_convex_type, domain_space_type::basis_type::nOrder, value_type> pset_type;
    typedef PointSetToMesh<domain_convex_type, value_type> p2m_type;

    typedef typename domain_mesh_type::gm_type gm_type;
    typedef typename domain_mesh_type::element_type element_type;
    typedef typename gm_type::template Context<vm::POINT, element_type> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
    typedef typename gm_type::precompute_ptrtype gmpc_ptrtype;

    //typedef typename mpl::at_c<functionspace_vector_type,SpaceIndex>::type functionspace_ptrtype;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * Construct a P1 Lagrange space (ie spanned by a P1 Lagrange
     * basis) out of the space \p space
     */
    OperatorLagrangeP1( domain_space_ptrtype const& space,
                        backend_ptrtype const& backend );

    /**
     * destructor. nothing really to be done here
     */
    ~OperatorLagrangeP1()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    /**
     * copy operator
     */
    OperatorLagrangeP1 const&  operator=( OperatorLagrangeP1 const& lp1 )
    {
        if ( this != &lp1 )
            {
                _M_el2el = lp1._M_el2el;
                _M_el2pt = lp1._M_el2pt;
                _M_pset = lp1._M_pset;
                _M_gmpc = lp1._M_gmpc;
                _M_p2m = lp1._M_p2m;
            }
        return *this;
    }


    //@}

    /** @name Accessors
     */
    //@{

    uint16_type
    localDof( typename domain_mesh_type::element_const_iterator it,
              uint16_type localptid ) const
    {
        //return localDof( it, localptid, mpl::bool_<is_modal>() );
        return localDof( it, localptid, mpl::bool_<false>() );
    }

    uint16_type
    localDof( typename domain_mesh_type::element_const_iterator /*it*/,
              uint16_type localptid,
              mpl::bool_<false> ) const
    {
        return localptid;
    }

    uint16_type
    localDof( typename domain_mesh_type::element_const_iterator it,
              uint16_type localptid,
              mpl::bool_<true> ) const;


    image_mesh_ptrtype
    mesh()
    {
        return _M_mesh;
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * check that the P1 mesh is properly constructed
     */
    void check() const;


    //@}

private:
    OperatorLagrangeP1();
    OperatorLagrangeP1( OperatorLagrangeP1 const & );

private:
    image_mesh_ptrtype _M_mesh;
    el2el_type _M_el2el;
    std::vector<std::vector<size_type> > _M_el2pt;

    pset_type _M_pset;
    gmpc_ptrtype _M_gmpc;

    mutable p2m_type _M_p2m;

};


//
// OperatorLagrangeP1
//
template<typename space_type>
OperatorLagrangeP1<space_type>::OperatorLagrangeP1( domain_space_ptrtype const& space,
                                                    backend_ptrtype const& backend )
    :
    super( space,
           dual_image_space_ptrtype( new dual_image_space_type( image_mesh_ptrtype( new image_mesh_type ) ) ),
           backend ),
    _M_mesh( new image_mesh_type ),
    _M_el2el(),
    _M_el2pt(),
    _M_pset( 0 ),
    _M_gmpc( new typename gm_type::precompute_type( this->domainSpace()->mesh()->gm(),
                                                    _M_pset.points() ) ),
    _M_p2m()
{
    // create a triangulation of the current convex
    // using equispaced points defined on the reference
    // element
    //_M_p2m.addBoundaryPoints( _M_pset.points() );
    _M_p2m.visit( &_M_pset );

    ExporterQuick<image_mesh_type> exp( "vtk", "ensight" );
    exp.save( _M_p2m.mesh() );

    // do not renumber the mesh entities
    _M_p2m.mesh()->components().clear ( MESH_RENUMBER );
    _M_p2m.mesh()->components().clear ( MESH_CHECK );
    _M_p2m.mesh()->updateForUse();

    Debug(5035) << "[P1 Lagrange] Pointset " << _M_pset.points() << "\n";

    typedef typename image_mesh_type::point_type point_type;
    typedef typename image_mesh_type::element_type element_type;

    //_M_el2el.resize( this->domainSpace()->mesh()->numElements() * _M_p2m.mesh()->numElements() );
    //_M_el2pt.resize( this->domainSpace()->mesh()->numElements() * _M_p2m.mesh()->numElements() );

    std::vector<bool> pts_done(  this->domainSpace()->nLocalDof()/domain_space_type::nComponents, false );
    size_type ne = this->domainSpace()->mesh()->numElements() * _M_p2m.mesh()->numElements();
    std::vector<boost::tuple<size_type, uint16_type, size_type> > dofindices( image_mesh_type::element_type::numVertices*ne, boost::make_tuple(0,0,0) );
    std::vector<bool> eltid_done(  this->domainSpace()->nLocalDof(), false );

    // construct the mesh
    typename domain_mesh_type::element_const_iterator it = this->domainSpace()->mesh()->beginElement();
    typename domain_mesh_type::element_const_iterator en = this->domainSpace()->mesh()->endElement();

    // use the geometric transformation to transform
    // the local equispace mesh to a submesh of the
    // current element
    gmc_ptrtype gmc( new gmc_type( this->domainSpace()->mesh()->gm(), *it, _M_gmpc ) );


    for ( size_type elid = 0, pt_image_id = 0; it != en; ++it )
        {
            Debug(5035) << "=========================================\n";
            Debug(5035) << "global element " << it->id() << " oriented ok ? : " << it->isAnticlockwiseOriented() << "\n";
            Debug(5035) << "global element G=" << it->G() << "\n";

            gmc->update( *it );

            // accumulate the local mesh in element *it in the new mesh
            typename image_mesh_type::element_const_iterator itl = _M_p2m.mesh()->beginElement();
            typename image_mesh_type::element_const_iterator enl = _M_p2m.mesh()->endElement();
            for( ; itl != enl; ++itl, ++elid )
                {
                    Debug(5035) << "************************************\n";
                    Debug(5035) << "local elt = " << elid << "\n";
                    Debug(5035) << "local element " << itl->id() << " oriented ok ? : " << itl->isAnticlockwiseOriented() << "\n";
                    Debug(5035) << "local element G=" << itl->G() << "\n";

                    //_M_el2pt[elid].resize( domain_mesh_type::element_type::numVertices );
                    element_type elt;

                    elt.setId( elid );
                    elt.setMarker( it->marker().value() );
                    // accumulate the points
                    for( int p = 0; p < image_mesh_type::element_type::numVertices; ++p )
                        {
                            Debug(5035) << "local In original element, vertex number " << itl->point( p ).id() << "\n";
                            uint16_type localptid = itl->point( p ).id();;//p;//it->point( p ).id(); //itl->point( p ).id();
                            uint16_type localptid_dof = localDof( it, localptid );
                            size_type ptid = boost::get<0>(this->domainSpace()->dof()->localToGlobal( it->id(),
                                                                                                      localptid_dof, 0 ));
                            FEELPP_ASSERT( ptid < this->domainSpace()->nLocalDof()/domain_space_type::nComponents )( ptid )( this->domainSpace()->nLocalDof()/domain_space_type::nComponents ).warn( "invalid domain dof index" );
                            //if ( pts_done[ptid] == false )
                            {
                                dofindices[pt_image_id] = boost::make_tuple( elid, p, ptid );
                                pt_image_id++;
                                //pts_done[ptid] = pt_image_id;
                            }

                            //point_type __pt( ptid, gmc->xReal(localptid) );
                            point_type __pt( ptid,
                                             boost::get<0>(this->domainSpace()->dof()->dofPoint( ptid ))  );

                            Debug(5035) << "[OperatorLagrangeP1] element id "
                                        << elid << "\n";

                            Debug(5035) << "[OperatorLagrangeP1] local point id "
                                        << localptid << " coord " << itl->point( p ).node() << "\n";
                            Debug(5035) << "[OperatorLagrangeP1] point local id "
                                        << localptid << " global id " << ptid << "\n"
                                        << " localptid_dof = " << localptid_dof << "\n";
                            Debug(5035) << "[OperatorLagrangeP1] adding point "
                                        << __pt.id() << " : " << __pt.node() << "\n";
                            Debug(5035) << "[OperatorLagrangeP1] domain point gid "
                                      << boost::get<1>(this->domainSpace()->dof()->dofPoint( ptid ))
                                //<< " lid: " << boost::get<2>(this->domainSpace()->dof()->dofPoint( ptid ))
                                      << " coords: " << boost::get<0>(this->domainSpace()->dof()->dofPoint( ptid )) << "\n";
                            FEELPP_ASSERT( ublas::norm_2( boost::get<0>(this->domainSpace()->dof()->dofPoint( ptid ))-__pt.node()) < 1e-10 )
                                (boost::get<0>(this->domainSpace()->dof()->dofPoint( ptid )))
                                (__pt.node())
                                (elid)
                                (itl->id())
                                (localptid)
                                (localptid_dof)
                                (ptid).warn( "inconsistent point coordinates" );


                            _M_mesh->addPoint( __pt );

                            elt.setPoint( p, _M_mesh->point( ptid ) );
                            //elt.setPoint( localptid, _M_mesh->point( ptid ) );
                            //_M_el2pt[elid][p] = localptid;
                            //FEELPP_ASSERT( ublas::norm_2( _M_mesh->point( ptid ).node(), __pt.node() ) ).error("invalid point");
#if 1
                            FEELPP_ASSERT( ublas::norm_2( boost::get<0>(this->domainSpace()->dof()->dofPoint( ptid ))-elt.point(p).node()) < 1e-10 )
                                (boost::get<0>(this->domainSpace()->dof()->dofPoint( ptid )))
                                (p)
                                (elt.point(p).node())
                                (elid)
                                (itl->id())
                                (localptid)
                                (localptid_dof)
                                (ptid).warn( "[after] inconsistent point coordinates" );

                            FEELPP_ASSERT( ublas::norm_2( elt.point(p).node()-ublas::column( elt.G(), p ) ) < 1e-10 )
                                (p)
                                (elt.point(p).node())
                                (elid)
                                (elt.G()).warn( "[after] inconsistent point coordinates" );
#endif
                        }
#if 0
                    static const int nDim = element_type::nDim;
                    ublas::vector<double> D( nDim + 1, 0 );
                    typename matrix_node<value_type>::type M( nDim+1,nDim+1);
                    typename matrix_node<value_type>::type P( nDim,nDim+1);
                    typename matrix_node<value_type>::type P0( nDim,nDim+1);
                    std::fill( M.data().begin(), M.data().end(), value_type( 1.0 ) );

                    for( int p = 0; p < element_type::numVertices; ++p )
                        {
                            ublas::column( P0, p ) = elt.vertex( p );
                        }
                    ublas::subrange( M, 0, nDim, 0, nDim+1 ) = P0;
                    double meas_times = details::det( M, mpl::int_<nDim+1>() );
                    FEELPP_ASSERT( meas_times > 0 )( elt.id() )( elt.G() ).warn( "negative area" );
#else
                    //FEELPP_ASSERT( elt.isAnticlockwiseOriented() )( elt.id() )( elt.G() ).warn( "negative area" );
                    //FEELPP_ASSERT( ublas::norm_inf( elt.G()- it->G() ) < 1e-10 )( it->G() )( elt.G() ).warn( "global: not same element" );
                    //FEELPP_ASSERT( ublas::norm_inf( elt.G()- itl->G() ) < 1e-10 )( itl->G() )( elt.G() ).warn( "local: not same element" );
#endif
                    //_M_el2el[it->id()].push_back( elt.id() );
                    _M_mesh->addElement ( elt );
                }
        }

    Debug(5035) << "[P1 Lagrange] Number of points in mesh: " << _M_mesh->numPoints() << "\n";

    _M_mesh->setNumVertices( _M_mesh->numPoints() );
    //_M_mesh->setRenumber( false );
    //_M_mesh->updateForUse( MESH_UPDATE_EDGES );
    _M_mesh->components().clear ( MESH_RENUMBER );

    for( size_type i = 0; i < this->domainSpace()->nLocalDof()/domain_space_type::nComponents; ++i )
        {
            //dofindices[i] = _M_mesh->point( i ).id();//boost::get<1>(this->domainSpace()->dof()->dofPoint( i ));
            FEELPP_ASSERT(ublas::norm_2( boost::get<0>(this->domainSpace()->dof()->dofPoint( i ))-
                                       _M_mesh->point( i ).node()
                                       ) < 1e-10 )
                (boost::get<0>(this->domainSpace()->dof()->dofPoint( i )))
                (_M_mesh->point( i ).node())
                (i).warn( "[mesh] check inconsistent point coordinates" );
            //Debug(5035) << "dofindices[" << i << "]=" << boost::get<0>( dofindices[i] ) << "," << boost::get<1>( dofindices[i] ) << "\n";
        }

    // construct the p1 space and set the operator
    this->init( this->domainSpace(),
                dual_image_space_ptrtype( new dual_image_space_type( _M_mesh, dofindices  ) ),
                backend );

    for( size_type i = 0; i < this->domainSpace()->nLocalDof()/domain_space_type::nComponents; ++i )
        {
            FEELPP_ASSERT(ublas::norm_2( boost::get<0>(this->domainSpace()->dof()->dofPoint( i ))-
                                       _M_mesh->point( i ).node()
                                       ) < 1e-10 )
                (boost::get<0>(this->domainSpace()->dof()->dofPoint( i )))
                (_M_mesh->point( i ).node())
                (i).warn( "[mesh2] check inconsistent point coordinates" );
        }
    //
    // Update generates the matrix associated with the interpolation operator
    // comment out for now
#if 0
    this->update();
#endif // 0

    for( size_type i = 0; i < this->domainSpace()->nLocalDof()/domain_space_type::nComponents; ++i )
        {
            FEELPP_ASSERT(ublas::norm_2( boost::get<0>(this->domainSpace()->dof()->dofPoint( i ))-
                                       _M_mesh->point( i ).node()
                                       ) < 1e-10 )
                (boost::get<0>(this->domainSpace()->dof()->dofPoint( i )))
                (_M_mesh->point( i ).node())
                (i).warn( "[mesh3] check inconsistent point coordinates" );
        }
#if 0
    for( size_type i = 0; i < this->domainSpace()->nLocalDof()/domain_space_type::nComponents; ++i )
        {
            FEELPP_ASSERT( boost::get<1>(this->domainSpace()->dof()->dofPoint( i )) ==
                         boost::get<1>(this->dualImageSpace()->dof()->dofPoint( i )) )
                (boost::get<0>(this->domainSpace()->dof()->dofPoint( i )))
                (boost::get<0>(this->dualImageSpace()->dof()->dofPoint( i )))
                (i).warn( "check inconsistent point id" );
            FEELPP_ASSERT( ublas::norm_2( boost::get<0>(this->domainSpace()->dof()->dofPoint( i ))-
                                        boost::get<0>(this->dualImageSpace()->dof()->dofPoint( i ))
                                        ) < 1e-10 )
                (boost::get<0>(this->domainSpace()->dof()->dofPoint( i )))
                (boost::get<0>(this->dualImageSpace()->dof()->dofPoint( i )))
                (i).warn( "check inconsistent point coordinates" );
        }
#endif
    this->check();

}
template<typename space_type>
void
OperatorLagrangeP1<space_type>::check() const
{


}
#if 0

template<typename space_type>
typename OperatorLagrangeP1<space_type>::image_space_type::element_type
OperatorLagrangeP1<space_type>::operator()( element_type const& u ) const
{
    typename image_space_type::element_type res( this->dualImageSpace(), u.name() );
    //res.resize( _M_mesh->numPoints() );
    res = ublas::scalar_vector<value_type>( res.size(), value_type(0) );
    // construct the mesh
    typename domain_mesh_type::element_const_iterator it = this->domainSpace()->mesh()->beginElement();
    typename domain_mesh_type::element_const_iterator en = this->domainSpace()->mesh()->endElement();
    for ( ; it != en; ++it )
        {
            gmc_ptrtype gmc( new gmc_type( this->domainSpace()->mesh()->gm(), *it, _M_gmpc ) );
            typename matrix_node<value_type>::type u_at_pts = u.id( *gmc );

            typename std::list<size_type>::const_iterator ite = _M_el2el[it->id()].begin();
            typename std::list<size_type>::const_iterator ene = _M_el2el[it->id()].end();
            for( ; ite != ene; ++ite )
                {
                    typename domain_mesh_type::element_type const& elt = _M_mesh->element( *ite );
                    for( int c = 0; c < nComponents; ++c )
                        for( int p = 0; p < domain_mesh_type::element_type::numVertices; ++p )
                            {
                                size_type ptid = boost::get<0>(this->dualImageSpace()->dof()->localToGlobal( elt.id(), p, c ));
                                res( ptid ) = ublas::column( u_at_pts, _M_el2pt[*ite][p] )( nComponents*c+c );
                            }
                }
        }

    return res;
}
#endif

template<typename space_type>
uint16_type
OperatorLagrangeP1<space_type>::localDof( typename domain_mesh_type::element_const_iterator it,
                                                   uint16_type localptid,
                                                   mpl::bool_<true> ) const
{
    typedef typename domain_mesh_type::element_type::edge_permutation_type edge_permutation_type;
    uint16_type localptid_dof = localptid;
    if ( domain_mesh_type::nDim == 2 &&
         _M_pset.pointToEntity( localptid ).first == 1 &&
         it->edgePermutation( _M_pset.pointToEntity( localptid ).second ) == edge_permutation_type::REVERSE_PERMUTATION )
        {
            //size_type edge_id = it->edge( _M_pset.pointToEntity( localptid ).second ).id();
            uint16_type edge_id = _M_pset.pointToEntity( localptid ).second;

            uint16_type l = localptid - ( domain_mesh_type::element_type::numVertices * domain_fe_type::nDofPerVertex +
                                          edge_id * domain_fe_type::nDofPerEdge );

            FEELPP_ASSERT( l != invalid_uint16_type_value && ( l < domain_fe_type::nDofPerEdge ) )
                ( l )
                ( domain_fe_type::nDofPerEdge )
                (localptid)
                (_M_pset.pointToEntity( localptid ).second)
                (domain_mesh_type::element_type::numVertices * domain_fe_type::nDofPerVertex +
                 edge_id * domain_fe_type::nDofPerEdge ).error( "invalid value" );

            localptid_dof = ( domain_mesh_type::element_type::numVertices * domain_fe_type::nDofPerVertex +
                              edge_id * domain_fe_type::nDofPerEdge +
                              domain_fe_type::nDofPerEdge - 1 - l );

            FEELPP_ASSERT( localptid_dof != invalid_uint16_type_value &&
                         ( localptid_dof < domain_fe_type::nLocalDof ) )
                ( l )
                (localptid)
                (localptid_dof)
                (domain_fe_type::nLocalDof)
                (_M_pset.pointToEntity( localptid ).second)
                (domain_mesh_type::element_type::numVertices * domain_fe_type::nDofPerVertex +
                 edge_id * domain_fe_type::nDofPerEdge ).error( "invalid value" );
        }
    return localptid_dof;
}

/**
 * \return the P1 Lagrange adaptor  associated to the space \p Xh
 */
template<typename space_type>
boost::shared_ptr<OperatorLagrangeP1<space_type> >
lagrangeP1( boost::shared_ptr<space_type> const& Xh )
{
    return new OperatorLagrangeP1<space_type>( Xh );
}

} // Feel
#endif /* __OperatorLagrangeP1_H */

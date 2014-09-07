/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@imag.fr>
       Date: 2008-01-03

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file quadptlocalization.hpp
   \author Vincent Chabannes <vincent.chabannes@imag.fr>
   \date 2011-06-10
 */

#ifndef __quadptlocalization_H
#define __quadptlocalization_H 1


#include <feel/feeldiscr/mesh.hpp>

#define FEELPP_EXPORT_QUADLOCALIZATION 0
#if FEELPP_EXPORT_QUADLOCALIZATION
#include <feel/feelfilters/exporter.hpp>
#endif


namespace Feel
{


template <typename IteratorRange,typename Im,typename Expr>
class QuadPtLocalization
{
public :

    static const size_type iDim = boost::tuples::template element<0, IteratorRange>::type::value;

    //static const size_type context = Expr::context|vm::POINT;
    static const size_type context = mpl::if_< boost::is_same< mpl::int_<iDim>, mpl::int_<MESH_FACES> >,
                           mpl::int_<Expr::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT>,
                           mpl::int_<Expr::context|vm::POINT> >::type::value;


    //expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT

    typedef IteratorRange range_iterator;


    typedef typename boost::tuples::template element<1, IteratorRange>::type element_iterator_type;
    typedef typename element_iterator_type::value_type::GeoShape GeoShape;
    //typedef Mesh<GeoShape> mesh_type;
#if 0
    typedef typename mesh_type::gm_type gm_type;
    typedef typename mesh_type::element_type geoelement_type;
#else

    //typedef typename element_iterator_type::value_type geoelement_type;
    //typedef typename geoelement_type::gm_type gm_type;

    typedef typename boost::remove_reference<typename element_iterator_type::reference>::type const_t;
    typedef typename boost::remove_const<const_t>::type the_face_element_type;
    typedef typename the_face_element_type::super2::template Element<the_face_element_type>::type the_element_type;
    typedef the_element_type geoelement_type;
    typedef typename geoelement_type::gm_type gm_type;



#endif


    typedef typename gm_type::template Context<context, geoelement_type> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
    typedef typename gm_type::precompute_type pc_type;
    typedef typename gm_type::precompute_ptrtype pc_ptrtype;
#if 0
    typedef typename mesh_type::value_type value_type;
    typedef typename mesh_type::node_type node_type;
    typedef typename matrix_node<value_type>::type matrix_node_type;
#else
    typedef typename geoelement_type::value_type value_type;
    typedef typename geoelement_type::node_type node_type;
    typedef typename geoelement_type::matrix_node_type matrix_node_type;
#endif
    //--------------------------------------------------------------------------------------//

    typedef typename mpl::if_<mpl::bool_<geoelement_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<geoelement_type::nDim, value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<geoelement_type::nDim, value_type, Hypercube>::type >
                              >::type::type im_type;

    typedef typename im_type::face_quadrature_type im_face_type;

    typedef typename QuadMapped<im_type>::permutation_type permutation_type;
    typedef typename QuadMapped<im_type>::permutation_points_type permutation_points_type;

    //--------------------------------------------------------------------------------------//

    // temporary container
    typedef std::map< size_type,std::list<boost::tuple<size_type,size_type,node_type,node_type> > > result_temp1_type;//idTrial->list( q, idq, ptRefTest, ptRefTrial )
    // sub container
    typedef std::map< size_type,boost::tuple< std::vector< boost::tuple<size_type,size_type > >,matrix_node_type,matrix_node_type> > result_temp2_type;//idTrial->( mapFor<q,idq>, ptRefsTest, ptsRefTrial)
    // result container
    typedef std::list< boost::tuple< size_type, result_temp2_type> > result_container_type;// list(idTest, idTrial->( mapForQ, ptRefsTest, ptsRefTrial) )

    // result container for linearform
    typedef std::list< boost::tuple<size_type,std::vector< boost::tuple<size_type,size_type > >,matrix_node_type > > result_container_linear_type;// list(idTest, mapForQ, ptRefsTest)

    //--------------------------------------------------------------------------------------//


    QuadPtLocalization( IteratorRange const& elts /*, im_type const& *//*__im*/ )
        :
        M_listRange(),
        M_im( ),
        M_qm( ),
        M_ppts( M_qm( this->im() ) ),
        M_hasPrecompute(false)
    {
        M_listRange.push_back( elts );
    }

    QuadPtLocalization( std::list<IteratorRange> const& elts )
        :
        M_listRange( elts ),
        M_im( ),
        M_qm( ),
        M_ppts( M_qm( this->im() ) ),
        M_hasPrecompute(false)
    {}


    /**
     * get the integration method
     */
    im_type const& im() const
    {
        return M_im;
    }

    /**
     * get the integration method on face f
     */
    im_face_type  im( uint16_type f ) const
    {
        return M_im.face( f );
    }

    /**
     * begin itrange
     */
    element_iterator_type beginElement() const
    {
        return M_listRange.front().template get<1>();// M_eltbegin;
    }
#if 0
    /**
     * end itrange
     */
    element_iterator_type endElement() const
    {
        return M_listRange.front().template get<2>();//M_eltend;
    }
#endif
    /**
     * container for linear form
     */
    result_container_linear_type const & resultLinear() const
    {
        return M_resLinear;
    }

    /**
     * container for bilinear form
     */
    result_container_type const & result() const
    {
        return M_resBilinear;
    }


    bool hasPrecompute() const { return M_hasPrecompute; }

    //--------------------------------------------------------------------------------------//


    size_type eltForThisQuadPt( size_type theIdq, mpl::int_<MESH_ELEMENTS> )
    {
        return theIdq/this->im().nPoints();
    }
    size_type eltForThisQuadPt( size_type theIdq, mpl::int_<MESH_FACES> )
    {
        return theIdq/this->im().nPointsOnFace();
    }

    size_type qIndexForThisGlobalQuadPt( size_type theIdq,mpl::int_<MESH_ELEMENTS> )
    {
        return theIdq%this->im().nPoints();
    }
    size_type qIndexForThisGlobalQuadPt( size_type theIdq,mpl::int_<MESH_FACES> )
    {
        return theIdq%this->im().nPointsOnFace();
    }

    //--------------------------------------------------------------------------------------//

    gmc_ptrtype gmcForThisElt( size_type theIdElt,
                               std::vector<boost::tuple<size_type,size_type> > const& indexLocalToQuad,
                               mpl::int_<MESH_ELEMENTS> )
    {
#if 1
        element_iterator_type elt_it = this->beginElement(), elt_en = this->beginElement();
        auto itListRange = M_listRange.begin();
        auto const enListRange = M_listRange.end();
        bool findElt=false;
        for ( size_type ide = 0 ; itListRange!=enListRange && !findElt; ++itListRange)
        {
            boost::tie( boost::tuples::ignore, elt_it, elt_en ) = *itListRange;
            const size_type distRange = std::distance(elt_it, elt_en);
            if ( (ide+distRange-1) < theIdElt) { ide+=distRange; continue; }

            for ( ; ide<theIdElt ; ++ide ) ++elt_it;
            findElt = true;
        }
#else
        //search element
        auto elt_it = this->beginElement();
        for ( size_type i=0; i<theIdElt; ++i ) ++elt_it;
#endif
        // get only usefull quad point and reoder
        uint16_type nContextPt = indexLocalToQuad.size();
        matrix_node_type const& quadPtsRef = this->im().points();
        matrix_node_type newquadPtsRef( quadPtsRef.size1() , nContextPt );

        for ( uint16_type i=0; i<nContextPt; ++i )
        {
            ublas::column( newquadPtsRef, i ) = ublas::column( quadPtsRef, indexLocalToQuad[i].template get<0>() );
        }

        // create context
        pc_ptrtype geopc( new pc_type( elt_it->gm(), newquadPtsRef ) );
        gmc_ptrtype gmc( new gmc_type( elt_it->gm(),*elt_it, geopc ) );

        return gmc;
    }

    //--------------------------------------------------------------------------------------//

    gmc_ptrtype gmcForThisElt( size_type theIdElt,
                               std::vector<boost::tuple<size_type,size_type> > const& indexLocalToQuad,
                               mpl::int_<MESH_FACES> )
    {
#if 1
        element_iterator_type elt_it=this->beginElement(), elt_en=this->beginElement();
        auto itListRange = M_listRange.begin();
        auto const enListRange = M_listRange.end();
        bool findElt=false;
        for ( size_type ide = 0 ; itListRange!=enListRange && !findElt; ++itListRange)
        {
            boost::tie( boost::tuples::ignore, elt_it, elt_en ) = *itListRange;
            const size_type distRange = std::distance(elt_it, elt_en);
            if ( (ide+distRange-1) < theIdElt) { ide+=distRange; continue; }

            for ( ; ide<theIdElt ; ++ide ) ++elt_it;
            findElt = true;
        }
#else
        //search element
        auto elt_it = this->beginElement();
        for ( size_type i=0; i<theIdElt; ++i ) ++elt_it;
#endif
        // get only usefull quad point and reoder
        const uint16_type nContextPt = indexLocalToQuad.size();
        const uint16_type __face_id_in_elt_0 = elt_it->pos_first();

        auto const __perm = elt_it->element( 0 ).permutation( __face_id_in_elt_0 );
        matrix_node_type const& quadPtsRef =  M_ppts[ __face_id_in_elt_0].find( __perm )->second;
        //matrix_node_type quadPtsRef = this->im().points();
        matrix_node_type newquadPtsRef( quadPtsRef.size1() , nContextPt );

        for ( uint16_type i=0; i<nContextPt; ++i )
        {
            ublas::column( newquadPtsRef, i ) = ublas::column( quadPtsRef, indexLocalToQuad[i].template get<0>() );
        }

        // create context
        std::vector<std::map<permutation_type, pc_ptrtype> > __geopc( this->im().nFaces() );

        for ( uint16_type __f = 0; __f < this->im().nFaces(); ++__f )
        {
            for ( permutation_type __p( permutation_type::IDENTITY );
                    __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
            {
                __geopc[__f][__p] = pc_ptrtype(  new pc_type( elt_it->element( 0 ).gm(),
                                                 newquadPtsRef ) );
            }
        }

        gmc_ptrtype gmc( new gmc_type( elt_it->element( 0 ).gm(),
                                       elt_it->element( 0 ),
                                       __geopc,
                                       __face_id_in_elt_0  ) );

        return gmc;
    }

    //--------------------------------------------------------------------------------------//

    std::vector< boost::tuple< std::vector<boost::tuple<size_type,size_type> >,gmc_ptrtype,matrix_node_type > >
    getUsableDataInFormContext( std::vector<boost::tuple<size_type,size_type> > const& indexLocalToQuad,
                                matrix_node_type const & ptsRefTest )
    {
        // compute the number of several elements
        uint16_type nContextPt = indexLocalToQuad.size();
        std::map<size_type,std::list<size_type> > mapEltId;

        for ( uint16_type i=0; i<nContextPt; ++i )
        {
            size_type eltId = this->eltForThisQuadPt( indexLocalToQuad[i].template get<1>(),mpl::int_<iDim>() );
            mapEltId[eltId].push_back( i );
        }

        // number of elements
        auto nEltInContext = mapEltId.size();

        // the vector result
        std::vector< boost::tuple< std::vector<boost::tuple<size_type,size_type> >, gmc_ptrtype, matrix_node_type > > vec_res( nEltInContext );

        auto map_it = mapEltId.begin();
        auto map_en = mapEltId.end();
        size_type cptId = 0;

        for ( ; map_it!= map_en ; ++map_it, ++cptId )
        {
            // subdivide into elements the map indexLocalToQuad
            std::vector<boost::tuple<size_type,size_type> > newindexLocalToQuad( map_it->second.size() );
            // subdivide into elements the ptsRef
            matrix_node_type newptsRefTest( ptsRefTest.size1(), map_it->second.size() );
            auto sublist_it = map_it->second.begin();
            auto sublist_en = map_it->second.end();
            size_type index = 0;

            for ( ; sublist_it != sublist_en ; ++sublist_it,++index )
            {
                newindexLocalToQuad[index] = indexLocalToQuad[*sublist_it];
                ublas::column( newptsRefTest, index ) = ublas::column( ptsRefTest,*sublist_it );
            }
            // get the corresponding gmc
            auto thegmc = gmcForThisElt( map_it->first,newindexLocalToQuad,mpl::int_<iDim>() );
            // add to result
            vec_res[cptId] = boost::make_tuple( newindexLocalToQuad, thegmc, newptsRefTest );
        }

        return vec_res;
    }

    //--------------------------------------------------------------------------------------//

    std::vector< boost::tuple< std::vector<boost::tuple<size_type,size_type> >,gmc_ptrtype,matrix_node_type, matrix_node_type > >
    getUsableDataInFormContext( std::vector<boost::tuple<size_type,size_type> > const& indexLocalToQuad,
                                matrix_node_type const & ptsRefTest,
                                matrix_node_type const & ptsRefTrial )
    {
        // compute the number of several elements
        uint16_type nContextPt = indexLocalToQuad.size();
        std::map<size_type,std::list<size_type> > mapEltId;

        for ( uint16_type i=0; i<nContextPt; ++i )
        {
            size_type eltId = this->eltForThisQuadPt( indexLocalToQuad[i].template get<1>(),mpl::int_<iDim>() );
            mapEltId[eltId].push_back( i );
        }

        // number of elements
        auto nEltInContext = mapEltId.size();

        // the vector result
        std::vector< boost::tuple< std::vector<boost::tuple<size_type,size_type> >, gmc_ptrtype, matrix_node_type, matrix_node_type > > vec_res( nEltInContext );

        auto map_it = mapEltId.begin();
        auto map_en = mapEltId.end();
        size_type cptId = 0;

        for ( ; map_it!= map_en ; ++map_it, ++cptId )
        {
            // subdivide into elements the map indexLocalToQuad
            std::vector<boost::tuple<size_type,size_type> > newindexLocalToQuad( map_it->second.size() );
            // subdivide into elements the ptsRef
            matrix_node_type newptsRefTest( ptsRefTest.size1(), map_it->second.size() );
            matrix_node_type newptsRefTrial( ptsRefTrial.size1(), map_it->second.size() );

            auto sublist_it = map_it->second.begin();
            auto sublist_en = map_it->second.end();
            size_type index = 0;

            for ( ; sublist_it != sublist_en ; ++sublist_it,++index )
            {
                newindexLocalToQuad[index] = indexLocalToQuad[*sublist_it];
                ublas::column( newptsRefTest, index ) = ublas::column( ptsRefTest,*sublist_it );
                ublas::column( newptsRefTrial, index ) = ublas::column( ptsRefTrial,*sublist_it );
            }

            // get the corresponding gmc
            auto thegmc = gmcForThisElt( map_it->first,newindexLocalToQuad,mpl::int_<iDim>() );
            // add to result
            vec_res[cptId] = boost::make_tuple( newindexLocalToQuad, thegmc, newptsRefTest,newptsRefTrial );
        }

        return vec_res;
    }

    //--------------------------------------------------------------------------------------//

    template <typename Mesh1Type>
    void
    localization( mpl::int_<MESH_FACES> /**/,
                  boost::shared_ptr<Mesh1Type> meshTest,
                  std::vector<std::list<boost::tuple< size_type,size_type,node_type> > > & testEltToPtsQuad )
    {

        auto begin_elt_it = this->beginElement();

        std::vector<std::map<permutation_type, pc_ptrtype> > __geopc( this->im().nFaces() );
        //typedef typename im_type::face_quadrature_type face_im_type;

        //std::vector<face_im_type> face_ims( this->im().nFaces() );

        for ( uint16_type __f = 0; __f < this->im().nFaces(); ++__f )
        {
            //face_ims[__f] = this->im( __f );

            for ( permutation_type __p( permutation_type::IDENTITY );
                    __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
            {
                //FEELPP_ASSERT( ppts[__f].find(__p)->second.size2() != 0 ).warn( "invalid quadrature type" );
                __geopc[__f][__p] = pc_ptrtype(  new pc_type( begin_elt_it->element( 0 ).gm(), M_ppts[__f].find( __p )->second ) );
            }
        }


        uint16_type __face_id_in_elt_0 = begin_elt_it->pos_first();

        gmc_ptrtype gmc( new gmc_type( begin_elt_it->element( 0 ).gm(),
                                       begin_elt_it->element( 0 ),
                                       __geopc,
                                       __face_id_in_elt_0 ) );

        auto meshTestLocalization = meshTest->tool_localization();
        meshTestLocalization->updateForUse();
        const auto nbNearNeighborAtStartTest = meshTestLocalization->kdtree()->nPtMaxNearNeighbor();
        bool notUseOptLocTest = Mesh1Type::nDim!=Mesh1Type::nRealDim;
        if (notUseOptLocTest) meshTestLocalization->kdtree()->nbNearNeighbor(Mesh1Type::element_type::numPoints);

        matrix_node_type ptsReal( begin_elt_it->vertices().size1(), 1 );
        size_type testIdElt=0;
        node_type testNodeRef;

        element_iterator_type elt_it, elt_en;

        auto itListRange = M_listRange.begin();
        auto const enListRange = M_listRange.end();
        for ( size_type ide = 0 ; itListRange!=enListRange ; ++itListRange)
        {
            boost::tie( boost::tuples::ignore, elt_it, elt_en ) = *itListRange;
            for ( ; elt_it != elt_en; ++elt_it, ++ide )
            {
                __face_id_in_elt_0 = elt_it->pos_first();
                //if ( elt_it->isConnectedTo1()) std::cout << "\n AIEEEEEE!!!!!!!!!!!\n";

                gmc->update( elt_it->element( 0 ), __face_id_in_elt_0 );

                //std::cout << "\n quad gmc "<< gmc->xReal();
                for ( int q = 0; q <  gmc->nPoints(); ++ q )
                {
                    size_type idq = gmc->nPoints()*ide+q;
#if 0
                    auto testAnalysis = meshTestLocalization->searchElement( gmc->xReal( q ) );
                    testIdElt = testAnalysis.template get<1>();
                    testNodeRef = testAnalysis.template get<2>();
#else
                    ublas::column(ptsReal,0 ) = gmc->xReal( q );
                    if (notUseOptLocTest) testIdElt=invalid_size_type_value;
                    auto resLocalisationTest = meshTestLocalization->run_analysis(ptsReal,testIdElt,elt_it->vertices(),mpl::int_<0>());
                    testIdElt = resLocalisationTest.template get<1>();
                    testNodeRef = meshTestLocalization->result_analysis().begin()->second.begin()->template get<1>();
#endif
                    testEltToPtsQuad[testIdElt].push_back( boost::make_tuple( idq,q,testNodeRef ) );
                }
            }
        }

        if (notUseOptLocTest) meshTestLocalization->kdtree()->nbNearNeighbor(nbNearNeighborAtStartTest);

    }

    //--------------------------------------------------------------------------------------//

    template <typename Mesh1Type>
    void
    localization( mpl::int_<MESH_ELEMENTS> /**/,
                  boost::shared_ptr<Mesh1Type> meshTest,
                  std::vector<std::list<boost::tuple< size_type,size_type,node_type> > > & testEltToPtsQuad )
    {

        auto begin_elt_it = this->beginElement();
        //auto elt_en = this->endElement();

        pc_ptrtype geopc( new pc_type( begin_elt_it->gm(), this->im().points() ) );
        gmc_ptrtype gmc( new gmc_type( begin_elt_it->gm(),*begin_elt_it, geopc ) );

        auto meshTestLocalization = meshTest->tool_localization();
        meshTestLocalization->updateForUse();
        const auto nbNearNeighborAtStartTest = meshTestLocalization->kdtree()->nPtMaxNearNeighbor();
        bool notUseOptLocTest = Mesh1Type::nDim!=Mesh1Type::nRealDim;
        if (notUseOptLocTest) meshTestLocalization->kdtree()->nbNearNeighbor(Mesh1Type::element_type::numPoints);

        matrix_node_type ptsReal( begin_elt_it->vertices().size1(), 1 );
        size_type testIdElt=0;
        node_type testNodeRef;

        element_iterator_type elt_it, elt_en;

        auto itListRange = M_listRange.begin();
        auto const enListRange = M_listRange.end();
        for ( size_type ide = 0 ; itListRange!=enListRange ; ++itListRange)
        {
            boost::tie( boost::tuples::ignore, elt_it, elt_en ) = *itListRange;
            for ( ; elt_it != elt_en; ++elt_it, ++ide )
            {
                gmc->update( *elt_it );

                for ( int q = 0; q <  gmc->nPoints(); ++ q )
                {
                    // cpt of quad pt
                    size_type idq = gmc->nPoints()*ide+q;
                    // search in test mesh
#if 0
                    auto testAnalysis = meshTestLocalization->searchElement( gmc->xReal( q ) );
                    testIdElt = testAnalysis.template get<1>();
                    testNodeRef = testAnalysis.template get<2>();
#else
                    ublas::column(ptsReal,0 ) = gmc->xReal( q );
                    if (notUseOptLocTest) testIdElt=invalid_size_type_value;
                    auto resLocalisationTest = meshTestLocalization->run_analysis(ptsReal,testIdElt,elt_it->vertices(),mpl::int_<0>());
                    testIdElt = resLocalisationTest.template get<1>();
                    testNodeRef = meshTestLocalization->result_analysis().begin()->second.begin()->template get<1>();
#endif

                    testEltToPtsQuad[testIdElt].push_back( boost::make_tuple( idq,q,testNodeRef ) );
                }
            }
        }

        if (notUseOptLocTest) meshTestLocalization->kdtree()->nbNearNeighbor(nbNearNeighborAtStartTest);

    }

    //--------------------------------------------------------------------------------------//

    template <typename Mesh1Type,typename Mesh2Type>
    void
    localization( mpl::int_<MESH_FACES> /**/,
                  boost::shared_ptr<Mesh1Type> meshTest,
                  boost::shared_ptr<Mesh2Type> meshTrial,
                  std::vector<std::list<boost::tuple< size_type,size_type,node_type> > > & testEltToPtsQuad,
                  std::vector<std::list<boost::tuple< size_type,size_type,node_type> > > & trialEltToPtsQuad,
                  std::vector<std::list<size_type> > & EltCoupled )
    {
        //std::cout << "[QuadPtLocalization] : localization<MESH_FACES>(bilinear form start" << std::endl;

        auto begin_elt_it = this->beginElement();

        std::vector<std::map<permutation_type, pc_ptrtype> > __geopc( this->im().nFaces() );

        for ( uint16_type __f = 0; __f < this->im().nFaces(); ++__f )
        {
            for ( permutation_type __p( permutation_type::IDENTITY );
                    __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
            {
                //FEELPP_ASSERT( ppts[__f].find(__p)->second.size2() != 0 ).warn( "invalid quadrature type" );
                __geopc[__f][__p] = pc_ptrtype(  new pc_type( begin_elt_it->element( 0 ).gm(), M_ppts[__f].find( __p )->second ) );
            }
        }

        uint16_type __face_id_in_elt_0 = begin_elt_it->pos_first();

        gmc_ptrtype gmc( new gmc_type( begin_elt_it->element( 0 ).gm(),
                                       begin_elt_it->element( 0 ),
                                       __geopc,
                                       __face_id_in_elt_0 ) );


        auto meshTrialLocalization = meshTrial->tool_localization();
        meshTrialLocalization->updateForUse();
        const auto nbNearNeighborAtStartTrial = meshTrialLocalization->kdtree()->nPtMaxNearNeighbor();
        bool notUseOptLocTrial = Mesh2Type::nDim!=Mesh2Type::nRealDim;
        if (notUseOptLocTrial) meshTrialLocalization->kdtree()->nbNearNeighbor(Mesh2Type::element_type::numPoints);

        auto meshTestLocalization = meshTest->tool_localization();
        meshTestLocalization->updateForUse();
        const auto nbNearNeighborAtStartTest = meshTestLocalization->kdtree()->nPtMaxNearNeighbor();
        bool notUseOptLocTest = Mesh1Type::nDim!=Mesh1Type::nRealDim;
        if (notUseOptLocTest) meshTestLocalization->kdtree()->nbNearNeighbor(Mesh1Type::element_type::numPoints);


        matrix_node_type ptsReal( begin_elt_it->vertices().size1(), 1 );
        size_type trialIdElt = 0, testIdElt=0;
        node_type trialNodeRef,testNodeRef;

        bool quadMeshIsSameThatTrialMesh=false,quadMeshIsSameThatTestMesh=false;
        if ( dynamic_cast<void*>( const_cast<MeshBase*>( begin_elt_it->mesh() ) ) == dynamic_cast<void*>( meshTrial.get() ) )
            quadMeshIsSameThatTrialMesh=true;
        if ( dynamic_cast<void*>( const_cast<MeshBase*>( begin_elt_it->mesh() ) ) == dynamic_cast<void*>( meshTest.get() ) )
            quadMeshIsSameThatTestMesh=true;

#if FEELPP_EXPORT_QUADLOCALIZATION
        std::map<size_type,std::list<size_type> > mapBetweenMeshes_test;
        std::map<size_type,std::list<size_type> > mapBetweenMeshes_trial;
#endif

        element_iterator_type elt_it, elt_en;

        auto itListRange = M_listRange.begin();
        auto const enListRange = M_listRange.end();
        for ( size_type ide = 0 ; itListRange!=enListRange ; ++itListRange)
        {
            boost::tie( boost::tuples::ignore, elt_it, elt_en ) = *itListRange;
            for ( ; elt_it != elt_en; ++elt_it, ++ide )
            {
                __face_id_in_elt_0 = elt_it->pos_first();
                //if ( elt_it->isConnectedTo1()) std::cout << "\n AIEEEEEE!!!!!!!!!!!\n";

                gmc->update( elt_it->element( 0 ), __face_id_in_elt_0 );

                //std::cout << "\n quad gmc "<< gmc->xReal();

                for ( int q = 0; q <  gmc->nPoints(); ++ q )
                {
                    // cpt of quad pt
                    size_type idq = gmc->nPoints()*ide+q;

                    // search in trial mesh
                    ublas::column(ptsReal,0 ) = gmc->xReal( q );
                    if (!quadMeshIsSameThatTrialMesh)
                    {
                        if (notUseOptLocTrial) trialIdElt=invalid_size_type_value;
                        auto resLocalisationTrial = meshTrialLocalization->run_analysis(ptsReal,trialIdElt,elt_it->vertices(),mpl::int_<0>());
                        trialIdElt = resLocalisationTrial.template get<1>();
                        trialNodeRef = meshTrialLocalization->result_analysis().begin()->second.begin()->template get<1>();
                    }
                    else
                    {
                        trialIdElt = gmc->id();
                        trialNodeRef = gmc->xRef(q);
                    }
                    trialEltToPtsQuad[trialIdElt].push_back( boost::make_tuple( idq,q,trialNodeRef ) );

                    // search in test mesh
                    if (!quadMeshIsSameThatTestMesh)
                    {
                        if (notUseOptLocTest) testIdElt=invalid_size_type_value;
                        auto resLocalisationTest = meshTestLocalization->run_analysis(ptsReal,testIdElt,elt_it->vertices(),mpl::int_<0>());
                        testIdElt = resLocalisationTest.template get<1>();
                        testNodeRef = meshTestLocalization->result_analysis().begin()->second.begin()->template get<1>();
                    }
                    else
                    {
                        testIdElt = gmc->id();
                        testNodeRef = gmc->xRef(q);
                    }
                    testEltToPtsQuad[testIdElt].push_back( boost::make_tuple( idq,q,testNodeRef ) );

                    // relation between test and trial
                    if ( std::find( EltCoupled[testIdElt].begin(),EltCoupled[testIdElt].end(),trialIdElt )==EltCoupled[testIdElt].end() )
                    {
                        EltCoupled[testIdElt].push_back( trialIdElt );
                        /*std::cout << "gmc->xRef(q)" << gmc->xReal(q)
                                  << " testIdElt " << testIdElt << " meshTest->element().G()" << meshTest->element(testIdElt).G()
                                  << " trialIdElt " << trialIdElt << " meshTrial->element().G()" << meshTrial->element(trialIdElt).G()
                                  << std::endl;*/
#if FEELPP_EXPORT_QUADLOCALIZATION
                        mapBetweenMeshes_test[testIdElt].push_back( trialIdElt );
                        mapBetweenMeshes_trial[trialIdElt].push_back( testIdElt );
#endif
                    }
                } // for ( int q = 0; q <  gmc->nPoints(); ++ q )
            } // for ( ; elt_it ...)
        } // end for( size_type ide ... )
#if FEELPP_EXPORT_QUADLOCALIZATION
        this->exportQuadLocalization(meshTest,meshTrial,mapBetweenMeshes_test,mapBetweenMeshes_trial);
#endif

        if (notUseOptLocTest) meshTestLocalization->kdtree()->nbNearNeighbor(nbNearNeighborAtStartTest);
        if (notUseOptLocTrial) meshTrialLocalization->kdtree()->nbNearNeighbor(nbNearNeighborAtStartTrial);

        //std::cout << "[QuadPtLocalization] : localization<MESH_FACES>(bilinear form finish" << std::endl;

    } // localization

    //--------------------------------------------------------------------------------------//

    template <typename Mesh1Type,typename Mesh2Type>
    void
    localization( mpl::int_<MESH_ELEMENTS> /**/,
                  boost::shared_ptr<Mesh1Type> meshTest,
                  boost::shared_ptr<Mesh2Type> meshTrial,
                  std::vector<std::list<boost::tuple< size_type,size_type,node_type> > > & testEltToPtsQuad,
                  std::vector<std::list<boost::tuple< size_type,size_type,node_type> > > & trialEltToPtsQuad,
                  std::vector<std::list<size_type> > & EltCoupled )
    {
        //std::cout << "[QuadPtLocalization] : localization<MESH_ELEMENTS>(bilinear form start" << std::endl;

        auto begin_elt_it = this->beginElement();
        //auto elt_it = this->beginElement();
        //auto elt_en = this->endElement();

        pc_ptrtype geopc( new pc_type( begin_elt_it->gm(), this->im().points() ) );
        gmc_ptrtype gmc( new gmc_type( begin_elt_it->gm(),*begin_elt_it, geopc ) );

        auto meshTrialLocalization = meshTrial->tool_localization();
        meshTrialLocalization->updateForUse();
        const auto nbNearNeighborAtStartTrial = meshTrialLocalization->kdtree()->nPtMaxNearNeighbor();
        bool notUseOptLocTrial = Mesh2Type::nDim!=Mesh2Type::nRealDim;
        if (notUseOptLocTrial) meshTrialLocalization->kdtree()->nbNearNeighbor(Mesh2Type::element_type::numPoints);

        auto meshTestLocalization = meshTest->tool_localization();
        meshTestLocalization->updateForUse();
        const auto nbNearNeighborAtStartTest = meshTestLocalization->kdtree()->nPtMaxNearNeighbor();
        bool notUseOptLocTest = Mesh1Type::nDim!=Mesh1Type::nRealDim;
        if (notUseOptLocTest) meshTestLocalization->kdtree()->nbNearNeighbor(Mesh1Type::element_type::numPoints);

#if FEELPP_EXPORT_QUADLOCALIZATION
        std::map<size_type,std::list<size_type> > mapBetweenMeshes_test;
        std::map<size_type,std::list<size_type> > mapBetweenMeshes_trial;
#endif


        //auto nQuadPtsInElt = this->im().nPoints();
        //auto nElts = std::distance(elt_it,elt_en);
        //auto nQuadPts = nElts*nQuadPtsInElt;
        // auto nEltTrial= meshTrial->numElements();
        //auto nEltTest= meshTest->numElements();

        matrix_node_type ptsReal( begin_elt_it->vertices().size1(), 1 );
        size_type trialIdElt = 0, testIdElt=0;
        node_type trialNodeRef,testNodeRef;

        bool quadMeshIsSameThatTrialMesh=false,quadMeshIsSameThatTestMesh=false;
        if ( dynamic_cast<void*>( const_cast<MeshBase*>( begin_elt_it->mesh() ) ) == dynamic_cast<void*>( meshTrial.get() ) )
            quadMeshIsSameThatTrialMesh=true;
        if ( dynamic_cast<void*>( const_cast<MeshBase*>( begin_elt_it->mesh() ) ) == dynamic_cast<void*>( meshTest.get() ) )
            quadMeshIsSameThatTestMesh=true;

        element_iterator_type elt_it, elt_en;

        auto itListRange = M_listRange.begin();
        auto const enListRange = M_listRange.end();
        for ( size_type ide = 0 ; itListRange!=enListRange ; ++itListRange)
        {
            boost::tie( boost::tuples::ignore, elt_it, elt_en ) = *itListRange;
            for ( ; elt_it != elt_en; ++elt_it, ++ide )
            {
                gmc->update( *elt_it );

                for ( int q = 0; q <  gmc->nPoints(); ++ q )
                {
                    // cpt of quad pt
                    size_type idq = gmc->nPoints()*ide+q;

                    // search in trial mesh
                    ublas::column(ptsReal,0 ) = gmc->xReal( q );
                    if (!quadMeshIsSameThatTrialMesh)
                    {
                        if (notUseOptLocTrial) trialIdElt=invalid_size_type_value;
                        auto resLocalisationTrial = meshTrialLocalization->run_analysis(ptsReal,trialIdElt,elt_it->vertices(),mpl::int_<0>());
                        trialIdElt = resLocalisationTrial.template get<1>();
                        trialNodeRef = meshTrialLocalization->result_analysis().begin()->second.begin()->template get<1>();
                    }
                    else
                    {
                        trialIdElt = gmc->id();
                        trialNodeRef = gmc->xRef(q);
                    }
                    //std::cout << "\n trialNodeRef " << trialNodeRef << std::endl;
                    trialEltToPtsQuad[trialIdElt].push_back( boost::make_tuple( idq,q,trialNodeRef ) );

                    // search in test mesh
                    if (!quadMeshIsSameThatTestMesh)
                    {
                        if (notUseOptLocTest) testIdElt=invalid_size_type_value;
                        auto resLocalisationTest = meshTestLocalization->run_analysis(ptsReal,testIdElt,elt_it->vertices(),mpl::int_<0>());
                        testIdElt = resLocalisationTest.template get<1>();
                        testNodeRef = meshTestLocalization->result_analysis().begin()->second.begin()->template get<1>();
                    }
                    else
                    {
                        testIdElt = gmc->id();
                        testNodeRef = gmc->xRef(q);
                    }
                    testEltToPtsQuad[testIdElt].push_back( boost::make_tuple( idq,q,testNodeRef ) );

                    //std::cout << " AAtestIdElt " << testIdElt << " AAtrialIdElt " << trialIdElt << std::endl;

                    // relation between test and trial
                    if ( std::find( EltCoupled[testIdElt].begin(),EltCoupled[testIdElt].end(),trialIdElt )==EltCoupled[testIdElt].end() )
                    {
                        EltCoupled[testIdElt].push_back( trialIdElt );
                        //std::cout << " AAtestIdElt " << testIdElt << " AAtrialIdElt " << trialIdElt << std::endl;
#if FEELPP_EXPORT_QUADLOCALIZATION
                        mapBetweenMeshes_test[testIdElt].push_back( trialIdElt );
                        mapBetweenMeshes_trial[trialIdElt].push_back( testIdElt );
#endif
                    }


#if 0 // try with neighbors
                    auto const& geoeltTrial = meshTrial->element( trialIdElt );
                    for ( uint16_type ms=0; ms < geoeltTrial.nNeighbors(); ms++ )
                    {
                        const size_type neighborTrial_id = geoeltTrial.neighbor( ms ).first;
                        if ( neighborTrial_id==invalid_size_type_value ) continue;

                        auto resNeighboor = meshTrialLocalization->isIn( neighborTrial_id,gmc->xReal( q ),elt_it->vertices(),mpl::int_<1>() );
                        if (resNeighboor.get<0>())
                        {
                            trialEltToPtsQuad[neighborTrial_id].push_back( boost::make_tuple( idq,q,resNeighboor.get<1>() ) );

                            if ( std::find( EltCoupled[testIdElt].begin(),EltCoupled[testIdElt].end(),neighborTrial_id )==EltCoupled[testIdElt].end() )
                            {
                                EltCoupled[testIdElt].push_back( neighborTrial_id );
                            }
#if FEELPP_EXPORT_QUADLOCALIZATION
                            mapBetweenMeshes_trial[neighborTrial_id].push_back( testIdElt );
#endif
                        }
                    } // for ( uint16_type ms=0;...

                    auto const& geoeltTest = meshTest->element( testIdElt );
                    for ( uint16_type ms=0; ms < geoeltTest.nNeighbors(); ms++ )
                    {
                        const size_type neighborTest_id = geoeltTest.neighbor( ms ).first;
                        if ( neighborTest_id==invalid_size_type_value ) continue;

                        auto resNeighboor = meshTestLocalization->isIn( neighborTest_id,gmc->xReal( q ),elt_it->vertices(),mpl::int_<1>() );
                        if (resNeighboor.get<0>())
                        {
                            testEltToPtsQuad[neighborTest_id].push_back( boost::make_tuple( idq,q,resNeighboor.get<1>() ) );

                            if ( std::find( EltCoupled[neighborTest_id].begin(),EltCoupled[neighborTest_id].end(),trialIdElt )==EltCoupled[neighborTest_id].end() )
                            {
                                EltCoupled[neighborTest_id].push_back( trialIdElt );
                            }
#if FEELPP_EXPORT_QUADLOCALIZATION
                            mapBetweenMeshes_test[neighborTest_id].push_back( trialIdElt );
#endif
                        }
                    } // for ( uint16_type ms=0;...
#endif
                } // for ( int q = 0; q <  gmc->nPoints(); ++ q )
            } // for ( ; elt_it != elt_en; ++elt_it, ++ide )
        } // end for( size_type ide ... )
#if FEELPP_EXPORT_QUADLOCALIZATION
        this->exportQuadLocalization(meshTest,meshTrial,mapBetweenMeshes_test,mapBetweenMeshes_trial);
#endif

        if (notUseOptLocTest) meshTestLocalization->kdtree()->nbNearNeighbor(nbNearNeighborAtStartTest);
        if (notUseOptLocTrial) meshTrialLocalization->kdtree()->nbNearNeighbor(nbNearNeighborAtStartTrial);

        //std::cout << "[QuadPtLocalization] : localization<MESH_ELEMENTS>(bilinear form finish" << std::endl;
    }

    //--------------------------------------------------------------------------------------//

    template <typename Mesh1Type>
    void
    update( boost::shared_ptr<Mesh1Type> meshTest )
    {
        typedef Mesh1Type mesh_test_type;

        //clean result
        M_resLinear.clear();

        // localize quad pts on meshTest -> storage in testEltToPtsQuad
        auto nEltTest= meshTest->numElements();
        std::vector<std::list<boost::tuple< size_type,size_type,node_type> > > testEltToPtsQuad( nEltTest );
        this->localization( mpl::int_<iDim>(),meshTest,testEltToPtsQuad );

        //build a efficent container : M_resLinear
        size_type theIdEltTest = 0;
        auto eltTest_it = testEltToPtsQuad.begin();
        auto eltTest_en = testEltToPtsQuad.end();

        for ( ; eltTest_it != eltTest_en ; ++eltTest_it, ++theIdEltTest )
        {
            if ( eltTest_it->size()>0 )
            {
                auto nPtsRef = eltTest_it->size();
                matrix_node_type ptsRefTest( mesh_test_type::nDim,nPtsRef );
                std::vector<boost::tuple<size_type,size_type> > indexLocalToQuad( nPtsRef );
                // build
                auto idq_it = eltTest_it->begin();
                auto idq_en = eltTest_it->end();
                uint16_type cptIdq = 0;

                for ( ; idq_it != idq_en ; ++idq_it, ++cptIdq )
                {
                    indexLocalToQuad[cptIdq] = boost::make_tuple( idq_it->template get<1>(),idq_it->template get<0>() );
                    ublas::column( ptsRefTest, cptIdq ) = idq_it->template get<2>();
                }

                M_resLinear.push_back( boost::make_tuple( theIdEltTest, indexLocalToQuad, ptsRefTest ) );
            }
        }

    }

    //--------------------------------------------------------------------------------------//

    template <typename Mesh1Type,typename Mesh2Type>
    void
    update( boost::shared_ptr<Mesh1Type> meshTest, boost::shared_ptr<Mesh2Type> meshTrial )
    {
        typedef Mesh1Type mesh_test_type;
        typedef Mesh2Type mesh_trial_type;

        //clean result
        M_resBilinear.clear();

        auto nEltTrial= meshTrial->numElements();
        auto nEltTest= meshTest->numElements();
        std::vector<std::list<boost::tuple< size_type,size_type,node_type> > > trialEltToPtsQuad( nEltTrial ) ;
        std::vector<std::list<boost::tuple< size_type,size_type,node_type> > > testEltToPtsQuad( nEltTest );
        std::vector<std::list<size_type> > EltCoupled( nEltTest );

        // localize quad pts on meshTest and meshTrial  -> storage in testEltToPtsQuad, trialEltToPtsQuad and EltCoupled
        this->localization( mpl::int_<iDim>(),
                            meshTest,
                            meshTrial,
                            testEltToPtsQuad,
                            trialEltToPtsQuad,
                            EltCoupled );

        //build a efficent container : M_resBilinear
        auto eltCoupled_it = EltCoupled.begin();
        auto eltCoupled_en = EltCoupled.end();
        size_type theIdEltTest = 0;

        for ( ; eltCoupled_it != eltCoupled_en ; ++eltCoupled_it, ++theIdEltTest )
        {
            if (eltCoupled_it->size()==0) continue;

            //init map
            result_temp1_type mapTrial2idq;

            // search
            auto quadPtTest_it = testEltToPtsQuad[theIdEltTest].begin();
            auto const quadPtTest_en = testEltToPtsQuad[theIdEltTest].end();
            for ( ; quadPtTest_it != quadPtTest_en ; ++quadPtTest_it )
            {
                // get test data
                auto theIdq = quadPtTest_it->template get<0>();
                auto theq = quadPtTest_it->template get<1>();
                //std::cout << "\nidq modulo : " << theIdq%nQuadPtsInElt << " q " << theq << " elt " << theIdq/nQuadPtsInElt;
                auto theNodeRefTest = quadPtTest_it->template get<2>();
                // search in trial data
                auto theRes = this->findQuadPt( *eltCoupled_it,trialEltToPtsQuad,theIdq );
                // get trial data
                auto const theIdEltTrial = theRes.template get<0>();
                auto const theNodeRefTrial = theRes.template get<1>();
                // add into the temporary map

                //std::cout << " theIdEltTrial " << theIdEltTrial <<"("<<meshTrial->numElements()<<")" << " theIdEltTest " << theIdEltTest <<"("<<meshTest->numElements()<<")"<<std::endl;
                mapTrial2idq[theIdEltTrial].push_back( boost::make_tuple( theq,theIdq,theNodeRefTest,theNodeRefTrial ) );
            }

            result_temp2_type mapiIdTrial2qAndPtRef;

            // build matrix_node of ptsRef (test and trial)
            auto map_it = mapTrial2idq.begin();
            auto map_en = mapTrial2idq.end();
            for ( ; map_it != map_en ; ++map_it )
            {
                // init matrix node
                auto nPtsRef = map_it->second.size();
                //std::cout << "nPtsRef " << nPtsRef << std::endl;
                matrix_node_type ptsRefTest( mesh_test_type::nDim,nPtsRef );
                matrix_node_type ptsRefTrial( mesh_trial_type::nDim,nPtsRef );
                // init map between index node and index quad point
                std::vector<boost::tuple<size_type,size_type> > indexLocalToQuad( nPtsRef );
                // build
                auto idq_it = map_it->second.begin();
                auto idq_en = map_it->second.end();
                uint16_type cptIdq = 0;

                for ( ; idq_it != idq_en ; ++idq_it, ++cptIdq )
                {
                    indexLocalToQuad[cptIdq] = boost::make_tuple( idq_it->template get<0>(),idq_it->template get<1>() );
                    ublas::column( ptsRefTest, cptIdq ) = idq_it->template get<2>();
                    ublas::column( ptsRefTrial, cptIdq ) = idq_it->template get<3>();
                }

                mapiIdTrial2qAndPtRef[map_it->first] = boost::make_tuple( indexLocalToQuad,ptsRefTest,ptsRefTrial );
            }

            // add to result container
            M_resBilinear.push_back( boost::make_tuple( theIdEltTest,mapiIdTrial2qAndPtRef ) );
        }

    } // update

private :

    boost::tuple<size_type,node_type >
    findQuadPt( std::list<size_type> const & eltCoupled_it,
                std::vector<std::list<boost::tuple< size_type,size_type,node_type> > > const & trialPtQuadToElt,
                size_type idq )
    {
        size_type resId=0;
        ublas::vector<double> resNodeRef;
        bool find=false;

        auto eltTrial_it = eltCoupled_it.begin();
        auto eltTrial_en = eltCoupled_it.end();

        while ( !find && eltTrial_it != eltTrial_en )
        {
            auto quadPtTrial_it = trialPtQuadToElt[*eltTrial_it].begin();
            auto quadPtTrial_en = trialPtQuadToElt[*eltTrial_it].end();

            while ( !find && quadPtTrial_it != quadPtTrial_en )
            {
                if ( quadPtTrial_it->template get<0>() == idq )
                {
                    find = true;
                    resNodeRef=quadPtTrial_it->template get<2>();
                }

                else ++quadPtTrial_it;
            }

            if ( !find ) ++eltTrial_it;

            else resId=*eltTrial_it;
        }

        if ( !find ) std::cout << "BUG in findQuadPt " << std::endl;

        return boost::make_tuple( resId,resNodeRef );
    }


#if FEELPP_EXPORT_QUADLOCALIZATION
    template <typename Mesh1Type,typename Mesh2Type>
    void
    exportQuadLocalization( boost::shared_ptr<Mesh1Type> meshTest,
                            boost::shared_ptr<Mesh2Type> meshTrial,
                            std::map<size_type,std::list<size_type> > const& mapBetweenMeshes_test,
                            std::map<size_type,std::list<size_type> > const& mapBetweenMeshes_trial ) const
    {
        typedef FunctionSpace<Mesh1Type, bases<Lagrange<0, Scalar,Discontinuous> > > space_test_disc_type;
        auto spaceGraphProjTest = space_test_disc_type::New( meshTest );
        auto elemTest_itt  = meshTest->beginElementWithProcessId( 0 );
        auto elemTest_ent  = meshTest->endElementWithProcessId( 0 );
        auto graphProjTest = spaceGraphProjTest->element();graphProjTest.zero();
        for ( ; elemTest_itt != elemTest_ent; ++elemTest_itt )
            {
                if ( mapBetweenMeshes_test.find( elemTest_itt->id() ) != mapBetweenMeshes_test.end() )
                    {
                        if ( mapBetweenMeshes_test.find( elemTest_itt->id() )->second.size()>0 )
                            {
                                auto element_dof1 = spaceGraphProjTest->dof()->getIndices( elemTest_itt->id() );
                                graphProjTest( element_dof1[0] ) = 1;
                            }
                    }
            }
        //auto exporterTest = boost::shared_ptr<Feel::Exporter<Mesh1Type> >( Feel::Exporter<Mesh1Type>::New( "ensight", "ExportQuadLocalizationTest" ) );
        auto exporterTest = Feel::Exporter<Mesh1Type>::New( "ensight", "ExportQuadLocalizationTest" );
        exporterTest->step( 0 )->setMesh( graphProjTest.mesh() );
        exporterTest->step( 0 )->add( "quadLocalizationProjTest", graphProjTest );
        exporterTest->save();

        typedef FunctionSpace<Mesh2Type, bases<Lagrange<0, Scalar,Discontinuous> > > space_trial_disc_type;
        auto spaceGraphProjTrial = space_trial_disc_type::New( meshTrial );
        auto elemTrial_itt  = meshTrial->beginElementWithProcessId( 0 );
        auto elemTrial_ent  = meshTrial->endElementWithProcessId( 0 );
        auto graphProjTrial = spaceGraphProjTrial->element();graphProjTrial.zero();
        for ( ; elemTrial_itt != elemTrial_ent; ++elemTrial_itt )
            {
                if ( mapBetweenMeshes_trial.find( elemTrial_itt->id() ) != mapBetweenMeshes_trial.end() )
                    {
                        if ( mapBetweenMeshes_trial.find( elemTrial_itt->id() )->second.size()>0 )
                            {
                                auto element_dof2 = spaceGraphProjTrial->dof()->getIndices( elemTrial_itt->id() );
                                graphProjTrial( element_dof2[0] ) = 1;
                            }
                    }
            }
        //auto exporterTrial = boost::shared_ptr<Feel::Exporter<Mesh2Type> >( Feel::Exporter<Mesh2Type>::New( "ensight", "ExportQuadLocalizationTrial" ) );
        auto exporterTrial = Feel::Exporter<Mesh2Type>::New( "ensight", "ExportQuadLocalizationTrial" );
        exporterTrial->step( 0 )->setMesh( graphProjTrial.mesh() );
        exporterTrial->step( 0 )->add( "quadLocalizationProjTrial", graphProjTrial );
        exporterTrial->save();

    }
#endif

public :



template<typename FE1,typename FE2,typename ElemContType>
struct bilinearformContext
{
    typedef Expr expression_type;

    // typedef on form (trial and test):
    typedef vf::detail::BilinearForm<FE1,FE2,ElemContType> FormType;
    // test
    typedef typename FormType::gm_1_type gm_formTest_type;
    typedef typename FormType::mesh_element_1_type geoelement_formTest_type;
    typedef typename gm_formTest_type::template Context<expression_type::context|vm::POINT,geoelement_formTest_type> gmc_formTest_type;
    typedef boost::shared_ptr<gmc_formTest_type> gmc_formTest_ptrtype;
    typedef typename gm_formTest_type::precompute_type pc_formTest_type;
    typedef typename gm_formTest_type::precompute_ptrtype pc_formTest_ptrtype;
    //typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTest_ptrtype> > map_gmc_formTest_type;
    // trial
    typedef typename FormType::gm_2_type gm_formTrial_type;
    typedef typename FormType::mesh_element_2_type geoelement_formTrial_type;
    typedef typename gm_formTrial_type::template Context<expression_type::context|vm::POINT,geoelement_formTrial_type> gmc_formTrial_type;
    typedef boost::shared_ptr<gmc_formTrial_type> gmc_formTrial_ptrtype;
    typedef typename gm_formTrial_type::precompute_type pc_formTrial_type;
    typedef typename gm_formTrial_type::precompute_ptrtype pc_formTrial_ptrtype;
    //typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTrial_ptrtype> > map_gmc_formTrial_type;

    typedef std::list<boost::tuple< std::vector<boost::tuple<size_type,size_type> >,
                                    gmc_ptrtype,
                                    gmc_formTest_ptrtype,
                                    gmc_formTrial_ptrtype > > return_loc_type;
    typedef std::list<return_loc_type> return_type;

};


template<typename FE1,typename FE2,typename ElemContType>
void
precompute(vf::detail::BilinearForm<FE1,FE2,ElemContType>const& __form)
{
    typedef typename bilinearformContext<FE1,FE2,ElemContType>::gmc_formTest_type gmc_formTest_type;
    typedef typename bilinearformContext<FE1,FE2,ElemContType>::gmc_formTest_ptrtype gmc_formTest_ptrtype;
    typedef typename bilinearformContext<FE1,FE2,ElemContType>::pc_formTest_type pc_formTest_type;
    typedef typename bilinearformContext<FE1,FE2,ElemContType>::pc_formTest_ptrtype pc_formTest_ptrtype;

    typedef typename bilinearformContext<FE1,FE2,ElemContType>::gmc_formTrial_type gmc_formTrial_type;
    typedef typename bilinearformContext<FE1,FE2,ElemContType>::gmc_formTrial_ptrtype gmc_formTrial_ptrtype;
    typedef typename bilinearformContext<FE1,FE2,ElemContType>::pc_formTrial_type pc_formTrial_type;
    typedef typename bilinearformContext<FE1,FE2,ElemContType>::pc_formTrial_ptrtype pc_formTrial_ptrtype;

    auto meshTrial = __form.trialSpace()->mesh();
    auto meshTest = __form.testSpace()->mesh();

    typename bilinearformContext<FE1,FE2,ElemContType>::return_type theres;

    auto res_it = this->result().begin();
    auto const res_en = this->result().end();
    for ( ; res_it != res_en ; ++res_it )
    {
        auto const idEltTest = res_it->template get<0>();
        auto const& map = res_it->template get<1>();
        auto map_it = map.begin();
        auto const map_en = map.end();
        for ( ; map_it != map_en ; ++map_it )
        {
            auto const idEltTrial = map_it->first;
            auto const& eltTrial = meshTrial->element( idEltTrial );
            auto const& eltTest = meshTest->element( idEltTest );

            auto const& ptRefTest = map_it->second.template get<1>();
            auto const& ptRefTrial = map_it->second.template get<2>();
            auto const& themapQuad = map_it->second.template get<0>();
            auto vec_gmcExpr = this->getUsableDataInFormContext( themapQuad,ptRefTest,ptRefTrial );
            auto gmcExpr_it = vec_gmcExpr.begin();
            auto const gmcExpr_en = vec_gmcExpr.end();
            typename bilinearformContext<FE1,FE2,ElemContType>::return_loc_type theresloc;
            for ( ; gmcExpr_it != gmcExpr_en ; ++gmcExpr_it )
            {
                pc_formTest_ptrtype geopcFormTest( new pc_formTest_type( __form.gm(), gmcExpr_it->template get<2>()/*__form.testSpace()->fe()->points()*/ ) );
                gmc_formTest_ptrtype gmcFormTest( new gmc_formTest_type( __form.gm(), eltTest /*__form.testSpace()->mesh()->element( 0 )*/, geopcFormTest ) );

                pc_formTrial_ptrtype geopcFormTrial( new pc_formTrial_type( __form.gmTrial(), gmcExpr_it->template get<3>() /* __form.trialSpace()->fe()->points()*/  ) );
                gmc_formTrial_ptrtype gmcFormTrial( new gmc_formTrial_type( __form.gmTrial(), eltTrial/*__form.trialSpace()->mesh()->element( 0 )*/, geopcFormTrial ) );

                theresloc.push_back(boost::make_tuple(gmcExpr_it->template get<0>(),gmcExpr_it->template get<1>(),gmcFormTest,gmcFormTrial));
            }
            theres.push_back(theresloc);
        }
    }

    M_precompute = theres;
    M_hasPrecompute=true;

}


template<typename FE,typename VectorType,typename ElemContType>
struct linearformContext
{
    typedef Expr expression_type;

    //typedef on form test :
    typedef vf::detail::LinearForm<FE,VectorType,ElemContType> FormType;
    typedef typename FormType::gm_type gm_form_type;
    typedef typename FormType::mesh_test_element_type geoelement_form_type;
    typedef typename gm_form_type::template Context<expression_type::context|vm::POINT,geoelement_form_type> gmc_form_type;
    typedef boost::shared_ptr<gmc_form_type> gmc_form_ptrtype;
    typedef typename gm_form_type::precompute_type pc_form_type;
    typedef typename gm_form_type::precompute_ptrtype pc_form_ptrtype;

    typedef std::list<boost::tuple< std::vector<boost::tuple<size_type,size_type> >,
                                    gmc_ptrtype,
                                    gmc_form_ptrtype > > return_loc_type;
    typedef std::list<return_loc_type> return_type;

};

template<typename FE,typename VectorType,typename ElemContType>
void
precompute(vf::detail::LinearForm<FE,VectorType,ElemContType> const& __form)
{

    typedef typename linearformContext<FE,VectorType,ElemContType>::gmc_form_type gmc_form_type;
    typedef typename linearformContext<FE,VectorType,ElemContType>::gmc_form_ptrtype gmc_form_ptrtype;
    typedef typename linearformContext<FE,VectorType,ElemContType>::pc_form_type pc_form_type;
    typedef typename linearformContext<FE,VectorType,ElemContType>::pc_form_ptrtype pc_form_ptrtype;

    auto meshTest = __form.testSpace()->mesh();

    typename linearformContext<FE,VectorType,ElemContType>::return_type theres;

    auto res_it = this->resultLinear().begin();
    auto const res_en = this->resultLinear().end();
    for ( ; res_it != res_en ; ++res_it )
    {
        auto const idEltTest = res_it->template get<0>();
        auto const& eltTest = meshTest->element( idEltTest );

        auto ptRefTest = res_it->template get<2>();
        auto themapQuad = res_it->template get<1>();

        auto vec_gmcExpr = this->getUsableDataInFormContext( themapQuad,ptRefTest );
        auto gmcExpr_it = vec_gmcExpr.begin();
        auto const gmcExpr_en = vec_gmcExpr.end();
        typename linearformContext<FE,VectorType,ElemContType>::return_loc_type theresloc;
        for ( ; gmcExpr_it != gmcExpr_en ; ++gmcExpr_it )
        {

            pc_form_ptrtype geopcForm( new pc_form_type( __form.gm(), gmcExpr_it->template get<2>() /*this->im().points()*/ ) );
            gmc_form_ptrtype gmcForm( new gmc_form_type( __form.gm(), eltTest, geopcForm ) );
            theresloc.push_back(boost::make_tuple(gmcExpr_it->template get<0>(),gmcExpr_it->template get<1>(),gmcForm));
        }
        theres.push_back(theresloc);

    }

    M_precompute = theres;
    M_hasPrecompute=true;

}


template<typename FE1,typename FE2,typename ElemContType>
typename bilinearformContext<FE1,FE2,ElemContType>::return_type const&
getPrecompute(vf::detail::BilinearForm<FE1,FE2,ElemContType>const& __form) const
{
    return boost::any_cast<typename bilinearformContext<FE1,FE2,ElemContType>::return_type const&>( M_precompute);
}

template<typename FE,typename VectorType,typename ElemContType>
typename linearformContext<FE,VectorType,ElemContType>::return_type const&
getPrecompute(vf::detail::LinearForm<FE,VectorType,ElemContType> const& __form) const
{
    return boost::any_cast<typename linearformContext<FE,VectorType,ElemContType>::return_type const&>( M_precompute);
}

private :

    std::list<range_iterator> M_listRange;

    element_iterator_type M_eltbegin;
    element_iterator_type M_eltend;
    mutable im_type M_im;
    QuadMapped<im_type> M_qm;
    permutation_points_type M_ppts;

    result_container_type M_resBilinear;

    result_container_linear_type M_resLinear;

    bool M_hasPrecompute;
    boost::any M_precompute;


}; // QuadPtLocalization




namespace detail
{
template <typename RangeType>
struct quadptlocrangetype
{
    typedef typename mpl::if_< boost::is_std_list<RangeType>,
                               mpl::identity<RangeType>,
                               mpl::identity<std::list<RangeType> > >::type::type::value_type type;
};
}


template<typename MeshTestType, typename MeshTrialType, typename IteratorRange,typename Im,typename Expr>
boost::shared_ptr<QuadPtLocalization<typename Feel::detail::quadptlocrangetype<IteratorRange>::type,Im,Expr> >
quadPtLocPtr( boost::shared_ptr<MeshTestType> meshTest, boost::shared_ptr<MeshTrialType> meshTrial, IteratorRange const& elts,Im const& im,Expr const& expr )
{
    typedef QuadPtLocalization<typename Feel::detail::quadptlocrangetype<IteratorRange>::type,Im,Expr> quadptloc_type;
    boost::shared_ptr<quadptloc_type> res(new quadptloc_type(elts) );
    res->update( meshTest,meshTrial );
    return res;
}

template<typename MeshTestType, typename IteratorRange,typename Im,typename Expr>
boost::shared_ptr<QuadPtLocalization<typename Feel::detail::quadptlocrangetype<IteratorRange>::type,Im,Expr> >
quadPtLocPtr( boost::shared_ptr<MeshTestType> meshTest, IteratorRange const& elts,Im const& im,Expr const& expr )
{
    typedef QuadPtLocalization<typename Feel::detail::quadptlocrangetype<IteratorRange>::type,Im,Expr> quadptloc_type;
    boost::shared_ptr<quadptloc_type> res(new quadptloc_type(elts) );
    res->update( meshTest );
    return res;
}



} // Feel
#endif

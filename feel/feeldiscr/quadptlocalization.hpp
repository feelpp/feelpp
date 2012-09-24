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
            mpl::identity<typename Im::template apply<geoelement_type::nDim, value_type, Simplex>::type >,
                mpl::identity<typename Im::template apply<geoelement_type::nDim, value_type, Hypercube>::type >
    >::type::type im_type;

    typedef typename im_type::face_quadrature_type im_face_type;

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


    QuadPtLocalization( IteratorRange const& elts, im_type const& /*__im*/ )
        :
        _M_eltbegin( elts.template get<1>() ),
        _M_eltend( elts.template get<2>() ),
        _M_im( )
    {}

    QuadPtLocalization( element_iterator_type  elts_it,
                        element_iterator_type  elts_en,
                        im_type const& /*__im*/ )
        :
        _M_eltbegin( elts_it ),
        _M_eltend( elts_en ),
        _M_im( )
    {}

    QuadPtLocalization( element_iterator_type  elts_it,
                        element_iterator_type  elts_en
                      )
        :
        _M_eltbegin( elts_it ),
        _M_eltend( elts_en ),
        _M_im( )
    {}

    /**
     * get the integration method
     */
    im_type const& im() const
    {
        return _M_im;
    }

    /**
     * get the integration method on face f
     */
    im_face_type  im( uint16_type f ) const
    {
        return _M_im.face( f );
    }

    /**
     * begin itrange
     */
    element_iterator_type beginElement() const
    {
        return _M_eltbegin;
    }

    /**
     * end itrange
     */
    element_iterator_type endElement() const
    {
        return _M_eltend;
    }

    /**
     * container for linear form
     */
    result_container_linear_type const & resultLinear() const
    {
        return _M_resLinear;
    }

    /**
     * container for bilinear form
     */
    result_container_type const & result() const
    {
        return _M_resBilinear;
    }

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
        //search element
        auto elt_it = this->beginElement();

        for ( size_type i=0; i<theIdElt; ++i ) ++elt_it;

        // get only usefull quad point and reoder
        uint16_type nContextPt = indexLocalToQuad.size();
        matrix_node_type quadPtsRef = this->im().points();
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
        //search element
        auto elt_it = this->beginElement();

        for ( size_type i=0; i<theIdElt; ++i ) ++elt_it;

        // get only usefull quad point and reoder
        uint16_type nContextPt = indexLocalToQuad.size();
        uint16_type __face_id_in_elt_0 = elt_it->pos_first();

        QuadMapped<im_type> qm;
        typedef typename QuadMapped<im_type>::permutation_type permutation_type;
        typename QuadMapped<im_type>::permutation_points_type ppts( qm( this->im() ) );
        auto __perm = elt_it->element( 0 ).permutation( __face_id_in_elt_0 );
        matrix_node_type quadPtsRef =  ppts[ __face_id_in_elt_0].find( __perm )->second;
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

        auto elt_it = this->beginElement();
        auto elt_en = this->endElement();


        QuadMapped<im_type> qm;
        typedef typename QuadMapped<im_type>::permutation_type permutation_type;
        typename QuadMapped<im_type>::permutation_points_type ppts( qm( this->im() ) );


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
                __geopc[__f][__p] = pc_ptrtype(  new pc_type( elt_it->element( 0 ).gm(), ppts[__f].find( __p )->second ) );
            }
        }


        uint16_type __face_id_in_elt_0 = elt_it->pos_first();

        gmc_ptrtype gmc( new gmc_type( elt_it->element( 0 ).gm(),
                                       elt_it->element( 0 ),
                                       __geopc,
                                       __face_id_in_elt_0 ) );

        auto meshTestLocalization = meshTest->tool_localization();
        meshTestLocalization->updateForUse();

        for ( size_type ide = 0; elt_it != elt_en; ++elt_it, ++ide )
        {
            __face_id_in_elt_0 = elt_it->pos_first();
            //if ( elt_it->isConnectedTo1()) std::cout << "\n AIEEEEEE!!!!!!!!!!!\n";

            gmc->update( elt_it->element( 0 ), __face_id_in_elt_0 );

            //std::cout << "\n quad gmc "<< gmc->xReal();
            for ( int q = 0; q <  gmc->nPoints(); ++ q )
            {
                size_type idq = gmc->nPoints()*ide+q;
                auto testAnalysis = meshTestLocalization->searchElement( gmc->xReal( q ) );
                auto testIdElt = testAnalysis.template get<1>();
                auto testNodeRef = testAnalysis.template get<2>();
                testEltToPtsQuad[testIdElt].push_back( boost::make_tuple( idq,q,testNodeRef ) );
            }
        }

    }

    //--------------------------------------------------------------------------------------//

    template <typename Mesh1Type>
    void
    localization( mpl::int_<MESH_ELEMENTS> /**/,
                  boost::shared_ptr<Mesh1Type> meshTest,
                  std::vector<std::list<boost::tuple< size_type,size_type,node_type> > > & testEltToPtsQuad )
    {

        auto elt_it = this->beginElement();
        auto elt_en = this->endElement();

        pc_ptrtype geopc( new pc_type( elt_it->gm(), this->im().points() ) );
        gmc_ptrtype gmc( new gmc_type( elt_it->gm(),*elt_it, geopc ) );

        auto meshTestLocalization = meshTest->tool_localization();
        meshTestLocalization->updateForUse();


        //auto nQuadPtsInElt = this->im().nPoints();
        //auto nElts = std::distance(elt_it,elt_en);
        //auto nQuadPts = nElts*nQuadPtsInElt;

        //auto nEltTest= meshTest->numElements();

        //for ( size_type theIdEltTest = 0; theIdEltTest < nEltTest ; ++theIdEltTest ) {testEltToPtsQuad[theIdEltTest].clear();}

        for ( size_type ide = 0; elt_it != elt_en; ++elt_it, ++ide )
        {
            gmc->update( *elt_it );

            for ( int q = 0; q <  gmc->nPoints(); ++ q )
            {
                // cpt of quad pt
                size_type idq = gmc->nPoints()*ide+q;
                // search in test mesh
                auto testAnalysis = meshTestLocalization->searchElement( gmc->xReal( q ) );
                auto testIdElt = testAnalysis.template get<1>();
                auto testNodeRef = testAnalysis.template get<2>();
                testEltToPtsQuad[testIdElt].push_back( boost::make_tuple( idq,q,testNodeRef ) );

            }
        }


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

        auto elt_it = this->beginElement();
        auto elt_en = this->endElement();

        QuadMapped<im_type> qm;
        typedef typename QuadMapped<im_type>::permutation_type permutation_type;
        typename QuadMapped<im_type>::permutation_points_type ppts( qm( this->im() ) );

        std::vector<std::map<permutation_type, pc_ptrtype> > __geopc( this->im().nFaces() );

        for ( uint16_type __f = 0; __f < this->im().nFaces(); ++__f )
        {
            for ( permutation_type __p( permutation_type::IDENTITY );
                    __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
            {
                //FEELPP_ASSERT( ppts[__f].find(__p)->second.size2() != 0 ).warn( "invalid quadrature type" );
                __geopc[__f][__p] = pc_ptrtype(  new pc_type( elt_it->element( 0 ).gm(), ppts[__f].find( __p )->second ) );
            }
        }

        uint16_type __face_id_in_elt_0 = elt_it->pos_first();

        gmc_ptrtype gmc( new gmc_type( elt_it->element( 0 ).gm(),
                                       elt_it->element( 0 ),
                                       __geopc,
                                       __face_id_in_elt_0 ) );


        auto meshTrialLocalization = meshTrial->tool_localization();
        meshTrialLocalization->updateForUse();

        auto meshTestLocalization = meshTest->tool_localization();
        meshTestLocalization->updateForUse();

        matrix_node_type ptsReal( elt_it->vertices().size1(), 1 );
        size_type trialIdElt = 0, testIdElt=0;
        node_type trialNodeRef,testNodeRef;

        bool quadMeshIsSameThatTrialMesh=false,quadMeshIsSameThatTestMesh=false;
        if ( dynamic_cast<void*>( const_cast<MeshBase*>( elt_it->mesh() ) ) == dynamic_cast<void*>( meshTrial.get() ) )
            quadMeshIsSameThatTrialMesh=true;
        if ( dynamic_cast<void*>( const_cast<MeshBase*>( elt_it->mesh() ) ) == dynamic_cast<void*>( meshTest.get() ) )
            quadMeshIsSameThatTestMesh=true;

#if FEELPP_EXPORT_QUADLOCALIZATION
        std::map<size_type,std::list<size_type> > mapBetweenMeshes_test;
        std::map<size_type,std::list<size_type> > mapBetweenMeshes_trial;
#endif

        for ( size_type ide = 0; elt_it != elt_en; ++elt_it, ++ide )
        {
            __face_id_in_elt_0 = elt_it->pos_first();
            //if ( elt_it->isConnectedTo1()) std::cout << "\n AIEEEEEE!!!!!!!!!!!\n";

            gmc->update( elt_it->element( 0 ), __face_id_in_elt_0 );

            for ( int q = 0; q <  gmc->nPoints(); ++ q )
            {
                // cpt of quad pt
                size_type idq = gmc->nPoints()*ide+q;

                // search in trial mesh
                ublas::column(ptsReal,0 ) = gmc->xReal( q );
                if (!quadMeshIsSameThatTrialMesh)
                    {
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
                    //std::cout << " testIdElt " << testIdElt << " trialIdElt " << trialIdElt << std::endl;
#if FEELPP_EXPORT_QUADLOCALIZATION
                    mapBetweenMeshes_test[testIdElt].push_back( trialIdElt );
                    mapBetweenMeshes_trial[trialIdElt].push_back( testIdElt );
#endif
                }
            } // for ( int q = 0; q <  gmc->nPoints(); ++ q )
        } // end for( size_type ide ... )

#if FEELPP_EXPORT_QUADLOCALIZATION
        this->exportQuadLocalization(meshTest,meshTrial,mapBetweenMeshes_test,mapBetweenMeshes_trial);
#endif
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

        auto elt_it = this->beginElement();
        auto elt_en = this->endElement();

        pc_ptrtype geopc( new pc_type( elt_it->gm(), this->im().points() ) );
        gmc_ptrtype gmc( new gmc_type( elt_it->gm(),*elt_it, geopc ) );

        auto meshTrialLocalization = meshTrial->tool_localization();
        meshTrialLocalization->updateForUse();

        auto meshTestLocalization = meshTest->tool_localization();
        meshTestLocalization->updateForUse();

#if FEELPP_EXPORT_QUADLOCALIZATION
        std::map<size_type,std::list<size_type> > mapBetweenMeshes_test;
        std::map<size_type,std::list<size_type> > mapBetweenMeshes_trial;
#endif


        //auto nQuadPtsInElt = this->im().nPoints();
        //auto nElts = std::distance(elt_it,elt_en);
        //auto nQuadPts = nElts*nQuadPtsInElt;
        // auto nEltTrial= meshTrial->numElements();
        //auto nEltTest= meshTest->numElements();

        matrix_node_type ptsReal( elt_it->vertices().size1(), 1 );
        size_type trialIdElt = 0, testIdElt=0;
        node_type trialNodeRef,testNodeRef;

        bool quadMeshIsSameThatTrialMesh=false,quadMeshIsSameThatTestMesh=false;
        if ( dynamic_cast<void*>( const_cast<MeshBase*>( elt_it->mesh() ) ) == dynamic_cast<void*>( meshTrial.get() ) )
            quadMeshIsSameThatTrialMesh=true;
        if ( dynamic_cast<void*>( const_cast<MeshBase*>( elt_it->mesh() ) ) == dynamic_cast<void*>( meshTest.get() ) )
            quadMeshIsSameThatTestMesh=true;

        for ( size_type ide = 0; elt_it != elt_en; ++elt_it, ++ide )
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
        } // end for( size_type ide ... )

#if FEELPP_EXPORT_QUADLOCALIZATION
        this->exportQuadLocalization(meshTest,meshTrial,mapBetweenMeshes_test,mapBetweenMeshes_trial);
#endif

        //std::cout << "[QuadPtLocalization] : localization<MESH_ELEMENTS>(bilinear form finish" << std::endl;
    }

    //--------------------------------------------------------------------------------------//

    template <typename Mesh1Type>
    void
    update( boost::shared_ptr<Mesh1Type> meshTest )
    {
        typedef Mesh1Type mesh_test_type;

        //clean result
        _M_resLinear.clear();

        // localize quad pts on meshTest -> storage in testEltToPtsQuad
        auto nEltTest= meshTest->numElements();
        std::vector<std::list<boost::tuple< size_type,size_type,node_type> > > testEltToPtsQuad( nEltTest );
        this->localization( mpl::int_<iDim>(),meshTest,testEltToPtsQuad );

        //build a efficent container : _M_resLinear
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

                _M_resLinear.push_back( boost::make_tuple( theIdEltTest, indexLocalToQuad, ptsRefTest ) );
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
        _M_resBilinear.clear();

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

        //build a efficent container : _M_resBilinear
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
            _M_resBilinear.push_back( boost::make_tuple( theIdEltTest,mapiIdTrial2qAndPtRef ) );
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






    element_iterator_type _M_eltbegin;
    element_iterator_type _M_eltend;
    mutable im_type _M_im;

    result_container_type _M_resBilinear;

    result_container_linear_type _M_resLinear;

}; // QuadPtLocalization


} // Feel
#endif

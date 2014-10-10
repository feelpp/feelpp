/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel++ library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date     : Mon Feb 24 16:03:52 2014

   Copyright (C) 2014 Feel++ Consortium

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
#ifndef FEELPP_MESH_LOCALIZATION_IMPL_HPP
#define FEELPP_MESH_LOCALIZATION_IMPL_HPP 1

namespace Feel {


template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::Localization::init()
{
    auto mesh = M_mesh.lock();
    if ( !mesh ) return;

    DLOG_IF( WARNING, this->isInit() == false ) << "You have already initialized the tool of localization\n";


    //clear data
    M_geoGlob_Elts.clear();
    M_kd_tree->clear();

    typename self_type::element_iterator el_it;
    typename self_type::element_iterator el_en;
    boost::tie( boost::tuples::ignore, el_it, el_en ) = Feel::elements( *mesh );

    for ( ; el_it != el_en; ++el_it )
    {
        for ( int i=0; i<el_it->nPoints(); ++i )
        {
            if ( boost::get<1>( M_geoGlob_Elts[el_it->point( i ).id()] ).size()==0 )
            {
                boost::get<0>( M_geoGlob_Elts[el_it->point( i ).id()] ) = el_it->point( i ).node();
                M_kd_tree->addPoint( el_it->point( i ).node(),el_it->point( i ).id() );
            }

            boost::get<1>( M_geoGlob_Elts[el_it->point( i ).id()] ).push_back( el_it->id() );
        }
    }

    this->computeBarycenter();

    M_isInit=true;
    M_isInitBoundaryFaces=false;

}

template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::Localization::initBoundaryFaces()
{
    auto mesh = M_mesh.lock();
    if ( !mesh ) return;

    DLOG_IF( WARNING, this->isInitBoundaryFaces() == false ) << "You have already initialized the tool of localization\n";

    //clear data
    M_geoGlob_Elts.clear();
    M_kd_tree->clear()
;
    typename self_type::location_face_iterator face_it;
    typename self_type::location_face_iterator face_en;
    boost::tie( boost::tuples::ignore, face_it, face_en ) = Feel::boundaryfaces( mesh );

    for ( ; face_it != face_en; ++face_it )
    {
        for ( int i=0; i<face_it->nPoints(); ++i )
        {
            if ( face_it->isConnectedTo0() )
                {
                    if ( boost::get<1>( M_geoGlob_Elts[face_it->point( i ).id()] ).size()==0 )
                        {
                            boost::get<0>( M_geoGlob_Elts[face_it->point( i ).id()] ) = face_it->point( i ).node();
                            M_kd_tree->addPoint( face_it->point( i ).node(),face_it->point( i ).id() );
                        }
                    boost::get<1>( M_geoGlob_Elts[face_it->point( i ).id()] ).push_back( face_it->element( 0 ).id() );
                }
        }
    }

    this->computeBarycenter();

    M_isInitBoundaryFaces=true;
    M_isInit=false;

}



template<typename Shape, typename T, int Tag>
boost::tuple<bool,typename Mesh<Shape, T, Tag>::node_type,double>
Mesh<Shape, T, Tag>::Localization::isIn( size_type _id, const node_type & _pt ) const
{
    bool isin=false;
    double dmin;
    node_type x_ref;

    auto mesh = M_mesh.lock();
    //get element with the id
    auto const& elt = mesh->element( _id );

    if ( elt.isOnBoundary() )
        {
            // get inverse geometric transformation
            gmc_inverse_type gic( mesh->gm(), elt, mesh->worldComm().subWorldCommSeq() );
            //apply the inverse geometric transformation for the point p
            gic.setXReal( _pt);
            x_ref=gic.xRef();
            // the point is in the reference element ?
            boost::tie( isin, dmin ) = M_refelem.isIn( gic.xRef() );
        }
    else
        {
            // get inverse geometric transformation
            gmc1_inverse_type gic( mesh->gm1(), elt, mpl::int_<1>(), mesh->worldComm().subWorldCommSeq() );
            //apply the inverse geometric transformation for the point p
            gic.setXReal( _pt);
            x_ref=gic.xRef();
            // the point is in the reference element ?
            boost::tie( isin, dmin ) = M_refelem1.isIn( gic.xRef() );
        }

    return boost::make_tuple(isin,x_ref,dmin);
}

template<typename Shape, typename T, int Tag>
boost::tuple<uint16_type,std::vector<bool> >
Mesh<Shape, T, Tag>::Localization::isIn( std::vector<size_type> _ids, const node_type & _pt )
{
    typedef typename self_type::gm_type::reference_convex_type ref_convex_type;
    typedef typename self_type::gm1_type::reference_convex_type ref_convex1_type;

    uint16_type nbId = _ids.size();
    std::vector<bool> isin( _ids.size(),false );
    bool isin2=false;
    double dmin;
    node_type __x_ref;

    auto mesh = M_mesh.lock();
    uint16_type nbIsIn=0;

    for ( uint16_type i = 0; i< nbId ; ++i )
    {
        //get element with the id
        auto const& elt = mesh->element( _ids[i] );

        if ( elt.isOnBoundary() )
            {
                // get inverse geometric transformation
                gmc_inverse_type gic( mesh->gm(), elt );
                //apply the inverse geometric transformation for the point p
                gic.setXReal( _pt);
                __x_ref=gic.xRef();
                // the point is in the reference element ?
                boost::tie( isin2, dmin ) = M_refelem.isIn( gic.xRef() );
                isin[i] = isin2;
            }
        else
            {
                // get inverse geometric transformation
                gmc1_inverse_type gic( mesh->gm1(), elt, mpl::int_<1>() );
                //apply the inverse geometric transformation for the point p
                gic.setXReal( _pt);
                __x_ref=gic.xRef();
                // the point is in the reference element ?
                boost::tie( isin2, dmin ) = M_refelem1.isIn( gic.xRef() );
                isin[i] = isin2;
            }
        if (isin[i]) ++nbIsIn;
    }

    return boost::make_tuple( nbIsIn,isin );
}


template<typename Shape, typename T, int Tag>
boost::tuple<bool, size_type, typename Mesh<Shape, T, Tag>::node_type>
Mesh<Shape, T, Tag>::Localization::searchElement( const node_type & p )
{

    DCHECK( this->isInit() || this->isInitBoundaryFaces() ) << "You don't have initialized the tool of localization\n";

    auto mesh = M_mesh.lock();
    bool isin=false;double dmin=0;
    node_type x_ref;
    size_type idEltFound = 0;//mesh->beginElementWithId(mesh->worldComm().localRank())->id();

    std::list< std::pair<size_type, uint> > ListTri;
    this->searchInKdTree(p,ListTri);

    //research the element which contains the point p
    auto itLT=ListTri.begin();
    auto const itLT_end=ListTri.end();

    DCHECK( std::distance( itLT,itLT_end )>0 ) << "problem in localization : listTri is empty\n";

    while ( itLT != itLT_end && !isin  )
    {
        //get element with the id
        //elt = M_mesh->element( itLT->first );

        // search point in this elt
        boost::tie(isin,x_ref,dmin) = this->isIn(itLT->first,p);
        // if not inside, continue the research with an other element
        if (!isin) ++itLT;
        else idEltFound=itLT->first;
    }


    if( !isin && this->doExtrapolation() )
        {
            // first elt
            idEltFound = ListTri.begin()->first;
            auto const& eltUsedForExtrapolation = mesh->element(idEltFound);
            DVLOG(1) << "localisation tool use extrapolation for the point" << p
                     << " with elt.id() " << eltUsedForExtrapolation.id()
                     << " and elt.G() " << eltUsedForExtrapolation.G()
                     << "\n";
            gmc_inverse_type gic( mesh->gm(), eltUsedForExtrapolation, mesh->worldComm().subWorldCommSeq() );
            //apply the inverse geometric transformation for the point p
            gic.setXReal( p);
            boost::tie(isin,idEltFound,x_ref) = boost::make_tuple(true,eltUsedForExtrapolation.id(),gic.xRef());
        }

    return boost::make_tuple( isin, idEltFound, x_ref);

}

template<typename Shape, typename T, int Tag>
boost::tuple<std::vector<bool>, size_type>
Mesh<Shape, T, Tag>::Localization::run_analysis( const matrix_node_type & m,
                                            const size_type & eltHypothetical )
{
    // if no init then init with all geo point
    if ( !( this->isInit() || this->isInitBoundaryFaces() ) )
        this->init();

    bool find_x;
    size_type cv_id=eltHypothetical;
    node_type x_ref;double dmin=0;
    std::vector<bool> hasFindPts(m.size2(),false);

    M_resultAnalysis.clear();

    bool doExtrapolationAtStart = this->doExtrapolation();
    auto nPtMaxNearNeighborAtStart = this->kdtree()->nPtMaxNearNeighbor();


    // first step : no extrapolation
    if ( doExtrapolationAtStart ) this->setExtrapolation( false );

    // first currentEltHypothetical
    auto currentEltHypothetical = eltHypothetical;//this->mesh()->beginElement()->id();//eltHypothetical;
    for ( size_type i=0; i< m.size2(); ++i )
        {
            bool testHypothetical_find = false;

            if ( eltHypothetical!=invalid_size_type_value )
                {
                    boost::tie( testHypothetical_find,x_ref,dmin ) = this->isIn( currentEltHypothetical,ublas::column( m, i ) );
                }
            if ( testHypothetical_find )
                {
                    cv_id = currentEltHypothetical;
                    M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
                    hasFindPts[i]=true;
                    DVLOG(1) << "localisation tool (hypthetical elt) : has found pt " << ublas::column( m, i )
                             << "in elt id " << cv_id
                             << "with G() " << M_mesh.lock()->element(cv_id).G()
                             << "\n";
                }
            else // search kdtree
                {
                    // kdtree call
                    boost::tie( find_x, cv_id, x_ref ) = this->searchElement(ublas::column( m, i ));
                    // traitement
                    if (find_x) // if find : OK
                        {
                            M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
                            currentEltHypothetical = cv_id;
                            hasFindPts[i]=true;
                            DVLOG(1)  << "localisation tool (first pass) : has found pt " << ublas::column( m, i )
                                      << " in elt id " << cv_id
                                      << " with G() " << M_mesh.lock()->element(cv_id).G()
                                      << " and ref point " << x_ref
                                      << "\n";
                        }
                    else// if (false) // try an other method (no efficient but maybe a solution)
                        {
                            // search in all element
                            this->kdtree()->nbNearNeighbor(2*nPtMaxNearNeighborAtStart /*15*/ /*this->mesh()->numElements()*/);

                            boost::tie( find_x, cv_id, x_ref ) = this->searchElement(ublas::column( m, i ));
                            //revert parameter
                            this->kdtree()->nbNearNeighbor(nPtMaxNearNeighborAtStart);
                            // if find : OK (but strange!)
                            if (find_x)
                                {
                                    M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
                                    currentEltHypothetical = cv_id;
                                    hasFindPts[i]=true;
                                    DVLOG(1)  << "localisation tool (second pass) : has found pt " << ublas::column( m, i )
                                              << " in elt id " << cv_id
                                              << " with G() " << M_mesh.lock()->element(cv_id).G()
                                              << " and ref point " << x_ref
                                              << "\n";
                                }
                            else if (doExtrapolationAtStart)
                                {
                                    this->setExtrapolation(true);

                                    boost::tie( find_x, cv_id, x_ref ) = this->searchElement(ublas::column( m, i ));

                                    CHECK( find_x ) << "localisation tool : invalid extrapolation \n";

                                    M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
                                    currentEltHypothetical = cv_id;
                                    hasFindPts[i]=true;

                                    this->setExtrapolation(false);
                                }
                        }
                } // search kdtree

            DLOG_IF(WARNING, !hasFindPts[i]) << "the localisation tool fails to find the point " << ublas::column( m, i ) << "\n";

        } // for (size_type i=0;i< m.size2();++i)

    //revert parameter
    this->setExtrapolation(doExtrapolationAtStart);

    return boost::make_tuple(hasFindPts,cv_id);

} //run_analysis



template<typename Shape, typename T, int Tag>
boost::tuple<bool, std::list<boost::tuple<size_type, typename Mesh<Shape, T, Tag>::node_type> > >
Mesh<Shape, T, Tag>::Localization::searchElements( const node_type & p )
{

    DCHECK( this->isInit() || this->isInitBoundaryFaces() ) << "You don't have initialized the tool of localization\n";

    //this->kdtree()->nbNearNeighbor(this->mesh()->numElements());

    std::list< std::pair<size_type, uint> > ListTri;
    searchInKdTree( p,ListTri );

    typename self_type::element_type elt;
    typename self_type::gm_type::reference_convex_type refelem;
    typename self_type::gm1_type::reference_convex_type refelem1;

    bool isin=false;
    double dmin;
    node_type x_ref;

    //research the element which contains the point p
    auto itLT=ListTri.begin();
    auto itLT_end=ListTri.end();

#if !defined( NDEBUG )
    //if(std::distance(itLT,itLT_end)==0) std::cout<<"\nListTri vide\n";
    FEELPP_ASSERT( std::distance( itLT,itLT_end )>0 ).error( " problem in list localization : is empty" );
#endif

    std::list<boost::tuple<size_type,node_type> > newlistelts;
    newlistelts.clear();
    bool find = false;
    bool finishSearch=false;
    while ( itLT != itLT_end && !finishSearch /*&& !isin*/  )
    {
#if 0
        //get element with the id
        elt= M_mesh->element( itLT->first );

        if ( elt.isOnBoundary() )
        {
            // get inverse geometric transformation
            typename self_type::Inverse::gic_type gic( M_mesh->gm(), elt );

            //apply the inverse geometric transformation for the point p
            gic.setXReal( p );
            __x_ref=gic.xRef();

            // the point is in the reference element ?
            boost::tie( isin, dmin ) = refelem.isIn( gic.xRef() );
        }

        else
        {
            // get inverse geometric transformation
            typename self_type::Inverse::gic1_type gic( M_mesh->gm1(), elt,mpl::int_<1>() );

            //apply the inverse geometric transformation for the point p
            gic.setXReal( p );
            __x_ref=gic.xRef();

            // the point is in the reference element ?
            boost::tie( isin, dmin ) = refelem1.isIn( gic.xRef() );
            //std::cout << "gic.xRef()" << gic.xRef() << std::endl;
        }
#else
        boost::tie(isin,x_ref, dmin) = this->isIn(itLT->first,p);
#endif
        if ( isin )
        {
            newlistelts.push_back( boost::make_tuple( itLT->first,x_ref ) );
            find = true;
            // if not on boundary -> finish for this point
            if (dmin>1e-7) finishSearch=true;
        }

        //if (find) std::cout << elt.G() << std::endl;

        //if not inside, continue the research with an other element
        //if (!isin) ++itLT;
        ++itLT;
    }

    if ( !find ) std::cout << "\n WARNING EXTRAPOLATION IN SEARCHELEMENTS!!!"<<std::endl;

    if ( find )
        return boost::make_tuple( true,newlistelts );

    else if ( !find && !M_doExtrapolation )
        return boost::make_tuple( false,newlistelts );

    else
    {
        auto mesh = M_mesh.lock();
        //std::cout << "\n WARNING EXTRAPOLATION \n";
        itLT=ListTri.begin();
        elt= mesh->element( itLT->first );
        typename self_type::Inverse::gic_type gic( mesh->gm(), elt );
        //apply the inverse geometric transformation for the point p
        //gic.setXReal(boost::get<0>(*ptsNN.begin()));
        gic.setXReal( p );
        x_ref=gic.xRef();
        //return boost::make_tuple( true, itLT->first, __x_ref);
        newlistelts.push_back( boost::make_tuple( itLT->first,x_ref ) );
        find = true;
        return boost::make_tuple( true,newlistelts );
    }
} // searchElements


template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::Localization::searchInKdTree( const node_type & p,
        std::list< std::pair<size_type, uint> > & ListTri )
{
    //search for nearest points
    M_kd_tree->search( p );

    //get the results of research
    typename KDTree::points_search_type ptsNN = M_kd_tree->pointsNearNeighbor();

    typename KDTree::points_search_const_iterator itNN = ptsNN.begin();
    typename KDTree::points_search_const_iterator itNN_end = ptsNN.end();

#if !defined( NDEBUG )
    FEELPP_ASSERT( std::distance( itNN,itNN_end )>0 ).error( "none Near Neighbor Points are find" );
#endif

    //iterator on a l(ist index element
    typename std::list<size_type>::iterator itL;
    typename std::list<size_type>::iterator itL_end;

    //ListTri will contain the indices of elements (size_type)
    //and the number of occurence(uint)
    //std::list< std::pair<size_type, uint> > ListTri;
    std::list< std::pair<size_type, uint> >::iterator itLT;
    std::list< std::pair<size_type, uint> >::iterator itLT_end;

    //create of ListTri : sort largest to smallest occurrences
    //In case of equality : if the point is closer than another then it will be before
    //                      if it is in the same point then  the lowest index will be before
    for ( ; itNN != itNN_end; ++itNN )
    {
        itL= boost::get<1>( M_geoGlob_Elts[boost::get<3>( *itNN )] ).begin();
        itL_end= boost::get<1>( M_geoGlob_Elts[boost::get<3>( *itNN )] ).end();

        for ( ; itL != itL_end; ++itL )
        {
            itLT=ListTri.begin();
            itLT_end=ListTri.end();
            bool find=false;

            while ( itLT != itLT_end && !find )
            {
                if ( itLT->first == *itL ) find=true;

                else ++itLT;
            }

            if ( find )
            {
                uint nb=itLT->second+1;
                size_type numEl=itLT->first;
                ListTri.remove( *itLT );
                itLT=ListTri.begin();
                itLT_end=ListTri.end();
                bool find=false;

                while ( itLT != itLT_end && !find )
                {
                    if ( itLT->second < nb ) find=true;

                    else ++itLT;
                }

                ListTri.insert( itLT,std::make_pair( numEl,nb ) );
            }

            else ListTri.push_back( std::make_pair( *itL,1 ) );
        }
    }

}

template<typename Shape, typename T, int Tag>
boost::tuple<bool,typename Mesh<Shape, T, Tag>::node_type,double>
Mesh<Shape, T, Tag>::Localization::isIn( size_type _id,
                                         const node_type & _pt,
                                         const matrix_node_type & setPoints,
                                         mpl::int_<1> /**/ ) const
{
    bool isin=true; // warning : start with true
    double dmin=0.;
    node_type x_ref;
    auto mesh = M_mesh.lock();
    //get element with the id
    auto const& elt= mesh->element( _id );
    auto const& eltG = elt.G();

    // check conformity between setPoints (given) and eltG (localize)
    std::vector<bool> find( setPoints.size2(),false );
    for ( size_type i=0; i< setPoints.size2() && isin; ++i )
    {
        auto const& thePt = ublas::column( setPoints,i );

        for ( size_type j=0; j<eltG.size2(); ++j )
        {
            auto ptjeltG = ublas::column( eltG,j );

            if ( ptjeltG.size()==1 )
            {
                if ( std::abs( thePt( 0 )-ptjeltG( 0 ) )<1e-8 )
                    find[i]=true;
            }
            else if ( ptjeltG.size()==2 )
            {
                if ( std::abs( thePt( 0 )-ptjeltG( 0 ) )<1e-8 &&
                     std::abs( thePt( 1 )-ptjeltG( 1 ) )<1e-8 )
                    find[i]=true;
            }
            else if ( ptjeltG.size()==3 )
            {
                if ( std::abs( thePt( 0 )-ptjeltG( 0 ) )<1e-8 &&
                     std::abs( thePt( 1 )-ptjeltG( 1 ) )<1e-8 &&
                     std::abs( thePt( 2 )-ptjeltG( 2 ) )<1e-8 )
                    find[i]=true;
            }
        }
        // up checking
        isin &= find[i];
    }

    // if find -> get ref point and check
    if ( isin )
    {
        bool isin2=false;
        boost::tie(isin2,x_ref,dmin) = this->isIn(_id,_pt);
        LOG_IF(ERROR, !isin2) << "Mesh::Localization::isIn<Conformal> : check fail -> maybe x_ref is not correct";
    }

    return boost::make_tuple(isin,x_ref,dmin);
}



template<typename Shape, typename T, int Tag>
boost::tuple<bool, size_type, typename Mesh<Shape, T, Tag>::node_type>
Mesh<Shape, T, Tag>::Localization::searchElement( const node_type & p,
        const matrix_node_type & setPoints,
        mpl::int_<1> /**/ )
{
    typename self_type::element_type elt;
    typename self_type::gm_type::reference_convex_type refelem;
    typename self_type::gm1_type::reference_convex_type refelem1;

    bool isin=false,isin2=false;
    double dmin=0.;
    node_type x_ref;
    auto mesh=this->mesh().lock();
    size_type idEltFound = mesh->beginElementWithId(mesh->worldComm().localRank())->id();

    std::list< std::pair<size_type, uint> > ListTri;
    searchInKdTree( p,ListTri );

    auto itLT=ListTri.begin();
    auto itLT_end=ListTri.end();

#if !defined( NDEBUG )
    //if(std::distance(itLT,itLT_end)==0) std::cout<<"\nListTri vide\n";
    FEELPP_ASSERT( std::distance( itLT,itLT_end )>0 ).error( " problem in list localization : is empty" );
#endif


    //research the element which contains the point p
    while ( itLT != itLT_end && !isin  )
    {
        boost::tie(isin,x_ref,dmin) = this->isIn(itLT->first,p,setPoints,mpl::int_<1>());

        //if not inside, continue the research with an other element
        if ( !isin ) ++itLT;
        else idEltFound=itLT->first;
    } //while ( itLT != itLT_end && !isin  )

    if (!isin)
        {
            if( this->doExtrapolation() )
                {
                    LOG(WARNING) << "WARNING EXTRAPOLATION for the point" << p;
                    //std::cout << "W";
                    auto const& eltUsedForExtrapolation = mesh->element(ListTri.begin()->first);
                    gmc_inverse_type gic( mesh->gm(), eltUsedForExtrapolation, mesh->worldComm().subWorldCommSeq() );
                    //apply the inverse geometric transformation for the point p
                    gic.setXReal( p);
                    boost::tie(isin,idEltFound,x_ref) = boost::make_tuple(true,eltUsedForExtrapolation.id(),gic.xRef());
                }
            else
                {
                    idEltFound = mesh->beginElementWithId(mesh->worldComm().localRank())->id();
                    isin = false;
                    //x_ref=?
                }
        }

    return boost::make_tuple( isin, idEltFound, x_ref);

} //searchElement


template<typename Shape, typename T, int Tag>
boost::tuple<std::vector<bool>, size_type>
Mesh<Shape, T, Tag>::Localization::run_analysis( const matrix_node_type & m,
        const size_type & eltHypothetical,
        const matrix_node_type & setPoints,
        mpl::int_<1> /**/ )
{

    DCHECK( this->isInit() || this->isInitBoundaryFaces() ) << "You don't have initialized the tool of localization\n";

    bool find_x=false;
    size_type cv_id=eltHypothetical;
    node_type x_ref;
    double dmin;
    std::vector<bool> hasFindPts(m.size2(),false);

    M_resultAnalysis.clear();
    auto currentEltHypothetical = eltHypothetical;
    for ( size_type i=0; i< m.size2(); ++i )
    {

        bool testHypothetical_find = false;
        if ( eltHypothetical!=invalid_size_type_value )
        {
            boost::tie( testHypothetical_find,x_ref,dmin ) = this->isIn( currentEltHypothetical,ublas::column( m, i ),setPoints,mpl::int_<1>() );
        }
        if ( testHypothetical_find )
        {
            cv_id = currentEltHypothetical;
            M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
            hasFindPts[i]=true;
        }
        else
        {
            boost::tie( find_x, cv_id, x_ref ) = this->searchElement( ublas::column( m, i ),setPoints,mpl::int_<1>() );

            if ( find_x )
            {
                M_resultAnalysis[cv_id].push_back( boost::make_tuple( i,x_ref ) );
                hasFindPts[i]=true;
                currentEltHypothetical = cv_id;
            }
        }
        //else std::cout<<"\nNew Probleme Localization\n" << std::endl;
    }

    return boost::make_tuple(hasFindPts,cv_id);

} // run_analysis


template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::Localization::computeBarycenter()
{
    M_barycenter = this->computeBarycenter(mpl::int_<nRealDim>());
}

template<typename Shape, typename T, int Tag>
typename Mesh<Shape, T, Tag>::node_type
Mesh<Shape, T, Tag>::Localization::computeBarycenter(mpl::int_<1> /**/) const
{
    node_type res(1);
    res(0)=0;
    if ( this->kdtree()->nPoints()==0 ) return res;
    for (size_type i = 0 ; i<this->kdtree()->nPoints() ; ++i)
        {
            auto const& pt = this->kdtree()->points()[i].template get<0>();
            res(0)+=pt(0);
        }
    res(0)/=this->kdtree()->nPoints();
    return res;
}
template<typename Shape, typename T, int Tag>
typename Mesh<Shape, T, Tag>::node_type
Mesh<Shape, T, Tag>::Localization::computeBarycenter(mpl::int_<2> /**/) const
{
    node_type res(2);
    res(0)=0;res(1)=0;
    if ( this->kdtree()->nPoints()==0 )
        return res;
    for (size_type i = 0 ; i<this->kdtree()->nPoints() ; ++i)
        {
            auto const& pt = this->kdtree()->points()[i].template get<0>();
            res(0)+=pt(0);res(1)+=pt(1);
        }
    res(0)/=this->kdtree()->nPoints();res(1)/=this->kdtree()->nPoints();
    return res;
}
template<typename Shape, typename T, int Tag>
typename Mesh<Shape, T, Tag>::node_type
Mesh<Shape, T, Tag>::Localization::computeBarycenter(mpl::int_<3> /**/) const
{
    node_type res(3);
    res(0)=0;res(1)=0;res(2)=0;
    if ( this->kdtree()->nPoints()==0 )
        return res;
    for (size_type i = 0 ; i<this->kdtree()->nPoints() ; ++i)
        {
            auto const& pt = this->kdtree()->points()[i].template get<0>();
            res(0)+=pt(0);res(1)+=pt(1);res(2)+=pt(2);
        }
    res(0)/=this->kdtree()->nPoints();res(1)/=this->kdtree()->nPoints();res(2)/=this->kdtree()->nPoints();
    return res;
}


template<typename Shape, typename T, int Tag>
void
Mesh<Shape, T, Tag>::Localization::computeBarycentersWorld()
{
    LOG(INFO) << "computeBarycentersWorld : run mpi::all_gather\n";
    auto mesh = M_mesh.lock();
    M_barycentersWorld = std::vector<boost::tuple<bool,node_type> >( mesh->worldComm().localSize() );
    mpi::all_gather( mesh->worldComm().localComm(),
                     boost::make_tuple( this->kdtree()->nPoints() > 0 , this->barycenter() ),
                     *M_barycentersWorld );
}



} // Feel

#endif // FEELPP_MESH_LOCALIZATION_IMPL_HPP

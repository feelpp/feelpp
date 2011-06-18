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

namespace Feel
{


template <typename IteratorRange,typename Im,typename Expr>
class QuadPtLocalization
{
public :

    //static const size_type context = vm::DIV|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE|vm::POINT; //Expr::context;
    static const size_type context = Expr::context|vm::POINT;


    typedef typename boost::tuples::template element<1, IteratorRange>::type element_iterator_type;
    typedef typename element_iterator_type::value_type::GeoShape GeoShape;
    typedef Mesh<GeoShape> mesh_type;
    typedef typename mesh_type::gm_type gm_type;
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename gm_type::template Context<context, geoelement_type> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
    typedef typename gm_type::precompute_type pc_type;
    typedef typename gm_type::precompute_ptrtype pc_ptrtype;

    typedef typename mesh_type::value_type value_type;
    typedef typename mesh_type::node_type node_type;
    typedef typename matrix_node<value_type>::type matrix_node_type;

    //--------------------------------------------------------------------------------------//

    typedef typename mpl::if_<mpl::bool_<geoelement_type::is_simplex>,
                              mpl::identity<typename Im::template apply<geoelement_type::nDim, value_type, Simplex>::type >,
                              mpl::identity<typename Im::template apply<geoelement_type::nDim, value_type, Hypercube>::type >
                              >::type::type im_type;

    //--------------------------------------------------------------------------------------//

    // temporary container
    typedef std::map< size_type,std::list<boost::tuple<size_type,size_type,node_type,node_type> > > result_temp1_type;//idTrial->list( q, idq, ptRefTest, ptRefTrial )
    // sub container
    typedef std::map< size_type,boost::tuple< std::vector< boost::tuple<size_type,size_type > >,matrix_node_type,matrix_node_type> > result_temp2_type;//idTrial->( mapFor<q,idq>, ptRefsTest, ptsRefTrial)
    // result container
    typedef std::list< boost::tuple< size_type, result_temp2_type> > result_container_type;// list(idTest, idTrial->( mapForQ, ptRefsTest, ptsRefTrial) )

    //--------------------------------------------------------------------------------------//


    QuadPtLocalization(IteratorRange const& elts, im_type const& /*__im*/ )
        :
        _M_eltbegin( elts.template get<1>() ),
        _M_eltend( elts.template get<2>() ),
        _M_im( )
    {}

    QuadPtLocalization(element_iterator_type  elts_it,
                       element_iterator_type  elts_en,
                       im_type const& /*__im*/ )
        :
        _M_eltbegin( elts_it ),
        _M_eltend( elts_en ),
        _M_im( )
    {}


    im_type const& im() const { return _M_im; }

    element_iterator_type beginElement() const { return _M_eltbegin; }

    element_iterator_type endElement() const { return _M_eltend; }

    result_container_type const & result() const { return _M_res; }

    //std::vector<double> const & exprValue() const { return M_exprValue; }

    size_type eltForThisQuadPt(size_type theIdq) { return theIdq/this->im().nPoints(); }

    size_type qIndexForThisGlobalQuadPt(size_type theIdq) { return theIdq%this->im().nPoints(); }

    gmc_ptrtype gmcForThisElt(size_type theIdElt)
    {
        auto elt_it = this->beginElement();
        for (uint i=0;i<theIdElt;++i) ++elt_it;

        pc_ptrtype geopc( new pc_type( elt_it->gm(), this->im().points() ) );
        gmc_ptrtype gmc( new gmc_type(elt_it->gm(),*elt_it, geopc ) );

        return gmc;
    }

    template <typename Mesh1Type,typename Mesh2Type>
    void
    update(boost::shared_ptr<Mesh1Type> meshTest, boost::shared_ptr<Mesh2Type> meshTrial )
    {
        typedef Mesh1Type mesh_test_type;
        typedef Mesh2Type mesh_trial_type;

        //clean result
        _M_res.clear();

        auto elt_it = this->beginElement();
        auto elt_en = this->endElement();

        pc_ptrtype geopc( new pc_type( elt_it->gm(), this->im().points() ) );
        gmc_ptrtype gmc( new gmc_type(elt_it->gm(),*elt_it, geopc ) );

        auto meshTrialLocalization = meshTrial->tool_localization();
        meshTrialLocalization->updateForUse();

        auto meshTestLocalization = meshTest->tool_localization();
        meshTestLocalization->updateForUse();

        auto nQuadPtsInElt = this->im().nPoints();
        auto nElts = std::distance(elt_it,elt_en);
        auto nQuadPts = nElts*nQuadPtsInElt;
        auto nEltTrial= meshTrial->numElements();
        auto nEltTest= meshTest->numElements();

        std::vector<std::list<boost::tuple< size_type,size_type,node_type> > > trialPtQuadToElt( nEltTrial) ;
        std::vector<std::list<boost::tuple< size_type,size_type,node_type> > > testPtQuadToElt( nEltTest );
        //std::vector<double> exprValue( nQuadPts );
        //M_vecGmc.resize( nElts );
        //M_exprValue.resize( nQuadPts );
        std::vector<bool> quadPtDone( nQuadPts );
        std::fill( quadPtDone.begin(), quadPtDone.end(), false );

        std::vector<std::list<size_type> > EltCoupled(nEltTest);
        for( size_type ide = 0; elt_it != elt_en; ++elt_it, ++ide )
            {
                gmc->update( *elt_it );

                for( int q = 0; q <  gmc->nPoints(); ++ q )
                    {
                        // cpt of quad pt
                        size_type idq = gmc->nPoints()*ide+q;
                        //M_exprValue[idq] = this->im().weight(q)*gmc->J(q);//inutile

                        // search in trial mesh
                        auto trialAnalysis= meshTrialLocalization->searchElement(gmc->xReal( q ));
                        auto trialIdElt = trialAnalysis.get<1>();
                        auto trialNodeRef = trialAnalysis.get<2>();
                        trialPtQuadToElt[trialIdElt].push_back( boost::make_tuple(idq,q,trialNodeRef) );
                        // search in test mesh
                        auto testAnalysis = meshTestLocalization->searchElement(gmc->xReal( q ));
                        auto testIdElt = testAnalysis.get<1>();
                        auto testNodeRef = testAnalysis.get<2>();
                        testPtQuadToElt[testIdElt].push_back( boost::make_tuple(idq,q,testNodeRef) );

                        // relation between test and trial
                        if (std::find( EltCoupled[testIdElt].begin(),EltCoupled[testIdElt].end(),trialIdElt )==EltCoupled[testIdElt].end())
                            EltCoupled[testIdElt].push_back(trialIdElt);
                    }
            } // end for( size_type ide ... )


        auto eltCoupled_it = EltCoupled.begin();
        auto eltCoupled_en = EltCoupled.end();
        size_type theIdEltTest = 0;
        for ( ; eltCoupled_it != eltCoupled_en ; ++eltCoupled_it, ++theIdEltTest )
            {
                //init map
                result_temp1_type mapTrial2idq;
                //auto eltTrial_it = eltCoupled_it->begin();
                //auto eltTrial_en = eltCoupled_it->end();
                //for ( ; eltTrial_it != eltTrial_en ; ++eltTrial_it ) mapTrial2idq[*eltTrial_it].clear();

                // search
                auto quadPtTest_it = testPtQuadToElt[theIdEltTest].begin();
                auto quadPtTest_en = testPtQuadToElt[theIdEltTest].end();
                for ( ; quadPtTest_it != quadPtTest_en ; ++quadPtTest_it )
                    {
                        // get test data
                        auto theIdq = quadPtTest_it->get<0>();
                        auto theq = quadPtTest_it->get<1>();
                        //std::cout << "\nidq modulo : " << theIdq%nQuadPtsInElt << " q " << theq << " elt " << theIdq/nQuadPtsInElt;
                        auto theNodeRefTest = quadPtTest_it->get<2>();
                        // search in trial data
                        auto theRes = this->findQuadPt(*eltCoupled_it,trialPtQuadToElt,theIdq );
                        // get trial data
                        auto theIdEltTrial = theRes.get<0>();
                        auto theNodeRefTrial = theRes.get<1>();
                        // add into the temporary map
                        mapTrial2idq[theIdEltTrial].push_back(boost::make_tuple(theq,theIdq,theNodeRefTest,theNodeRefTrial));
                    }

                result_temp2_type mapiIdTrial2qAndPtRef;
                // build matrix_node of ptsRef (test and trial)
                auto map_it = mapTrial2idq.begin();
                auto map_en = mapTrial2idq.end();
                for ( ; map_it != map_en ; ++map_it)
                    {
                        // init matrix node
                        auto nPtsRef = map_it->second.size();
                        matrix_node_type ptsRefTest(mesh_test_type::nDim,nPtsRef);
                        matrix_node_type ptsRefTrial(mesh_trial_type::nDim,nPtsRef);
                        // init map between index node and index quad point
                        std::vector<boost::tuple<size_type,size_type> > indexLocalToQuad(nPtsRef);
                        // build
                        auto idq_it = map_it->second.begin();
                        auto idq_en = map_it->second.end();
                        uint cptIdq = 0;
                        for ( ; idq_it != idq_en ; ++idq_it, ++cptIdq )
                            {
                                indexLocalToQuad[cptIdq] = boost::make_tuple(idq_it->get<0>(),idq_it->get<1>());
                                ublas::column( ptsRefTest, cptIdq ) = idq_it->get<2>();
                                ublas::column( ptsRefTrial, cptIdq ) = idq_it->get<3>();
                            }
                        mapiIdTrial2qAndPtRef[map_it->first] = boost::make_tuple(indexLocalToQuad,ptsRefTest,ptsRefTrial);

                    }

                // add to result container
                _M_res.push_back(boost::make_tuple(theIdEltTest,mapiIdTrial2qAndPtRef));
            }




    } // update

private :

    boost::tuple<size_type,node_type >
    findQuadPt(std::list<size_type> const & eltCoupled_it,
               std::vector<std::list<boost::tuple< size_type,size_type,node_type> > > const & trialPtQuadToElt,
               size_type idq)
{
    size_type resId;
    ublas::vector<double> resNodeRef;
    bool find=false;

    auto eltTrial_it = eltCoupled_it.begin();
    auto eltTrial_en = eltCoupled_it.end();
    while (!find && eltTrial_it != eltTrial_en )
        {
            auto quadPtTrial_it = trialPtQuadToElt[*eltTrial_it].begin();
            auto quadPtTrial_en = trialPtQuadToElt[*eltTrial_it].end();
            while (!find && quadPtTrial_it != quadPtTrial_en )
                {
                    if ( quadPtTrial_it->get<0>() == idq ) { find = true;resNodeRef=quadPtTrial_it->get<2>(); }
                    else ++quadPtTrial_it;
                }
            if (!find) ++eltTrial_it;
            else resId=*eltTrial_it;
        }

    return boost::make_tuple(resId,resNodeRef);
}




    element_iterator_type _M_eltbegin;
    element_iterator_type _M_eltend;
    mutable im_type _M_im;

    result_container_type _M_res;

    //std::vector<gmc_ptrtype> M_vecGmc;

    //std::vector<double> M_exprValue;//temporaire!!!

}; // QuadPtLocalization


} // Feel
#endif

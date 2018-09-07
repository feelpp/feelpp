#ifndef __FUNCTIONSUP_H
#define __FUNCTIONSUP_H 1

#include <feel/feeldiscr/functionspace.hpp>
//#include <boost/mpl/int.hpp>
//#include <feel/feelvf/vf.hpp>

namespace Feel
{

    template<typename ElementRange,typename eltType,typename ExprType,typename vectorType >
    void
    modifVec(std::list<ElementRange> const& __r, eltType const& u,vectorType & UnVec,ExprType const& expr,
             size_type rowstart, int ComponentShiftFactor,
             mpl::int_<MESH_FACES> /**/ )
    {
        //using namespace Feel::vf;
        typedef typename eltType::functionspace_type::mesh_type mesh_type;
        typedef typename mesh_type::face_iterator face_iterator;
        typedef typename mesh_type::face_const_iterator face_const_iterator;

        typedef typename mesh_type::element_type geoelement_type;
        typedef typename geoelement_type::face_type face_type;
        // basis
        typedef typename eltType::functionspace_type::fe_type fe_type;

        typedef typename eltType::functionspace_type::dof_type dof_type;

        const size_type context = ExprType::context|vm::POINT;

        // geometric mapping context
        typedef typename mesh_type::gm_type gm_type;
        typedef std::shared_ptr<gm_type> gm_ptrtype;
        typedef typename gm_type::template Context<context, geoelement_type> gmc_type;
        typedef std::shared_ptr<gmc_type> gmc_ptrtype;
        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;

        typedef typename ExprType::template tensor<map_gmc_type> t_expr_type;

        if ( __r.size() == 0 ) return;
        auto __face_it =  __r.begin()->template get<1>();
        auto __face_en =  __r.begin()->template get<2>();
        //if ( __face_it == __face_en ) return;
        bool findAFace = false;
        for( auto lit = __r.begin(), len = __r.end(); lit != len; ++lit )
        {
            __face_it = lit->template get<1>();
            __face_en = lit->template get<2>();
            if ( __face_it != __face_en )
            {
                findAFace=true;
                break;
            }
        }
        if ( !findAFace ) return;

        // get the first face properly connected
        bool findAFaceToInit=false;
        for( auto lit = __r.begin(), len = __r.end(); lit != len; ++lit )
        {
            __face_it = lit->template get<1>();
            __face_en = lit->template get<2>();
            for( ; __face_it != __face_en; ++__face_it )
            {
                if ( boost::unwrap_ref(*__face_it).isConnectedTo0() )
                {
                    findAFaceToInit=true;
                    break;
                }
            }
            if ( findAFaceToInit ) break;
        }
        CHECK( findAFaceToInit ) << "not find a face to init\n";


        size_type nbFaceDof = invalid_size_type_value;
        if ( !fe_type::is_modal )
            nbFaceDof = ( face_type::numVertices * fe_type::nDofPerVertex +
                          face_type::numEdges * fe_type::nDofPerEdge +
                          face_type::numFaces * fe_type::nDofPerFace );
        else
            nbFaceDof = face_type::numVertices * fe_type::nDofPerVertex;



        //dof_type const* __dof = u.functionSpace()->dof().get();
        fe_type const* __fe = u.functionSpace()->fe().get();
        gm_ptrtype __gm( new gm_type );

        //
        // Precompute some data in the reference element for
        // geometric mapping and reference finite element
        //
        typedef typename geoelement_type::permutation_type permutation_type;
        typedef typename gm_type::precompute_ptrtype geopc_ptrtype;
        typedef typename gm_type::precompute_type geopc_type;
        std::vector<std::map<permutation_type, geopc_ptrtype> > __geopc( geoelement_type::numTopologicalFaces );
        for ( uint16_type __f = 0; __f < geoelement_type::numTopologicalFaces; ++__f )
        {
            permutation_type __p( permutation_type::IDENTITY );
            __geopc[__f][__p] = geopc_ptrtype(  new geopc_type( __gm, __fe->points( __f ) ) );
        }
        auto const& initface = boost::unwrap_ref(*__face_it);
        uint16_type __face_id = initface.pos_first();
        gmc_ptrtype __c( new gmc_type( __gm, initface.element(0), __geopc, __face_id ) );

        map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
        t_expr_type LExpr( expr, mapgmc );

        std::vector<bool> dofdone( u.functionSpace()->dof()->nLocalDofWithGhost(), false );

        auto const& dmVec = UnVec->map();
        int basisIndexDmInVec = dmVec.databaseIndexFromContainerId( rowstart );
        int subBasisIndexDm = (u.start()>0)? 1 : 0;
        basisIndexDmInVec+=subBasisIndexDm;//TODO!!!!!!!!!
        //face_const_iterator __face_it, __face_en;

        for( auto lit = __r.begin(), len = __r.end(); lit != len; ++lit )
        {
            __face_it = lit->template get<1>();
            __face_en = lit->template get<2>();
            for ( ; __face_it != __face_en; ++__face_it )
            {
                auto const& curface = boost::unwrap_ref( *__face_it );
                uint16_type __face_id = curface.pos_first();
                __c->update( curface.element(0), __face_id );
                map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
                LExpr.update( mapgmc );

                for (uint c1=0;c1<eltType::nComponents1;c1++)
                    for (uint c2=0;c2<eltType::nComponents2;c2++)
                    {
                        for ( uint16_type l = 0; l < nbFaceDof; ++l )
                        {
                            size_type index = boost::get<0>(u.functionSpace()->dof()->faceLocalToGlobal( curface.id(), l, c1 ));
                            if ( dofdone[index] ) continue;
                            double __value=LExpr.evalq( c1, c2, l );
#if 0
                            size_type thedof =  u.start() + ComponentShiftFactor*index;
                            //u( thedof ) =  __value;
                            UnVec->set(rowstart+thedof,__value);
#else
                            size_type thedof = /*(is_comp_space?Elem1::nComponents:1)*/ComponentShiftFactor*index;
                            thedof = dmVec.dofIdToContainerId( basisIndexDmInVec ,thedof );
                            UnVec->set(thedof,__value);
#endif
                            dofdone[index] = true;
                        }
                    }
            }
        }
        //UnVec->close();
    } // modifVec


    template<typename ElementRange,typename eltType,typename ExprType,typename vectorType >
    void
    modifVec(std::list<ElementRange> const& __r, eltType const& u,vectorType & UnVec,ExprType const& expr,
             size_type rowstart, int ComponentShiftFactor,
             mpl::int_<MESH_EDGES> /**/ )
    {
        const size_type context = ExprType::context|vm::POINT;

        auto mesh = u.functionSpace()->mesh().get();
        auto const* dof = u.functionSpace()->dof().get();
        auto const* fe = u.functionSpace()->fe().get();


        if ( __r.size() == 0 ) return;
        auto edge_it =  __r.begin()->template get<1>();
        auto edge_en =  __r.begin()->template get<2>();

        bool findAEdge = false;
        for( auto lit = __r.begin(), len = __r.end(); lit != len; ++lit )
        {
            edge_it = lit->template get<1>();
            edge_en = lit->template get<2>();
            if ( edge_it != edge_en )
            {
                findAEdge=true;
                break;
            }
        }
        if ( !findAEdge ) return;

        auto const& edgeForInit = boost::unwrap_ref( *edge_it );

        auto gm = mesh->gm();
        //auto const& firstEntity = *entity_it;
        size_type eid = edgeForInit.elements().begin()->first;
        size_type ptid_in_element = edgeForInit.elements().begin()->second;
        auto const& elt = mesh->element( eid );
        auto geopc = gm->preCompute( fe->edgePoints(ptid_in_element) );
        auto ctx = gm->template context<context>( elt, geopc );
        auto expr_evaluator = expr.evaluator( mapgmc(ctx) );
        auto IhLoc = fe->edgeLocalInterpolant();

        std::vector<bool> dofdone( dof->nLocalDofWithGhost(), false );

        for( auto const& lit : __r )
        {
            edge_it = lit.template get<1>();
            edge_en = lit.template get<2>();
            for ( ; edge_it != edge_en;++edge_it )
            {
                auto const& theedge = boost::unwrap_ref( *edge_it );

                if ( theedge.isGhostCell() )
                {
                    LOG(WARNING) << "edge id : " << theedge.id() << " is a ghost edge";
                    continue;
                }

                size_type eid = theedge.elements().begin()->first;
                size_type edgeid_in_element = theedge.elements().begin()->second;
                auto const& elt = mesh->element( eid );
                geopc = gm->preCompute( fe->edgePoints(edgeid_in_element) );
                ctx->update( elt, geopc );
                expr_evaluator.update( mapgmc( ctx ) );
                fe->edgeInterpolate( expr_evaluator, IhLoc );

                for( auto const& ldof : u.functionSpace()->dof()->edgeLocalDof( eid, edgeid_in_element ) )
                {
                    size_type index = ldof.index();
                    if ( dofdone[index] ) continue;
                    //size_type thedof = u.start()+ (is_comp_space?Elem1::nComponents:1)*ldof.index();
                    size_type thedof = u.start() + ComponentShiftFactor*index;
                    double value = ldof.sign()*IhLoc( ldof.localDofInFace() );
                    UnVec->set(rowstart+thedof,value);
                    dofdone[index] = true;
                }
            } // edge_it
        } // lit
    }

    template<typename ElementRange,typename eltType,typename ExprType,typename vectorType >
    void
    modifVec(std::list<ElementRange> const& __r, eltType const& u,vectorType & UnVec,ExprType const& expr,
             size_type rowstart, int ComponentShiftFactor,
             mpl::int_<MESH_POINTS> /**/ )
    {
        const size_type context = ExprType::context|vm::POINT;

        auto mesh = u.functionSpace()->mesh().get();
        auto const* dof = u.functionSpace()->dof().get();
        auto const* fe = u.functionSpace()->fe().get();

        if ( __r.size() == 0 ) return;
        auto point_it =  __r.begin()->template get<1>();
        auto point_en =  __r.begin()->template get<2>();

        bool findAPoint = false;
        for( auto lit = __r.begin(), len = __r.end(); lit != len; ++lit )
        {
            point_it = lit->template get<1>();
            point_en = lit->template get<2>();
            if ( point_it != point_en )
            {
                findAPoint=true;
                break;
            }
        }
        if ( !findAPoint ) return;

        auto const& pointForInit = boost::unwrap_ref( *point_it );

        size_type eid = pointForInit.elements().begin()->first;
        size_type ptid_in_element = pointForInit.elements().begin()->second;

        auto const& elt = mesh->element( eid );
        auto gm = mesh->gm();
        auto geopc = gm->preCompute( fe->vertexPoints(ptid_in_element) );
        auto ctx = gm->template context<context>( elt, geopc );
        auto expr_evaluator = expr.evaluator( mapgmc(ctx) );
        auto IhLoc = fe->vertexLocalInterpolant();

        std::vector<bool> dofdone( dof->nLocalDofWithGhost(), false );

        for( auto const& lit : __r )
        {
            point_it = lit.template get<1>();
            point_en = lit.template get<2>();
            DVLOG(2) << "point " << point_it->id() << " with marker " << point_it->marker() << " nb: " << std::distance(point_it,point_en);

            if ( point_it == point_en )
                continue;

            for ( ; point_it != point_en;++point_it )
            {
                auto const& thept = boost::unwrap_ref( *point_it );

                size_type eid = thept.elements().begin()->first;
                size_type ptid_in_element = thept.elements().begin()->second;
                auto const& elt = mesh->element( eid );
                geopc = gm->preCompute( fe->vertexPoints(ptid_in_element) );
                ctx->update( elt, ptid_in_element, geopc, mpl::int_<0>() );
                expr_evaluator.update( mapgmc( ctx ) );
                fe->vertexInterpolate( expr_evaluator, IhLoc );

                for (int c1=0;c1<eltType::nComponents1;c1++)
                    //for( int c = 0; c < (is_product?nComponents:1); ++c )
                {
                    size_type index = dof->localToGlobal( eid, ptid_in_element, c1 ).index();
                    //size_type thedof = u.start()+ (is_comp_space?Elem1::nComponents:1)*index; // global dof
                    size_type thedof = u.start() + ComponentShiftFactor*index;
                    if ( dofdone[index] ) continue;
                    double value = IhLoc( c1 );
                    UnVec->set(rowstart+thedof,value);
                    dofdone[index] = true;
                }
            }
        }

    }



    template<typename RangeType,typename eltType,typename ExprType,typename vectorType >
    void
    modifVec(std::list<RangeType> const& __r, eltType const& u,vectorType & UnVec,ExprType const& expr,size_type rowstart=0, int ComponentShiftFactor = 1)
    {
        const int iDim = boost::tuples::template element<0, RangeType>::type::value;
        modifVec( __r, u,UnVec,expr,rowstart,ComponentShiftFactor,mpl::int_<iDim>() );
    }

    template<typename RangeType,typename eltType,typename ExprType,typename vectorType>
    void
    modifVec(RangeType const& __r, eltType const& u,vectorType & UnVec,ExprType const& expr,size_type rowstart=0,int ComponentShiftFactor = 1)
    {
        std::list<RangeType> listRange = { __r };
        modifVec( listRange,u,UnVec,expr,rowstart,ComponentShiftFactor);
    }


} // namespace Feel

#endif // FUNCTIONSUP

#ifndef __FUNCTIONSUP_H
#define __FUNCTIONSUP_H 1

#include <feel/feeldiscr/functionspace.hpp>
//#include <feel/feelvf/vf.hpp>

namespace Feel
{

    template<typename ElementRange,typename eltType,typename ExprType,typename vectorType >
    void
    modifVec(std::list<ElementRange> const& __r, eltType const& u,vectorType & UnVec,ExprType const& expr,size_type rowstart=0, int ComponentShiftFactor = 1)
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
        typedef boost::shared_ptr<gm_type> gm_ptrtype;
        typedef typename gm_type::template Context<context, geoelement_type> gmc_type;
        typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
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

        uint16_type __face_id = __face_it->pos_first();
        gmc_ptrtype __c( new gmc_type( __gm, __face_it->element(0), __geopc, __face_id ) );

        map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
        t_expr_type LExpr( expr, mapgmc );




        //face_const_iterator __face_it, __face_en;

        for( auto lit = __r.begin(), len = __r.end(); lit != len; ++lit )
        {
            __face_it = lit->template get<1>();
            __face_en = lit->template get<2>();
        for ( ; __face_it != __face_en; ++__face_it )
            {

                uint16_type __face_id = __face_it->pos_first();
                __c->update( __face_it->element(0), __face_id );

                map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );

                LExpr.update( mapgmc );



                for (uint c1=0;c1<eltType::nComponents1;c1++)
                    for (uint c2=0;c2<eltType::nComponents2;c2++)
                        {
                            for ( uint16_type l = 0; l < nbFaceDof; ++l )
                                {
                                    double __value=LExpr.evalq( c1, c2, l );
                                    size_type thedof =  u.start() + ComponentShiftFactor*boost::get<0>(u.functionSpace()->dof()->faceLocalToGlobal( __face_it->id(), l, c1 ));
                                    //u( thedof ) =  __value;
                                    UnVec->set(rowstart+thedof,__value);
                                }
                        }
            }
        }
        //UnVec->close();
    } // modifVec

    template<typename ElementRange,typename eltType,typename ExprType,typename vectorType>
    void
    modifVec(ElementRange const& __r, eltType const& u,vectorType & UnVec,ExprType const& expr,size_type rowstart=0,int ComponentShiftFactor = 1)
    {
        std::list<ElementRange> listRange = { __r };
        modifVec( listRange,u,UnVec,expr,rowstart,ComponentShiftFactor);
    }


} // namespace Feel

#endif // FUNCTIONSUP

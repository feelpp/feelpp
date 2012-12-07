#ifndef __FUNCTIONSUP_H
#define __FUNCTIONSUP_H 1

namespace Feel
{

//<int Dim, int Order, <uint16_type,uint16_type,uint16_type> class Entity>
template <typename eltType,typename ExprType,typename vectorType>
void
modifVec222( eltType/*element_fluid_velocity_type*/ & u,vectorType/*vector_ptrtype*/ & UnVec,ExprType expr,std::string marker )
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


    size_type nbFaceDof = invalid_size_type_value;

    if ( !fe_type::is_modal )
        nbFaceDof = ( face_type::numVertices * fe_type::nDofPerVertex +
                      face_type::numEdges * fe_type::nDofPerEdge +
                      face_type::numFaces * fe_type::nDofPerFace );

    else
        nbFaceDof = face_type::numVertices * fe_type::nDofPerVertex;



    dof_type const* __dof = u.functionSpace()->dof().get();

    fe_type const* __fe = u.functionSpace()->fe().get();


    gm_ptrtype __gm( new gm_type );



    typename MeshTraits<mesh_type>::marker_face_const_iterator __face_it, __face_en;
    boost::tie( boost::tuples::ignore, __face_it, __face_en ) =
        markedfaces( u.functionSpace()->mesh(), u.functionSpace()->mesh()->markerName( marker/*"ParoiH"*/ ) );

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
    gmc_ptrtype __c( new gmc_type( __gm, __face_it->element( 0 ), __geopc, __face_id ) );

    map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
    t_expr_type LExpr( expr, mapgmc );




    //face_const_iterator __face_it, __face_en;

    for ( ; __face_it != __face_en; ++__face_it )
    {

        uint16_type __face_id = __face_it->pos_first();
        __c->update( __face_it->element( 0 ), __face_id );

        map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );

        LExpr.update( mapgmc );



        for ( uint c1=0; c1<eltType::nComponents1; c1++ )
            for ( uint c2=0; c2<eltType::nComponents2; c2++ )
            {
                for ( uint16_type l = 0; l < nbFaceDof; ++l )
                {
                    double __value=LExpr.evalq( c1, c2, l );
                    size_type thedof =  u.start() + boost::get<0>( u.functionSpace()->dof()->faceLocalToGlobal( __face_it->id(), l, c1 ) );
                    //u( thedof ) =  __value;
                    UnVec->set( thedof,__value );
                }
            }
    }
}




template <typename ElementRange,typename eltType,typename ExprType,typename vectorType>
void
modifVec( ElementRange const& __r, eltType & u,vectorType & UnVec,ExprType expr )
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

    size_type nbFaceDof = invalid_size_type_value;

    if ( !fe_type::is_modal )
        nbFaceDof = ( face_type::numVertices * fe_type::nDofPerVertex +
                      face_type::numEdges * fe_type::nDofPerEdge +
                      face_type::numFaces * fe_type::nDofPerFace );

    else
        nbFaceDof = face_type::numVertices * fe_type::nDofPerVertex;



    dof_type const* __dof = u.functionSpace()->dof().get();

    fe_type const* __fe = u.functionSpace()->fe().get();


    gm_ptrtype __gm( new gm_type );


#if 0
    typename MeshTraits<mesh_type>::marker_face_const_iterator __face_it, __face_en;
    boost::tie( boost::tuples::ignore, __face_it, __face_en ) =
        markedfaces( u.functionSpace()->mesh(), u.functionSpace()->mesh()->markerName( marker/*"ParoiH"*/ ) );
#endif

    auto __face_it =  __r.template get<1>();
    auto __face_en =  __r.template get<2>();


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
    gmc_ptrtype __c( new gmc_type( __gm, __face_it->element( 0 ), __geopc, __face_id ) );

    map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
    t_expr_type LExpr( expr, mapgmc );




    //face_const_iterator __face_it, __face_en;

    for ( ; __face_it != __face_en; ++__face_it )
    {

        uint16_type __face_id = __face_it->pos_first();
        __c->update( __face_it->element( 0 ), __face_id );
        //std::cout << "marker: " << __face_it->marker() << "\n";
        map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );

        LExpr.update( mapgmc );



        for ( uint c1=0; c1<eltType::nComponents1; c1++ )
            for ( uint c2=0; c2<eltType::nComponents2; c2++ )
            {
                for ( uint16_type l = 0; l < nbFaceDof; ++l )
                {
                    double __value=LExpr.evalq( c1, c2, l );
                    size_type thedof =  u.start() + boost::get<0>( u.functionSpace()->dof()->faceLocalToGlobal( __face_it->id(), l, c1 ) );
                    //u( thedof ) =  __value;
                    UnVec->set( thedof,__value );
                }
            }
    }

    UnVec->close();
} // modifVec



} // namespace Feel

#endif // FUNCTIONSUP

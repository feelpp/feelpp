/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_CORE_MEASURE_POINTS_EVALUATION_HPP
#define FEELPP_TOOLBOXES_CORE_MEASURE_POINTS_EVALUATION_HPP 1

#include <feel/feelmodels/modelpostprocess.hpp>
#include <feel/feelmodels/modelcore/traits.hpp>
#include <feel/feelcore/tuple_utils.hpp>

namespace Feel
{
namespace FeelModels
{

template <typename TupleMeshType,typename... SpacesType>
class MeasurePointsEvaluation
{
    struct TransformSpaceToFieldNamesAndSpace
    {
        template <typename T>
        struct apply {
            using type = std::pair< std::set<std::string>,std::shared_ptr<T>>;
        };

        template <typename T>
        constexpr auto operator()(T const& t) const
            {
                return typename TransformSpaceToFieldNamesAndSpace::template apply<T>::type{};
            }
    };

    struct TransformSpaceToContext
    {
        template <typename T>
        struct apply {
            using context_type = typename T::second_type::element_type::Context;
            using context_ptrtype = std::shared_ptr<context_type>;
            using map_fieldname_to_nodeid_type = std::map<std::string, std::map<int,std::string>>; // fieldname -> ( (nodeIdInCtx1 -> pointName1) , (nodeIdInCtx2 -> pointName2 ), ... )
            using type = std::tuple<context_ptrtype, map_fieldname_to_nodeid_type >;
        };

        template <typename T>
        constexpr auto operator()(T const& t) const
            {
                using context_type = typename TransformSpaceToContext::template apply<T>::context_type;
                using context_ptrtype = typename TransformSpaceToContext::template apply<T>::context_ptrtype;
                using map_fieldname_to_nodeid_type = typename TransformSpaceToContext::template apply<T>::map_fieldname_to_nodeid_type;
                map_fieldname_to_nodeid_type mapFieldToExport;
                for ( std::string const& s : t.first )
                    mapFieldToExport[s];
                return typename TransformSpaceToContext::template apply<T>::type( std::make_shared<context_type>( t.second->context() ), mapFieldToExport );
            }
    };
#if 1
    struct TransformGeometricSpaceToContext
    {
        template <typename T>
        struct apply {
            using context_type = typename T::second_type::element_type::Context;
            using context_ptrtype = std::shared_ptr<context_type>;
            //using map_fieldname_to_nodeid_type = std::map<std::string, std::map<int,std::string>>; // fieldname -> ( (nodeIdInCtx1 -> pointName1) , (nodeIdInCtx2 -> pointName2 ), ... )
            using map_fieldname_to_nodeid_type = std::map<std::string, std::set<std::string> >; // fieldname -> ( pointName1,pointName2,...)
            using map_pointsname_to_expr_nodeid_type = std::map<std::string, std::tuple< std::map<std::string,ModelExpression> /*,std::set<index_type>*/ > >; // pointName1 -> ( ((exprname1->expr1),(exprname2->expr2),...), (nodeIdInCtx1,nodeIdInCtx2,...) )
            // exprname -> ( (nodeIdInCtx1 -> pointName1) , (nodeIdInCtx2 -> pointName2 ), ... )
            // exprname -> ( ( pointName1 -> ( nodeIdInCtx1,nodeIdInCtx1b,...) ), ( pointName2 -> (nodeIdInCtx2,nodeIdInCtx2b,...) ) )
            // pointName1 -> ( (exprname1,exprname2,...), (nodeIdInCtx1,nodeIdInCtx2,...) ),
            // pointName1 -> ( ((expr,exprname1),(expr,exprname2),...), (nodeIdInCtx1,nodeIdInCtx2,...) )
            using type = std::tuple<context_ptrtype, std::map<std::string,std::set<index_type>>, map_fieldname_to_nodeid_type, map_pointsname_to_expr_nodeid_type >;
        };

        template <typename T>
        constexpr auto operator()(T const& t) const
            {
                using apply_type = typename TransformGeometricSpaceToContext::template apply<T>;
                using context_type = typename apply_type::context_type;
                using context_ptrtype = typename apply_type::context_ptrtype;
                using map_fieldname_to_nodeid_type = typename apply_type::map_fieldname_to_nodeid_type;
                map_fieldname_to_nodeid_type mapFieldToExport;
                for ( std::string const& s : t.first )
                    mapFieldToExport[s];
                return typename apply_type::type( std::make_shared<context_type>( t.second->context() ), {}, mapFieldToExport, {} );
            }
    };
#endif

public :

    using tuple_mesh_type = TupleMeshType;
    using spaces_tuple_type = hana::tuple<SpacesType...>;
    using fieldnames_space_tuple_type = std::decay_t<decltype( hana::transform( spaces_tuple_type{}, TransformSpaceToFieldNamesAndSpace{} ) ) >;
    using context_fields_tuple_type = std::decay_t<decltype( hana::transform( fieldnames_space_tuple_type{}, TransformSpaceToContext{} ) ) >;

    using tuple_geospaces_with_names_type = std::decay_t<decltype( hana::transform( tuple_mesh_type{}, TransformSpaceToFieldNamesAndSpace{} ) ) >;
    using tuple_geocontext_type = std::decay_t<decltype( hana::transform( tuple_geospaces_with_names_type{}, TransformGeometricSpaceToContext{} ) ) >;

    MeasurePointsEvaluation( tuple_geospaces_with_names_type const& geospaces, fieldnames_space_tuple_type const& fieldNamesSpaces )
        :
        M_geoContexts( hana::transform( geospaces, TransformGeometricSpaceToContext{} ) ),
        M_contextFields( hana::transform( fieldNamesSpaces, TransformSpaceToContext{} ) )
        {}

    void
    init( ModelPostprocessPointPosition const& evalPoints )
        {
            auto const& ptPos = evalPoints.pointPosition();
            node_type ptCoord(3);
            for ( int c=0;c<3;++c )
                ptCoord[c]=ptPos.value()(c);

            auto const& fields = evalPoints.fields();

#if 0
            std::map<int,int> nodeAddedInCtx;
            for ( std::string const& field : fields )
            {

                int ctxId = 0;
                hana::for_each( M_contextFields,
                                [&field,&ptPos,&ptCoord,&nodeAddedInCtx,&ctxId](auto & x) {
                                    auto itFindField = std::get<1>( x ).find( field );
                                    if ( itFindField != std::get<1>( x ).end() )
                                    {
                                        if ( nodeAddedInCtx.find( ctxId ) == nodeAddedInCtx.end() )
                                        {
                                            auto const& fectx = std::get<0>( x );
                                            int nodeIdInCtx = fectx->nPoints();
                                            fectx->add( ptCoord );
                                            nodeAddedInCtx[ ctxId ] = nodeIdInCtx;
                                        }
                                        itFindField->second[nodeAddedInCtx.find( ctxId )->second] = ptPos.name();
                                    }
                                    ++ctxId;
                                });
            }
#endif


            std::map<int,int> nodeAddedInGeoCtx;
            int ctxId = 0;
            hana::for_each( M_geoContexts,
                            [&evalPoints,&fields,&ptPos,&ptCoord,&nodeAddedInGeoCtx,&ctxId](auto & x) {
                                auto & geoctx = std::get<0>( x );
                                auto & ptPosNameToNodeIds = std::get<1>( x );
                                auto & fieldsInCtx = std::get<2>( x );

                                std::set<std::string> fieldsUsed;
                                for ( std::string const& field : fields )
                                {
                                    if ( true ) // TODO check compatibility between field/geoctx
                                    {
#if 0
                                        if ( nodeAddedInGeoCtx.find( ctxId ) == nodeAddedInGeoCtx.end() )
                                        {
                                            int nodeIdInCtx = geoctx->nPoints();
                                            geoctx->add( ptCoord );
                                            nodeAddedInGeoCtx[ ctxId ] = nodeIdInCtx;
                                        }
                                        //nodeAddedInCtx[ ctxId ] = nodeIdInCtx;

                                        fieldsInCtx[field][nodeAddedInGeoCtx.find( ctxId )->second] = ptPos.name();
#else
                                        fieldsUsed.insert( field );
#endif
                                    }
                                }

                                auto & exprInCtx = std::get<3>( x );
                                std::map<std::string,ModelExpression> mapExprNameToExpr;
                                for ( auto const& [name,exprData] : evalPoints.expressions() )
                                {
                                    auto const& [mexpr,tag] = exprData;
                                    if ( true ) // TODO check compatibility between expr tag/geoctx
                                    {
                                        mapExprNameToExpr.emplace( name, mexpr );
                                    }
                                }
                                if ( !mapExprNameToExpr.empty() || !fieldsUsed.empty() )
                                {
                                    std::set<index_type> & nodeIds = ptPosNameToNodeIds[ptPos.name()];
                                    //std::set<index_type> nodeIds;
                                    if ( true )//nodeAddedInGeoCtx.find( ctxId ) == nodeAddedInGeoCtx.end() )
                                    {
                                        int nodeIdInCtx = geoctx->nPoints();
                                        geoctx->add( ptCoord );
                                        nodeIds.insert( nodeIdInCtx );
                                        nodeAddedInGeoCtx[ ctxId ] = nodeIdInCtx;
                                    }
                                    for ( std::string const& field : fieldsUsed )
                                        fieldsInCtx[field].insert( ptPos.name() );
                                    if ( !mapExprNameToExpr.empty() )
                                        exprInCtx[ptPos.name()] = std::make_tuple( std::move(mapExprNameToExpr)/*, std::move( nodeIds )*/ );
                                }
                                ++ctxId;
                            });




        }



    template <typename SymbolsExprType,typename FieldTupleType>
    void
    eval( std::vector<ModelPostprocessPointPosition> const& allEvalPoints, std::map<std::string,double> & res, SymbolsExprType const& se, FieldTupleType const& fieldTuple )
        {
            for ( auto const& evalPoints : allEvalPoints )
            {
                // TODO : update point position if node has moved
            }

            // update gmc ctx
            hana::for_each( M_geoContexts,
                            [this,&se](auto & x) {
                                
                                std::get<0>( x )->updateForUse(); // TODO VINCENT!!!!!!!!!!!!!!!!!
                                
                                size_type dynctx=0;
                                auto const& exprInCtx = std::get<3>( x );
                                for ( auto const& [ptPosName, exprDataAndNodeIds ] : exprInCtx )
                                {
                                    for ( auto const& exprData : std::get<0>( exprDataAndNodeIds ) )
                                    {
                                        //auto const& exprName = std::get<0>( exprData );
                                        auto const& mexpr = std::get<1>( exprData );
                                        dynctx = dynctx | mexpr.dynamicContext( se );
                                    }
                                }
                                std::get<0>( x )->template updateGmcContext<vm::DYNAMIC>( dynctx );
                            });
#if 0
            hana::for_each( M_contextFields,
                            [this,&fieldTuple,&res](auto const& x) {
                                hana::for_each( fieldTuple.tuple(),
                                                [this,&x,&res](auto const& y) {
                                                    if constexpr ( is_iterable_v<decltype(y)> )
                                                    {
                                                        for ( auto const& mfield : y )
                                                        {
                                                            this->evalFieldImpl( x, mfield, res );
                                                        }
                                                    }
                                                    else
                                                    {
                                                        //this->evalFieldImpl( x,y.first,y.second,res );
                                                    }
                                                }); // for_each fieldTuple
                            }); // for_each M_contextFields
#endif
            hana::for_each( M_geoContexts,
                            [this,&se,&fieldTuple,&res](auto const& x) {
                                auto const& geoctx = std::get<0>( x );
                                // fields
                                auto const& fieldsInCtx = std::get<2>( x );
                                for ( auto const& fieldInCtx : fieldsInCtx )
                                {
                                    auto const& fieldName = fieldInCtx.first;
                                    auto const& ptPosSet = fieldInCtx.second;
                                    hana::for_each( fieldTuple.tuple(), [this,&x,&geoctx,&fieldName,&ptPosSet,&res](auto const& y)
                                                    {
                                                    for ( auto const& mfield : y )
                                                    {
#if 0
                                                        this->evalFieldImpl( x, mfield, res );
#else
                                                        if ( fieldName != mfield.nameWithPrefix() )
                                                            continue;
                                                        //std::string fieldName = mfield.nameWithPrefix();
                                                        // auto itFindField = fieldInCtx.find( fieldName );
                                                        // if ( itFindField == fieldsInCtx.end() )
                                                        //     continue;

                                                        auto const& fieldFunc = mfield.field();

                                                        using FieldType = std::decay_t<decltype(fieldFunc)>;
                                                        if constexpr ( is_shared_ptr< FieldType >::value )
                                                        {
                                                            if ( !fieldFunc )
                                                                continue;
                                                        }
                                                        if constexpr( std::is_same_v< typename Feel::decay_type<FieldType>::functionspace_type::mesh_type,
                                                                      typename std::decay_t<decltype(*geoctx)>::functionspace_type::mesh_type > )
                                                        {
                                                            //std::set<index_type> nodeIds;
                                                            for ( std::string const& ptPosName : ptPosSet )
                                                            {
                                                                auto itFindNodeIds = std::get<1>( x ).find( ptPosName );
                                                                CHECK( itFindNodeIds != std::get<1>( x ).end() ) << "ptPosName not find";
                                                                auto const& nodeIds = itFindNodeIds->second;

                                                                std::string outputNameBase = (boost::format("Points_%1%_field_%2%")%ptPosName %fieldName).str();
                                                                this->evalImpl( geoctx, nodeIds, idv(fieldFunc), outputNameBase, res );
                                                            }
                                                        }
#endif
                                                    }
                                                    }); // for_each fieldTuple
                                }
                                // expressions
                                auto const& exprInCtx = std::get<3>( x );
                                for ( auto const& [ptPosName, exprDataAndNodeIds ] : exprInCtx )
                                {
                                    //auto const& nodeIds = std::get<1>( exprDataAndNodeIds );
                                    auto itFindNodeIds = std::get<1>( x ).find( ptPosName );
                                    CHECK( itFindNodeIds != std::get<1>( x ).end() ) << "ptPosName not find";
                                    auto const& nodeIds = itFindNodeIds->second;
                                    for ( auto const& exprData : std::get<0>( exprDataAndNodeIds ) )
                                    {
                                        auto const& exprName = std::get<0>( exprData );
                                        auto const& mexpr = std::get<1>( exprData );
                                        std::string outputNameBase = (boost::format("Points_%1%_expr_%2%")%ptPosName %exprName).str();
                                        hana::for_each( ModelExpression::expr_shapes, [this,&se,&x,&nodeIds,&outputNameBase,&mexpr,&res]( auto const& e_ij )
                                                        {
                                                            constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                                                            constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;
                                                            if ( mexpr.template hasExpr<ni,nj>() )
                                                            {
                                                                auto theExpr = expr( mexpr.template expr<ni,nj>(), se );
                                                                this->evalImpl( std::get<0>( x ), nodeIds, theExpr, outputNameBase, res );
                                                            }
                                                        });
                                    }
                                }
                            }); // for_each M_geoContexts

        }

private :
#if 0
    template <typename ContextDataType,typename MFieldType>
    void
    evalFieldImpl( ContextDataType const& x, MFieldType const& mfield, std::map<std::string,double> & res )
        {
            std::string fieldName = mfield.nameWithPrefix();
            auto const& fieldFunc = mfield.field();

            using FieldType = std::decay_t<decltype(fieldFunc)>;
            if constexpr ( is_shared_ptr< FieldType >::value )
                {
                    if ( !fieldFunc )
                        return;
                }

            auto const& fectx = std::get<0>( x );

            if constexpr( std::is_same_v< typename Feel::decay_type<FieldType>::functionspace_type::Context, std::decay_t<decltype(*fectx)> > ||
                          std::is_same_v< typename Feel::decay_type<FieldType>::functionspace_type::mesh_type, typename std::decay_t<decltype(*fectx)>::functionspace_type::mesh_type >
                          )
                {
                    auto itFindField = std::get<2>( x ).find( fieldName );
                    if ( itFindField != std::get<2>( x ).end() && !itFindField->second.empty() )
                    {
                        mfield.applyUpdateFunction();
                        auto expr = idv( fieldFunc );
                        auto evalAtNodes = evaluateFromContext( _context=*fectx,
                                                                _expr=expr );

                        typedef typename ExprTraitsFromContext<std::decay_t<decltype(*fectx)>,std::decay_t<decltype(expr)>>::shape shape_type;

                        for ( int nodeId=0;nodeId<fectx->nPoints();++nodeId )
                        {
                            auto itFindNodeIdInCtx = itFindField->second.find( nodeId );
                            if ( itFindNodeIdInCtx != itFindField->second.end() )
                            {
                                std::string ptPosName = itFindNodeIdInCtx->second;
                                std::string pointNameOutputBase = (boost::format("Points_%1%_%2%")%ptPosName %fieldName).str();
                                for ( int i=0;i<shape_type::M ;++i )
                                    for ( int j=0;j<shape_type::N ;++j )
                                    {
                                        double val = evalAtNodes( nodeId*shape_type::M*shape_type::N+i+j*shape_type::M );
                                        std::string pointNameOutput = pointNameOutputBase;
                                        if ( shape_type::M > 1 && shape_type::N > 1 )
                                            pointNameOutput = (boost::format("%1%_%2%_%3%")%pointNameOutputBase %i %j).str();
                                        else if ( shape_type::M > 1 )
                                            pointNameOutput = (boost::format("%1%_%2%")%pointNameOutputBase %i).str();
                                        else if ( shape_type::N > 1 )
                                            pointNameOutput = (boost::format("%1%_%2%")%pointNameOutputBase %j).str();
                                        res[pointNameOutput] = val;
                                    }
                            }
                        } //
                    }
                }

        }
#endif
    template <typename ContextDataType,typename ExprType>
    void
    evalImpl( ContextDataType const& ctx, std::set<index_type> const& nodeIds, ExprType const& theExpr, std::string const& outputNameBase, std::map<std::string,double> & res )
        {
            std::cout << "ModelMeasurePointEval evalImpl " << outputNameBase << " nodeIds " << nodeIds << std::endl;
            auto evalAtNodes = evaluateFromContext( _context=*ctx,
                                                    _expr=theExpr,
                                                    _points_used=nodeIds );

            typedef typename ExprTraitsFromContext<std::decay_t<decltype(*ctx)>,std::decay_t<decltype(theExpr)>>::shape shape_type;

            for ( int nodeId=0;nodeId<nodeIds.size()/*ctx->nPoints()*/;++nodeId )
            {
                for ( int i=0;i<shape_type::M ;++i )
                {
                    for ( int j=0;j<shape_type::N ;++j )
                    {
                        double val = evalAtNodes( nodeId*shape_type::M*shape_type::N+i+j*shape_type::M );
                        std::string outputName = outputNameBase;
                        if ( shape_type::M > 1 && shape_type::N > 1 )
                            outputName = (boost::format("%1%_%2%_%3%")%outputNameBase %i %j).str();
                        else if ( shape_type::M > 1 )
                            outputName = (boost::format("%1%_%2%")%outputNameBase %i).str();
                        else if ( shape_type::N > 1 )
                            outputName = (boost::format("%1%_%2%")%outputNameBase %j).str();
                        res[outputName] = val;
                    }
                }
            }
        }


private :
    tuple_geocontext_type M_geoContexts;
    context_fields_tuple_type M_contextFields;
};

} // namespace FeelModels
} // namespace Feel


#endif

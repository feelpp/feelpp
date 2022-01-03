/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_CORE_MEASURE_POINTS_EVALUATION_HPP
#define FEELPP_TOOLBOXES_CORE_MEASURE_POINTS_EVALUATION_HPP 1

#include <feel/feelmodels/modelpostprocess.hpp>
#include <feel/feelmodels/modelcore/traits.hpp>
#include <feel/feelcore/tuple_utils.hpp>
#include <feel/feelmodels/modelcore/modelmeasures.hpp>

namespace Feel
{
namespace FeelModels
{

class MeasurePointsEvaluationBase {
public :
    virtual ~MeasurePointsEvaluationBase() {}
};

template <typename MeshType>
class MeasurePointsEvaluation : public MeasurePointsEvaluationBase
{
#if 0
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
#endif

    struct PointOnContext
    {
        PointOnContext( std::string const& outputName, std::string const& outputType, bool includeCoordinates )
            :
            M_outputName( outputName ), M_outputType( outputType),
            M_includeCoordinates( includeCoordinates )
            {}
        PointOnContext( PointOnContext const& ) = default;//delete;
        PointOnContext( PointOnContext && ) = default;

        std::string const& outputName() const { return M_outputName; }
        std::string const& outputType() const { return M_outputType; }
        bool includeCoordinates() const { return M_includeCoordinates; }
        std::set<index_type> & nodeIds() { return M_nodeIds; }
        std::set<index_type> const& nodeIds() const { return M_nodeIds; }

    private :
        std::string M_outputName, M_outputType;
        bool M_includeCoordinates;
        std::set<index_type> M_nodeIds;
    };
#if 0
    struct TransformGeometricSpaceToContext
    {
        template <typename T>
        struct apply {
            using context_type = typename T::second_type::element_type::Context;
            using context_ptrtype = std::shared_ptr<context_type>;
            using map_fieldname_to_nodeid_type = std::map<std::string, std::set<std::string> >; // fieldname -> ( pointName1,pointName2,...)
            using map_pointsname_to_expr_nodeid_type = std::map<std::string, std::tuple< std::map<std::string,ModelExpression> > >; // pointName1 -> ( (exprname1->expr1),(exprname2->expr2), ... )
            using type = std::tuple<context_ptrtype, std::map<std::string,PointOnContext>, map_fieldname_to_nodeid_type, map_pointsname_to_expr_nodeid_type >;
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
    using mesh_type = MeshType;
    using geometric_space_type = GeometricSpace<mesh_type>;
    using geometric_context_type = typename geometric_space_type::Context;
    using geometric_context_ptrtype = std::shared_ptr<geometric_context_type>;

    using map_fieldname_to_nodeid_type = std::map<std::string, std::set<std::string> >; // fieldname -> ( pointName1,pointName2,...)
    using map_pointsname_to_expr_nodeid_type = std::map<std::string, std::tuple< std::map<std::string,ModelExpression> > >; // pointName1 -> ( (exprname1->expr1),(exprname2->expr2), ... )

#if 0
    using internal_container_type = std::tuple<geometric_context_ptrtype, std::map<std::string,PointOnContext>, map_fieldname_to_nodeid_type, map_pointsname_to_expr_nodeid_type >;
#endif
    using points_eval_desc_type = std::map<std::string, std::tuple<std::map<std::string,PointOnContext>, map_fieldname_to_nodeid_type, map_pointsname_to_expr_nodeid_type > >;
public :

#if 0
    using tuple_mesh_type = TupleMeshType;
    using tuple_geospaces_with_names_type = std::decay_t<decltype( hana::transform( tuple_mesh_type{}, TransformSpaceToFieldNamesAndSpace{} ) ) >;
    using tuple_geocontext_type = std::decay_t<decltype( hana::transform( tuple_geospaces_with_names_type{}, TransformGeometricSpaceToContext{} ) ) >;

    MeasurePointsEvaluation( tuple_geospaces_with_names_type const& geospaces )
        :
        M_geoContexts( hana::transform( geospaces, TransformGeometricSpaceToContext{} ) )
        {}
#endif
    MeasurePointsEvaluation( std::shared_ptr<geometric_space_type> geospace )
        :
        M_geoContext( std::make_shared<geometric_context_type>( geospace ) )
    {
        //std::get<0>( M_internalContainer ) = std::make_shared<geometric_context_type>( geospace ); //geospace->context();
    }

    void
    init( std::vector<ModelPostprocessPointPosition> const& allEvalPoints, std::string const& groupName = "" )
        {
            for ( auto const& evalPoints : allEvalPoints )
                this->init( evalPoints, groupName );
        }

    void
    init( ModelPostprocessPointPosition const& evalPoints, std::string const& groupName = "" )
        {
            auto const& fields = evalPoints.fields();

            auto & M_internalContainer = M_pointsEvalDesc[groupName];

            auto & geoctx = M_geoContext;//std::get<0>( M_internalContainer );
            auto & ptPosNameToNodeIds = std::get<0>( M_internalContainer );
            auto & fieldsInCtx = std::get<1>( M_internalContainer );

            std::set<std::string> fieldsUsed;
            for ( std::string const& field : fields )
            {
                if ( true ) // TODO check compatibility between field/geoctx
                {
                    fieldsUsed.insert( field );
                }
            }

            auto & exprInCtx = std::get<2>( M_internalContainer );
            std::map<std::string,ModelExpression> mapExprNameToExpr;
            for ( auto const& [name,exprData] : evalPoints.expressions() )
            {
                auto const& [mexpr,tag] = exprData;
                if ( true ) // TODO check compatibility between expr tag/geoctx
                {
                    mapExprNameToExpr.emplace( name, mexpr );
                }
            }

            bool includeCoord = evalPoints.includeCoordinates();

            if ( !mapExprNameToExpr.empty() || !fieldsUsed.empty() || includeCoord )
            {
                auto [itPointOnCtx,isInserted] = ptPosNameToNodeIds.emplace( std::make_pair( evalPoints.name(), PointOnContext( evalPoints.measuresOutput().name(),
                                                                                                                                evalPoints.measuresOutput().type(),
                                                                                                                                includeCoord ) ) );
                std::set<index_type> & nodeIds = itPointOnCtx->second.nodeIds();
                for ( auto const& ptOverGeometry : evalPoints.pointsOverAllGeometry() )
                {
                    for ( auto const& ptCoordEig : ptOverGeometry->coordinates() )
                    {
                        int nodeIdInCtx = geoctx->nPoints();
                        node_type ptCoord(3);
                        for ( int c=0;c</*3*/ptCoordEig.size();++c )
                            ptCoord[c]=ptCoordEig(c);
                        geoctx->add( ptCoord,false );
                        nodeIds.insert( nodeIdInCtx );
                    }
                }

                for ( std::string const& field : fieldsUsed )
                    fieldsInCtx[field].insert( evalPoints.name() );
                if ( !mapExprNameToExpr.empty() )
                    exprInCtx[evalPoints.name()] = std::make_tuple( std::move(mapExprNameToExpr) );
            }
        }



    template <typename SymbolsExprType,typename ModelFieldsType>
    void
    eval( std::string const& groupName, ModelMeasuresStorage & res, SymbolsExprType const& se, ModelFieldsType const& mfields )
        {
            //#if 0 // TODO

            auto itFindGroup = M_pointsEvalDesc.find( groupName );
            if ( itFindGroup == M_pointsEvalDesc.end() )
                return;
            auto & M_internalContainer = itFindGroup->second;

            auto mfieldIsCompatible = []( auto const& geoctx, auto const& mfield, std::string const& fieldName) {
                                          if ( fieldName != mfield.nameWithPrefix() )
                                              return false;

                                          auto const& fieldFunc = mfield.field();
                                          using FieldType = std::decay_t<decltype(fieldFunc)>;
                                          if constexpr ( is_shared_ptr< FieldType >::value )
                                                       {
                                                           if ( !fieldFunc )
                                                               return false;
                                                       }
                                          if constexpr ( !std::is_same_v< typename Feel::decay_type<FieldType>::functionspace_type::mesh_type,
                                                         typename std::decay_t<decltype(*geoctx)>::functionspace_type::mesh_type > )
                                                           return false;
                                          return true;
                                      };

            // update gmc ctx
            //hana::for_each( M_geoContexts,[this,&se,&mfields,&mfieldIsCompatible](auto & x)
            if ( true )
            {
                auto & geoctx = M_geoContext;//std::get<0>( M_internalContainer );
                                
                geoctx->updateForUse(); // TODO VINCENT!!!!!!!!!!!!!!!!!
                                

                size_type dynctx=0;

                auto const& fieldsInCtx = std::get<1>( M_internalContainer );
                for ( auto const& fieldInCtx : fieldsInCtx )
                {
                    auto const& fieldName = fieldInCtx.first;
                    hana::for_each( mfields.tuple(), [&geoctx,&fieldName,&dynctx,&mfieldIsCompatible](auto const& y)
                    {
                        for ( auto const& mfield : y )
                        {
                            if ( !mfieldIsCompatible(geoctx,mfield,fieldName) )
                                continue;

                            mfield.applyUpdateFunction();
                            auto idexpr = idv(mfield.field());
                            dynctx = dynctx | Feel::vf::dynamicContext( idexpr );
                        }
                    });
                }

                auto const& exprInCtx = std::get<2>( M_internalContainer );
                for ( auto const& [ptPosName, exprDataAndNodeIds ] : exprInCtx )
                {
                    for ( auto const& exprData : std::get<0>( exprDataAndNodeIds ) )
                    {
                        //auto const& exprName = std::get<0>( exprData );
                        auto const& mexpr = std::get<1>( exprData );
                        dynctx = dynctx | mexpr.dynamicContext( se );
                    }
                }
                geoctx->template updateGmcContext<vm::DYNAMIC>( dynctx );
            }
                //});


            std::map<std::string, std::map<std::string, std::vector<std::pair<std::string,std::vector<double>>>>> measuresStored; //type -> ( outputName -> ( ( measureName1, measureValue1), ... ) )

            //hana::for_each( M_geoContexts,[this,&se,&mfields,&measuresStored,&mfieldIsCompatible](auto const& x)
            //{
            auto const& geoctx = M_geoContext;//std::get<0>( M_internalContainer );
                auto const& pointsOnCtx = std::get<0>( M_internalContainer );
                // fields
                auto const& fieldsInCtx = std::get<1>( M_internalContainer );
                for ( auto const& fieldInCtx : fieldsInCtx )
                {
                    auto const& fieldName = fieldInCtx.first;
                    auto const& ptPosSet = fieldInCtx.second;
                    hana::for_each( mfields.tuple(), [this,&geoctx,&pointsOnCtx,&fieldName,&ptPosSet,&measuresStored,&mfieldIsCompatible](auto const& y)
                    {
                        for ( auto const& mfield : y )
                        {
                            if ( !mfieldIsCompatible(geoctx,mfield,fieldName) )
                                continue;

                            auto const& fieldFunc = mfield.field();
                            using FieldType = std::decay_t<decltype(fieldFunc)>;


                            if constexpr( std::is_same_v< typename Feel::decay_type<FieldType>::functionspace_type::mesh_type,
                                          typename std::decay_t<decltype(*geoctx)>::functionspace_type::mesh_type > )
                            {
                                //mfield.applyUpdateFunction();
                                for ( std::string const& ptPosName : ptPosSet )
                                {
                                    auto itFindPointOnCtx = pointsOnCtx.find( ptPosName );
                                    CHECK( itFindPointOnCtx != pointsOnCtx.end() ) << "ptPosName not find";
                                    auto const& pointOnCtx = itFindPointOnCtx->second;

                                    std::string measurePrefix = fmt::format("Points_{}",ptPosName);
                                    std::string measureEvaluatedFrom = fmt::format("field_{}",fieldName);
                                    this->updateMeasureImpl( geoctx, pointOnCtx, idv(fieldFunc), measurePrefix, measureEvaluatedFrom, measuresStored );
                                }
                             }
                        }
                    }); // for_each mfields
                }
                // expressions
                auto const& exprInCtx = std::get<2>( M_internalContainer );
                for ( auto const& [ptPosName, exprDataAndNodeIds ] : exprInCtx )
                {
                    auto itFindPointOnCtx = pointsOnCtx.find( ptPosName );
                    CHECK( itFindPointOnCtx != pointsOnCtx.end() ) << "ptPosName not find";
                    auto const& pointOnCtx = itFindPointOnCtx->second;
                    std::string measurePrefix = fmt::format("Points_{}",ptPosName);
                    for ( auto const& exprData : std::get<0>( exprDataAndNodeIds ) )
                    {
                        auto const& exprName = std::get<0>( exprData );
                        auto const& mexpr = std::get<1>( exprData );
                        std::string measureEvaluatedFrom = fmt::format("expr_{}",exprName);
                        hana::for_each( ModelExpression::expr_shapes, [this,&se,&geoctx,&pointOnCtx,&measurePrefix,&measureEvaluatedFrom,&mexpr,&measuresStored]( auto const& e_ij )
                        {
                            constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                            constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;
                            if ( mexpr.template hasExpr<ni,nj>() )
                            {
                                auto theExpr = expr( mexpr.template expr<ni,nj>(), se );
                                this->updateMeasureImpl( geoctx, pointOnCtx, theExpr, measurePrefix, measureEvaluatedFrom, measuresStored );
                            }
                        });
                    }
                }

                // coordinates
                for ( auto const& [ptPosName,pointOnCtx] : pointsOnCtx )
                {
                    if ( !pointOnCtx.includeCoordinates() )
                        continue;
                    auto const& nodeIds = pointOnCtx.nodeIds();
                    static const uint16_type geoctxRealDim = std::decay_t<decltype(unwrap_ptr(geoctx))>::functionspace_type::mesh_type::nRealDim;
                    std::vector<eigen_vector_type<geoctxRealDim>> coords;
                    coords.reserve( nodeIds.size() );
                    auto const& ptsFromGeoCtx = geoctx->points();
                    eigen_vector_type<geoctxRealDim> coordEigen = eigen_vector_type<geoctxRealDim>::Zero();
                    for ( index_type nId : nodeIds )
                    {
                        auto const& ptFromGeoCtx = ptsFromGeoCtx[nId];
                        for (int d=0;d<geoctxRealDim;++d)
                            coordEigen(d) = ptFromGeoCtx[d];
                        coords.push_back( std::move( coordEigen ) );
                    }
                    std::string measurePrefix = fmt::format("Points_{}",ptPosName);
                    this->updateMeasureImpl( coords,  pointOnCtx.outputType(), pointOnCtx.outputName(), measurePrefix, "coordinates", measuresStored );
                }
                //}); // for_each M_geoContexts

            // update measures in ModelMeasuresStorage
            for ( auto & [type,byType] : measuresStored )
            {
                if ( type == "values" )
                {
                    for ( auto & [oname,byos] : byType )
                    {
                        for ( auto & byo : byos )
                        {
                            std::string const& measureName = byo.first;
                            for ( auto & val : byo.second )
                                res.setValue( oname, measureName, val );
                        }
                    }
                }
                else if ( type == "table" )
                {
                    for ( auto & [oname,byos] : byType )
                    {
                        if ( byos.empty() || byos.front().second.empty() )
                            continue;
                        int nRow = byos.front().second.size() + 1;
                        int nCol = byos.size();
                        Feel::Table table(nRow,nCol);
                        int j=0;
                        for ( auto & byo : byos )
                        {
                            std::string const& measureName = byo.first;
                            table(0,j) = measureName;
                            table.set_col(j,byo.second,1);
                            ++j;
                        }
                        res.setTable( oname, std::move( table ) );
                    }
                }
            }
            //#endif // TODO
        }

private :

    template <typename ContextDataType,typename ExprType>
    void
    updateMeasureImpl( ContextDataType const& ctx, PointOnContext const& pointOnCtx, ExprType const& theExpr,
                       std::string const& measurePrefix,std::string const& measureEvaluatedFrom,
                       std::map<std::string, std::map<std::string, std::vector<std::pair<std::string,std::vector<double>>>>> & measuresStored )
        {
            auto const& nodeIds = pointOnCtx.nodeIds();
            if ( nodeIds.empty() )
                return;
            //std::cout << "ModelMeasurePointEval evalImpl " << outputNameBase << " nodeIds " << nodeIds << std::endl;
            auto evalAtNodes = evaluateFromContext( _context=*ctx,
                                                    _expr=theExpr,
                                                    _points_used=nodeIds );

            //std::cout << "evalAtNodes=\n" << evalAtNodes << std::endl;

            typedef typename ExprTraitsFromContext<std::decay_t<decltype(*ctx)>,std::decay_t<decltype(theExpr)>>::shape shape_type;

            // copy result of evaluateFromContext in vector<Eigen::Matrix> : TODO evaluateFromContext should return this container
            index_type numberOfNodes = nodeIds.size();
            std::vector<eigen_matrix_type<shape_type::M,shape_type::N>> measureValues;
            measureValues.reserve( numberOfNodes );
            eigen_matrix_type<shape_type::M,shape_type::N> valMat = eigen_matrix_type<shape_type::M,shape_type::N>::Zero();
            for ( int nodeId=0;nodeId<numberOfNodes;++nodeId )
            {
                for ( int i=0;i<shape_type::M ;++i )
                    for ( int j=0;j<shape_type::N ;++j )
                        valMat(i,j) = evalAtNodes( nodeId*shape_type::M*shape_type::N+i+j*shape_type::M );
                measureValues.push_back( std::move( valMat ) );
            }

            this->updateMeasureImpl( measureValues, pointOnCtx.outputType(), pointOnCtx.outputName(), measurePrefix, measureEvaluatedFrom, measuresStored );
        }


    template <typename T>
    void
    updateMeasureImpl( std::vector<T> const& measures, std::string const& outputType, std::string const& outputName,
                       std::string const& measurePrefix, std::string const& measureEvaluatedFrom,
                       std::map<std::string, std::map<std::string, std::vector<std::pair<std::string,std::vector<double>>>>> & measuresStored )
        {
            if ( measures.empty() )
                return;

            int M = measures.front().rows();
            int N = measures.front().cols();
            if ( outputType == "values" )
            {
                auto & currentMeasure = measuresStored["values"][outputName];
                for ( int k=0;k<measures.size();++k )
                {
                    std::string measureNameBase = measures.size() > 1 ? fmt::format("{}_{}_{}",measurePrefix,k,measureEvaluatedFrom) : fmt::format("{}_{}",measurePrefix,measureEvaluatedFrom);
                    auto const& valEig = measures[k];
                    for ( int i=0;i<M ;++i )
                    {
                        for ( int j=0;j<N ;++j )
                        {
                            std::string measureName = measureNameBase;
                            if ( M > 1 && N > 1 )
                                measureName = fmt::format("{}_{}_{}",measureNameBase,i,j);
                            else if ( M > 1 )
                                measureName = fmt::format("{}_{}",measureNameBase,i);
                            else if ( N > 1 )
                                measureName = fmt::format("{}_{}",measureNameBase,j);

                            double val = valEig( i,j );
                            currentMeasure.push_back( std::make_pair(measureName, std::vector<double>(1,val) ) );
                        }
                    }
                }

            }
            else if ( outputType == "table" )
            {
                std::string measureNameBase = fmt::format("{}_{}",measurePrefix,measureEvaluatedFrom);
                for ( int i=0;i<M ;++i )
                {
                    for ( int j=0;j<N ;++j )
                    {
                        std::string measureName = measureNameBase;
                        if ( M > 1 && N > 1 )
                            measureName = fmt::format("{}_{}_{}",measureNameBase,i,j);
                        else if ( M > 1 )
                            measureName = fmt::format("{}_{}",measureNameBase,i);
                        else if ( N > 1 )
                            measureName = fmt::format("{}_{}",measureNameBase,j);

                        std::vector<double> currentMeasureValues;
                        currentMeasureValues.reserve( measures.size() );
                        for ( auto const& valEig : measures )
                        {
                            double val = valEig( i,j );
                            currentMeasureValues.push_back( val );
                        }
                        measuresStored["table"][outputName].push_back( std::make_pair( measureName, std::move(currentMeasureValues) ) );
                    }
                }
            }

        }

private :
    //tuple_geocontext_type M_geoContexts;
    //internal_container_type M_internalContainer;
    geometric_context_ptrtype M_geoContext;
    points_eval_desc_type M_pointsEvalDesc;
};

} // namespace FeelModels
} // namespace Feel


#endif

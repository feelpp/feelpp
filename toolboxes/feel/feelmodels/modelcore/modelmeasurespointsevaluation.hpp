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
    struct PointOnContext
    {
        PointOnContext() = default;
        PointOnContext( PointOnContext const& ) = default;//delete;
        PointOnContext( PointOnContext && ) = default;

        std::set<index_type> & nodeIds() { return M_nodeIds; }
        std::set<index_type> const& nodeIds() const { return M_nodeIds; }

    private :
        std::set<index_type> M_nodeIds;
    };

    using mesh_type = MeshType;
    using geometric_space_type = GeometricSpace<mesh_type>;
    using geometric_context_type = typename geometric_space_type::Context;
    using geometric_context_ptrtype = std::shared_ptr<geometric_context_type>;

    //using map_fieldname_to_nodeid_type = std::map<std::string, std::set<std::string> >; // fieldname -> ( pointName1,pointName2,...)
    //using map_pointsname_to_expr_nodeid_type = std::map<std::string, std::tuple< std::map<std::string,ModelExpression> > >; // pointName1 -> ( (exprname1->expr1),(exprname2->expr2), ... )
    //using points_eval_desc_type = std::map<std::string, std::tuple<std::map<std::string,PointOnContext>, map_fieldname_to_nodeid_type, map_pointsname_to_expr_nodeid_type > >;

    using points_eval_desc_type = std::map<std::string, std::map<std::string,PointOnContext> >; // groupName -> ( ( ptPosName -> PointOnContext ), ... )
public :

    MeasurePointsEvaluation( std::shared_ptr<geometric_space_type> geospace )
        :
        M_geoContext( std::make_shared<geometric_context_type>( geospace ) )
    {}

    void
    init( std::vector<ModelPostprocessPointPosition> const& allEvalPoints, std::string const& groupName = "" )
        {
            for ( auto const& evalPoints : allEvalPoints )
                this->init( evalPoints, groupName );
        }

    void
    init( ModelPostprocessPointPosition const& evalPoints, std::string const& groupName = "" )
        {
            // TODO check compatibility between tag and geoctx
            if ( !evalPoints.fields().empty() || !evalPoints.expressions().empty() || evalPoints.includeCoordinates() )
            {
                auto & ptPosNameToNodeIds = M_pointsEvalDesc[groupName];
                auto [itPointOnCtx,isInserted] = ptPosNameToNodeIds.emplace( std::make_pair( evalPoints.name(), PointOnContext{} ) );

                std::set<index_type> & nodeIds = itPointOnCtx->second.nodeIds();
                for ( auto const& ptOverGeometry : evalPoints.pointsOverAllGeometry() )
                {
                    for ( auto const& ptCoordEig : ptOverGeometry->coordinates() )
                    {
                        int nodeIdInCtx = M_geoContext->nPoints();
                        node_type ptCoord(3);
                        for ( int c=0;c</*3*/ptCoordEig.size();++c )
                            ptCoord[c]=ptCoordEig(c);
                        M_geoContext->add( ptCoord,false ); // TODO detect if the same point already registered
                        nodeIds.insert( nodeIdInCtx );
                    }
                }
            }
        }

    template <typename SymbolsExprType,typename ModelFieldsType>
    void
    apply( std::string const& groupName, ModelPostprocessPointPosition const& evalPoints, ModelMeasuresStorage & res, SymbolsExprType const& se, ModelFieldsType const& mfields )
        {
            auto itFindGroup = M_pointsEvalDesc.find( groupName );
            if ( itFindGroup == M_pointsEvalDesc.end() )
                return;
            auto const& pointsOnCtx = itFindGroup->second;

            std::string const& ptPosName = evalPoints.name();
            auto itFindPointOnCtx = pointsOnCtx.find( ptPosName );
            //CHECK( itFindPointOnCtx != pointsOnCtx.end() ) << "ptPosName not find";
            if ( itFindPointOnCtx == pointsOnCtx.end() )
                return;
            auto const& pointOnCtx = itFindPointOnCtx->second;


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
            if ( true )
            {
                // update geo ctx for use (necessary after adding a point in geoctx)
                M_geoContext->updateForUse();

                size_type dynctx=0;
                for ( auto const& fieldName : evalPoints.fields() )
                {
                    hana::for_each( mfields.tuple(), [this,&fieldName,&dynctx,&mfieldIsCompatible](auto const& y)
                    {
                        for ( auto const& mfield : y )
                        {
                            if ( !mfieldIsCompatible(M_geoContext,mfield,fieldName) )
                                continue;

                            mfield.applyUpdateFunction();
                            auto idexpr = idv(mfield.field());
                            dynctx = dynctx | Feel::vf::dynamicContext( idexpr );
                        }
                    });
                }

                for ( auto const& [exprName,mexpr] : evalPoints.expressions() )
                    dynctx = dynctx | mexpr.dynamicContext( se );

                // TODO : not apply this update if already done (i.e. same dynctx and number of points)
                M_geoContext->template updateGmcContext<vm::DYNAMIC>( dynctx );
            }


            std::string const& outputType = evalPoints.measuresOutput().type();
            std::string const& outputName = evalPoints.measuresOutput().name();
            std::string measurePrefix = fmt::format("Points_{}",ptPosName);

            std::vector<std::pair<std::string,std::vector<double>>> measuresStored; // ( ( measureName1, (measureValue1a,measureValue1b,...) ), ... )

            // fields
            for ( auto const& fieldName : evalPoints.fields() )
            {
                hana::for_each( mfields.tuple(), [this,&outputType,&pointOnCtx,&fieldName,&ptPosName,&measuresStored,&mfieldIsCompatible,&measurePrefix](auto const& y)
                {
                    for ( auto const& mfield : y )
                    {
                        if ( !mfieldIsCompatible(M_geoContext,mfield,fieldName) )
                            continue;

                        auto const& fieldFunc = mfield.field();
                        using FieldType = std::decay_t<decltype(fieldFunc)>;

                        if constexpr( std::is_same_v< typename Feel::decay_type<FieldType>::functionspace_type::mesh_type,
                                      typename std::decay_t<decltype(*M_geoContext)>::functionspace_type::mesh_type > )
                        {
                            std::string measureEvaluatedFrom = fmt::format("field_{}",fieldName);
                            this->updateMeasureImpl( M_geoContext, outputType, pointOnCtx, idv(fieldFunc), measurePrefix, measureEvaluatedFrom, measuresStored );
                         }
                    }
                }); // for_each mfields
            }


            // expressions
            for ( auto const& exprData : evalPoints.expressions() )
            {
                auto const& exprName = exprData.first;
                auto const& mexpr = exprData.second;
                std::string measureEvaluatedFrom = fmt::format("expr_{}",exprName);
                hana::for_each( ModelExpression::expr_shapes, [this,&se,&outputType,&pointOnCtx,&measurePrefix,&measureEvaluatedFrom,&mexpr,&measuresStored]( auto const& e_ij )
                {
                    constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                    constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;
                    if ( mexpr.template hasExpr<ni,nj>() )
                    {
                        auto theExpr = expr( mexpr.template expr<ni,nj>(), se );
                        this->updateMeasureImpl( M_geoContext, outputType, pointOnCtx, theExpr, measurePrefix, measureEvaluatedFrom, measuresStored );
                    }
                });
            }

            // coordinates
            if ( evalPoints.includeCoordinates() )
            {
                auto const& nodeIds = pointOnCtx.nodeIds();
                constexpr uint16_type geoctxRealDim = std::decay_t<decltype(unwrap_ptr(M_geoContext))>::functionspace_type::mesh_type::nRealDim;
                std::vector<eigen_vector_type<geoctxRealDim>> coords;
                coords.reserve( nodeIds.size() );
                auto const& ptsFromGeoCtx = M_geoContext->points();
                eigen_vector_type<geoctxRealDim> coordEigen = eigen_vector_type<geoctxRealDim>::Zero();
                for ( index_type nId : nodeIds )
                {
                    auto const& ptFromGeoCtx = ptsFromGeoCtx[nId];
                    for (int d=0;d<geoctxRealDim;++d)
                        coordEigen(d) = ptFromGeoCtx[d];
                    coords.push_back( std::move( coordEigen ) );
                }
                this->updateMeasureImpl( coords, outputType, measurePrefix, "coordinates", measuresStored );
            }


            // update measures in ModelMeasuresStorage
            if ( outputType == "values" )
            {
                for ( auto & byo : measuresStored )
                {
                    std::string const& measureName = byo.first;
                    for ( auto & val : byo.second )
                        res.setValue( outputName, measureName, val );
                }
            }
            else if ( outputType == "table" )
            {
                if ( !measuresStored.empty() && !measuresStored.front().second.empty() )
                {
                    int nRow = measuresStored.front().second.size() + 1;
                    int nCol = measuresStored.size();
                    Feel::Table table(nRow,nCol);
                    int j=0;
                    for ( auto & byo : measuresStored )
                    {
                        std::string const& measureName = byo.first;
                        table(0,j) = measureName;
                        table.set_col(j,byo.second,1);
                        ++j;
                    }
                    res.setTable( outputName, std::move( table ) );
                }
            }
        }
private :

    template <typename ContextDataType,typename ExprType>
    void
    updateMeasureImpl( ContextDataType const& ctx, std::string const& outputType, PointOnContext const& pointOnCtx, ExprType const& theExpr,
                       std::string const& measurePrefix,std::string const& measureEvaluatedFrom,
                       std::vector<std::pair<std::string,std::vector<double>>> & measuresStored )
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

            this->updateMeasureImpl( measureValues, outputType, measurePrefix, measureEvaluatedFrom, measuresStored );
        }


    template <typename T>
    void
    updateMeasureImpl( std::vector<T> const& measures, std::string const& outputType,
                       std::string const& measurePrefix, std::string const& measureEvaluatedFrom,
                       std::vector<std::pair<std::string,std::vector<double>>> & measuresStored )
        {
            if ( measures.empty() )
                return;

            int M = measures.front().rows();
            int N = measures.front().cols();
            if ( outputType == "values" )
            {
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
                            measuresStored.push_back( std::make_pair(measureName, std::vector<double>(1,val) ) );
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
                        measuresStored.push_back( std::make_pair( measureName, std::move(currentMeasureValues) ) );
                    }
                }
            }

        }

private :
    geometric_context_ptrtype M_geoContext;
    points_eval_desc_type M_pointsEvalDesc;
};

} // namespace FeelModels
} // namespace Feel


#endif

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

template <typename... SpacesType>
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


public :

    using spaces_tuple_type = hana::tuple<SpacesType...>;
    using fieldnames_space_tuple_type = std::decay_t<decltype( hana::transform( spaces_tuple_type{}, TransformSpaceToFieldNamesAndSpace{} ) ) >;
    using context_fields_tuple_type =  std::decay_t<decltype( hana::transform( fieldnames_space_tuple_type{}, TransformSpaceToContext{} ) ) >;

    MeasurePointsEvaluation( fieldnames_space_tuple_type const& fieldNamesSpaces )
        :
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

        }



    template <typename FieldTupleType>
    void
    eval( std::vector<ModelPostprocessPointPosition> const& allEvalPoints, std::map<std::string,double> & res, FieldTupleType const& fieldTuple )
        {
            for ( auto const& evalPoints : allEvalPoints )
            {
                // TODO : update point position if node has moved
            }

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

        }

private :
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

            if constexpr( std::is_same_v< typename Feel::decay_type<FieldType>::functionspace_type::Context, std::decay_t<decltype(*fectx)> > )
                {
                    auto itFindField = std::get<1>( x ).find( fieldName );
                    if ( itFindField != std::get<1>( x ).end() && !itFindField->second.empty() )
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

private :
    context_fields_tuple_type M_contextFields;
};

} // namespace FeelModels
} // namespace Feel


#endif

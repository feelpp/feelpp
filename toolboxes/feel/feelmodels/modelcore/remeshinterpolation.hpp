/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */
#pragma once


namespace Feel
{
namespace FeelModels
{

/**
 * @brief Interpolation Helper class for remeshing
 * @ingroup Fluid
 */
class RemeshInterpolation
{
public:
    using functionspace_base_ptrtype = std::shared_ptr<FunctionSpaceBase>;
    using sparse_matrix_ptrtype = std::shared_ptr<MatrixSparse<double>>;
    using vector_type = Vector<double>;
    using vector_ptrtype = std::shared_ptr<vector_type>;
    using datamap_ptrtype = datamap_ptr_t<>;
private:
    using key_functionspaces_type = std::pair<functionspace_base_ptrtype,functionspace_base_ptrtype>;
    using key_datamaps_type = std::tuple<std::string,datamap_ptrtype,datamap_ptrtype>;
    using matrix_interpolation_key_type = std::variant<key_functionspaces_type,key_datamaps_type>;
public:

    sparse_matrix_ptrtype matrixInterpolation( functionspace_base_ptrtype domainSpace, functionspace_base_ptrtype imageSpace ) const
        {
            auto itFind = M_matrixInterpolations.find( matrix_interpolation_key_type{std::make_pair(domainSpace, imageSpace)} );
            if ( itFind != M_matrixInterpolations.end() )
                return itFind->second;
            return sparse_matrix_ptrtype{};
        }

    sparse_matrix_ptrtype matrixInterpolation( std::string const& name, datamap_ptrtype const& domainDataMap, datamap_ptrtype const& imageDataMap ) const
        {
            return this->matrixInterpolation( std::make_tuple( name,domainDataMap, imageDataMap ) );
        }

    sparse_matrix_ptrtype matrixInterpolation( std::string const& name ) const
        {
            auto itFind = std::find_if(M_matrixInterpolations.begin(),M_matrixInterpolations.end(),
                                       [&name]( auto const& keyRelated ) {
                                           if ( const key_datamaps_type* dmRelated = std::get_if<key_datamaps_type>(&keyRelated.first) )
                                               return std::get<0>( *dmRelated ) == name;
                                           return false;
                                       });
            if ( itFind != M_matrixInterpolations.end() )
                return itFind->second;
            return sparse_matrix_ptrtype{};
        }

    template <typename DomainSpaceType, typename ImageSpaceType, typename RangeType>
    sparse_matrix_ptrtype computeMatrixInterpolation( std::shared_ptr<DomainSpaceType> domainSpace, std::shared_ptr<ImageSpaceType> imageSpace, RangeType const& range )
        {
            auto opI = opInterpolation(_domainSpace=domainSpace,
                                       _imageSpace=imageSpace,
                                       //_range=elements(support( imageSpace ) )
                                       _range=range );
            sparse_matrix_ptrtype matInterp = opI->matPtr();
            this->setMatrixInterpolation( domainSpace, imageSpace, matInterp );
            return matInterp;
        }
    template <typename DomainSpaceType, typename ImageSpaceType>
    sparse_matrix_ptrtype computeMatrixInterpolation( std::shared_ptr<DomainSpaceType> domainSpace, std::shared_ptr<ImageSpaceType> imageSpace )
        {
            return this->computeMatrixInterpolation( domainSpace,imageSpace,elements(support( imageSpace ) ) );
        }

    void setMatrixInterpolation( functionspace_base_ptrtype domainSpace, functionspace_base_ptrtype imageSpace, sparse_matrix_ptrtype mat )
        {
            M_matrixInterpolations[matrix_interpolation_key_type{std::make_pair(domainSpace,imageSpace)}] = mat;
        }

    void setMatrixInterpolation( std::string const& name, datamap_ptrtype domainDataMap, datamap_ptrtype imageDataMap, sparse_matrix_ptrtype mat )
        {
            M_matrixInterpolations[matrix_interpolation_key_type{std::make_tuple(name,domainDataMap,imageDataMap)}] = mat;
        }

    template <typename DomainElementType, typename ImageElementType>
    bool interpolate( DomainElementType const& u, ImageElementType & v ) const
        {
            auto domainSpace = unwrap_ptr( u ).functionSpace();
            auto imageSpace = unwrap_ptr( v ).functionSpace();
            auto matInterp = this->matrixInterpolation( domainSpace, imageSpace );
            if ( !matInterp )
                return false;
            matInterp->multVector( unwrap_ptr( u ),  unwrap_ptr( v ) );
            return true;
        }
    bool interpolate( std::string const& name, vector_ptrtype const& u, vector_ptrtype & v ) const
        {
            auto matInterp = this->matrixInterpolation( name );
            if ( !matInterp )
                return false;
            matInterp->multVector( unwrap_ptr( u ),  unwrap_ptr( v ) );
            return true;
        }

    bool interpolateBlockVector( std::string const& name, vector_ptrtype oldVecMonolithic, vector_ptrtype newVecMonolithic,  BlocksBaseVector<double> const& bvs, bool close = true ) const
        {
            auto thebackend = Feel::backend();
            auto itFindName = M_blockIndexToSpace.find( name );
            if ( itFindName == M_blockIndexToSpace.end() )
                return false;
            auto const& blockIndexToSpace = itFindName->second;

            auto const& old_dm = oldVecMonolithic->map();
            auto const& new_dm = newVecMonolithic->map();
            CHECK( old_dm.numberOfDofIdToContainerId() == new_dm.numberOfDofIdToContainerId() ) << "vectors not compatible";
            for ( size_type tag=0 ; tag<old_dm.numberOfDofIdToContainerId() ; ++tag )
            {
                auto matInterp = this->matrixInterpolation( blockIndexToSpace, tag );
                CHECK( matInterp ) << "missing block index or interpolation matrix";
#if 0
                auto oldBlockField = thebackend->newVector( matInterp->mapColPtr() );
                auto newBlockField = thebackend->newVector( matInterp->mapRowPtr() );
#else
                auto [mapDomainPtr,mapImagePtr] = this->dataMap( blockIndexToSpace, tag );
                auto oldBlockField = thebackend->newVector( mapDomainPtr );
                auto newBlockField = thebackend->newVector( mapImagePtr );
#endif

                bvs.setSubVector( *oldBlockField, *oldVecMonolithic, tag );

                matInterp->multVector( unwrap_ptr( oldBlockField ),  unwrap_ptr( newBlockField ) );

                bvs.setVector( *newVecMonolithic, *newBlockField, tag, false );
            }
            if ( close )
                newVecMonolithic->close();
            return true;
        }

    //! registering a the link between a \blockUndex (of block vector called \nameOfBlockVector) and interpolation matrix (represented by \domainSpace and \nameOfBlockVector)
    void registeringBlockIndex( std::string const& nameOfBlockVector, size_type blockIndex, functionspace_base_ptrtype domainSpace, functionspace_base_ptrtype imageSpace )
        {
            M_blockIndexToSpace[nameOfBlockVector].insert( { blockIndex, matrix_interpolation_key_type{std::make_pair(domainSpace,imageSpace)} } );
        }
    //! registering a the link between a \blockUndex (of block vector called \nameOfBlockVector) and interpolation matrix (called matInterpName)
    void registeringBlockIndex( std::string const& nameOfBlockVector, size_type blockIndex, std::string const& matInterpName )
        {
            auto itFind = std::find_if(M_matrixInterpolations.begin(),M_matrixInterpolations.end(),
                                       [&matInterpName]( auto const& keyRelated ) {
                                           if ( const key_datamaps_type* dmRelated = std::get_if<key_datamaps_type>(&keyRelated.first) )
                                               return std::get<0>( *dmRelated ) == matInterpName;
                                           return false;
                                       });
            CHECK( itFind != M_matrixInterpolations.end() ) << "no matrix interpolation with name" << matInterpName;
            auto const& datmapsRelated = std::get<key_datamaps_type>( itFind->first );
            M_blockIndexToSpace[nameOfBlockVector].insert( { blockIndex, matrix_interpolation_key_type{std::make_tuple(matInterpName,std::get<1>(datmapsRelated),std::get<2>(datmapsRelated))} } );
        }
private :
    sparse_matrix_ptrtype matrixInterpolation( key_datamaps_type const& datamaps ) const
        {
            auto itFind = M_matrixInterpolations.find( matrix_interpolation_key_type{datamaps} );
            if ( itFind != M_matrixInterpolations.end() )
                return itFind->second;
            return sparse_matrix_ptrtype{};
        }
    sparse_matrix_ptrtype matrixInterpolation( std::map<size_type,matrix_interpolation_key_type> const& blockIndexToSpace, size_type blockIndex ) const
        {
            auto itFindBlockIndex = blockIndexToSpace.find( blockIndex );
            if ( itFindBlockIndex == blockIndexToSpace.end() )
                return sparse_matrix_ptrtype{};
            auto const& keyRelated = itFindBlockIndex->second;
            if ( const key_functionspaces_type* spacesRelated = std::get_if<key_functionspaces_type>(&keyRelated) )
                return this->matrixInterpolation( spacesRelated->first, spacesRelated->second );
            else if ( const key_datamaps_type* datamapsRelated = std::get_if<key_datamaps_type>(&keyRelated) )
                return this->matrixInterpolation( *datamapsRelated );
            return sparse_matrix_ptrtype{};
        }

    std::pair<datamap_ptr_t<>,datamap_ptr_t<>> dataMap( std::map<size_type,matrix_interpolation_key_type> const& blockIndexToSpace, size_type blockIndex ) const
        {
            auto itFindBlockIndex = blockIndexToSpace.find( blockIndex );
            if ( itFindBlockIndex == blockIndexToSpace.end() )
                return {};
            auto const& keyRelated = itFindBlockIndex->second;
            if ( const key_functionspaces_type* spacesRelated = std::get_if<key_functionspaces_type>(&keyRelated) )
                return std::make_pair( spacesRelated->first->mapPtr(), spacesRelated->second->mapPtr() );
            else if ( const key_datamaps_type* datamapsRelated = std::get_if<key_datamaps_type>(&keyRelated) )
                return std::make_pair( std::get<1>( *datamapsRelated ), std::get<2>( *datamapsRelated ) );
            return {};
        }

private:
    std::map<matrix_interpolation_key_type, sparse_matrix_ptrtype > M_matrixInterpolations;
    std::map<std::string,std::map<size_type,matrix_interpolation_key_type> > M_blockIndexToSpace;
};


} // namespace FeelModels
} // namespace Feel

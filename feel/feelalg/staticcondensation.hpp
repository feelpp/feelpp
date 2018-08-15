/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 02 Aug 2016

 Copyright (C) 2016 Feel++ Consortium

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
#ifndef FEELPP_STATICCONDENSATION_HPP
#define FEELPP_STATICCONDENSATION_HPP 1

#include <future>


#include <boost/functional/hash.hpp>

#include <unordered_map>
#include <Eigen/Core>
#include <boost/hana/equal.hpp>
#include <boost/hana/integral_constant.hpp>
#include <boost/hana/length.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/feelio.hpp>
#include <feel/feeldiscr/traits.hpp>
#include <feel/feeltiming/tic.hpp>
#include <feel/feelmesh/iterator.hpp>

namespace Feel {

namespace detail {

//! this data structure stores element faces ids which can possibly split onto
//! two sets of face ids depending on the type of approximation used on the faces
class ElementFaces
{
public:
    ElementFaces() : M_has_same_trace(true), M_space_index(-1) {}
    ElementFaces( bool t ) : M_has_same_trace(t), M_space_index(-1) {}
    ElementFaces( std::vector<size_type>& f1 )
        : M_has_same_trace(true), M_space_index(-1), M_f1(f1) {}
    ElementFaces( std::vector<size_type>&& f1 )
        : M_has_same_trace(true), M_space_index(-1), M_f1(f1) {}
    ElementFaces( std::vector<size_type>& f1, std::vector<size_type>& f2, int si )
        : M_has_same_trace(false), M_space_index(si), M_f1(f1), M_f2(f2) {}
    ElementFaces( std::vector<size_type>&& f1, std::vector<size_type>&& f2, int si )
        :
        M_has_same_trace(false), M_space_index(si), M_f1(f1), M_f2(f2) {}

    bool hasSameTrace() const { return M_has_same_trace; }
    int spaceIndex() const { return M_space_index; }
    
    //! @return list of faces of first type
    std::vector<size_type> const& faces1() const { return M_f1; }
    //! @return list of faces of second type
    std::vector<size_type> const& faces2() const { return M_f2; }
    
private:
    
    bool M_has_same_trace;
    int M_space_index;
    std::vector<size_type> M_f1;
    std::vector<size_type> M_f2;
};

}

template<typename T>
class StaticCondensation
{
public:
    using block_index_t = std::pair<int,int>;
    using block_element_t = std::pair<size_type,size_type>;
    using value_type = T;

    StaticCondensation() = default;
    template<typename E, typename M_ptrtype, typename V_ptrtype>
    void condense( std::shared_ptr<StaticCondensation<T>> const& rhs, E& e, M_ptrtype& S, V_ptrtype& V,
                   std::enable_if_t<std::decay_t<E>::nspaces == 3>* = nullptr );
    
    template<typename E, typename M_ptrtype, typename V_ptrtype>
    void condense( std::shared_ptr<StaticCondensation<T>> const& rhs, E& e, M_ptrtype& S, V_ptrtype& V,
                   std::enable_if_t<std::decay_t<E>::nspaces == 4>* = nullptr );

    template<typename E>
    void localSolve( std::shared_ptr<StaticCondensation<T>> const& rhs, E& e,
                     std::enable_if_t<std::decay_t<E>::nspaces == 3>* = nullptr );
    
    template<typename E>
    void localSolve( std::shared_ptr<StaticCondensation<T>> const& rhs, E& e,
                     std::enable_if_t<std::decay_t<E>::nspaces == 4>* = nullptr );

    void addLocalMatrix( int* rows, int nrows,
                         int* cols, int ncols,
                         value_type* data,
                         size_type K = 0,
                         size_type K2 = invalid_size_type_value );

    void addLocalVector( int* rows, int nrows,
                         value_type* data,
                         size_type K = 0,
                         size_type K2 = invalid_size_type_value );


    template <typename Space1,typename Space2>
    void syncLocalMatrix( std::shared_ptr<Space1> const& rowSpace, std::shared_ptr<Space2> const& colSpace )
        {
            if ( rowSpace->worldComm().localSize() == 1 )
                return;

            std::map<size_type,rank_type> allghostIdsRow;
            for( auto const& gw : elements( rowSpace->mesh(), EntityProcessType::GHOST_ONLY ) )
            {
                //for ( auto it=rowSpace->mesh()->beginGhostElement(),en=rowSpace->mesh()->endGhostElement();it!=en;++it )
                auto const& g = unwrap_ref(gw);
                allghostIdsRow[ g.id() ] = g.processId();
            }
            std::map<size_type,rank_type> allghostIdsCol;
            //for ( auto it=colSpace->mesh()->beginGhostElement(),en=colSpace->mesh()->endGhostElement();it!=en;++it )
            for( auto const& gw : elements( colSpace->mesh(), EntityProcessType::GHOST_ONLY ) )
            {
                auto const& g = unwrap_ref(gw);
                allghostIdsCol[ g.id() ] = g.processId();
            }

            std::map< rank_type, std::set< std::tuple<size_type,size_type,size_type,size_type> > > localMatKeysToSynchronize;

            for ( auto const& localMat : this->M_local_matrices[this->M_block_rowcol] )
            {
                size_type K = localMat.first.first;
                size_type K2 = localMat.first.second;
                std::tuple<size_type,size_type,size_type,size_type> keySync = std::make_tuple( K,K2,K,K2 );

                //if ( rowSpace->mesh()->dimension() == rowSpace->mesh()->realDimension() )
                // {
                auto itFindId = allghostIdsRow.find( K );
                if ( itFindId != allghostIdsRow.end() )
                {
                    auto const& eltRow = rowSpace->mesh()->element( itFindId->first );
                    rank_type otherpid = eltRow.processId();
                    auto itFindIdCol = allghostIdsCol.find( K2 );
                    auto const& eltCol = (itFindIdCol != allghostIdsCol.end())?
                        colSpace->mesh()->element( itFindIdCol->first ) :
                        colSpace->mesh()->element( K2 );
                    CHECK( eltRow.idInOthersPartitions().find( otherpid ) != eltRow.idInOthersPartitions().end() ) << "aie no idInOtherPartition = " << eltRow.G();
                    CHECK( eltCol.idInOthersPartitions().find( otherpid ) != eltCol.idInOthersPartitions().end() ) << "aie no idInOtherPartition = " << eltCol.G();
                    std::get<2>( keySync ) = eltRow.idInOthersPartitions( otherpid );
                    std::get<3>( keySync ) = eltCol.idInOthersPartitions( otherpid );
                    localMatKeysToSynchronize[otherpid].insert( keySync );
                }
                else if ( rowSpace->mesh()->dimension() != rowSpace->mesh()->realDimension() )
                {
                    CHECK( rowSpace->mesh()->hasElement( K ) ) << "mesh doesnt have an element with id " << K;
                    auto const& eltRow = rowSpace->mesh()->element( K );
                    for ( auto const& otherProcess : eltRow.idInOthersPartitions() )
                    {
                        rank_type otherpid = otherProcess.first;

                        auto itFindIdCol = allghostIdsCol.find( K2 );
                        auto const& eltCol = ( itFindIdCol != allghostIdsCol.end() )?
                            colSpace->mesh()->element( itFindIdCol->first ) :
                            colSpace->mesh()->element( K2 );
                        auto itFindOtherPartitionCol = eltCol.idInOthersPartitions().find( otherpid );
                        if ( itFindOtherPartitionCol != eltCol.idInOthersPartitions().end() )
                        {
                            std::get<2>( keySync ) = otherProcess.second;
                            std::get<3>( keySync ) = itFindOtherPartitionCol->second;
                            localMatKeysToSynchronize[otherpid].insert( keySync );
                        }
                    }
                }

                // }
                //else if ( colSpace->mesh()->dimension() == colSpace->mesh()->realDimension() )
                // else
                // {
                auto itFindIdCol = allghostIdsCol.find( K2 );
                if ( itFindIdCol != allghostIdsCol.end() )
                {
                    auto const& eltCol = colSpace->mesh()->element( itFindIdCol->first );
                    rank_type otherpid = eltCol.processId();

                    auto itFindIdRow = allghostIdsRow.find( K );
                    auto const& eltRow = ( itFindIdRow != allghostIdsRow.end() )?
                        rowSpace->mesh()->element( itFindIdRow->first ) :
                        rowSpace->mesh()->element( K );
                    CHECK( eltCol.idInOthersPartitions().find( otherpid ) != eltCol.idInOthersPartitions().end() ) << "no idInOtherPartition("<<otherpid<<") with " << eltRow.G();
                    CHECK( eltRow.idInOthersPartitions().find( otherpid ) != eltRow.idInOthersPartitions().end() ) << "no idInOtherPartition("<<otherpid<<") with " << eltRow.G();
                    std::get<2>( keySync ) = eltRow.idInOthersPartitions( otherpid );
                    std::get<3>( keySync ) = eltCol.idInOthersPartitions( otherpid );
                    localMatKeysToSynchronize[otherpid].insert( keySync );
                }
                else if ( colSpace->mesh()->dimension() != colSpace->mesh()->realDimension() )
                {
                    CHECK( colSpace->mesh()->hasElement( K2 ) ) << "mesh doesnt have an element with id " << K;
                    auto const& eltCol = colSpace->mesh()->element( K2 );
                    for ( auto const& otherProcess : eltCol.idInOthersPartitions() )
                    {
                        rank_type otherpid = otherProcess.first;
                        auto itFindIdRow = allghostIdsRow.find( K );
                        auto const& eltRow = ( itFindIdRow != allghostIdsRow.end() )?
                            rowSpace->mesh()->element( itFindIdRow->first ) :
                            rowSpace->mesh()->element( K );
                        auto itFindOtherPartitionRow = eltRow.idInOthersPartitions().find( otherpid );
                        if ( itFindOtherPartitionRow != eltRow.idInOthersPartitions().end() )
                        {
                            std::get<2>( keySync ) = itFindOtherPartitionRow->second;
                            std::get<3>( keySync ) = otherProcess.second;
                            localMatKeysToSynchronize[otherpid].insert( keySync );
                        }
                    }
                }

                // }
            }


#undef FEELPP_STATICCONDENSATION_DETECT_LOCALDOF_PERMUTATIONS
#define FEELPP_STATICCONDENSATION_DETECT_LOCALDOF_PERMUTATIONS 0

#if FEELPP_STATICCONDENSATION_DETECT_LOCALDOF_PERMUTATIONS
            std::map<rank_type,std::vector<boost::tuple<size_type,size_type,local_matrix_t,std::vector<size_type>,std::vector<size_type> > > > dataToSend;
            std::map<rank_type,std::vector<boost::tuple<size_type,size_type,local_matrix_t,std::vector<size_type>,std::vector<size_type> > > > dataToRecv;
#else
            std::map<rank_type,std::vector<boost::tuple<size_type,size_type,local_matrix_t> > > dataToSend;
            std::map<rank_type,std::vector<boost::tuple<size_type,size_type,local_matrix_t> > > dataToRecv;
#endif

            std::set<rank_type> neighborSubdomainsRowCol;
            for ( rank_type neighborPid : rowSpace->mesh()->neighborSubdomains() )
                neighborSubdomainsRowCol.insert( neighborPid );
            for ( rank_type neighborPid : colSpace->mesh()->neighborSubdomains() )
                neighborSubdomainsRowCol.insert( neighborPid );

            // update dataToSend
            for ( auto const& dataByProcess : localMatKeysToSynchronize )
            {
                rank_type pid = dataByProcess.first;
                for ( auto const& dataSync : dataByProcess.second )
                {
                    auto key = std::make_pair(std::get<0>(dataSync),std::get<1>(dataSync));
                    local_matrix_t localMat = this->M_local_matrices.find(this->M_block_rowcol)->second.find(key)->second;

#if FEELPP_STATICCONDENSATION_DETECT_LOCALDOF_PERMUTATIONS
                    auto localDofInEltRow = rowSpace->dof()->localDof( key.first/*K*/ );
                    auto localDofInEltCol = colSpace->dof()->localDof( key.second/*K2*/ );
                    if ( localDofInEltRow.first == localDofInEltRow.second || localDofInEltCol.first == localDofInEltCol.second )
                        continue;
                    dataToSend[pid].push_back( boost::make_tuple( std::get<2>( dataSync ),std::get<3>( dataSync ), localMat,
                                                                  rowSpace->dof()->getIndicesOnGlobalCluster( key.first ),
                                                                  colSpace->dof()->getIndicesOnGlobalCluster( key.second ) ) );
#else
                    dataToSend[pid].push_back( boost::make_tuple( std::get<2>( dataSync ),std::get<3>( dataSync ), localMat ) );
#endif
                }
            }
            // send/recv
            int neighborSubdomains = neighborSubdomainsRowCol.size();
            int nbRequest = 2*neighborSubdomains;
            mpi::request * reqs = new mpi::request[nbRequest];
            int cptRequest=0;
            WorldComm const& worldComm = rowSpace->mesh()->worldComm();
            for ( rank_type neighborPid : neighborSubdomainsRowCol )
            {
                reqs[cptRequest++] = worldComm.localComm().isend( neighborPid , 0, dataToSend[neighborPid] );
                reqs[cptRequest++] = worldComm.localComm().irecv( neighborPid , 0, dataToRecv[neighborPid] );
            }
            mpi::wait_all(reqs, reqs + nbRequest);

            // update local matrix
            for ( auto const& dataRecvByProcess : dataToRecv )
            {
                rank_type pid = dataRecvByProcess.first;
                for ( auto const& dataSync : dataRecvByProcess.second )
                {
                    // std::cout << "recv mat\n " << boost::get<2>(dataSync) << "\n";
                    auto key = std::make_pair(boost::get<0>(dataSync),boost::get<1>(dataSync));

#if FEELPP_STATICCONDENSATION_DETECT_LOCALDOF_PERMUTATIONS
                    auto localDofInEltRow = rowSpace->dof()->localDof( key.first/*K*/ );
                    auto localDofInEltCol = colSpace->dof()->localDof( key.second/*K2*/ );
                    if ( localDofInEltRow.first == localDofInEltRow.second || localDofInEltCol.first == localDofInEltCol.second )
                        continue;

                    auto const& locIndicesFromCurrentProcRow = rowSpace->dof()->getIndicesOnGlobalCluster( key.first );
                    auto const& locIndicesFromCurrentProcCol = colSpace->dof()->getIndicesOnGlobalCluster( key.second );
                    auto const& locIndicesFromOtherProcRow = boost::get<3>(dataSync);
                    auto const& locIndicesFromOtherProcCol = boost::get<4>(dataSync);
                    for (int locDofCurrentProc = 0; locDofCurrentProc<locIndicesFromCurrentProcRow.size(); ++locDofCurrentProc )
                    {
                        size_type gcdof = locIndicesFromCurrentProcRow[locDofCurrentProc];
                        for (int locDofOtherProc = 0; locDofOtherProc<locIndicesFromOtherProcRow.size(); ++locDofOtherProc )
                        {
                            if ( locIndicesFromOtherProcRow[locDofOtherProc] == gcdof )
                            {
                                if ( locDofCurrentProc != locDofOtherProc )
                                    std::cout << "["<<Environment::worldComm().rank() << "] FIND permutation in local matrix (row) : "
                                              << locDofCurrentProc << " and " << locDofOtherProc << "\n";
                                break;
                            }
                        }
                    }
                    // BK.col(4).swap(BK.col(5));
                    // CK.row(4).swap(CK.row(5));
#endif

                    auto entry = this->M_local_matrices[this->M_block_rowcol].find(key);
                    if ( entry == this->M_local_matrices[this->M_block_rowcol].end() )
                        this->M_local_matrices[this->M_block_rowcol][key] = boost::get<2>(dataSync);
                    else
                        this->M_local_matrices[this->M_block_rowcol][key] += boost::get<2>(dataSync);
                }
            }

        } // syncLocalMatrix
    /**
     * get current block
     */
    void block( int row, int col )
        {
            //M_block_nrow = nrow;
            //M_block_ncol = ncol;
            M_block_rowcol = std::make_pair(row,col);
        }
    void block( int row )
        {
            M_block_row = row;
        }
    //!
    //! zero out the local vectors and matrices
    //! @note the allocated memory is preserved
    //!
    void zero()
        {
            this->zeroMatrix();
            this->zeroVector();
        }
    //!
    //! zero out the set of local matrices
    //! @note the allocated memory is preserved    
    //!
    void zeroMatrix()
        {
            std::for_each( M_local_matrices.begin(), M_local_matrices.end(),
                           []( auto& m )
                           {
                               for( auto& e : m.second )
                               {
                                   e.second.setZero();
                               }
                           } );
                           
        }
    //!
    //! zero out the block matrix @arg n1,n2
    //! @note the allocated memory is preserved    
    //!
    void zero( int n1, int n2 )
        {
            for( auto& e: M_local_matrices[std::make_pair(n1,n2)] )
            {
                e.second.setZero();
            }
        }
    //!
    //! zero out the set of local vectors
    //! @note the allocated memory is preserved    
    //!
    void zeroVector()
        {
            std::for_each( M_local_vectors.begin(), M_local_vectors.end(),
                           []( auto& m )
                           {
                               for( auto& e : m.second )
                               {
                                   e.second.setZero();
                               }
                           } );
        }
    //!
    //! zero out the block vector @arg n1
    //!
    void zero( int n1 )
        {
            for( auto& e: M_local_vectors[n1] )
            {
                e.second.setZero();
            }
        }
    void setDim4( int dim4 ) { M_dim4 = dim4; }
    
    using local_vector_t = Eigen::Matrix<value_type,Eigen::Dynamic,1>;
    using local_matrix_t = Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    using local_index_t = Eigen::Matrix<int,Eigen::Dynamic,1>;
    using raw_vector_map_t = Eigen::Map<local_vector_t>;
    using raw_matrix_map_t = Eigen::Map<local_matrix_t>;
    using raw_index_map_t = Eigen::Map<local_index_t>;

#if 0
    local_vector_t& localVector( size_type K ) { return M_local_vectors[this->M_block_row][K]; }
    local_vector_t const& localVector( size_type K ) const { return M_local_vectors.at(this->M_block_row).at(K); }

    local_index_t& localVRow( size_type K ) { return M_local_vrows[this->M_block_row][K]; }
    local_index_t const& localVRow( size_type K ) const { return M_local_vrows.at(this->M_block_row).at(K); }

    local_matrix_t& localMatrix( block_element_t const& K ) { return M_local_matrices[this->M_block_rowcol][K]; }
    local_matrix_t const& localMatrix( block_element_t const& K ) const { return M_local_matrices.at(this->M_block_rowcol).at(K); }

    local_index_t& localRow( block_element_t const& K ) { return M_local_rows[this->M_block_rowcol][K]; }
    local_index_t const& localRow( block_element_t const& K ) const { return M_local_rows.at(this->M_block_rowcol).at(K); }
    local_index_t& localCol( block_element_t const& K ) { return M_local_cols[this->M_block_rowcol][K]; }
    local_index_t const& localCol( block_element_t const& K ) const { return M_local_cols.at(this->M_block_rowcol).at(K); }
#endif
private:
    /**
     * unassembled view of the matrix
     */
    std::unordered_map<int,std::unordered_map<size_type,local_vector_t>> M_local_vectors;
    std::unordered_map<int,std::unordered_map<size_type,local_index_t>> M_local_vrows;
    std::unordered_map<block_index_t,std::unordered_map<block_element_t,local_matrix_t,boost::hash<block_element_t>>,boost::hash<block_index_t>> M_local_matrices;
    std::unordered_map<block_index_t,std::unordered_map<block_element_t,local_index_t,boost::hash<block_element_t>>,boost::hash<block_index_t>> M_local_rows;
    std::unordered_map<block_index_t,std::unordered_map<block_element_t,local_index_t,boost::hash<block_element_t>>,boost::hash<block_index_t>> M_local_cols;

    std::unordered_map<int,local_matrix_t> M_AinvB;
    std::unordered_map<int,local_vector_t> M_AinvF;

    using dK_iterator_type = std::unordered_map<int,Feel::detail::ElementFaces>::iterator;
    std::unordered_map<int,Feel::detail::ElementFaces> M_dK;
    block_index_t M_block_rowcol;
    int M_block_row;
    int M_dim4;
    std::mutex mutex_add_v, mutex_add_m;

};


template<typename A00_t,typename A01_t,typename A10_t, typename M2, typename E1, typename E2>
void
extractBlock( A00_t const& a00, A01_t const& a01, A10_t const& a10,
              M2& ak, E1 const& e1, E2 const& e2, tensor2symm_true )
{
    uint16_type N0 = e1.dof()->nRealLocalDof( false );
    uint16_type N0c = e1.dof()->nLocalDof( true );
    uint16_type N1 = e2.dof()->nLocalDof();
    // A00
    for( int c1e1 = 0; c1e1 < e1.nComponents1; ++c1e1 )
    {
        for( int c2e1 = c1e1+1; c2e1 < e1.nComponents1; ++c2e1 )
        {
            const int k1 = Feel::detail::symmetricIndex(c1e1,c2e1,e1.nComponents1);
            for( int c1e2 = 0; c1e2 < e1.nComponents1; ++c1e2 )
            {
                for( int c2e2 = c1e2+1; c2e2 < e1.nComponents1; ++c2e2 )
                {
                    const int k2 = Feel::detail::symmetricIndex(c1e2,c2e2,e1.nComponents1);
                    ak.block( N0c*k1, N0c*k2, N0c, N0c ) = a00.block( N0c*(c1e1*e1.nComponents1+c2e1),
                                                                      N0c*(c1e2*e2.nComponents1+c2e2), N0c, N0c );
                    ak.block( N0c*k1, N0c*k2, N0c, N0c ) += a00.block( N0c*(c2e1*e1.nComponents1+c1e1),
                                                                       N0c*(c1e2*e2.nComponents1+c2e2), N0c, N0c );
                    ak.block( N0c*k1, N0c*k2, N0c, N0c ) += a00.block( N0c*(c2e1*e1.nComponents1+c1e1),
                                                                       N0c*(c2e2*e2.nComponents1+c1e2), N0c, N0c );
                    ak.block( N0c*k1, N0c*k2, N0c, N0c ) += a00.block( N0c*(c1e1*e1.nComponents1+c2e1),
                                                                       N0c*(c2e2*e2.nComponents1+c1e2), N0c, N0c );

                }
                const int k2 = Feel::detail::symmetricIndex(c1e2,c1e2,e2.nComponents1);
                ak.block( N0c*k1, N0c*k2, N0c, N0c ) = a00.block( N0c*(c1e1*e1.nComponents1+c2e1),
                                                                  N0c*(c1e2*e2.nComponents1+c1e2), N0c, N0c );
                ak.block( N0c*k1, N0c*k2, N0c, N0c ) += a00.block( N0c*(c2e1*e1.nComponents1+c1e1),
                                                                   N0c*(c1e2*e2.nComponents1+c1e2), N0c, N0c );
            }
        }
        const int k1 = Feel::detail::symmetricIndex(c1e1,c1e1,e1.nComponents1);
        for( int c1e2 = 0; c1e2 < e2.nComponents1; ++c1e2 )
        {
                for( int c2e2 = c1e2+1; c2e2 < e2.nComponents1; ++c2e2 )
                {
                    const int k2 = Feel::detail::symmetricIndex(c1e2,c2e2,e2.nComponents1);
                    ak.block( N0c*k1, N0c*k2, N0c, N0c ) = a00.block( N0c*(c1e1*e1.nComponents1+c1e1),
                                                                      N0c*(c1e2*e2.nComponents1+c2e2), N0c, N0c );
                    ak.block( N0c*k1, N0c*k2, N0c, N0c ) += a00.block( N0c*(c1e1*e1.nComponents1+c1e1),
                                                                       N0c*(c2e2*e2.nComponents1+c1e2), N0c, N0c );

                }
                const int k2 = Feel::detail::symmetricIndex(c1e2,c1e2,e2.nComponents1);
                ak.block( N0c*k1, N0c*k2, N0c, N0c ) += a00.block( N0c*(c1e1*e1.nComponents1+c1e1),
                                                                   N0c*(c1e2*e2.nComponents1+c1e2), N0c, N0c );
        }

    }
    // A01
    for( int c1e1 = 0; c1e1 < e1.nComponents1; ++c1e1 )
    {
        for( int c2e1 = c1e1+1; c2e1 < e1.nComponents1; ++c2e1 )
        {
            const int k1 = Feel::detail::symmetricIndex(c1e1,c2e1,e1.nComponents1);
            ak.block( N0c*k1, N0, N0c, N1 ) = a01.block( N0c*(c1e1*e1.nComponents1+c2e1), 0, N0c, N1 );
            ak.block( N0c*k1, N0, N0c, N1 ) += a01.block( N0c*(c2e1*e1.nComponents1+c1e1), 0, N0c, N1 );

            ak.block( N0, N0c*k1, N1, N0c ) = a10.block( 0, N0c*(c1e1*e1.nComponents1+c2e1), N1, N0c );
            ak.block( N0, N0c*k1, N1, N0c ) += a10.block( 0, N0c*(c2e1*e1.nComponents1+c1e1), N1, N0c );
        }
        const int k1 = Feel::detail::symmetricIndex(c1e1,c1e1,e1.nComponents1);

        ak.block( N0c*k1, N0, N0c, N1 ) = a01.block( N0c*(c1e1*e1.nComponents1+c1e1), 0, N0c, N1 );
        ak.block( N0, N0c*k1, N1, N0c ) += a10.block( 0, N0c*(c1e1*e1.nComponents1+c1e1), N1, N0c );
    }
}

template<typename A00_t,typename A01_t,typename A10_t, typename M2, typename E1, typename E2>
void
extractBlock( A00_t const& a00, A01_t const& a01, A10_t const& a10,
              M2& ak, E1 const& e1, E2 const& e2, tensor2symm_false )
{
    int N0 = e1.dof()->nLocalDof();
    int N1 = e2.dof()->nLocalDof();
    ak.topLeftCorner(N0, N0 ) = a00;
    ak.bottomLeftCorner(N1, N0 ) = a10;
    ak.topRightCorner(N0, N1 ) = a01;
}


template<typename A00_t,typename A01_t,typename A10_t, typename M2, typename E1, typename E2>
void
extractBlock( A00_t const& a00, A01_t const& a01, A10_t const& a10,
              M2& ak, E1 const& e1, E2 const& e2 )
{
    return extractBlock( a00, a01, a10, ak, e1, e2, is_tensor2symm_field<E1>{} );
}

template<typename A02_t, typename A20_t, typename Key1_t, typename Key2_t, typename BK_t, typename CK_t, typename E1, typename E3>
void
extractBlock( A02_t const& A02K, Key1_t const& key2,
              A20_t const& A20K, Key2_t const& key3,
              BK_t& BK, CK_t& CK,
              int n,
              E1 const& e1, E3 const& e3, int start = 0,
              std::enable_if_t<is_tensor2symm_field_v<E1>>* = nullptr )
{
    uint16_type N0 = e1.dof()->nRealLocalDof( false );
    uint16_type N0c = e1.dof()->nLocalDof( true );
    int N2 = e3.dof()->nLocalDof();

    if ( A02K.count(key2) )
    {
        for( int c1e1 = 0; c1e1 < e1.nComponents1; ++c1e1 )
        {
            for( int c2e1 = c1e1+1; c2e1 < e1.nComponents1; ++c2e1 )
            {
                const int k1 = Feel::detail::symmetricIndex(c1e1,c2e1,e1.nComponents1);
                BK.block(N0c*k1, start+n*N2, N0c, N2 ) += A02K.at(key2).block( N0c*(c1e1*e1.nComponents1+c2e1), 0, N0c, N2 );
                BK.block(N0c*k1, start+n*N2, N0c, N2 ) += A02K.at(key2).block( N0c*(c2e1*e1.nComponents1+c1e1), 0, N0c, N2 );
            }
            const int k1 = Feel::detail::symmetricIndex(c1e1,c1e1,e1.nComponents1);
            BK.block(N0c*k1, start+n*N2, N0c, N2 ) += A02K.at(key2).block( N0c*(c1e1*e1.nComponents1+c1e1), 0, N0c, N2 );
        }
    }

    if ( A20K.count(key3) )
    {
        for( int c1e1 = 0; c1e1 < e1.nComponents1; ++c1e1 )
        {
            for( int c2e1 = c1e1+1; c2e1 < e1.nComponents1; ++c2e1 )
            {
                const int k1 = Feel::detail::symmetricIndex(c1e1,c2e1,e1.nComponents1);
                CK.block(start+n*N2, N0c*k1, N2, N0c ) += A20K.at(key3).block(0, N0c*(c1e1*e1.nComponents1+c2e1), N2, N0c );
                CK.block(start+n*N2, N0c*k1, N2, N0c ) += A20K.at(key3).block(0, N0c*(c2e1*e1.nComponents1+c1e1), N2, N0c );
            }
            const int k1 = Feel::detail::symmetricIndex(c1e1,c1e1,e1.nComponents1);
            CK.block(start+n*N2, N0c*k1, N2, N0c ) += A20K.at(key3).block(0, N0c*(c1e1*e1.nComponents1+c1e1), N2, N0c );
        }
    }

}



template<typename A02_t, typename A20_t, typename Key1_t, typename Key2_t, typename BK_t, typename CK_t, typename E1, typename E3>
void
extractBlock( A02_t const& A02K, Key1_t const& key2,
              A20_t const& A20K, Key2_t const& key3,
              BK_t& BK, CK_t& CK, int n,
              E1 const& e1, E3 const& e3, int start = 0,
              std::enable_if_t<!is_tensor2symm_field_v<E1>>* = nullptr )
{
    int N0 = e1.dof()->nLocalDof();
    int N2 = e3.dof()->nLocalDof();
    if ( A02K.count(key2) )
        BK.block(0, start+n*N2, N0, N2 ) = A02K.at(key2).block(0,0,N0,N2);

    if ( A20K.count(key3) )
        CK.block(start+n*N2, 0, N2, N0 ) = A20K.at(key3).block(0,0,N2,N0);

}

template<typename F0K_t, typename FK_t, typename E1>
void
extractBlock( F0K_t const& F0K, size_type K,
              FK_t& FK, E1 const& e1,
              std::enable_if_t<is_tensor2symm_field_v<E1>>* = nullptr )
{
    int N0 = e1.dof()->nRealLocalDof();
    int N0c = e1.dof()->nRealLocalDof(true);
    if ( F0K.count(K) )
    {
        for( int c1e1 = 0; c1e1 < e1.nComponents1; ++c1e1 )
        {
            for( int c2e1 = c1e1+1; c2e1 < e1.nComponents1; ++c2e1 )
            {
                const int k1 = Feel::detail::symmetricIndex(c1e1,c2e1,e1.nComponents1);
                FK.segment(N0c*k1, N0c) += F0K.at(K).segment( N0c*(c1e1*e1.nComponents1+c2e1), N0c );
                FK.segment(N0c*k1, N0c) += F0K.at(K).segment( N0c*(c2e1*e1.nComponents1+c1e1), N0c );
            }
            const int k1 = Feel::detail::symmetricIndex(c1e1,c1e1,e1.nComponents1);
            FK.segment(N0c*k1, N0c) += F0K.at(K).segment( N0c*(c1e1*e1.nComponents1+c1e1), N0c );
        }
    }
}
template<typename F0K_t, typename FK_t, typename E1>
void
extractBlock( F0K_t const& F0K, size_type K,
              FK_t& FK, E1 const& e1,
              std::enable_if_t<!is_tensor2symm_field_v<E1>>* = nullptr )
{
    int N0 = e1.dof()->nLocalDof();
    if ( F0K.count(K) )
        FK.head(N0) = F0K.at(K);
}

template<typename T>
template<typename E, typename M_ptrtype, typename V_ptrtype>
void
StaticCondensation<T>::condense( std::shared_ptr<StaticCondensation<T>> const& rhs, E& e, M_ptrtype& S, V_ptrtype& V,
                                 std::enable_if_t<std::decay_t<E>::nspaces == 3>* )
{
    using Feel::cout;

    auto& e1 = e(0_c);
    auto& e2 = e(1_c);
    auto& e3 = e(2_c);
    auto const& A00K = M_local_matrices[std::make_pair(0,0)];
    LOG(INFO) << "A00K.size=" << A00K.size();
    auto const& A10K = M_local_matrices[std::make_pair(1,0)];
    LOG(INFO) << "A10K.size=" << A10K.size();
    auto const& A01K = M_local_matrices[std::make_pair(0,1)];
    LOG(INFO) << "A01K.size=" << A01K.size();
    auto const& A11K = M_local_matrices[std::make_pair(1,1)];
    LOG(INFO) << "A11K.size=" << A11K.size();
    auto const& A02K = M_local_matrices[std::make_pair(0,2)];
    LOG(INFO) << "A02K.size=" << A02K.size();
    auto const& A12K = M_local_matrices[std::make_pair(1,2)];
    LOG(INFO) << "A12K.size=" << A12K.size();
    auto const& A20K = M_local_matrices[std::make_pair(2,0)];
    LOG(INFO) << "A20K.size=" << A20K.size();
    auto const& A21K = M_local_matrices[std::make_pair(2,1)];
    LOG(INFO) << "A21K.size=" << A21K.size();
    auto const& A22K = M_local_matrices[std::make_pair(2,2)];
    LOG(INFO) << "A22K.size=" << A22K.size();

    auto const& F0K = rhs->M_local_vectors[0];
    LOG(INFO) << "F0K.size=" << F0K.size();
    auto const& F1K = rhs->M_local_vectors[1];
    LOG(INFO) << "F1K.size=" << F1K.size();
    auto const& F2K = rhs->M_local_vectors[2];
    LOG(INFO) << "F2K.size=" << F2K.size();

    auto it = A00K.begin();
    auto en = A00K.end();

    int N00 = e1.dof()->nLocalDof();
    int N0 = e1.dof()->nRealLocalDof();
    int N1 = e2.dof()->nLocalDof();
    int N = N0+N1;
    int N2 = e3.dof()->nLocalDof();
    int N3 = N2*e1.mesh()->numLocalTopologicalFaces();
    cout << "[staticcondensation] N=" << N << " N0=" << N0 << " N1=" << N1 << " N2=" << N2 << " N3=" << N3 << " ntf=" << e1.mesh()->numLocalTopologicalFaces()<< std::endl;
    local_matrix_t AK( N, N ),A00(N0,N0),A01(N0,N1),A10(N1,N0), A11(N1,N1), A20(N3,N0), A21(N3,N1), A22(N3,N3);
    local_matrix_t BK( N, N3 );
    local_matrix_t CK( N3, N );
    local_matrix_t DK( N3, N3 );
    local_vector_t FK( N );
    local_vector_t DKF( N3 );
    local_vector_t F2( N3 );
    //ocal_matrix_t Aldlt( N, N );
    for( ; it != en ; ++it )
    {
        tic();
        auto key = it->first;
        size_type K = key.first;

        DVLOG(2) << "======= Key=" << key ;

        DVLOG(2) << "A00K=" << A00K.at(key);
        DVLOG(2) << "A01K=" << A01K.at(key);
        DVLOG(2) << "A10K=" << A10K.at(key);

        AK = local_matrix_t::Zero( N, N );
        BK = local_matrix_t::Zero( N, N3 );
        CK = local_matrix_t::Zero( N3, N );


        extractBlock( A00K.at(key), A01K.at(key), A10K.at(key), AK, e1, e2 );
        AK.bottomRightCorner(N1, N1 ) = A11K.at(key);

        A22 = local_matrix_t::Zero(N3,N3);
        A20 = local_matrix_t::Zero(N3,N0);
        A21 = local_matrix_t::Zero(N3,N1);
        CK = local_matrix_t::Zero(N3,N);
        FK = local_vector_t::Zero(N);
        F2 = local_vector_t::Zero(N3);

        //std::cout << "AK=\n" << AK << std::endl;
        //std::cout << "AK-AK^T=\n" << (AK-AK.transpose()) << std::endl;
        // dK contains the set of faces ids in the submesh associated to the boundary of K
        auto dK = M_dK.emplace(key.first, Feel::detail::ElementFaces(e3.mesh()->meshToSubMesh( e1.mesh()->element(key.first).facesId()).first) );
        
        int n = 0;
        std::for_each( dK.first->second.faces1().begin(), dK.first->second.faces1().end(), [&]( auto dKi )
                       {
                           auto key2 = std::make_pair(key.first, dKi );
                           auto key3 = std::make_pair(dKi,key.first);
                           auto key4 = std::make_pair(dKi, dKi);

                           DVLOG(2) << "A02.count(" << key2 << ")"  << A02K.count(key2);
                           DVLOG(2) << "A20.count(" << key3 << ")"  << A20K.count(key3);
                           DVLOG(2) << "A12.count(" << key2 << ")"  << A12K.count(key2);
                           DVLOG(2) << "A21.count(" << key3 << ")"  << A21K.count(key3);
                           DVLOG(2) << "A22.count(" << key4 << ")"  << A22K.count(key4);
                           DVLOG(2) << "F2.count(" << dKi << ")"  << F2K.count(dKi);

                           extractBlock( A02K, key2, A20K, key3, BK, CK, n, e1, e3 );

                           if ( A12K.count(key2) )
                               BK.block(N0,n*N2, N1, N2 ) = A12K.at(key2);

                           if ( A21K.count(key3) )
                               CK.block(n*N2, N0, N2, N1 ) = A21K.at(key3);

                           if ( A22K.count(key4) )
                               A22.block(n*N2, n*N2, N2, N2 ) = A22K.at(key4);

                           if ( F2K.count(dKi) )
                           {
                               F2.segment(n*N2,N2)=F2K.at(dKi);
                           }

                           ++n;
                       } );

        extractBlock( F0K, K, FK, e1 );
        FK.tail(N1) = F1K.at(K);

#if 1
        auto Aldlt = AK.ldlt();
        ////LOG(INFO) << "Aldlt=" << Aldlt;
#else
        auto Aldlt = AK.lu();
#endif
        auto AinvB = Aldlt.solve( BK );
        auto AinvF = Aldlt.solve( FK );

        // local assemble DK and DKF
        DK=-CK*AinvB+A22 ;
        DKF=-CK*AinvF+F2;

        M_AinvB.emplace( K, AinvB );
        M_AinvF.emplace( K, AinvF );
        toc("sc.condense.localassembly", FLAGS_v>0);
        tic();
        auto dofs = e3.dofs(dK.first->second.faces1());

        S.addMatrix( dofs.data(), dofs.size(), dofs.data(), dofs.size(), DK.data(), invalid_size_type_value, invalid_size_type_value );
        V.addVector( dofs.data(), dofs.size(), DKF.data(), invalid_size_type_value, invalid_size_type_value );
        toc("sc.condense.globalassembly", FLAGS_v>0);
    }
}

template<typename T>
template<typename E, typename M_ptrtype, typename V_ptrtype>
void
StaticCondensation<T>::condense( std::shared_ptr<StaticCondensation<T>> const& rhs, E &e, M_ptrtype& S, V_ptrtype& V,
                                 std::enable_if_t<std::decay_t<E>::nspaces == 4>* ) 
{
    using Feel::cout;
    auto& e1 = e(0_c);
    auto& e2 = e(1_c);
    auto& e3 = e(2_c);
    auto const& A00K = M_local_matrices[std::make_pair(0,0)];
    //LOG(INFO) << "A00K.size=" << A00K.size();
    auto const& A10K = M_local_matrices[std::make_pair(1,0)];
    //LOG(INFO) << "A10K.size=" << A10K.size();
    auto const& A01K = M_local_matrices[std::make_pair(0,1)];
    //LOG(INFO) << "A01K.size=" << A01K.size();
    auto const& A11K = M_local_matrices[std::make_pair(1,1)];
    //LOG(INFO) << "A11K.size=" << A11K.size();
    auto const& A02K = M_local_matrices[std::make_pair(0,2)];
    //LOG(INFO) << "A02K.size=" << A02K.size();
    auto const& A12K = M_local_matrices[std::make_pair(1,2)];
    //LOG(INFO) << "A12K.size=" << A12K.size();
    auto const& A20K = M_local_matrices[std::make_pair(2,0)];
    //LOG(INFO) << "A20K.size=" << A20K.size();
    auto const& A21K = M_local_matrices[std::make_pair(2,1)];
    //LOG(INFO) << "A21K.size=" << A21K.size();
    auto const& A22K = M_local_matrices[std::make_pair(2,2)];
    //LOG(INFO) << "A22K.size=" << A22K.size();
    auto const& A30K = M_local_matrices[std::make_pair(3,0)];
    auto const& A03K = M_local_matrices[std::make_pair(0,3)];
    //LOG(INFO) << "A20K.size=" << A20K.size();
    auto const& A31K = M_local_matrices[std::make_pair(3,1)];
    auto const& A13K = M_local_matrices[std::make_pair(1,3)];
    //LOG(INFO) << "A21K.size=" << A21K.size();
    auto const& A33K = M_local_matrices[std::make_pair(3,3)];
    //LOG(INFO) << "A22K.size=" << A22K.size();
    auto A4K = [this]( int i, int j ) { return this->M_local_matrices[std::make_pair(i,j)]; };
    auto const& F0K = rhs->M_local_vectors[0];
    //LOG(INFO) << "F0K.size=" << F0K.size();
    auto const& F1K = rhs->M_local_vectors[1];
    //LOG(INFO) << "F1K.size=" << F1K.size();
    auto const& F2K = rhs->M_local_vectors[2];
    auto const& F3K = rhs->M_local_vectors[3];
    //LOG(INFO) << "F2K.size=" << F2K.size();
    auto F4K = [&rhs]( int i ) { return rhs->M_local_vectors[i]; };

    auto it = A00K.begin();
    auto en = A00K.end();

    int N00 = e1.dof()->nLocalDof();
    int N0 = e1.dof()->nRealLocalDof();
    int N1 = e2.dof()->nLocalDof();
    int N = N0+N1;
    int N2 = e3.dof()->nLocalDof();
    int N3 = N2*e1.mesh()->numLocalTopologicalFaces();

    int N41 = N2*(e1.mesh()->numLocalTopologicalFaces());
    int N42 = 0;
    
    if ( e.functionSpace(3_c)->numberOfSpaces() > 0 )
    {
        N41 = N2*(e1.mesh()->numLocalTopologicalFaces()-1);
        N42 = e(3_c,0).dof()->nLocalDof();
    }
    int N4 = N41+N42;
    cout << "[staticcondensation] N=" << N << " N0=" << N0 << " N1=" << N1 << " N2=" << N2 << " N3=" << N3 << " N4=" << N4 << "(" << N41 << "+" << N42 << ")" << " ntf=" << e1.mesh()->numLocalTopologicalFaces()<< std::endl;
    local_matrix_t AK( N, N ),A00(N0,N0),A01(N0,N1),A10(N1,N0), A11(N1,N1), A20(N3,N0), A21(N3,N1), A22(N3,N3);
    local_matrix_t A03(N0,N42),A30(N42,N0), A13(N1,N42),A31(N42,N1), A23(N2,N42),A32(N42,N2), A33(N42,N42);
    local_matrix_t BK( N, N3 ),BK3( N,N4);
    local_matrix_t CK( N3, N ),CK3( N4,N);
    local_matrix_t DK( N3, N3 ), DK3( N4,N4);
    local_vector_t FK( N );
    local_vector_t DKF( N3 );
    local_vector_t DKF3(N4);
    local_vector_t F2( N3 ), F3(N4);
    for( ; it != en ; ++it )
    {
        auto key = it->first;
        size_type K = key.first;
        DVLOG(2) << "======= Key=" << key ;

        DVLOG(2) << "A00K=" << A00K.at(key);
        DVLOG(2) << "A01K=" << A01K.at(key);
        DVLOG(2) << "A10K=" << A10K.at(key);

        AK = local_matrix_t::Zero( N, N );
        BK = local_matrix_t::Zero( N, N3 );
        CK = local_matrix_t::Zero( N3, N );
        
        tic();
        extractBlock( A00K.at(key), A01K.at(key), A10K.at(key), AK, e1, e2 );
        AK.bottomRightCorner(N1, N1 ) = A11K.at(key);

        dK_iterator_type dK_it;
        // index of the function space in e4
        int index = 0;
        {
            auto dK = e3.mesh()->meshToSubMesh( e1.mesh()->element(key.first).facesId());
            decltype(dK) dK1;
            dK_it = M_dK.end();
            if ( dK.second ) // there are missing
            {
                LOG(INFO) << "missing faces" << key.first << ":" << dK;
                for ( int i = 0; i < 1/*M_dim4*/; ++i )
                {
                    
                    dK1 = e(3_c,i).mesh()->meshToSubMesh( e1.mesh()->element(key.first).facesId());
                    DLOG(INFO) << "possibly adding faces from space " << i << " element " << key.first << ":" << dK1 << " fIds" << e1.mesh()->element(key.first).facesId();
                    if ( !dK1.first.empty() )
                    {
                        DLOG(INFO) << " . added face from space index " << i;
                        index = i;
                        dK_it = M_dK.emplace( K, Feel::detail::ElementFaces( std::move(dK.first), std::move(dK1.first), index ) ).first;
                        break;
                    }
                    
                }
                DCHECK( dK_it != M_dK.end() ) << "could not find missing faces";
            }
            else
                dK_it = M_dK.emplace( K, Feel::detail::ElementFaces( std::move(dK.first) ) ).first;
        }

        auto & dK = dK_it->second;
        if ( dK.hasSameTrace() )
        {
            A22 = local_matrix_t::Zero(N3,N3);
            A20 = local_matrix_t::Zero(N3,N0);
            A21 = local_matrix_t::Zero(N3,N1);
            CK = local_matrix_t::Zero(N3,N);
            FK = local_vector_t::Zero(N);
            F2 = local_vector_t::Zero(N3);
            //std::cout << "AK=\n" << AK << std::endl;
            //std::cout << "AK-AK^T=\n" << (AK-AK.transpose()) << std::endl;
            // dK contains the set of faces ids in the submesh associated to the boundary of K
        
            int n = 0;
            std::for_each( dK.faces1().begin(), dK.faces1().end(), [&]( auto dKi )
                           {
                               auto key2 = std::make_pair(key.first, dKi );
                               auto key3 = std::make_pair(dKi,key.first);
                               auto key4 = std::make_pair(dKi, dKi);

                               DVLOG(2) << "A02.count(" << key2 << ")"  << A02K.count(key2);
                               DVLOG(2) << "A20.count(" << key3 << ")"  << A20K.count(key3);
                               DVLOG(2) << "A12.count(" << key2 << ")"  << A12K.count(key2);
                               DVLOG(2) << "A21.count(" << key3 << ")"  << A21K.count(key3);
                               DVLOG(2) << "A22.count(" << key4 << ")"  << A22K.count(key4);
                               DVLOG(2) << "F2.count(" << dKi << ")"  << F2K.count(dKi);

                               extractBlock( A02K, key2, A20K, key3, BK, CK, n, e1, e3 );

                               if ( A12K.count(key2) )
                                   BK.block(N0,n*N2, N1, N2 ) = A12K.at(key2);

                               if ( A21K.count(key3) )
                                   CK.block(n*N2, N0, N2, N1 ) = A21K.at(key3);

                               if ( A22K.count(key4) )
                                   A22.block(n*N2, n*N2, N2, N2 ) = A22K.at(key4);

                               if ( F2K.count(dKi) )
                               {
                                   F2.segment(n*N2,N2)=F2K.at(dKi);
                               }

                               ++n;
                           } );
            extractBlock( F0K, K, FK, e1 );
            FK.tail(N1) = F1K.at(K);
            if ( VLOG_IS_ON(2) )
            {
                cout<< "A22=" << A22 << std::endl;
                cout<< "BK=" << BK << std::endl;
                cout<< "CK=" << CK << std::endl;
                cout<< "CK.T=" << CK.transpose() << std::endl;
                cout<< "FK=" << FK << std::endl;
                cout<< "F2=" << F2 << std::endl;
            }


#if 1
            auto Aldlt = AK.ldlt();
            ////LOG(INFO) << "Aldlt=" << Aldlt;
#else
            auto Aldlt = AK.lu();
#endif
            auto AinvB = Aldlt.solve( BK );
            auto AinvF = Aldlt.solve( FK );

            // local assemble DK and DKF
            DK=-CK*AinvB+A22 ;
            DKF=-CK*AinvF+F2;

            M_AinvB.emplace( K, AinvB );
            M_AinvF.emplace( K, AinvF );
            toc("sc.condense.localassembly", FLAGS_v>0);
            tic();
            auto dofs = e3.dofs(dK.faces1());

            S(0_c,0_c).addMatrix( dofs.data(), dofs.size(), dofs.data(), dofs.size(), DK.data(), invalid_size_type_value, invalid_size_type_value );
            V(0_c).addVector( dofs.data(), dofs.size(), DKF.data(), invalid_size_type_value, invalid_size_type_value );
            toc("sc.condense.globalassembly", FLAGS_v>0);

        }
        else
        {
            A22 = local_matrix_t::Zero(N4,N4);
            A20 = local_matrix_t::Zero(N4,N0);
            A21 = local_matrix_t::Zero(N4,N1);
            CK = local_matrix_t::Zero(N4,N);
            FK = local_vector_t::Zero(N);
            F2 = local_vector_t::Zero(N4);
            BK3 = local_matrix_t::Zero(N, N4);

            int n = 0;
            std::for_each( dK.faces1().begin(), dK.faces1().end(), [&]( auto dKi )
                           {
                               //std::cout << "face1 key.first=" << key.first << " dKi=" << dKi << std::endl;
                               auto key2 = std::make_pair(key.first, dKi );
                               auto key3 = std::make_pair(dKi,key.first);
                               auto key4 = std::make_pair(dKi, dKi);

                               DVLOG(2) << "A02.count(" << key2 << ")"  << A02K.count(key2);
                               DVLOG(2) << "A20.count(" << key3 << ")"  << A20K.count(key3);
                               DVLOG(2) << "A12.count(" << key2 << ")"  << A12K.count(key2);
                               DVLOG(2) << "A21.count(" << key3 << ")"  << A21K.count(key3);
                               DVLOG(2) << "A22.count(" << key4 << ")"  << A22K.count(key4);
                               DVLOG(2) << "F2.count(" << dKi << ")"  << F2K.count(dKi);

                               extractBlock( A02K, key2, A20K, key3, BK3, CK, n, e1, e3 );
                               //std::cout << "a CK=" << CK << std::endl;
                               if ( A12K.count(key2) )
                                   BK3.block(N0,n*N2, N1, N2 ) = A12K.at(key2);

                               if ( A21K.count(key3) )
                                   CK.block(n*N2, N0, N2, N1 ) = A21K.at(key3);
                               //std::cout << "d CK=" << CK << std::endl;
                               if ( A22K.count(key4) )
                                   A22.block(n*N2, n*N2, N2, N2 ) = A22K.at(key4);

                               if ( F2K.count(dKi) )
                               {
                                   F2.segment(n*N2,N2)=F2K.at(dKi);
                               }

                               ++n;
                           } );
            if ( VLOG_IS_ON(2) )
            {
                cout<< "A22=" << A22 << std::endl;
                cout<< "BK3=" << BK3 << std::endl;
                cout<< "CK=" << CK << std::endl;
                cout<< "CK.T=" << CK.transpose() << std::endl;
                cout<< "FK=" << FK << std::endl;
                cout<< "F2=" << F2 << std::endl;
                std::cout << "AK=" << AK << std::endl;
            }
            int n2 = 0;
            std::for_each( dK.faces2().begin(), dK.faces2().end(), [&]( auto dKi )
                           {
                               //std::cout << "face2 key.first=" << key.first << " dKi=" << dKi << std::endl;
                               auto key2 = std::make_pair(key.first, dKi );
                               auto key3 = std::make_pair(dKi,key.first);
                               auto key4 = std::make_pair(dKi, dKi);

                               DVLOG(2) << "A03.count(" << key2 << ")"  << A02K.count(key2);
                               DVLOG(2) << "A30.count(" << key3 << ")"  << A20K.count(key3);
                               DVLOG(2) << "A13.count(" << key2 << ")"  << A12K.count(key2);
                               DVLOG(2) << "A31.count(" << key3 << ")"  << A21K.count(key3);
                               DVLOG(2) << "A33.count(" << key4 << ")"  << A22K.count(key4);
                               DVLOG(2) << "F3.count(" << dKi << ")"  << F3K.count(dKi);

                               //std::cout << "d BK=" << BK3 << std::endl;
                               extractBlock( A03K, key2, A30K, key3, BK3, CK, n2, e1, e(3_c,index), n*e3.dof()->nLocalDof() );
                               //std::cout << "e CK=" << CK << std::endl;
                               //std::cout << "e BK=" << BK3 << std::endl;
                               if ( A13K.count(key2) )
                                   BK3.block(N0,n*N2+n2*N42, N1, N42 ) = A13K.at(key2);
                               //std::cout << "g BK=" << BK3 << std::endl;
                               if ( A31K.count(key3) )
                                   CK.block(n*N2+n2*N42, N0, N42, N1 ) = A31K.at(key3);
                               //std::cout << "f CK=" << CK << std::endl;
                               if ( A33K.count(key4) )
                                   A22.block(n*N2+n2*N42, n*N2, N42, N42 ) = A33K.at(key4);

                               if ( F3K.count(dKi) )
                               {
                                   F2.segment(n*N2+n2*N42,N42)=F3K.at(dKi);
                               }

                               ++n2;
                           } );
            extractBlock( F0K, K, FK, e1 );
            FK.tail(N1) = F1K.at(K);
            if ( VLOG_IS_ON(2) )
            {
                cout<< "A22=" << A22 << std::endl;
                cout<< "BK3=" << BK3 << std::endl;
                cout<< "CK=" << CK << std::endl;
                cout<< "CK.T=" << CK.transpose() << std::endl;
                cout<< "FK=" << FK << std::endl;
                cout<< "F2=" << F2 << std::endl;
                std::cout << "AK=" << AK << std::endl;
            }


#if 0
            auto Aldlt = AK.ldlt();
#else
            auto Aldlt = AK.lu();
#endif
#if 0
            local_vector_t pK(N4);
            pK=local_vector_t::Constant(N4,1);
            Eigen::VectorXd upK( N0 + N1 );
            upK.segment(0,N0)=local_vector_t::Zero(N0);
            upK.segment(N0,N1)=local_vector_t::Constant(N1,1.);
            
            std::cout << "AK*U=" << AK*upK << std::endl;
            std::cout << "BK*K=" << BK3*pK << std::endl;
            std::cout << "FK=" << FK << std::endl;
#endif
            auto AinvB = Aldlt.solve( BK3 );
            auto AinvF = Aldlt.solve( FK );

            // local assemble DK and DKF
            DK3=-CK*AinvB+A22 ;
            DKF3=-CK*AinvF+F2;

#if 0
            std::cout << "DK3 = " << DK3 << std::endl;
            std::cout << "DKF3 = " << DKF3 << std::endl;
            std::cout << "DK3*pK=" << DK3*pK << std::endl;

            std::cout << "AinvB = " << AinvB << std::endl;
            std::cout << "AinvF = " << AinvF << std::endl;
            auto Dldlt = DK3.lu();
            auto p  = Dldlt.solve( DKF3);
            std::cout << "pdK=" << p  << std::endl;
            std::cout << "uK = " << -AinvB*p + AinvF << std::endl;
#endif            
            M_AinvB.emplace( K, AinvB );
            M_AinvF.emplace( K, AinvF );
            toc("sc.condense.localassembly", FLAGS_v>0);
            tic();
            auto dofs1 = e3.dofs(dK.faces1(),S.matrixPtr()->mapRow(),0);

            auto dofs2 = e(3_c,dK.spaceIndex()).dofs(dK.faces2(),S.matrixPtr()->mapRow(),1);

            decltype(dofs1) dofs( dofs1.size()+dofs2.size() );
            std::copy(dofs1.begin(), dofs1.end(), dofs.begin() );
            std::copy(dofs2.begin(), dofs2.end(), dofs.begin()+dofs1.size() );

            local_matrix_t DK00( DK3.block( 0, 0, N41, N41 ) );
            local_matrix_t DK01( DK3.block( 0, N41, N41, N42 ) );
            local_matrix_t DK10( DK3.block( N41, 0, N42, N41 ) );
            local_matrix_t DK11( DK3.block( N41, N41, N42, N42 ) );
            local_vector_t DKF1( DKF3.head( N41 ) );
            local_vector_t DKF2( DKF3.tail( N42 ) );
            S(0_c,0_c).addMatrix( dofs1.data(), dofs1.size(), dofs1.data(), dofs1.size(), DK00.data(), invalid_size_type_value, invalid_size_type_value );
            //S.matrixPtr()->printMatlab("S1.m");
            S(0_c,1_c,0,dK.spaceIndex()).addMatrix( dofs1.data(), dofs1.size(), dofs2.data(), dofs2.size(), DK01.data(), invalid_size_type_value, invalid_size_type_value );
            //S.matrixPtr()->printMatlab("S2.m");
            S(1_c,0_c,dK.spaceIndex(),0).addMatrix( dofs2.data(), dofs2.size(), dofs1.data(), dofs1.size(), DK10.data(), invalid_size_type_value, invalid_size_type_value );
            //S.matrixPtr()->printMatlab("S3.m");
            S(1_c,1_c,dK.spaceIndex(),dK.spaceIndex()).addMatrix( dofs2.data(), dofs2.size(), dofs2.data(), dofs2.size(), DK11.data(), invalid_size_type_value, invalid_size_type_value );
            //S.matrixPtr()->printMatlab("S4.m");
            V(0_c).addVector( dofs1.data(), dofs1.size(), DKF1.data(), invalid_size_type_value, invalid_size_type_value );
            //V.vectorPtr()->printMatlab("g1.m");
            V(1_c,dK.spaceIndex()).addVector( dofs2.data(), dofs2.size(), DKF2.data(), invalid_size_type_value, invalid_size_type_value );
            //V.vectorPtr()->printMatlab("g2.m");
            toc("sc.condense.globalassembly", FLAGS_v>0);

        }

    }
}

template<typename T>
template<typename E>
void
StaticCondensation<T>::localSolve( std::shared_ptr<StaticCondensation<T>> const& rhs, E& e, std::enable_if_t<std::decay_t<E>::nspaces == 3>* )
{
    using Feel::cout;
    auto& e1 = e(0_c);
    auto& e2 = e(1_c);
    auto& e3 = e(2_c);
    
    int N0 = e1.dof()->nRealLocalDof();
    int N1 = e2.dof()->nLocalDof();
    int N2 = e3.dof()->nLocalDof();
    int N3 = N2*e1.mesh()->numLocalTopologicalFaces();
    LOG(INFO) << "[staticcondensation::localSolve]  N(flux)=" << N0 << " N(potential)=" << N1;
    
    using trace_interpolant_t  = typename std::decay_t<decltype(e3)>::local_interpolant_type;
    std::vector<std::future<void>> futs;
    futs.reserve( M_AinvB.size() );
    for( auto const& e : M_AinvB )
    {
        auto f = [N0,N1,Nup=N0+N1,Nt=N3,&e,&e1,&e2,&e3,this]()
            {
                auto K = e.first;
                auto const& A = e.second;
                auto const& F = M_AinvF.at( K );
                Eigen::VectorXd upK( Nup );
                trace_interpolant_t pdK( Nt );
                e3.element( M_dK[K].faces1(), pdK );
                upK.noalias() = -A*pdK + F;

                e1.assignE( K, upK.head( N0 ) );
                e2.assignE( K, upK.tail( N1 ) );
            };
        f();
        //futs.push_back( std::async(std::launch::async, f ) );
    }
#if 0
    for( auto& fut : futs )
    {
        fut.get();
    }
#endif
}

template<typename T>
template<typename E>
void
StaticCondensation<T>::localSolve( std::shared_ptr<StaticCondensation<T>> const& rhs, E& e, std::enable_if_t<std::decay_t<E>::nspaces == 4>*  )
{
    using Feel::cout;
    auto& e1 = e(0_c);
    auto& e2 = e(1_c);
    auto& e3 = e(2_c);
    int N0 = e1.dof()->nRealLocalDof();
    int N1 = e2.dof()->nLocalDof();
    cout << "[staticcondensation::localSolve]  N0=" << N0 << " N1=" << N1 << std::endl;
    int N2 = e3.dof()->nLocalDof();
    int N4 = 0;
    int N4d1 = N2*(e1.mesh()->numLocalTopologicalFaces());
    int N4d11 = N2*(e1.mesh()->numLocalTopologicalFaces());
    if ( e.functionSpace( 3_c )->numberOfSpaces() > 0 )
    {
        N4 = e(3_c,0).dof()->nLocalDof();
        N4d1 = N2*(e1.mesh()->numLocalTopologicalFaces()-1);
    }
    int N3 = N2*e1.mesh()->numLocalTopologicalFaces();
    using trace_interpolant_t  = typename std::decay_t<decltype(e(2_c))>::local_interpolant_type;
    std::vector<std::future<void>> futs;
    futs.reserve( M_AinvB.size() );
    for( auto const& m : M_AinvB )
    {
        //f(m);

        auto f = [N0,N1,N4d1,N4d11,N4,&m,this,&e,&e1,&e2,&e3](  )
            {
                auto K = m.first;
                auto const& A = m.second;
                auto const& F = M_AinvF.at( K );
                Eigen::VectorXd upK( N0 + N1 );
                auto & dK = M_dK[K];
                if ( dK.hasSameTrace() )
                {
                    trace_interpolant_t pdK1( N4d11 );

                    e3.element( dK.faces1(), pdK1 );
                    upK.noalias() = -A*pdK1 + F;
                    e1.assignE( K, upK.head( N0 ) );
                    e2.assignE( K, upK.tail( N1 ) );
                }
                else
                {
                    trace_interpolant_t pdK2( N4d1+N4 );
                    e3.element( dK.faces1(), pdK2.topLeftCorner(N4d1,1) );
                    e(3_c,dK.spaceIndex()).element( dK.faces2(), pdK2.bottomLeftCorner(N4,1) );
                    //std::cout << "pdK2=" << pdK2 << std::endl;
                    //std::cout << "A=" << A << std::endl;
                    //std::cout << "F=" << F << std::endl;
                    upK.noalias() = -A*pdK2 + F;
                    //std::cout << "upK = " << upK << std::endl;
                    //pdK2 = local_vector_t::Constant( pdK2.size(), 1. );
                    //pdK2( pdK2.size()-1)=1;
                    //std::cout << "pdK2=" << pdK2 << std::endl;
                    //upK.noalias() = -A*pdK2 + F;
                    //std::cout << "upK2 = " << upK  << std::endl;
                    e1.assignE( K, upK.head( N0 ) );
                    e2.assignE( K, upK.tail( N1 ) );
                }
            };

        //futs.push_back( std::async(std::launch::async, f ) );
        f();
    }
#if 0
    for( auto& fut : futs )
    {
        fut.get();
    }
#endif
}


}
#endif




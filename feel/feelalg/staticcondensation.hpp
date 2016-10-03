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

#include <boost/functional/hash.hpp>

#include <unordered_map>
#include <Eigen/Core>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/feelio.hpp>
#include <feel/feeldiscr/traits.hpp>


namespace Feel {

template<typename T>
class StaticCondensation
{
public:
    using block_index_t = std::pair<int,int>;
    using block_element_t = std::pair<size_type,size_type>;
    using value_type = T;

    StaticCondensation() = default;
    template<typename E1, typename E2, typename E3, typename M_ptrtype, typename V_ptrtype>
    void condense( boost::shared_ptr<StaticCondensation<T>> const& rhs, E1& e1, E2& e2, E3 & e3, M_ptrtype& S, V_ptrtype& V );
    template<typename E1, typename E2, typename E3, typename E4, typename M_ptrtype, typename V_ptrtype>
    void condense( boost::shared_ptr<StaticCondensation<T>> const& rhs, E1& e1, E2& e2, E3 & e3, E4& e4, M_ptrtype& S, V_ptrtype& V );

    template<typename E1, typename E2, typename E3>
    void localSolve( boost::shared_ptr<StaticCondensation<T>> const& rhs, E1& e1, E2& e2, E3 & e3 );
    template<typename E1, typename E2, typename E3, typename E4>
    void localSolve( boost::shared_ptr<StaticCondensation<T>> const& rhs, E1& e1, E2& e2, E3 & e3, E4& e4 );

    void addLocalMatrix( int* rows, int nrows,
                         int* cols, int ncols,
                         value_type* data,
                         size_type K = 0,
                         size_type K2 = invalid_size_type_value );

    void addLocalVector( int* rows, int nrows,
                         value_type* data,
                         size_type K = 0,
                         size_type K2 = invalid_size_type_value );

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
    block_index_t M_block_rowcol;
    int M_block_row;
};


template<typename A00_t,typename A01_t,typename A10_t, typename M2, typename E1, typename E2>
void
extractBlock( A00_t const& a00, A01_t const& a01, A10_t const& a10,
              M2& ak, E1 const& e1, E2 const& e2, tensor2symm_true )
{
    uint16_type N0 = e1.dof()->nRealLocalDof( false );
    uint16_type N0c = e1.dof()->nLocalDof( true );
    uint16_type N1 = e2.dof()->nLocalDof();
    ak.Zero( N0+N1, N0+N1 );
    //cout << "A00.rows=" << a00.rows() <<  " A00.cols=" << a00.cols() << std::endl;
    //cout << "AK.rows=" << ak.rows() <<  " AK.cols=" << ak.cols() << std::endl;
    //cout << "AK.0 = " << std::endl << ak << std::endl;
    // A00
    for( int c1e1 = 0; c1e1 < e1.nComponents1; ++c1e1 )
    {
        for( int c2e1 = 0; c2e1 < c1e1; ++c2e1 )
        {
            const int k1 = Feel::detail::symmetricIndex(c1e1,c2e1,e1.nComponents1);
            for( int c1e2 = 0; c1e2 < e2.nComponents1; ++c1e2 )
            {
                for( int c2e2 = 0; c2e2 < c1e2; ++c2e2 )
                {
                    const int k2 = Feel::detail::symmetricIndex(c1e2,c2e2,e2.nComponents1);
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
                for( int c2e2 = 0; c2e2 < c1e2; ++c2e2 )
                {
                    const int k2 = Feel::detail::symmetricIndex(c1e2,c2e2,e2.nComponents1);
                    ak.block( N0c*k1, N0c*k2, N0c, N0c ) = a00.block( N0c*(c1e1*e1.nComponents1+c1e1),
                                                                      N0c*(c1e2*e2.nComponents1+c2e2), N0c, N0c );
                    ak.block( N0c*k1, N0c*k2, N0c, N0c ) += a00.block( N0c*(c1e1*e1.nComponents1+c1e1),
                                                                       N0c*(c2e2*e2.nComponents1+c1e2), N0c, N0c );

                }
                const int k2 = Feel::detail::symmetricIndex(c1e2,c1e2,e2.nComponents1);
                ak.block( N0c*k1, N0c*k2, N0c, N0c ) = a00.block( N0c*(c1e1*e1.nComponents1+c1e1),
                                                                  N0c*(c1e2*e2.nComponents1+c1e2), N0c, N0c );
        }

    }


    //cout << "AK.1 = " << std::endl << ak << std::endl;
    // A01
    for( int c1e1 = 0; c1e1 < e1.nComponents1; ++c1e1 )
    {
        for( int c2e1 = 0; c2e1 < c1e1; ++c2e1 )
        {
            const int k1 = Feel::detail::symmetricIndex(c1e1,c2e1,e1.nComponents1);
            ak.block( N0c*k1, N0, N0c, N1 ) = a01.block( N0c*(c1e1*e1.nComponents1+c2e1), 0, N0c, N1 );
            ak.block( N0c*k1, N0, N0c, N1 ) += a01.block( N0c*(c2e1*e1.nComponents1+c1e1), 0, N0c, N1 );

            ak.block( N0, N0c*k1, N1, N0c ) = a10.block( 0, N0c*(c1e1*e1.nComponents1+c2e1), N1, N0c );
            ak.block( N0, N0c*k1, N1, N0c ) += a10.block( 0, N0c*(c2e1*e1.nComponents1+c1e1), N1, N0c );
        }
        const int k1 = Feel::detail::symmetricIndex(c1e1,c1e1,e1.nComponents1);
        ak.block( N0c*k1, N0, N0c, N1 ) = a01.block( N0c*(c1e1*e1.nComponents1+c1e1), 0, N0c, N1 );
        ak.block( N0, N0c*k1, N1, N0c ) = a10.block( 0, N0c*(c1e1*e1.nComponents1+c1e1), N1, N0c );
    }
    //cout << "AK.2 = " << std::endl << ak << std::endl;
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
              E1 const& e1, E3 const& e3,
              std::enable_if_t<is_tensor2symm_field_v<E1>>* = nullptr )
{
    uint16_type N0 = e1.dof()->nRealLocalDof( false );
    uint16_type N0c = e1.dof()->nLocalDof( true );
    int N2 = e3.dof()->nLocalDof();

    if ( A02K.count(key2) )
    {
        for( int c1e1 = 0; c1e1 < e1.nComponents1; ++c1e1 )
        {
            for( int c2e1 = 0; c2e1 < c1e1; ++c2e1 )
            {
                const int k1 = Feel::detail::symmetricIndex(c1e1,c2e1,e1.nComponents1);
                BK.block(N0c*k1, n*N2, N0c, N2 ) = A02K.at(key2).block( N0c*(c1e1*e1.nComponents1+c2e1), 0, N0c, N2 );
                BK.block(N0c*k1, n*N2, N0c, N2 ) += A02K.at(key2).block( N0c*(c2e1*e1.nComponents1+c1e1), 0, N0c, N2 );
            }
            const int k1 = Feel::detail::symmetricIndex(c1e1,c1e1,e1.nComponents1);
            BK.block(N0c*k1, n*N2, N0c, N2 ) = A02K.at(key2).block( N0c*(c1e1*e1.nComponents1+c1e1), 0, N0c, N2 );
        }
    }

    if ( A20K.count(key3) )
    {
        for( int c1e1 = 0; c1e1 < e1.nComponents1; ++c1e1 )
        {
            for( int c2e1 = 0; c2e1 < c1e1; ++c2e1 )
            {
                const int k1 = Feel::detail::symmetricIndex(c1e1,c2e1,e1.nComponents1);
                CK.block(n*N2, N0c*k1, N2, N0c ) = A20K.at(key3).block(0, N0c*(c1e1*e1.nComponents1+c2e1), N2, N0c );
                CK.block(n*N2, N0c*k1, N2, N0c ) += A20K.at(key3).block(0, N0c*(c2e1*e1.nComponents1+c1e1), N2, N0c );
            }
            const int k1 = Feel::detail::symmetricIndex(c1e1,c1e1,e1.nComponents1);
            CK.block(n*N2, N0c*k1, N2, N0c ) = A20K.at(key3).block(0, N0c*(c1e1*e1.nComponents1+c1e1), N2, N0c );
        }
    }

}



template<typename A02_t, typename A20_t, typename Key1_t, typename Key2_t, typename BK_t, typename CK_t, typename E1, typename E3>
void
extractBlock( A02_t const& A02K, Key1_t const& key2,
              A20_t const& A20K, Key2_t const& key3,
              BK_t& BK, CK_t& CK, int n,
              E1 const& e1, E3 const& e3,
              std::enable_if_t<!is_tensor2symm_field_v<E1>>* = nullptr )
{
    int N0 = e1.dof()->nLocalDof();
    int N2 = e3.dof()->nLocalDof();
    if ( A02K.count(key2) )
        BK.block(0, n*N2, N0, N2 ) = A02K.at(key2);

    if ( A20K.count(key3) )
        CK.block(n*N2, 0, N2, N0 ) = A20K.at(key3);

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
            for( int c2e1 = 0; c2e1 < c1e1; ++c2e1 )
            {
                const int k1 = Feel::detail::symmetricIndex(c1e1,c2e1,e1.nComponents1);
                FK.segment(N0c*k1, N0c) = F0K.at(K).segment( N0c*(c1e1*e1.nComponents1+c2e1), N0c );
                FK.segment(N0c*k1, N0c) += F0K.at(K).segment( N0c*(c2e1*e1.nComponents1+c1e1), N0c );
            }
            const int k1 = Feel::detail::symmetricIndex(c1e1,c1e1,e1.nComponents1);
            FK.segment(N0c*k1, N0c) = F0K.at(K).segment( N0c*(c1e1*e1.nComponents1+c1e1), N0c );
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
template<typename E1, typename E2, typename E3, typename M_ptrtype, typename V_ptrtype>
void
StaticCondensation<T>::condense( boost::shared_ptr<StaticCondensation<T>> const& rhs, E1 &e1, E2& e2, E3 & e3, M_ptrtype& S, V_ptrtype& V )
{
    using Feel::cout;

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

    auto const& F0K = rhs->M_local_vectors[0];
    //LOG(INFO) << "F0K.size=" << F0K.size();
    auto const& F1K = rhs->M_local_vectors[1];
    //LOG(INFO) << "F1K.size=" << F1K.size();
    auto const& F2K = rhs->M_local_vectors[2];
    //LOG(INFO) << "F2K.size=" << F2K.size();

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
        auto dK = e3.mesh()->meshToSubMesh( e1.mesh()->element(key.first).facesId());

        int n = 0;
        std::for_each( dK.begin(), dK.end(), [&]( auto dKi )
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
        auto dofs = e3.dofs(dK);

        S->addMatrix( dofs.data(), dofs.size(), dofs.data(), dofs.size(), DK.data(), invalid_size_type_value, invalid_size_type_value );
        V->addVector( dofs.data(), dofs.size(), DKF.data(), invalid_size_type_value, invalid_size_type_value );
        toc("sc.condense.globalassembly", FLAGS_v>0);
    }
}

template<typename T>
template<typename E1, typename E2, typename E3, typename E4, typename M_ptrtype, typename V_ptrtype>
void
StaticCondensation<T>::condense( boost::shared_ptr<StaticCondensation<T>> const& rhs, E1 &e1, E2& e2, E3 & e3, E4& e4, M_ptrtype& S, V_ptrtype& V )
{
    using Feel::cout;

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

    auto const& F0K = rhs->M_local_vectors[0];
    //LOG(INFO) << "F0K.size=" << F0K.size();
    auto const& F1K = rhs->M_local_vectors[1];
    //LOG(INFO) << "F1K.size=" << F1K.size();
    auto const& F2K = rhs->M_local_vectors[2];
    //LOG(INFO) << "F2K.size=" << F2K.size();

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
    for( ; it != en ; ++it )
    {
        auto key = it->first;
        size_type K = key.first;
        DVLOG(2) << "======= Key=" << key ;

        DVLOG(2) << "A00K=" << A00K.at(key);
        DVLOG(2) << "A01K=" << A01K.at(key);
        DVLOG(2) << "A10K=" << A10K.at(key);


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
        auto dK = e3.mesh()->meshToSubMesh( e1.mesh()->element(key.first).facesId());
        int n = 0;
        std::for_each( dK.begin(), dK.end(), [&]( auto dKi )
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




#if 0
        auto Aldlt = AK.ldlt();
        ////LOG(INFO) << "Aldlt=" << Aldlt;
#else
        auto Aldlt = AK.lu();
#endif
        auto AinvB = Aldlt.solve( BK );
        auto AinvF = Aldlt.solve( FK );

        auto pdK = e3.element( dK );

        Eigen::VectorXd upK = -AinvB*pdK + AinvF;

        e1.assignE( K, upK.head( N0 ) );
        e2.assignE( K, upK.tail( N1 ) );

        // local assemble DK and DKF
        DK=-CK*AinvB+A22 ;
        DKF=-CK*AinvF+F2;

        M_AinvB.emplace( K, AinvB );
        M_AinvF.emplace( K, AinvF );

        auto dofs = e3.dofs(dK);

        S->addMatrix( dofs.data(), dofs.size(), dofs.data(), dofs.size(), DK.data(), invalid_size_type_value, invalid_size_type_value );
        V->addVector( dofs.data(), dofs.size(), DKF.data(), invalid_size_type_value, invalid_size_type_value );
    }
}

template<typename T>
template<typename E1, typename E2, typename E3>
void
StaticCondensation<T>::localSolve( boost::shared_ptr<StaticCondensation<T>> const& rhs, E1& e1, E2& e2, E3 & e3 )
{
    using Feel::cout;

    int N0 = e1.dof()->nRealLocalDof();
    int N1 = e2.dof()->nLocalDof();
    cout << "[staticcondensation::localSolve]  N0=" << N0 << " N1=" << N1 << std::endl;
    Eigen::VectorXd upK( N0 + N1 );
    for( auto const& e : M_AinvB )
    {
        auto K = e.first;
        auto const& A = e.second;
        auto const& F = M_AinvF.at( K );

        auto dK = e3.mesh()->meshToSubMesh( e1.mesh()->element(K).facesId());
        auto pdK = e3.element( dK );
        upK.noalias() = -A*pdK + F;

        e1.assignE( K, upK.head( N0 ) );
        e2.assignE( K, upK.tail( N1 ) );
    }

}

template<typename T>
template<typename E1, typename E2, typename E3, typename E4>
void
StaticCondensation<T>::localSolve( boost::shared_ptr<StaticCondensation<T>> const& rhs, E1& e1, E2& e2, E3 & e3, E4& e4 )
{
    using Feel::cout;

    int N0 = e1.dof()->nRealLocalDof();
    int N1 = e2.dof()->nLocalDof();
    cout << "[staticcondensation::localSolve]  N0=" << N0 << " N1=" << N1 << std::endl;

    Eigen::VectorXd upK( N0 + N1 );
    for( auto const& e : M_AinvB )
    {
        auto K = e.first;
        auto const& A = e.second;
        auto const& F = M_AinvF.at( K );

        auto dK = e3.mesh()->meshToSubMesh( e1.mesh()->element(K).facesId());
        auto pdK = get_trace( dK, e3, e4 );

        upK.noalias() = -A*pdK + F;
        e1.assignE( K, upK.head( N0 ) );
        e2.assignE( K, upK.tail( N1 ) );
    }

}


}
#endif

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


namespace Feel {

template<typename T>
class StaticCondensation
{
public:
    using block_index_t = std::pair<int,int>;
    using block_element_t = std::pair<size_type,size_type>;
    using value_type = T;

    StaticCondensation() = default;
    template<typename E1, typename E2, typename E3>
    void localSolve( boost::shared_ptr<StaticCondensation<T>> const& rhs, E1& e1, E2& e2, E3 const& e3 );

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

private:
    /**
     * unassembled view of the matrix
     */
    std::unordered_map<int,std::unordered_map<size_type,local_vector_t>> M_local_vectors;
    std::unordered_map<int,std::unordered_map<size_type,local_index_t>> M_local_vrows;
    std::unordered_map<block_index_t,std::unordered_map<block_element_t,local_matrix_t,boost::hash<block_element_t>>,boost::hash<block_index_t>> M_local_matrices;
    std::unordered_map<block_index_t,std::unordered_map<block_element_t,local_index_t,boost::hash<block_element_t>>,boost::hash<block_index_t>> M_local_rows;
    std::unordered_map<block_index_t,std::unordered_map<block_element_t,local_index_t,boost::hash<block_element_t>>,boost::hash<block_index_t>> M_local_cols;
    block_index_t M_block_rowcol;
    int M_block_row;

};

template<typename T>
template<typename E1, typename E2, typename E3>
void
StaticCondensation<T>::localSolve( boost::shared_ptr<StaticCondensation<T>> const& rhs, E1& e1, E2& e2, E3 const& e3 )
{
    using Feel::cout;

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

    auto it = A00K.begin();
    auto en = A00K.end();

    int N = M_local_rows[std::make_pair(0,0)].begin()->second.size() + M_local_rows[std::make_pair(1,0)].begin()->second.size();
    int N0 = M_local_rows[std::make_pair(0,0)].begin()->second.size();
    int N1 = M_local_rows[std::make_pair(1,0)].begin()->second.size();
    int N2 = M_local_cols[std::make_pair(0,2)].begin()->second.size();
    int N3 = N2*e1.mesh()->numLocalTopologicalFaces();
    cout << "N=" << N << " N0=" << N0 << " N1=" << N1 << " N2=" << N2 << std::endl;
    local_matrix_t AK( N, N ),A00(N0,N0),A01(N0,N1),A10(N1,N0), A11(N1,N1), A20(N3,N0), A21(N3,N1), A22(N3,N3);
    local_matrix_t BK( N, N3 );
    local_matrix_t CK( N3, N );
    local_matrix_t DK( N3, N3 );
    local_vector_t FK( N );
    local_vector_t F3K( N3 );
    for( ; it != en ; ++it )
    {
        auto key = it->first;
        size_type K = key.first;
        LOG(INFO) << "======= Key=" << key ;
        AK.topLeftCorner(N0, N0 ) = A00K.at(key);
        A00=A00K.at(key);
        LOG(INFO) << "A00K=" << A00K.at(key);
        AK.bottomLeftCorner(N1, N0 ) = A10K.at(key);
        LOG(INFO) << "A10K=" << A10K.at(key);
        AK.topRightCorner(N0, N1 ) = A01K.at(key);
        A01=A01K.at(key);
        LOG(INFO) << "A01K=" << A01K.at(key);
        LOG(INFO) << "A01K.t=" << A01K.at(key).transpose();
        AK.bottomRightCorner(N1, N1 ) = A11K.at(key);
        LOG(INFO) << "A11K=" << A11K.at(key);
        A22K = A22K.at(key);
        LOG(INFO) << "A22K=" << A22K.at(key);
        LOG(INFO) << "AK=" << AK;
        local_matrix_t Z=(AK-AK.transpose());
        cout << "Z=" << Z << std::endl;
        cout << "AK=" << AK << std::endl;
        cout << "AK.T="  << AK.transpose() << std::endl;
        // dK contains the set of faces ids in the submesh associated to the boundary of K
        auto dK = e3.mesh()->meshToSubMesh( e1.mesh()->element(key.first).facesId());
        LOG(INFO) << "element K faceids: " << e1.mesh()->element(key.first).facesId();
        LOG(INFO) << "dK=" << dK;
        int n = 0;
        local_matrix_t B0K(N0,N3);
        local_matrix_t B1K(N1,N3);
        std::for_each( dK.begin(), dK.end(), [&]( auto dKi )
                       {
                           auto key2 = std::make_pair(key.first, dKi );
                           LOG(INFO) << "dKi=" << key2;
                           if ( A02K.find(key2) == A02K.end() )
                           {
                               LOG(INFO) << "A02K not ok";
                               return;
                           }

                           if ( A12K.find(key2) == A12K.end() )
                           {
                               LOG(INFO) << "A12K not ok";
                               return;
                           }
                           BK.block(0, n*N2, N0, N2 ) = A02K.at(key2);
                           CK.block(n*N2, 0, N0, N2 ) = A20K.at(key2);
                           //F3K.segment(n*N2)=F1K.at(K);
                           B0K.block(0, n*N2, N0, N2 ) = A02K.at(key2);
                           cout << "BK.block(0, " << n*N2 << "," << N0 <<  N2 << " ) = " << BK.block(0, n*N2, N0, N2 ) << std::endl;
                           A02K.at(key2);
                           BK.block(N0,n*N2, N1, N2 ) = A12K.at(key2);
                           CK.block(n*N2, N0, N1, N2 ) = A21K.at(key2);
                           B1K.block(0,n*N2, N1, N2 ) = A12K.at(key2);
                           cout << "BK.block(" << N0 << ", " << n*N2 << "," << N1 <<  N2 << " ) = " << BK.block(N0, n*N2, N1, N2 ) << std::endl;
                           ++n;
                       } );
        LOG(INFO) << "BK=" << BK;
        cout<< "BK=" << BK << std::endl;
        cout<< "B0K=" << B0K << std::endl;
        cout<< "B1K=" << B1K << std::endl;
        cout<< "CK=" << CK << std::endl;
        FK.head(N0) = F0K.size()?F0K.at(K):local_vector_t::Zero( N0 );
        FK.tail(N1) = F1K.at(K);

        LOG(INFO) << "FK=" << FK;
        cout<< "FK=" << FK << std::endl;

#if 0
        auto Aldlt = AK.ldlt();
        //LOG(INFO) << "Aldlt=" << Aldlt;
#else
        auto Aldlt = AK.lu();
#endif
        auto AinvB = Aldlt.solve( BK );
        LOG(INFO) << "AinvB=" << AinvB;
        auto AinvF = Aldlt.solve( FK );
        LOG(INFO) << "AinvF=" << AinvF;

        // loop over each face on boundary of K
        auto pdK = e3.element( dK );
        auto uK = e1.element( { K } );
        cout << "uK=" << uK << std::endl;
        auto pK = e2.element( { K } );
        cout << "pK=" << pK << std::endl;
        LOG(INFO) << "pdK=" << pdK;
        CHECK(pdK.size() == N2*e1.mesh()->numLocalTopologicalFaces() ) << "Invalid sizes " << pdK.size() << "," << N2*e1.mesh()->numLocalTopologicalFaces() << "," << N2 << "," << e1.mesh()->numLocalTopologicalFaces() ;
        cout << "A00*uK=" << A00*uK << std::endl;
        cout << "A01*pK=" << A01*pK << std::endl;
        cout << "B0K*pdK=" << B0K*pdK << std::endl;
        cout << "A10*uK=" << A10*uK << std::endl;
        cout << "A11*pK=" << A11*pK << std::endl;
        cout << "B1K*pdK=" << B1K*pdK << std::endl;
        Eigen::VectorXd upK = -AinvB*pdK + AinvF;
        LOG(INFO) << "upK=" << upK;
        e1.assignE( K, upK.head( N0 ) );
        e2.assignE( K, upK.tail( N1 ) );

        // now calculate  pdK
        // local assemble DK and DKF
        DK=CK*Ainv+A22K;
        DKF=CK*AinvF;
        // global assemble contribution to element K
#if 0
        auto dofs = e3.dof(dK);
        M_mat->addMatrix( dofs.data(), dofs.size(), dofs.data(), dofs.size(), DK.data(), invalid_size_type_value, invalid_size_type_value );
         M_vec->addVector( dofs.data(), dofs.size(), DFK.data(), invalid_size_type_value, invalid_size_type_value );
#endif

        LOG(INFO) << "======= done";

    }
    M_mat->close();
    M_vec->close();

}


}
#endif

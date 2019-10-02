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
#include <future>

#include <feel/feelalg/staticcondensation.hpp>


namespace Feel {

template<typename T, typename IndexT>
void StaticCondensation<T,IndexT>::addLocalMatrix ( int* rows, int nrows,
                                             int* cols, int ncols,
                                             value_type* data,
                                             size_type K, size_type K2  )
{
    std::lock_guard<std::mutex> guard(mutex_add_m);
    //tic();
    if ( K == invalid_v<size_type> ) return;
    if ( K2 == invalid_v<size_type> ) return;
    auto key = std::make_pair(K,K2);
    auto entry = this->M_local_matrices[this->M_block_rowcol].find(key);
    if ( entry == this->M_local_matrices[this->M_block_rowcol].end() )
    {
        this->M_local_matrices[this->M_block_rowcol][key].noalias() = raw_matrix_map_t( data, nrows, ncols );
        //cout << "SC inserting matrix entry [" << this->M_block_rowcol << "][" << key << "]=" << this->M_local_matrices[this->M_block_rowcol][key] << std::endl;
    }
    else
    {
        this->M_local_matrices[this->M_block_rowcol][key].noalias()+=raw_matrix_map_t( data, nrows, ncols );
        //cout << "SC adding matrix entry [" << this->M_block_rowcol << "][" << key << "]=" << this->M_local_matrices[this->M_block_rowcol][key] << std::endl;
    }
    //LOG(INFO) << "[" << this->M_block_rowcol << "][" << key << "]=" << this->M_local_matrices[this->M_block_rowcol][key];
#if 0
    this->M_local_rows[this->M_block_rowcol][key] = raw_index_map_t( rows, nrows );
    this->M_local_cols[this->M_block_rowcol][key] = raw_index_map_t( cols, ncols );
#else
    //LOG(INFO) << "ROWS=" << raw_index_map_t( rows, nrows );
    //LOG(INFO) << "COLS=" << raw_index_map_t( cols, ncols );
#endif
    //toc("sc.addLocalMatrix",FLAGS_v>0);
}

template<typename T, typename IndexT>
void StaticCondensation<T,IndexT>::addLocalVector ( int* rows, int nrows,
                                             value_type* data,
                                             size_type K, size_type K2  )
{
    std::lock_guard<std::mutex> guard(mutex_add_v);
    //tic();
    if ( K == invalid_v<size_type> ) return;
    auto entry = this->M_local_vectors[this->M_block_row].find(K);
    if ( entry == this->M_local_vectors[this->M_block_row].end() )
    {
        this->M_local_vectors[this->M_block_row][K].noalias() = raw_vector_map_t( data, nrows );
        //cout << "SC vec inserting F entry " << this->M_block_row << "," << K << " =" << this->M_local_vectors[this->M_block_row][K] << std::endl;
    }
    else
    {
        this->M_local_vectors[this->M_block_row][K].noalias()+=raw_vector_map_t( data, nrows );
        //cout << "SC vec add F entry " << this->M_block_row << "," << K << " =" << this->M_local_vectors[this->M_block_row][K] << std::endl;
    }
#if 0
    //LOG(INFO) << "F add entry " << this->M_block_row << "," << K << " =" << this->M_local_vectors[this->M_block_row][K];
    // cout << "F add entry " << this->M_block_row << "," << K << " =" << this->M_local_vectors[this->M_block_row][K] << std::endl;
    this->M_local_vrows[this->M_block_row][K] = raw_index_map_t( rows, nrows );
#endif
    //toc("sc.addLocalVector",FLAGS_v>0);
}

template class StaticCondensation<double>;
template class StaticCondensation<std::complex<double>>;

}

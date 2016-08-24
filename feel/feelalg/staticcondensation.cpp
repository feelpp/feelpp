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
#include <feel/feelalg/staticcondensation.hpp>


namespace Feel {

template<typename T>
void StaticCondensation<T>::addLocalMatrix ( int* rows, int nrows,
                                          int* cols, int ncols,
                                          value_type* data,
                                             size_type K, size_type K2  )
{
    auto key = std::make_pair(K,K2);
    auto entry = this->M_local_matrices[this->M_block_rowcol].find(key);
    if ( entry == this->M_local_matrices[this->M_block_rowcol].end() )
        this->M_local_matrices[this->M_block_rowcol][key] = raw_matrix_map_t( data, nrows, ncols );
    else
        this->M_local_matrices[this->M_block_rowcol][key]+=raw_matrix_map_t( data, nrows, ncols );
    LOG(INFO) << "[" << this->M_block_rowcol << "][" << key << "]=" << this->M_local_matrices[this->M_block_rowcol][key];
    // cout << "[" << this->M_block_rowcol << "][" << key << "]=" << this->M_local_matrices[this->M_block_rowcol][key] << std::endl;
    this->M_local_rows[this->M_block_rowcol][key] = raw_index_map_t( rows, nrows );
    this->M_local_cols[this->M_block_rowcol][key] = raw_index_map_t( cols, ncols );

}

template<typename T>
void StaticCondensation<T>::addLocalVector ( int* rows, int nrows,
                                             value_type* data,
                                             size_type K, size_type K2  )
{
    auto entry = this->M_local_vectors[this->M_block_row].find(K);
    if ( entry == this->M_local_vectors[this->M_block_row].end() )
        this->M_local_vectors[this->M_block_row][K] = raw_vector_map_t( data, nrows );
    else
        this->M_local_vectors[this->M_block_row][K]+=raw_vector_map_t( data, nrows );
    LOG(INFO) << "F add entry " << this->M_block_row << "," << K << " =" << this->M_local_vectors[this->M_block_row][K];
    // cout << "F add entry " << this->M_block_row << "," << K << " =" << this->M_local_vectors[this->M_block_row][K] << std::endl;
    this->M_local_vrows[this->M_block_row][K] = raw_index_map_t( rows, nrows );
}

template class StaticCondensation<double>;
template class StaticCondensation<std::complex<double>>;

}

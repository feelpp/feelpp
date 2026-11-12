/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2016-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
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
    this->M_local_rows[this->M_block_rowcol][key] = raw_index_map_t( rows, nrows );
    this->M_local_cols[this->M_block_rowcol][key] = raw_index_map_t( cols, ncols );
    //toc("sc.addLocalMatrix",Environment::logVerbosityLevel()>0);
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
    this->M_local_vrows[this->M_block_row][K] = raw_index_map_t( rows, nrows );
    //toc("sc.addLocalVector",Environment::logVerbosityLevel()>0);
}

template class StaticCondensation<double>;
template class StaticCondensation<std::complex<double>>;

}

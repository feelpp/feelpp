//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 24 Sep 2017
//! @copyright 2017 Feel++ Consortium
//!
#define FEELPP_BILINEARFORMBASE_NOEXTERN 1
#include <feel/feelvf/bilinearformbase.hpp>



namespace Feel {




template<typename T>
BilinearFormBase<T>::BilinearFormBase( BilinearFormBase const& __vf )
    :
    super( __vf ),
    M_name( __vf.M_name ),
    M_pattern( __vf.M_pattern ),
    M_matrix( __vf.M_matrix ),
    M_lb( __vf.M_lb ),
    M_row_startInMatrix( __vf.M_row_startInMatrix ),
    M_col_startInMatrix( __vf.M_col_startInMatrix ),
    M_do_build( __vf.M_do_build ),
    M_do_threshold( __vf.M_do_threshold ),
    M_threshold( __vf.M_threshold ),
    b_mutex()
{
    auto dmTest = M_matrix->mapRowPtr();
    auto dmTrial = M_matrix->mapColPtr();
    this->setDofIdToContainerIdTest( dmTest->dofIdToContainerId( M_row_startInMatrix ) );
    this->setDofIdToContainerIdTrial( dmTrial->dofIdToContainerId( M_col_startInMatrix ) );
}


template<typename T>
BilinearFormBase<T>&
BilinearFormBase<T>::operator=( BilinearFormBase const& form )
{
    if ( this != &form )
    {
        super::operator=( form );
        M_name = form.M_name;
        M_pattern = form.M_pattern;
        M_matrix->zero();
        M_matrix->addMatrix( 1.0, form.M_matrix );
        M_row_startInMatrix = form.M_row_startInMatrix;
        M_col_startInMatrix = form.M_col_startInMatrix;
        M_lb = form.M_lb;
        auto dmTest = M_matrix->mapRowPtr();
        auto dmTrial = M_matrix->mapColPtr();
        this->setDofIdToContainerIdTest( dmTest->dofIdToContainerId( M_row_startInMatrix ) );
        this->setDofIdToContainerIdTrial( dmTrial->dofIdToContainerId( M_col_startInMatrix ) );
    }

    return *this;
}
template<typename T>
void
BilinearFormBase<T>::addMatrix( int* rows, int nrows,
                                int* cols, int ncols,
                                value_type* data,
                                size_type K,
                                size_type K2)
{
    std::lock_guard<std::mutex> guard(b_mutex);
    //std::cout << "add matrix\n";
    M_matrix->addMatrix( rows, nrows, cols, ncols, data, K, K2 );
}

template<typename T>
void
BilinearFormBase<T>::zeroRows( std::vector<int> const& __dofs,
                               Vector<value_type> const&__values,
                               Vector<value_type>& rhs,
                               Feel::Context const& on_context,
                               double value_on_diagonal )
{
    M_matrix->zeroRows( __dofs, __values, rhs, on_context, value_on_diagonal );
}


template class BilinearFormBase<double>;
} // Feel

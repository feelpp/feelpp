/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@imag.fr>
       Date: 2012-01-18

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file vectorblock.cpp
   \author Vincent Chabannes <vincent.chabannes@imag.fr>
   \date 2012-01-18
 */

#include <feel/feelalg/vectorblock.hpp>



namespace Feel
{


template <typename T>
VectorBlockBase<T>::VectorBlockBase( vf::BlocksBase<vector_ptrtype> const & blockVec,
                                   backend_type &backend,
                                   bool copy_values )
    :
    M_vec()
{
    auto NR = blockVec.nRow();

    size_type _size = 0;
    for (uint i=0;i<NR;++i)
        _size += blockVec(i,0)->size();

    M_vec = backend.newVector(_size,_size);
    M_vec->zero();

    if (copy_values)
        {
            size_type start_i=0;
            for (uint i=0;i<NR;++i)
                {
                    this->updateBlockVec(blockVec(i,0),start_i);
                    start_i += blockVec(i,0)->size();
                }
        }
}

template <typename T>
void
VectorBlockBase<T>::updateBlockVec(vector_ptrtype const& m, size_type start_i)
{
    const size_type size = m->size();
    for (uint i=0;i<size;++i)
        M_vec->set(start_i+i,m->operator()(i));
}


template class VectorBlockBase<double>;

} // Feel



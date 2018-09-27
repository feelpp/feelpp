/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 17 May 2015

 Copyright (C) 2015 Feel++ Consortium

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
#ifndef FEELPP_CSRGRAPHMETIS_HPP
#define FEELPP_CSRGRAPHMETIS_HPP 1

#include <feel/feelcore/feel.hpp>


namespace Feel {


/**
 * This utility class provides a convenient implementation for
 * building the compressed-row-storage graph required for the METIS/ParMETIS
 * graph partitioning schemes.
 */
template <class IndexType>
class CSRGraphMetis
{
public:
    std::vector<IndexType> offsets, vals;

    void prepareNumberNonZeros(const dof_id_type row,
                               const dof_id_type n_nonzeros_in)
        {
            CHECK( row+1 < offsets.size() ) << "Invalid row id " << row+1 << " greater than offset size " << offsets.size();
            offsets[row+1] = n_nonzeros_in;
        }



    dof_id_type
    numberNonZero (const dof_id_type row) const
        {
            CHECK( row+1 < offsets.size() ) << "Invalid row id " << row+1 << " greater than offset size " << offsets.size();
            return (offsets[row+1] - offsets[row]);
        }


    void prepareForUse()
        {
            std::partial_sum (offsets.begin(), offsets.end(), offsets.begin());

            CHECK( !offsets.empty() ) << "Invalid offsets in CSR Graph for Metis";

            vals.resize(offsets.back());

            if (vals.empty())
                vals.push_back(0);
        }



    IndexType& operator()( const dof_id_type row, const dof_id_type nonzero )
        {
            CHECK( vals.size() > offsets[row]+nonzero ) << "Invalid row and associated offset " << row << " (offset:" << offsets[row] << ")";

            return vals[offsets[row]+nonzero];
        }



    const IndexType& operator()(const dof_id_type row, const dof_id_type nonzero) const
        {
            CHECK( vals.size() > offsets[row]+nonzero ) << "Invalid row and associated offset " << row << " (offset:" << offsets[row] << ")";
            return vals[offsets[row]+nonzero];
        }

};


}
#endif

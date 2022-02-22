/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-07-20

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
/**
   \file opusmodelfactory.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-07-20
 */
#include "opusmodelfactory.hpp"
#include "opusmodel.hpp"

#define ALL_MODELS 0
namespace Feel
{
OpusModelFactory::opusmodel_ptrtype
OpusModelFactory::New( Feel::po::variables_map const& vm )
{
    using namespace Feel;
    int order_temp = vm["order-temp"].as<int>();

    switch ( order_temp  )
    {
    case 1:
#if ALL_MODELS
        return opusmodel_ptrtype( new OpusModel<2,1,1>( vm ) );
#endif
#if OPUS_WITH_THERMAL_ORDER >= 2

    case 2:
        return opusmodel_ptrtype( new OpusModel<2,1,2>( vm ) );
#endif
#if ALL_MODELS
#if OPUS_WITH_THERMAL_ORDER >= 3

    case 3:
        return opusmodel_ptrtype( new OpusModel<2,1,3>( vm ) );
#endif
#if OPUS_WITH_THERMAL_ORDER >= 4

    case 4:
        return opusmodel_ptrtype( new OpusModel<2,1,4>( vm ) );
#endif
#endif

    default:
        std::cout << "Approximation order for temperature: " << order_temp << "\n";
        std::cout << "fallback to order 1\n";
        //return opusmodel_ptrtype( new OpusModel<2,1,1>() );


    }
    LOG(FATAL) << "model for temperature approximation order " << order_temp << " does not exist\n";
    return opusmodel_ptrtype();
}
}

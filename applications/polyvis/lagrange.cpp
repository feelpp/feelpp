/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2006-02-01

  Copyright (C) 2006 EPFL

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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file lagrange.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2006-02-01
 */
#include <boost/any.hpp>
#include <map>
#include <utility>

#include <boost/plugin/export_plugin.hpp>

#include <polyvis.hpp>

#include <life/lifepoly/polynomialset.hpp>
#include <life/lifepoly/lagrange.hpp>

class Lagrange
    :
    public Polyvis
{
public:
    Lagrange(std::string const& s, int d = 2) : Polyvis(),name(s)
    {
        std::cout << "Created lagrange '" << s << "'\n";
        std::cout << "dimension is " << d << "\n";

        using namespace Life;

        if ( d == 2 )
            {
                this->save( "lag_2_1", fem::Lagrange<2,1,Scalar>().functionShape() );
            }
        if ( d == 3 )
            {
                this->save( "lag_3_1", fem::Lagrange<3,1,Scalar>().functionShape() );
            }
    }
    Lagrange(std::string const& s, int d, int o) : Polyvis(),name(s)
    {
        std::cout << "Created lagrange '" << s << "'\n";
        std::cout << "dimension is " << d << "\n";
        std::cout << "order is " << o << "\n";

        using namespace Life;

        if ( d == 2 )
            {
                switch( o )
                    {
                    case 1:
                        { this->save( "lag_2_1", fem::Lagrange<2,1,Scalar>().functionShape() ); }
                        break;
                    case 2:
                        { this->save( "lag_2_2", fem::Lagrange<2,2,Scalar>().functionShape() ); }
                        break;
                    case 3:
                        { this->save( "lag_2_3", fem::Lagrange<2,3,Scalar>().functionShape() ); }
                        break;
                    case 4:
                        { this->save( "lag_2_4", fem::Lagrange<2,4,Scalar>().functionShape() ); }
                        break;
                    default:
                        { this->save( "lag_2_1", fem::Lagrange<2,1,Scalar>().functionShape() ); }
                    }

            }
        if ( d == 3 )
            {
                switch( o )
                    {
                    case 1:
                        { this->save( "lag_3_1", fem::Lagrange<3,1,Scalar>().functionShape() ); }
                        break;
                    case 2:
                        { this->save( "lag_3_2", fem::Lagrange<3,2,Scalar>().functionShape() ); }
                        break;
                    case 3:
                        { this->save( "lag_3_3", fem::Lagrange<3,3,Scalar>().functionShape() ); }
                        break;
                    case 4:
                        { this->save( "lag_3_4", fem::Lagrange<3,4,Scalar>().functionShape() ); }
                        break;
                    default:
                        { this->save( "lag_3_1", fem::Lagrange<3,1,Scalar>().functionShape() ); }
                    }

            }
    }

    ~Lagrange()
    {}

    std::string name;
};

BOOST_PLUGIN_EXPORT(Polyvis, Lagrange, "Lagrange");



/* -*- mode: c++ -*-

  This file is part of the LifeV library

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
   \file raviartthomas.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2006-02-01
 */
#include <boost/any.hpp>
#include <map>
#include <utility>

#include <boost/plugin/export_plugin.hpp>

#include <polyvis.hpp>


#include <life/lifepoly/raviartthomas.hpp>

class RaviartThomas
    :
    public Polyvis
{
public:
    RaviartThomas(std::string const& s, int d = 2) : Polyvis(),name(s)
    {
        std::cout << "Created raviartthomas '" << s << "'\n";
        std::cout << "dimension is " << d << "\n";

        using namespace LifeV;

        if ( d == 2 )
            {
                this->save( "rt_2_1", fem::RaviartThomas<2,1>().functionShape() );
            }
        if ( d == 3 )
            {
                this->save( "rt_3_1", fem::RaviartThomas<3,1>().functionShape() );
            }
    }
    RaviartThomas(std::string const& s, int d, int o) : Polyvis(),name(s)
    {
        std::cout << "Created raviart-thomas '" << s << "'\n";
        std::cout << "dimension is " << d << "\n";
        std::cout << "order is " << o << "\n";

        using namespace LifeV;

        if ( d == 2 )
            {
                switch( o )
                    {
                    case 1:
                        { this->save( "rt_2_1", fem::RaviartThomas<2,1>().functionShape() ); }
                        break;
                    case 2:
                        { this->save( "rt_2_2", fem::RaviartThomas<2,2>().functionShape() ); }
                        break;
                    case 3:
                        { this->save( "rt_2_3", fem::RaviartThomas<2,3>().functionShape() ); }
                        break;
                    case 4:
                        { this->save( "rt_2_4", fem::RaviartThomas<2,4>().functionShape() ); }
                        break;
                    default:
                        { this->save( "rt_2_1", fem::RaviartThomas<2,1>().functionShape() ); }
                    }

            }
        if ( d == 3 )
            {
                switch( o )
                    {
                    case 1:
                        { this->save( "rt_3_1", fem::RaviartThomas<3,1>().functionShape() ); }
                        break;
                    case 2:
                        { this->save( "rt_3_2", fem::RaviartThomas<3,2>().functionShape() ); }
                        break;
                    case 3:
                        { this->save( "rt_3_3", fem::RaviartThomas<3,3>().functionShape() ); }
                        break;
                    case 4:
                        { this->save( "rt_3_4", fem::RaviartThomas<3,4>().functionShape() ); }
                        break;
                    default:
                        { this->save( "rt_3_1", fem::RaviartThomas<3,1>().functionShape() ); }
                    }

            }
    }

    ~RaviartThomas()
    {}

    std::string name;
};

BOOST_PLUGIN_EXPORT(Polyvis, RaviartThomas, "RaviartThomas");



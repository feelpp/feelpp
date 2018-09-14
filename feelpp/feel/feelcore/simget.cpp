/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-05-24

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file simget.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-05-24
 */

#include <feel/feelcore/application.hpp>
#include <feel/feelcore/simget.hpp>
#include <feel/feelcore/environment.hpp>

namespace Feel
{
Simget::Simget()
    :
    M_vm( Environment::vm() ),
    M_about( Environment::about() )
{}

Simget&
Simget::changeRepository( boost::format fmt )
{
    Environment::changeRepository( fmt );
    return  *this;
}

void
Simget::print( std::ostream& out, std::vector<ptree::ptree> & stats )
{
    const std::string key = M_about.appName();

    if ( stats.front().count( "e" ) )
    {
        out << std::setw( 10 ) << std::right << "levels"
            << std::setw( 10 ) << std::right << "h";
        BOOST_FOREACH( ptree::ptree::value_type &v,
                       stats.front().get_child( "e" ) )
        {
            out << std::setw( 15 ) << std::right << v.first
                << std::setw( 15 ) << std::right << "ROC";
        }
        out << "\n" ;
        int l=1;

        for ( auto it = stats.begin(), en =stats.end(); it!=en; ++it,++l )
        {
            //std::for_each( it->begin(),it->end(), []( std::pair<std::string,boost::any> const& o ) { std::cout << o.first << "\n"; } );
            //std::map<std::string,boost::any> data = *it;
            //std::map<std::string,boost::any> datap;
            double h  = it->get<double>( "h" );
#if 0
            double rocu = 1, rocp=1;

            double u  = boost::any_cast<double>( data.find( "||u_error||_L2" )->second );
            double p  =  boost::any_cast<double>( data.find( "||p_error||_L2" )->second );

            if ( l > 1 )
            {
                datap = *boost::prior( it );

                double hp  = boost::any_cast<double>( datap.find( "h" )->second );
                double up  = boost::any_cast<double>( datap.find( "||u_error||_L2" )->second );
                double pp  = boost::any_cast<double>( datap.find( "||p_error||_L2" )->second );
                rocu = std::log10( up/u )/std::log10( hp/h );
                rocp = std::log10( pp/p )/std::log10( hp/h );
            }

#endif
            out << std::right << std::setw( 10 ) << l
                << std::right << std::setw( 10 ) << std::fixed  << std::setprecision( 4 ) << h
                //<< std::right << std::setw(15) << std::scientific << std::setprecision( 2 ) << u
                //<< std::right << std::setw(15) << std::fixed << std::setprecision( 2 ) << rocu
                //<< std::right << std::setw(15) << std::scientific << std::setprecision( 2 ) << p
                //<< std::right << std::setw(15) << std::fixed << std::setprecision( 2 ) << rocp
                << "\n";
        }
    }

#if 0
    out << std::setw( 10 ) << std::right << "levels"
        << std::setw( 10 ) << std::right << "h"
        << std::setw( 10 ) << std::right << "nElts"
        << std::setw( 10 ) << std::right << "nDof"
        << std::setw( 10 ) << std::right << "nDofu"
        << std::setw( 10 ) << std::right << "nDofp"
        << std::setw( 10 ) << std::right << "nNz"
        << std::setw( 10 ) << std::right << "T_space"
        << std::setw( 10 ) << std::right << "T_matrix"
        << std::setw( 15 ) << std::right << "T_m_assembly"
        << std::setw( 15 ) << std::right << "T_v_assembly"
        << std::setw( 15 ) << std::right << "T_assembly"
        << std::setw( 10 ) << std::right << "T_solve"
        << "\n" ;
    l=1;

    for ( auto it = stats.begin(), en =stats.end(); it!=en; ++it,++l )
    {
        //std::for_each( it->begin(),it->end(), []( std::pair<std::string,boost::any> const& o ) { std::cout << o.first << "\n"; } );
        std::map<std::string,boost::any> data = *it;
        std::map<std::string,boost::any> datap;
        double rocu = 1, rocp=1;
        double h  = boost::any_cast<double>( data.find( "h" )->second ) ;
        using namespace Feel;
        out << std::right << std::setw( 10 ) << l
            << std::right << std::setw( 10 ) << std::fixed  << std::setprecision( 4 ) << h
            << std::right << std::setw( 10 ) << boost::any_cast<size_type>( data.find( "space.nelts" )->second )
            << std::right << std::setw( 10 ) << boost::any_cast<size_type>( data.find( "space.ndof" )->second )
            << std::right << std::setw( 10 ) << boost::any_cast<size_type>( data.find( "space.ndof.u" )->second )
            << std::right << std::setw( 10 ) << boost::any_cast<size_type>( data.find( "space.ndof.p" )->second )
            << std::right << std::setw( 10 ) << boost::any_cast<size_type>( data.find( "matrix.nnz" )->second )
            << std::right << std::setw( 10 ) << std::scientific << std::setprecision( 2 ) << boost::any_cast<double>( data.find( "space.time" )->second )
            << std::right << std::setw( 10 ) << std::scientific << std::setprecision( 2 ) << boost::any_cast<double>( data.find( "matrix.init" )->second )
            << std::right << std::setw( 15 ) << std::scientific << std::setprecision( 2 ) << boost::any_cast<double>( data.find( "matrix.assembly" )->second )
            << std::right << std::setw( 15 ) << std::scientific << std::setprecision( 2 ) << boost::any_cast<double>( data.find( "vector.assembly" )->second )
            << std::right << std::setw( 15 ) << std::scientific << std::setprecision( 2 ) << boost::any_cast<double>( data.find( "vector.assembly" )->second )+boost::any_cast<double>( data.find( "matrix.assembly" )->second )
            << std::right << std::setw( 10 ) << std::scientific << std::setprecision( 2 ) << boost::any_cast<double>( data.find( "solver.time" )->second )
            << "\n";
    }

#endif // 0
}

}

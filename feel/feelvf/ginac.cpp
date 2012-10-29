/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2012-10-24

  Copyright (C) 2012 Université de Strasbourg

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
   \file ginac.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2012-10-24
 */
#include <ginac/ginac.h>

namespace GiNaC
{
matrix
grad( ex const& f, std::vector<symbol> const& l )
{
	lst g;
	std::for_each( l.begin(), l.end(), [&] ( symbol const& x ) { g.append( f.diff( x ) ); } );
	return matrix( 1, l.size(), g );
}
matrix
grad( matrix const& f, std::vector<symbol> const& l )
{
	matrix ff = f;
	if ( f.rows() == 1 )
		ff = f.transpose();

	matrix g( ff.rows(), l.size() );
	for( int i = 0; i < f.rows(); ++i )
	{
		for( int j = 0; j < l.size(); ++j )
		{
			g.set( i, j, f(i,0).diff( l[j] ) );
		}
	}
	return g;
}

matrix
div( matrix const& f, std::vector<symbol> const& l )
{
	matrix ff = f;
	if ( f.cols() == 1 )
		ff = f.transpose();

	matrix g(ff.rows(), 1 );
	for( int i = 0; i < ff.rows(); i++ )
	{
		ex e;
		for( int j = 0; j < ff.cols(); j++ )
			e += ff(i,j).diff( l[j] );
		g(i,0) = e;
	}
	return g;
}

ex
laplacian( ex const& f, std::vector<symbol> const& l )
{
	ex g;
	std::for_each( l.begin(), l.end(),
		       [&] ( symbol const& x )
		       {
			       g += f.diff( x,2 );
		       } );
	return g;
}

matrix
laplacian( matrix const& f, std::vector<symbol> const& l )
{
	matrix g(f.rows(),1);
	for(int i = 0; i < f.rows(); ++i )
	{
		ex lap;
		std::for_each( l.begin(), l.end(),
			       [&] ( symbol const& x )
			       {
				       lap += f(i,0).diff( x,2 );
			       } );
		g(i,0) = lap;
	}
	return g;
}

} // GiNaC


/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Gilles Steiner <gilles.steiner@epfl.ch>
       Date: 2005-11-10

  Copyright (C) 2005,2006 EPFL

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
   @file pointset.hpp
   @author Gilles Steiner
   @date 2005-11-10
 */

#ifndef __PointSet_H
#define __PointSet_H 1


#include <feel/feelmesh/refentity.hpp>
#include <feel/feelcore/visitor.hpp>


namespace Feel
{
namespace ublas = boost::numeric::ublas;

/** @class PointSet
 *  @brief  Class of all PointSet on a Convex.
 *  @author Gilles Steiner
 *
 *   This is a template class to represent every type of Point Set we may encounter,
 *   for example : Equidistant points, Fekete points, quadrature points (Gauss, Gauss-Lobatto, etc).
 *   The class will then be specialised into QuadraturePointSet.
 *
 */


template<class Convex, typename T>
class PointSet : public VisitableBase<>
{
    typedef VisitableBase<> super;

public:

    typedef Convex convex_type;
    typedef T value_type;
    typedef PointSet<convex_type, value_type> self_type;

    typedef Reference<Convex, Convex::nDim, Convex::nOrder, Convex::nDim, value_type> RefElem;

    typedef typename node<value_type>::type node_type;
    typedef ublas::matrix<value_type, ublas::column_major> nodes_type;

    PointSet()
        :
        super(),
        M_npoints( 0 ),
        M_points(),
        M_points_face( Convex::numTopologicalFaces )
    {}

    PointSet( const self_type& P )
        :
        super(),
        M_npoints( P.nPoints() ),
        M_points( P.points() ),
        M_points_face( P.M_points_face )
    {}

    PointSet( uint32_type Npoints )
        :
        super(),
        M_npoints( Npoints ),
        M_points( Convex::nDim, Npoints ),
        M_points_face( Convex::numTopologicalFaces )
    {}

    PointSet( uint32_type Npoints, uint16_type Dim )
        :
        super(),
        M_npoints( Npoints ),
        M_points( Dim, Npoints ),
        M_points_face( Convex::numTopologicalFaces )
    {}

    PointSet( nodes_type const& SomePoints )
        :
        super(),
        M_npoints( SomePoints.size2() ),
        M_points( SomePoints ),
        M_points_face( Convex::numTopologicalFaces )
    {}

    PointSet( std::vector<std::map<uint16_type,nodes_type> > const& SomePoints, uint16_type perm=1)
        :
        super(),
        M_npoints( 0 ),
        M_points(),
        M_points_face( Convex::numTopologicalFaces )
    {
        if (SomePoints.size()==0) return;

        auto it=SomePoints.begin();
        auto const en=SomePoints.end();
        uint32_type size1=it->find( perm )->second.size1(), size2=0;
        for ( ; it!=en ; ++it)
            size2+=it->find( perm )->second.size2();
        M_npoints = size2;
        M_points.resize(size1,size2);

        it=SomePoints.begin();
        uint32_type start=0;
        for ( ; it!=en ; ++it)
        {
            auto const& thenodes = it->find( perm )->second;
            for ( uint32_type i=0 ; i< thenodes.size2();++i )
                ublas::column( M_points, start + i ) = ublas::column( thenodes, i );
            start+=thenodes.size2();
        }
    }

    virtual ~PointSet()
    {}

    self_type& operator=( self_type const& p )
    {
        if ( this != &p )
        {
            M_npoints = p.M_npoints;
            M_points = p.M_points;
            M_points_face = p.M_points_face;
        }

        return *this;
    }


    uint32_type nPoints() const
    {
        return M_npoints;
    }
    nodes_type const& points() const
    {
        return M_points;
    }
    ublas::matrix_column<nodes_type const> point( uint32_type __i ) const
    {
        return ublas::column( M_points, __i );
    }
    ublas::matrix_column<nodes_type> point( uint32_type __i )
    {
        return ublas::column( M_points, __i );
    }


    nodes_type const& points( uint16_type f ) const
    {
        return M_points_face[f];
    }
    ublas::matrix_column<nodes_type const> point( uint16_type f, uint32_type __i ) const
    {
        return ublas::column( M_points_face[f], __i );
    }
    ublas::matrix_column<nodes_type> point( uint16_type f, uint32_type __i )
    {
        return ublas::column( M_points_face[f], __i );
    }

    void setName( std::string name, uint32_type order )
    {
        std::ostringstream ostr;

        ostr << Convex::nDim
             << "_" << order
             << "_" << Convex::nRealDim;

        pointsInfo = ostr.str();

        pointsName = name;
    }

    std::string getName()
    {
        return pointsName;
    }

    std::string getPointsInfo()
    {
        return pointsInfo;
    }

    void toPython()
    {
        RefElem RefConv;

        std::ostringstream ostr;

        ostr << getName() << "_" << getPointsInfo() << ".py";

        std::ofstream ofs( ostr.str().c_str() );

        ofs << "from pyx import *\n";

        ofs << "p=path.path(";

        for ( uint16_type i = 0; i < Convex::numEdges; ++i )
        {
            for ( uint16_type  j = 0; j < 2; ++j )
            {
                node_type x( 2 );

                if ( Convex::nRealDim == 1 )
                {
                    x( 0 ) = RefConv.edgeVertex( i,j )( 0 );
                    x( 1 ) = value_type( 0 );
                }

                if ( Convex::nRealDim == 2 )
                {
                    x = RefConv.edgeVertex( i, j );
                }

                if ( Convex::nRealDim == 3 )
                {
                    x( 0 ) = RefConv.edgeVertex( i, j )( 0 )+RefConv.edgeVertex( i, j )( 1 )*std::cos( M_PI/4 );
                    x( 1 ) = RefConv.edgeVertex( i, j )( 2 )+RefConv.edgeVertex( i, j )( 1 )*std::sin( M_PI/4 );
                }

                if ( j == 0 )
                    ofs << "path.moveto(" << double( x( 0 ) )<< "," << double( x( 1 ) ) << "),\n";

                else if ( j == 1 )
                    ofs << "path.lineto(" << double( x( 0 ) )<< "," << double( x( 1 ) ) << "),\n";
            }
        }

        ofs << "path.closepath() )\n";

        ofs << "text.set(mode=\"latex\")\n"
            << "c = canvas.canvas()\n"
            << "c.stroke(p, [style.linewidth.Thin])\n";

        for ( uint16_type i = 0; i < nPoints(); ++i )
        {
            node_type x( 2 );

            if ( Convex::nRealDim == 1 )
            {
                x( 0 ) = this->point( i )( 0 );
                x( 1 ) = value_type( 0 );
            }

            if ( Convex::nRealDim == 2 )
            {
                x = this->point( i );
            }

            if ( Convex::nRealDim == 3 )
            {
                x( 0 ) = this->point( i )( 0 ) + this->point( i )( 1 )*std::cos( M_PI/4 );
                x( 1 ) = this->point( i )( 2 ) + this->point( i )( 1 )*std::sin( M_PI/4 );
            }


            ofs << "c.fill ( path.circle(" << double( x( 0 ) ) << "," << double( x( 1 ) )<< ", 0.02 ),[deco.filled([color.rgb.red])])\n";
            ofs << "c.text( " << double( x( 0 )+0.025 ) << "," << double( x( 1 )+0.025 )<< ", r\"{\\scriptsize " << i << "}\")\n";
        }

        ofs << "c.writeEPSfile(\"" << getName() << "_" << getPointsInfo()
            << "\", document.paperformat.A4)\n";
    }

    FEELPP_DEFINE_VISITABLE();

protected:

    /**
     * set the points of the pointset
     */
    void setPoints( nodes_type const& pts )
    {
        M_points = pts;
        M_npoints = pts.size2();
    }

    /**
     * set the pointset at face \c f using nodes \c n
     */
    void setPoints( uint16_type f, nodes_type const& n )
    {
        M_points_face[f] = n;
    }
protected:

    uint32_type M_npoints;
    nodes_type M_points;
    std::vector<nodes_type> M_points_face;

    //Identifies if points are equispaced, warpblend, fekete
    std::string pointsName;

    //Identifies the Order, dim and realdim
    std::string pointsInfo;

};

} // Feel

#endif /* _PointSet_H */

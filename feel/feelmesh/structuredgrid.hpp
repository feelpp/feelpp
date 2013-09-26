/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel Imaging project.

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-02-28

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
   \file structuredgrid.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-02-28
 */
#ifndef __StructuredGrid_H
#define __StructuredGrid_H 1
#include <iostream>
#include <fstream>
#include <iomanip>

#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include <feel/feelalg/glas.hpp>

namespace Feel
{
/// \cond disabled
/**
 * @class bbox
 * @brief bounding box for the mesh
 *
 * @author Christophe Prud'homme
 * @sa StructuredGrid
 */
template<uint Dim>
struct bbox
{
    /**
     * default bounding box contructor
     */
    bbox()
        :
        min( Dim ),
        max( Dim )
    {
        min = ublas::scalar_vector<double>( Dim, -1 );
        max = ublas::scalar_vector<double>( Dim, 1 );
    }

    /**
     * constructs the bounding box of a cloud of nodes
     *
     * @param __cp cloud of nodes
     */
    bbox( ublas::vector<node_type> const& __cp )
        :
        min( Dim ),
        max( Dim ),
        eps( ublas::scalar_vector<double>( min.size(), 1e-10 ) )
    {
        min = __cp( 0 );
        max = __cp( 0 );

        for ( int __i = 1; __i < __cp.size(); ++__i )
        {
            for ( int __d = 0; __d < Dim; ++__d )
            {
                min[__d] = ( __cp( __i )[__d] > min[__d] )?min[__d]: __cp( __i )[__d];
                max[__d] = ( __cp( __i )[__d] < max[__d] )?max[__d]: __cp( __i )[__d];
            }
        }

        // enlarge a tad the bounding box
        min -= eps;
        max += eps;
    }

    /**
     * constructs the bounding box of a cloud of nodes
     *
     * @param bboxes set of bounding boxes
     */
    bbox( std::vector<boost::shared_ptr<bbox<Dim> > > const& bboxes )
        :
        min( Dim ),
        max( Dim ),
        eps( ublas::scalar_vector<double>( min.size(), 1e-10 ) )
    {
        min = bboxes[0]->min;
        max = bboxes[0]->max;

        for ( int __i = 1; __i < bboxes.size(); ++__i )
        {
            for ( int __d = 0; __d < Dim; ++__d )
            {
                min[__d] = ( bboxes[__i]->min[__d] > min[__d] )?min[__d]: bboxes[ __i ]->min[__d];
                max[__d] = ( bboxes[__i]->max[__d] < max[__d] )?max[__d]: bboxes[ __i ]->max[__d];
            }
        }

        // enlarge a tad the bounding box
        min -= eps;
        max += eps;
    }

    /**
     * enlarge the bounding box
     *
     * @param e size of enlargment of the bounding box
     */
    void enlarge( double e )
    {
        min -= ublas::scalar_vector<double>( min.size(), e );
        max += ublas::scalar_vector<double>( min.size(), e );
    }
    node_type min;
    node_type max;
    node_type eps;
};
/// \endcond

/**
 * \class StructuredGrid
 * \brief class to represent a Structured Grid
 *
 * @author Christophe Prud'homme
 */
template<uint16_type Dim, uint16_type Order = 1>
class StructuredGrid
{
public:


    /** @name Typedefs
     */
    //@{

    typedef boost::multi_array<node_type, Dim> coord_type;
    typedef boost::multi_array<double, Dim> f_type;
    typedef typename coord_type::index index_type;
    typedef boost::array<index_type, Dim> shape_type;
    typedef boost::array<index_type, Dim> idx_type;

    typedef typename mpl::if_<mpl::equal<mpl::int_<Dim>, mpl::int_<1> >,
            mpl::identity<GeoElement1D<Dim, LinearLine, DefMarkerCommon> >,
            mpl::if_<mpl::equal<mpl::int_<Dim>, mpl::int_<2> >,
            mpl::identity<GeoElement2D<Dim, LinearQuad, DefMarkerCommon> >,
            mpl::if_<mpl::equal<mpl::int_<Dim>, mpl::int_<3> >,
            mpl::identity<GeoElement3D<Dim, LinearHexa, DefMarkerCommon> > >::type::type element_type;
    typedef typename mpl::if_<mpl::equal<mpl::int_<Dim>, mpl::int_<1> >,
            mpl::identity<GeoElement0D<Dim, DefMarkerCommon> >,
            mpl::if_<mpl::equal<mpl::int_<Dim>, mpl::int_<2> >,
            mpl::identity<GeoElement1D<Dim, LinearLine, DefMarkerCommon> >,
            mpl::if_<mpl::equal<mpl::int_<Dim>, mpl::int_<3> >,
            mpl::identity<GeoElement2D<Dim, LinearQuad, DefMarkerCommon> > >::type::type face_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    StructuredGrid( bbox<Dim> const& __bbox,  shape_type const& __shape )
        :
        M_shape( __shape ),
        M_coord( __shape )
    {
#if 0
        std::for_each( M_coord.begin(), M_coord.end(), lambda::bind( &node_type::resize,
                       lambda::_1,
                       Dim ) );
#else
        idx_type __idx;
        node_type __pt( Dim );

        for ( index_type i = 0; i < M_shape[0]; ++i )
        {
            __idx[0] = i;
            __pt[0] = __bbox.min[0]+i*( __bbox.max[0]-__bbox.min[0] )/( M_shape[0]-1 );

            for ( index_type j = 0; j < M_shape[1]; ++j )
            {
                __idx[1] = j;
                __pt[1] = __bbox.min[1]+j*( __bbox.max[1]-__bbox.min[1] )/( M_shape[1]-1 );

                for ( index_type k = 0; k < M_shape[2]; ++k )
                {
                    __idx[2] = k;
                    __pt[2] = __bbox.min[2]+k*( __bbox.max[2]-__bbox.min[2] )/( M_shape[2]-1 );
                    M_coord( __idx ).resize( Dim );
                    M_coord( __idx ) = __pt;
                }
            }
        }

#endif
    }
    ~StructuredGrid()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     */
    size_type numElements() const
    {
        return M_coord.num_elements();
    }

    /**
     *
     *
     *
     * @return
     */
    shape_type const& shape() const
    {
        return M_shape;
    }

    /**
     *
     *
     * @param __i
     *
     * @return
     */
    uint extent( uint __i ) const
    {
        return M_shape[__i];
    }

    /**
     *
     *
     *
     * @return
     */
    node_type& operator()( idx_type const& __idx )
    {
        return M_coord( __idx );
    }

    /**
     *
     *
     *
     * @return
     */
    node_type operator()( idx_type const& __idx ) const
    {
        return M_coord( __idx );
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     *
     *
     * @param __f
     */
    void save( std::string const& prefix,
               std::vector<boost::shared_ptr<f_type> > const& __f );

    //@}

protected:

    shape_type M_shape;

    coord_type M_coord;

private:

    StructuredGrid();
    StructuredGrid( StructuredGrid const & );

};
template<uint16_type Dim,  uint16_type Order>
void
StructuredGrid<Dim, Order>::save( std::string const& prefix,
                                  std::vector<boost::shared_ptr<f_type> > const& __f )
{
    std::ofstream __ofs( ( prefix + ".geo" ).c_str() );
    __ofs <<
          "Structured grid\n"
          "Structured grid\n"
          "node id off\n"
          "element id off\n"
          "coordinates\n"
          << std::setw( 8 ) << 0 << "\n";

    __ofs <<
          "part 1\n"
          "structured block\n"
          "block\n"
          << std::setw( 8 ) << M_shape[0]
          << std::setw( 8 ) << M_shape[1]
          << std::setw( 8 ) << M_shape[2];

    idx_type __idx;
    size_t __count = 0;

    for ( index_type i = 0; i < M_shape[0]; ++i )
    {
        __idx[0] = i;

        for ( index_type j = 0; j < M_shape[1]; ++j )
        {
            __idx[1] = j;

            for ( index_type k = 0; k < M_shape[2]; ++k )
            {
                __idx[2] = k;

                if ( __count ++ % 6 == 0 )
                    __ofs << "\n";

                __ofs.precision( 5 );
                __ofs.setf( std::ios::scientific );
                __ofs << std::setw( 12 ) << M_coord( __idx )[0];
            }
        }
    }

    //if ( __count-1 % 6 != 0 )
    //    __ofs << "\n";
    __count = 0;

    for ( index_type i = 0; i < M_shape[0]; ++i )
    {
        __idx[0] = i;

        for ( index_type j = 0; j < M_shape[1]; ++j )
        {
            __idx[1] = j;

            for ( index_type k = 0; k < M_shape[2]; ++k )
            {
                __idx[2] = k;

                if ( __count ++ % 6 == 0 )
                    __ofs << "\n";

                __ofs.precision( 5 );
                __ofs.setf( std::ios::scientific );
                __ofs << std::setw( 12 ) << M_coord( __idx )[1];
            }
        }
    }

    //if ( __count-1 % 6 != 0 )
    //    __ofs << "\n";
    __count = 0;

    for ( index_type i = 0; i < M_shape[0]; ++i )
    {
        __idx[0] = i;

        for ( index_type j = 0; j < M_shape[1]; ++j )
        {
            __idx[1] = j;

            for ( index_type k = 0; k < M_shape[2]; ++k )
            {
                __idx[2] = k;

                if ( __count ++ % 6 == 0 )
                    __ofs << "\n";

                __ofs.precision( 5 );
                __ofs.setf( std::ios::scientific );
                __ofs << std::setw( 12 ) << M_coord( __idx )[2];
            }
        }
    }

    __ofs << "\n";

    __ofs.close();
    __ofs.open( ( prefix + ".case" ).c_str() );
    __ofs <<
          "FORMAT\n"
          "type: ensight\n"
          "GEOMETRY\n"
          "model: toto.geo\n"
          "VARIABLE\n";

    for ( int i = 0; i < __f.size(); ++i )
    {
        __ofs << "scalar per node: f" << i << " f" << i << ".001\n";
    }

    __ofs.close();

    for ( int l = 0; l < __f.size(); ++l )
    {
        std::ostringstream __ostr;
        __ostr << "f" << l << ".001";
        __ofs.open( __ostr.str().c_str() );

        __ofs.precision( 5 );
        __ofs.setf( std::ios::scientific );
        __ofs << "function\n"
              << "part 1\n"
              << "block";
        __count = 0;

        for ( index_type i = 0; i < M_shape[0]; ++i )
        {
            __idx[0] = i;

            for ( index_type j = 0; j < M_shape[1]; ++j )
            {
                __idx[1] = j;

                for ( index_type k = 0; k < M_shape[2]; ++k )
                {
                    __idx[2] = k;

                    if ( __count ++ % 6 == 0 )
                        __ofs << "\n";

                    __ofs << std::setw( 12 ) << __f[l]->operator()( __idx );
                }
            }
        }

        __ofs.close();
    }
} // end StructuredGrid<>::save

}
#endif /* __StructuredGrid_H */

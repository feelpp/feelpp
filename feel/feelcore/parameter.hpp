/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-04-03

  Copyright (C) 2009 Université de Grenoble 1

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
   \file parameter.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-04-03
 */
#ifndef __feelcore_parameter_H
#define __feelcore_parameter_H 1


#if !defined(BOOST_PARAMETER_MAX_ARITY)
#define BOOST_PARAMETER_MAX_ARITY 10
#endif

#include <boost/parameter.hpp>
#include <boost/type_traits.hpp>
#if 0
#include <boost/parameter/keyword.hpp>
#include <boost/parameter/aux_/maybe.hpp>
#include <boost/parameter/name.hpp>
#include <boost/parameter/preprocessor.hpp>
#endif

namespace Feel
{
namespace parameter = boost::parameter;

BOOST_PARAMETER_NAME( vm )  // Note: no semicolon
BOOST_PARAMETER_NAME( options )
    BOOST_PARAMETER_NAME( about )
    BOOST_PARAMETER_NAME( prefix )
    BOOST_PARAMETER_NAME( sub )
    BOOST_PARAMETER_NAME( opt )
    BOOST_PARAMETER_NAME( path )
    BOOST_PARAMETER_NAME( suffix )
    BOOST_PARAMETER_NAME( filename )
    BOOST_PARAMETER_NAME( sep )
    BOOST_PARAMETER_NAME( directory )
    BOOST_PARAMETER_NAME( subdir )
    BOOST_PARAMETER_NAME( format )
    BOOST_PARAMETER_NAME( argc )
    BOOST_PARAMETER_NAME( argv )

    BOOST_PARAMETER_NAME( verbose )
    BOOST_PARAMETER_NAME( threading )


    BOOST_PARAMETER_NAME( matrix )
    BOOST_PARAMETER_NAME( buildGraphWithTranspose )
    BOOST_PARAMETER_NAME( matrixA )
    BOOST_PARAMETER_NAME( matrixB )
    BOOST_PARAMETER_NAME( formA )
    BOOST_PARAMETER_NAME( formB )
    BOOST_PARAMETER_NAME( rhs )
    BOOST_PARAMETER_NAME( solution )
    BOOST_PARAMETER_NAME( prec )
    BOOST_PARAMETER_NAME( transpose )
    BOOST_PARAMETER_NAME( reuse_prec )
    BOOST_PARAMETER_NAME( reuse_jac )
    BOOST_PARAMETER_NAME( maxit )
    BOOST_PARAMETER_NAME( tolerance )
    BOOST_PARAMETER_NAME( rtolerance )
    BOOST_PARAMETER_NAME( atolerance )
    BOOST_PARAMETER_NAME( dtolerance )
    BOOST_PARAMETER_NAME( stolerance )
    BOOST_PARAMETER_NAME( ksp )
    BOOST_PARAMETER_NAME( pc )
    BOOST_PARAMETER_NAME( pcfactormatsolverpackage )
    BOOST_PARAMETER_NAME( constant_null_space )
    BOOST_PARAMETER_NAME( null_space )
    BOOST_PARAMETER_NAME( near_null_space )
    BOOST_PARAMETER_NAME( test )
    BOOST_PARAMETER_NAME( trial )
    BOOST_PARAMETER_NAME( vector )
    BOOST_PARAMETER_NAME( pattern )
    BOOST_PARAMETER_NAME( pattern_block )
    BOOST_PARAMETER_NAME( diag_is_nonzero )
    BOOST_PARAMETER_NAME( block )
    BOOST_PARAMETER_NAME( copy_values )
    BOOST_PARAMETER_NAME( properties )
    BOOST_PARAMETER_NAME( do_threshold )
    BOOST_PARAMETER_NAME( threshold )
    BOOST_PARAMETER_NAME( init )
    BOOST_PARAMETER_NAME( rowstart )
    BOOST_PARAMETER_NAME( colstart )
    BOOST_PARAMETER_NAME( name )
    BOOST_PARAMETER_NAME( nev )
    BOOST_PARAMETER_NAME( ncv )
    BOOST_PARAMETER_NAME( backend )
    BOOST_PARAMETER_NAME( problem )
    BOOST_PARAMETER_NAME( solver )
    BOOST_PARAMETER_NAME( spectrum )
    BOOST_PARAMETER_NAME( transform )
// parameter for exporter
    BOOST_PARAMETER_NAME( geo )
    BOOST_PARAMETER_NAME( fileset )
// parameter for description of geometries
    BOOST_PARAMETER_NAME( h )
    BOOST_PARAMETER_NAME( dim )
    BOOST_PARAMETER_NAME( order )
    BOOST_PARAMETER_NAME( geo_parameters )
    BOOST_PARAMETER_NAME( in_memory )
    BOOST_PARAMETER_NAME( addmidpoint )
    BOOST_PARAMETER_NAME( usenames )
    BOOST_PARAMETER_NAME( xmin )
    BOOST_PARAMETER_NAME( xmax )
    BOOST_PARAMETER_NAME( ymin )
    BOOST_PARAMETER_NAME( ymax )
    BOOST_PARAMETER_NAME( zmin )
    BOOST_PARAMETER_NAME( zmax )
    BOOST_PARAMETER_NAME( nx )
    BOOST_PARAMETER_NAME( ny )
    BOOST_PARAMETER_NAME( nz )
    BOOST_PARAMETER_NAME( refine )
    BOOST_PARAMETER_NAME( update )
    BOOST_PARAMETER_NAME( physical_are_elementary_regions )
    BOOST_PARAMETER_NAME( parametricnodes )
    BOOST_PARAMETER_NAME( force_rebuild )
    BOOST_PARAMETER_NAME( rebuild )
    BOOST_PARAMETER_NAME( shear )
    BOOST_PARAMETER_NAME( recombine )
    BOOST_PARAMETER_NAME( files_path )
    BOOST_PARAMETER_NAME( depends )
    BOOST_PARAMETER_NAME( optimize3d_netgen )
// parameter for adapt
    BOOST_PARAMETER_NAME( model )
    BOOST_PARAMETER_NAME( geotracking )
    BOOST_PARAMETER_NAME( snapthickness )
    BOOST_PARAMETER_NAME( statistics )
    BOOST_PARAMETER_NAME( hmin )
    BOOST_PARAMETER_NAME( hmax )
    BOOST_PARAMETER_NAME( collapseOnBoundary )
    BOOST_PARAMETER_NAME( collapseOnBoundaryTolerance )
// parameter for xmlParse
    BOOST_PARAMETER_NAME( kind )
    BOOST_PARAMETER_NAME( type )
    BOOST_PARAMETER_NAME( latex )
    BOOST_PARAMETER_NAME( cmdName )
    BOOST_PARAMETER_NAME( values )
    BOOST_PARAMETER_NAME( dependencies )
    BOOST_PARAMETER_NAME( funcs )
    BOOST_PARAMETER_NAME( mesh )
    BOOST_PARAMETER_NAME( geoentity )
    BOOST_PARAMETER_NAME( pointset )
    BOOST_PARAMETER_NAME( desc )
    BOOST_PARAMETER_NAME( desc_lib )
    BOOST_PARAMETER_NAME( shape )
    BOOST_PARAMETER_NAME( convex )
// project and integrate
    BOOST_PARAMETER_NAME( sum )
    BOOST_PARAMETER_NAME( accumulate )
    BOOST_PARAMETER_NAME( geomap )
    BOOST_PARAMETER_NAME( straighten )
    BOOST_PARAMETER_NAME( expr )
    BOOST_PARAMETER_NAME( grad_expr)
    BOOST_PARAMETER_NAME( div_expr)
    BOOST_PARAMETER_NAME( curl_expr)
    BOOST_PARAMETER_NAME( pset )
    BOOST_PARAMETER_NAME( quad )
    BOOST_PARAMETER_NAME( quad1 )
    BOOST_PARAMETER_NAME( arg )

    BOOST_PARAMETER_NAME( quadptloc )

    BOOST_PARAMETER_NAME( extended_doftable )



// orders
    BOOST_PARAMETER_NAME( order_u )
    BOOST_PARAMETER_NAME( order_p )

    BOOST_PARAMETER_NAME( initial_time )
    BOOST_PARAMETER_NAME( final_time )
    BOOST_PARAMETER_NAME( time_step )
    BOOST_PARAMETER_NAME( strategy )
    BOOST_PARAMETER_NAME( steady )
    BOOST_PARAMETER_NAME( restart )
    BOOST_PARAMETER_NAME( restart_path )
    BOOST_PARAMETER_NAME( restart_at_last_save )
    BOOST_PARAMETER_NAME( rank_proc_in_files_name )
    BOOST_PARAMETER_NAME( freq )

    BOOST_PARAMETER_NAME( markerName )
    BOOST_PARAMETER_NAME( markerAll )
    BOOST_PARAMETER_NAME( marker1 )
    BOOST_PARAMETER_NAME( marker2 )
    BOOST_PARAMETER_NAME( marker3 )
    BOOST_PARAMETER_NAME( marker4 )
    BOOST_PARAMETER_NAME( marker5 )
    BOOST_PARAMETER_NAME( marker6 )
    BOOST_PARAMETER_NAME( marker7 )
    BOOST_PARAMETER_NAME( marker8 )
    BOOST_PARAMETER_NAME( marker9 )
    BOOST_PARAMETER_NAME( marker10 )
    BOOST_PARAMETER_NAME( marker11 )
    BOOST_PARAMETER_NAME( marker12 )

    BOOST_PARAMETER_NAME( domain )
    BOOST_PARAMETER_NAME( image )
    BOOST_PARAMETER_NAME( domainSpace )
    BOOST_PARAMETER_NAME( imageSpace )
    BOOST_PARAMETER_NAME( range )
    BOOST_PARAMETER_NAME( range_extended )
    BOOST_PARAMETER_NAME( element )
    BOOST_PARAMETER_NAME( parameter )
    BOOST_PARAMETER_NAME( sampling )
    BOOST_PARAMETER_NAME( context )
    BOOST_PARAMETER_NAME( mpi_communications )

    BOOST_PARAMETER_NAME( components )
    BOOST_PARAMETER_NAME( periodicity )
    BOOST_PARAMETER_NAME( periodic )

    BOOST_PARAMETER_NAME( collect_garbage )

    BOOST_PARAMETER_NAME( partitions )
    BOOST_PARAMETER_NAME( partition_file )
    BOOST_PARAMETER_NAME( respect_partition )
    BOOST_PARAMETER_NAME( rebuild_partitions )
    BOOST_PARAMETER_NAME( rebuild_partitions_filename )
    BOOST_PARAMETER_NAME( worldcomm )
    BOOST_PARAMETER_NAME( worldscomm )
    BOOST_PARAMETER_NAME( parallel )
    BOOST_PARAMETER_NAME( substructuring )
    BOOST_PARAMETER_NAME( structured )

    BOOST_PARAMETER_NAME( jacobian )
    BOOST_PARAMETER_NAME( residual )
    BOOST_PARAMETER_NAME( currentElt )
    BOOST_PARAMETER_NAME( newElt )
    BOOST_PARAMETER_NAME( space )
    BOOST_PARAMETER_NAME( initial_theta )
    BOOST_PARAMETER_NAME( min_theta )
    BOOST_PARAMETER_NAME( forceRelaxation )

    BOOST_PARAMETER_NAME( use_tbb )
    BOOST_PARAMETER_NAME( use_harts )
    BOOST_PARAMETER_NAME( grainsize )
    BOOST_PARAMETER_NAME( partitioner )

    BOOST_PARAMETER_NAME( save )
    BOOST_PARAMETER_NAME( ddmethod )
    BOOST_PARAMETER_NAME( penaldir )

    BOOST_PARAMETER_NAME( close )

    BOOST_PARAMETER_NAME( author )
    BOOST_PARAMETER_NAME( task )
    BOOST_PARAMETER_NAME( email )
    BOOST_PARAMETER_NAME( license )
    BOOST_PARAMETER_NAME( copyright )
    BOOST_PARAMETER_NAME( home )
    BOOST_PARAMETER_NAME( bugs )
    BOOST_PARAMETER_NAME( version )

    BOOST_PARAMETER_NAME( max_points_used )
    BOOST_PARAMETER_NAME( projection )

    BOOST_PARAMETER_NAME( bc )
    BOOST_PARAMETER_NAME( nu )
    BOOST_PARAMETER_NAME( alpha )
} // Feel


namespace Feel
{
namespace detail
{
template<typename TheArgs, typename Tag>
struct remove_pointer_const_reference_type
{
    typedef typename boost::remove_pointer<
    typename boost::remove_const<
    typename boost::remove_reference<
    typename parameter::binding<TheArgs, Tag>::type
    >::type
    >::type
    >::type type;
};
template<typename TheArgs, typename Tag, typename Default>
struct remove_pointer_const_reference_default_type
{
    typedef typename boost::remove_pointer<
    typename boost::remove_const<
    typename boost::remove_reference<
    typename parameter::binding<TheArgs, Tag, Default>::type
    >::type
    >::type
    >::type type;
};
} // detail
} // Feel

#endif /* __feelcore_parameter_H */

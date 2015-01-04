/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@imag.fr>
       Date: 2011-03-03

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
   \file geotool.hpp
   \author Vincent Chabannes <vincent.chabannes@imag.fr>
   \date 2011-03-03
 */


#ifndef FEELPP_GEOTOOL_HPP
#define FEELPP_GEOTOOL_HPP 1

#include <iostream>
#include <string>
#include <sstream>
#include <list>
#include <map>

#include <boost/preprocessor/tuple/elem.hpp>

#include <feel/feelalg/glas.hpp>
//#include <boost/parameter/keyword.hpp>
//#include <boost/parameter/preprocessor.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/straightenmesh.hpp>
#include <feel/feelfilters/importergmsh.hpp>
#include <feel/feelfilters/detail/mesh.hpp>





/*_________________________________________________*/
/*_________________________________________________*/
/*_________________________________________________*/

# define GEOTOOL_SHAPE                                                  \
    ( 19, ( ( Line          , 1, 0, 0, "line"         , 2, LINE        ), \
            ( Triangle      , 2, 1, 0, "triangle"     , 3, TRIANGLE    ), \
            ( Rectangle     , 2, 1, 0, "rectangle"    , 2, RECTANGLE   ), \
            ( Quadrangle    , 2, 1, 0, "quadrangle"   , 4, QUADRANGLE  ), \
            ( Pentagon      , 2, 1, 0, "pentagon"     , 5, PENTAGON    ), \
            ( Hexagon       , 2, 1, 0, "hexagon"      , 6, HEXAGON     ), \
            ( Circle        , 2, 1, 0, "circle"       , 2, CIRCLE      ), \
            ( Ellipse       , 2, 1, 0, "ellipse"      , 3, ELLIPSE     ), \
            ( Pie           , 2, 1, 0, "pie"          , 3, PIE         ), \
            ( Special_1a    , 2, 2, 0, "special_1a"   , 1, SPECIAL_1A  ), \
            ( Special_1b    , 2, 1, 0, "special_1b"   , 1, SPECIAL_1B  ), \
            ( Peanut        , 2, 1, 0, "peanut"       , 4, PEANUT      ), \
            ( Tetrahedron   , 3, 4, 1, "tetrahedron"  , 4, TETRAHEDRON ), \
            ( Hexahedron    , 3, 6, 1, "hexahedron"   , 8, HEXAHEDRON  ), \
            ( Cube          , 3, 6, 1, "cube"         , 2, CUBE        ), \
            ( Cylindre      , 3, 6, 1, "cylindre"     , 4, CYLINDRE    ), \
            ( Sphere        , 3, 8, 1, "sphere"       , 2, SPHERE      ), \
            ( Tube          , 3,20, 4, "tube"         , 5, TUBE        ), \
            ( Special3D_1   , 3, 9, 1, "special3D_1"  , 1, SPECIAL3D_1 ) \
            )                                                           \
      )                                                                 \
    /**/


// Accessors

# define GEOTOOL_SHAPE_NAME_CLASS(i) BOOST_PP_TUPLE_ELEM(7, 0, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_DIM(i) BOOST_PP_TUPLE_ELEM(7, 1, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_NBSURFACE(i) BOOST_PP_TUPLE_ELEM(7, 2, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_NBVOLUME(i) BOOST_PP_TUPLE_ELEM(7, 3, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_NAME_STR(i) BOOST_PP_TUPLE_ELEM(7, 4, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_NBPARAM(i) BOOST_PP_TUPLE_ELEM(7, 5, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))
# define GEOTOOL_SHAPE_NAME_MACRO(i) BOOST_PP_TUPLE_ELEM(7, 6, BOOST_PP_ARRAY_ELEM(i,GEOTOOL_SHAPE))


/*_________________________________________________*/
/*_________________________________________________*/
/*_________________________________________________*/


namespace Feel
{

namespace GeoTool
{

typedef node<double>::type node_type;

class GeoGMSHTool;
typedef boost::shared_ptr< GeoGMSHTool> GeoGMSHTool_ptrtype;

typedef std::map<uint16_type,uint16_type> map_data_type;
typedef std::vector<map_data_type> vec_map_data_type;
typedef boost::shared_ptr<vec_map_data_type> vec_map_data_ptrtype;

//if bool=true => surface stoker dans un tableau gmsh
typedef std::vector<std::map<uint16_type,bool> > vec_map_data_surf1_type;
typedef boost::shared_ptr<vec_map_data_surf1_type> vec_map_data_surf1_ptrtype;
//=> la string est le nom de ce tableau
typedef std::vector<std::map<uint16_type,std::string> > vec_map_data_surf2_type;
typedef boost::shared_ptr<vec_map_data_surf2_type> vec_map_data_surf2_ptrtype;
// list of pt define in more in the surface
typedef std::vector<std::map<uint16_type,std::list<uint16_type> > > vec_map_data_ptsinsurf_type;
typedef boost::shared_ptr<vec_map_data_ptsinsurf_type> vec_map_data_ptsinsurf_ptrtype;

typedef std::map<int,std::list<int> > map_surfaceLoop_type;
//typedef boost::shared_ptr<map_surfaceLoop_type> map_surfaceLoop_ptrtype;


typedef boost::tuple< GeoGMSHTool_ptrtype,
        vec_map_data_ptrtype,
        std::string,
        std::string,
        vec_map_data_surf1_ptrtype,
        vec_map_data_surf2_ptrtype,
        vec_map_data_surf1_ptrtype,
        vec_map_data_ptsinsurf_ptrtype,
        map_surfaceLoop_type,
        double /*meshSize*/ > data_geo_type;
typedef boost::shared_ptr<data_geo_type> data_geo_ptrtype;


void run( data_geo_ptrtype __dg );



#define GEOTOOL_INSTANTIATES_FOR_COMP(r, state)                         \
        BOOST_PP_NOT_EQUAL( BOOST_PP_TUPLE_ELEM(2, 0, state),           \
                            BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, 1, state)) \
                            )                                           \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_INSTANTIATES_FOR_INCR(r, state)             \
        (                                                   \
         BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, 0, state)),	\
         BOOST_PP_TUPLE_ELEM(2, 1, state) )                 \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_INSTANTIATES_FOR(r,state)                               \
        void BOOST_PP_CAT(run,GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state))) (data_geo_ptrtype dg); \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
BOOST_PP_FOR( ( 0, BOOST_PP_SUB( BOOST_PP_ARRAY_SIZE( GEOTOOL_SHAPE ),1 ) ),
              GEOTOOL_INSTANTIATES_FOR_COMP,
              GEOTOOL_INSTANTIATES_FOR_INCR,
              GEOTOOL_INSTANTIATES_FOR )




} // namespace GeoTool

namespace GeoTool
{

/*_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*
 * GeoGMSHTool :                                   *
 *_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*/


class Node
{
public :

    Node()
        :
        M_node( new node_type() )
    {}

    Node( double __x ) :
        M_node( new node_type( 1 ) )
    {
        ( *M_node )( 0 )=__x;
    }

    Node( double __x, double __y ) :
        M_node( new node_type( 2 ) )
    {
        ( *M_node )( 0 )=__x;
        ( *M_node )( 1 )=__y;
    }

    Node( double __x, double __y, double __z ) :
        M_node( new node_type( 3 ) )
    {
        ( *M_node )( 0 )=__x;
        ( *M_node )( 1 )=__y;
        ( *M_node )( 2 )=__z;
    }

    Node( Node const & m )
        :
        M_node( m.M_node )
    {}

    Node operator=( Node const & m )
    {
        M_node.reset( new node_type( *( m.M_node ) ) );
        return *this;
    }

    double operator()( uint16_type n ) const
    {
        return this->getNode()( n );
    }

    double & operator()( uint16_type n )
    {
        return ( *M_node )( n );
    }

    node_type
    getNode() const
    {
        return *M_node;
    }

    node_type &
    getNode()
    {
        return *M_node;
    }

    boost::shared_ptr<node_type> M_node;
};

/*_________________________________________________*/

class Loop
{
public :

    Loop( Loop const & L ) : M_loop( L.M_loop ) {}

    Loop()
    {
        M_loop.clear();
    }

    void  operator=( Loop m )
    {
        this->M_loop=m.M_loop;
    }
    Loop  operator>>( int __n )
    {
        M_loop.push_back( __n );
        return *this;
    }

    uint16_type size()
    {
        return M_loop.size();
    }

    std::list<int>::const_iterator begin() const
    {
        return M_loop.begin();
    }
    std::list<int>::const_iterator end() const
    {
        return M_loop.end();
    }

    std::list<int> M_loop;
};

/*_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*
 * GeoGMSHTool data build :                        *
 *_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*/

namespace detail {

class GeoToolPoint
{
public :
    GeoToolPoint()
        :
        M_localId( invalid_size_type_value ),
        M_globalId( invalid_size_type_value ),
        M_node(0),
        M_hSize(0)
        {}
    GeoToolPoint(double x,double y,double z,double hsize, size_type localId,size_type globalId)
        :
        M_localId( localId ),
        M_globalId( globalId ),
        M_node(3),
        M_hSize( hsize )
        {
            M_node[0] = x;
            M_node[1] = y;
            M_node[2] = z;
        }
    GeoToolPoint( GeoToolPoint const& p )
        :
        M_localId( p.M_localId ),
        M_globalId( p.M_globalId ),
        M_node( p.M_node ),
        M_hSize( p.M_hSize ),
        M_physicalMarker( p.M_physicalMarker )
        {}

    size_type localId() const { return M_localId; }
    size_type globalId() const { return M_globalId; }
    node_type const& node() const { return M_node; }
    double hSize() const { return M_hSize; }
    std::string physicalMarker() const { return M_physicalMarker; }

    void setPhysicalMarker( std::string s ) { M_physicalMarker=s; }

    bool representSameNode( GeoToolPoint const& p ) const
    {
        return
            ( math::abs( this->node()[0]-p.node()[0] ) < 1e-9 ) &&
            ( math::abs( this->node()[1]-p.node()[1] ) < 1e-9 ) &&
            ( math::abs( this->node()[2]-p.node()[2] ) < 1e-9 );
    }
    void showMe() const
    {
        std::cout << "POINT -> "
                  << "localId  : " << this->localId() << " ; "
                  << "globalId : " << this->globalId() << " ; "
                  << "node     : " << M_node[0] << " " << M_node[1] << " " << M_node[2] << " ; "
                  << "hSize : " << M_hSize  << " ; "
                  << "physicalMarker : " << this->physicalMarker()
                  << "\n";
    }

private :
    size_type M_localId,M_globalId;
    node_type M_node;
    double M_hSize;
    std::string M_physicalMarker;

};

class GeoToolLine
{
public :
    GeoToolLine()
        :
        M_localId( invalid_size_type_value ),
        M_globalId( invalid_size_type_value )
        {}
    GeoToolLine(std::string name,size_type localId,size_type globalId,std::string type)
        :
        M_name( name ),
        M_localId( localId ),
        M_globalId( globalId ),
        M_lineType( type )
        {
            CHECK( type == "line" ||
                   type == "circle" || type == "ellipse" ||
                   type == "spline" || type == "bspline" ) << "invalid line type " << type;
        }
    GeoToolLine( GeoToolLine const& l )
        :
        M_name( l.M_name ),
        M_localId( l.M_localId ),
        M_globalId( l.M_globalId ),
        M_lineType( l.M_lineType ),
        M_listPt( l.M_listPt ),
        M_physicalMarker( l.M_physicalMarker )
        {}
    std::string name() const { return M_name; }
    size_type localId() const { return M_localId; }
    size_type globalId() const { return M_globalId; }
    std::string lineType() const { return M_lineType; }
    std::list<size_type> const& listPt() const { return M_listPt; }
    std::string physicalMarker() const { return M_physicalMarker; }
    size_type firstPointIdConnection() const
    {
        // warning must depend on type
        return M_listPt.front();
    }
    size_type secondPointIdConnection() const
    {
        // warning must depend on type
        return M_listPt.back();
    }

    void setPhysicalMarker( std::string s ) { M_physicalMarker=s; }

    //void setPoints( std::initializer_list<size_type> const& list )
    void setPoints( std::list<size_type> const& listPt )
    {
        M_listPt = listPt;
    }
    void replacePointId( int idIn, int idOut )
    {
        std::replace( M_listPt.begin(), M_listPt.end(), idIn, idOut );
    }

    void showMe() const
    {
        std::cout << "LINE -> "
                  << "name : " << this->name() << " ; "
                  << "localId  : " << this->localId() << " ; "
                  << "globalId : " << this->globalId() << " ; "
                  << "lineType : " << this->lineType() << " ; "
                  << "physicalMarker : " << this->physicalMarker() << " ; "
                  << "ptsId    : ";
        for ( size_type ptId : M_listPt )
            std::cout << " " << ptId;
        std::cout << "\n";

    }

private :
    std::string M_name;
    size_type M_localId,M_globalId;
    std::string M_lineType;
    std::list<size_type> M_listPt;
    std::string M_physicalMarker;
};

class GeoToolLineLoop
{
public :
    GeoToolLineLoop()
        :
        M_localId( invalid_size_type_value ),
        M_globalId( invalid_size_type_value )
        {}
    GeoToolLineLoop(std::string name,size_type localId,size_type globalId)
        :
        M_name( name ),
        M_localId( localId ),
        M_globalId( globalId )
        {}
    GeoToolLineLoop( GeoToolLineLoop const& l )
        :
        M_name( l.M_name ),
        M_localId( l.M_localId ),
        M_globalId( l.M_globalId ),
        M_listLine( l.M_listLine )
        {}
    std::string name() const { return M_name; }
    size_type localId() const { return M_localId; }
    size_type globalId() const { return M_globalId; }
    std::list<int> const& listLine() const { return M_listLine; }

    void replaceLineId( int idIn, int idOut )
    {
        std::replace( M_listLine.begin(), M_listLine.end(), idIn, idOut );
    }
    void setLines( std::list<int> const& listLine )
    {
        M_listLine = listLine;
    }
    void showMe() const
    {
        std::cout << "LINELOOP -> "
                  << "name : " << this->name() << " ; "
                  << "localId  : " << this->localId() << " ; "
                  << "globalId : " << this->globalId() << " ; "
                  << "linesId  : ";
        for ( int lineId : M_listLine )
            std::cout << " " << lineId;
        std::cout << "\n";
    }

private :
    std::string M_name;
    size_type M_localId,M_globalId;
    std::list<int> M_listLine;
};

class GeoToolSurface
{
public :

    GeoToolSurface()
        :
        M_localId( invalid_size_type_value ),
        M_globalId( invalid_size_type_value )
        {}
    GeoToolSurface(std::string name,size_type localId,size_type globalId,std::string type)
        :
        M_name( name ),
        M_localId( localId ),
        M_globalId( globalId ),
        M_surfaceType( type )
        {}
    GeoToolSurface( GeoToolSurface const& l )
        :
        M_name( l.M_name ),
        M_localId( l.M_localId ),
        M_globalId( l.M_globalId ),
        M_surfaceType( l.M_surfaceType ),
        M_lineLoopId( l.M_lineLoopId ),
        M_physicalMarker( l.M_physicalMarker ),
        M_ptsInSurface( l.M_ptsInSurface )
        {}
    std::string name() const { return M_name; }
    size_type localId() const { return M_localId; }
    size_type globalId() const { return M_globalId; }
    std::string surfaceType() const { return M_surfaceType; }
    std::list<int> const& listLineLoop() const { return M_lineLoopId; }
    std::string physicalMarker() const { return M_physicalMarker; }
    std::set<int> const& ptsInSurface() const { return M_ptsInSurface; }

    void setPhysicalMarker( std::string s ) { M_physicalMarker=s; }

    void setLineLoop( int lineLoopId )
    {
        M_lineLoopId = { lineLoopId };
    }
    void setLineLoop( std::list<int> const& lineLoopId )
    {
        M_lineLoopId = lineLoopId;
    }
    void addLineLoop( int lineLoopId )
    {
        M_lineLoopId.push_back( lineLoopId );
    }
    void addLineLoop( std::list<int> const& lineLoopIds )
    {
        for ( auto lineLoopId : lineLoopIds )
            M_lineLoopId.push_back( lineLoopId );
    }

    void setPtInSurface( std::set<int> const& ptIds )
    {
        M_ptsInSurface = ptIds;
    }
    void addPtInSurface( int ptId )
    {
        M_ptsInSurface.insert( ptId );
    }

    void showMe() const
    {
        std::cout << "SURFACE -> "
                  << "localId : " << this->localId() << " ; "
                  << "globalId : " << this->globalId() << " ; "
                  << "surfaceType" << this->surfaceType() << " ; "
                  << "physicalMarker : " << this->physicalMarker() << " ; "
                  << "lineLoopId : ";
        for ( auto llid : this->listLineLoop() )
            std::cout << llid << " ";
        std::cout << "\n";
    }


private :
    std::string M_name;
    size_type M_localId,M_globalId;
    std::string M_surfaceType;
    std::list<int> M_lineLoopId;
    std::string M_physicalMarker;
    std::set<int> M_ptsInSurface;

};

class GeoToolSurfaceLoop
{
public :

    GeoToolSurfaceLoop()
        :
        M_localId( invalid_size_type_value ),
        M_globalId( invalid_size_type_value )
        {}
    GeoToolSurfaceLoop(std::string name,size_type localId,size_type globalId)
        :
        M_name( name ),
        M_localId( localId ),
        M_globalId( globalId )
        {}
    GeoToolSurfaceLoop( GeoToolSurfaceLoop const& l )
        :
        M_name( l.M_name ),
        M_localId( l.M_localId ),
        M_globalId( l.M_globalId ),
        M_listSurface( l.M_listSurface )
        {}

    std::string name() const { return M_name; }
    size_type localId() const { return M_localId; }
    size_type globalId() const { return M_globalId; }
    std::list<int> const& listSurface() const { return M_listSurface; }

    void replaceSurfaceId( int idIn, int idOut )
    {
        std::replace( M_listSurface.begin(), M_listSurface.end(), idIn, idOut );
    }
    void setSurfaces( std::list<int> const& listSurface )
    {
        M_listSurface = listSurface;
    }
    void showMe() const
    {
        std::cout << "SURFACELOOP -> "
                  << "localId  : " << this->localId() << " ; "
                  << "globalId : " << this->globalId() << " ; "
                  << "surfacesId  : ";
        for ( int surfId : M_listSurface )
            std::cout << " " << surfId;
        std::cout << "\n";
    }

private :
    std::string M_name;
    size_type M_localId,M_globalId;
    std::list<int> M_listSurface;

};


class GeoToolVolume
{
public :

    GeoToolVolume()
        :
        M_localId( invalid_size_type_value ),
        M_globalId( invalid_size_type_value )
        {}
    GeoToolVolume(size_type localId,size_type globalId)
        :
        M_localId( localId ),
        M_globalId( globalId )
        {}
    GeoToolVolume( GeoToolVolume const& l )
        :
        M_localId( l.M_localId ),
        M_globalId( l.M_globalId ),
        M_surfaceLoopId( l.M_surfaceLoopId ),
        M_physicalMarker( l.M_physicalMarker )
        {}
    //std::string name() const { return M_name; }
    size_type localId() const { return M_localId; }
    size_type globalId() const { return M_globalId; }
    std::list<int> const& listSurfaceLoop() const { return M_surfaceLoopId; }
    std::string physicalMarker() const { return M_physicalMarker; }

    void setPhysicalMarker( std::string s ) { M_physicalMarker=s; }

    void setSurfaceLoop( int surfLoopId )
    {
        M_surfaceLoopId = { surfLoopId };
    }
    void setSurfaceLoop( std::list<int> const& surfLoopIds )
    {
        M_surfaceLoopId = surfLoopIds;
    }
    void addSurfaceLoop( int surfLoopId )
    {
        M_surfaceLoopId.push_back( surfLoopId );
    }

    void showMe() const
    {
        std::cout << "VOLUME -> "
                  << "localId : " << this->localId() << " ; "
                  << "globalId : " << this->globalId() << " ; "
                  << "physicalMarker : " << this->physicalMarker() << " ; "
                  << "surfaceLoopId : ";
        for ( auto llid : this->listSurfaceLoop() )
            std::cout << llid << " ";
        std::cout << "\n";
    }

private :
    //std::string M_name;
    size_type M_localId,M_globalId;
    std::string M_surfaceType;
    std::list<int> M_surfaceLoopId;
    std::string M_physicalMarker;
};



class ApplyFusionMarkersLine
{
public :

    ApplyFusionMarkersLine()
        :
        M_globalId1Line( 0 ),
        M_globalId2Line( 0 ),
        M_globalId1LineLoop( 0 ),
        M_globalId2LineLoop( 0 )
        {}

    ApplyFusionMarkersLine( ApplyFusionMarkersLine const& a )
        :
        M_globalId1Line( a.M_globalId1Line ),
        M_globalId2Line( a.M_globalId2Line ),
        M_globalId1LineLoop( a.M_globalId1LineLoop ),
        M_globalId2LineLoop( a.M_globalId2LineLoop )
        {}

    void setGlobalId1( int globalId1Line, int globalId1LineLoop )
    {
        M_globalId1Line = globalId1Line;
        M_globalId1LineLoop = globalId1LineLoop;
    }
    void setGlobalId2( int globalId2Line, int globalId2LineLoop )
    {
        M_globalId2Line = globalId2Line;
        M_globalId2LineLoop = globalId2LineLoop;
    }

    int globalId1Line() const { return M_globalId1Line; }
    int globalId2Line() const { return M_globalId2Line; }
    int globalId1LineLoop() const { return M_globalId1LineLoop; }
    int globalId2LineLoop() const { return M_globalId2LineLoop; }

private :

    int M_globalId1Line, M_globalId2Line;
    int M_globalId1LineLoop, M_globalId2LineLoop;
};

} // namespace detail

//class GeoGMSHTool;

class FusionMarkers
{
public :
    FusionMarkers(GeoGMSHTool const& gt1, int marker1,GeoGMSHTool const& gt2, int marker2);

    FusionMarkers( FusionMarkers const& f )
        :
        M_nameGT1( f.M_nameGT1 ), M_nameGT2( f.M_nameGT2 ),
        M_marker1( f.M_marker1 ),M_marker2( f.M_marker2 )
        {}

    std::string nameGT1() const { return M_nameGT1; }
    std::string nameGT2() const { return M_nameGT2; }
    int marker1() const { return M_marker1; }
    int marker2() const { return M_marker2; }

private :
    std::string M_nameGT1, M_nameGT2;
    int M_marker1, M_marker2;
};

/*_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*
 * GeoGMSHTool :                                   *
 *_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*/


class GeoGMSHTool
{
public:

    typedef node<double>::type node_type;

    /*            // list de < nameMesh, meshSize >
    typedef boost::tuple<std::string,double> names_base_type;
    typedef std::list< names_base_type > names_type;
    typedef std::map< std::string, names_type > map_shape_names_type;
    typedef names_type::const_iterator names_const_iterator_type;
    typedef map_shape_names_type::const_iterator map_shape_names_const_iterator_type;
    */
    typedef boost::tuple<std::string,std::string,uint16_type> marker_base_type;
    //typedef std::map<std::string,std::list<marker_base_type > > marker_markerName_type;
    typedef std::map<std::string,std::vector<marker_base_type > > marker_markerName_type;
    typedef std::map< std::string, marker_markerName_type > marker_type_type;
    typedef std::map< std::string, marker_type_type > marker_name_type;
    typedef std::map< std::string, marker_type_type > marker_shape_type;

    typedef marker_markerName_type::const_iterator marker_markerName_const_iterator_type;
    typedef marker_type_type::const_iterator marker_type_const_iterator_type;
    typedef marker_name_type::const_iterator marker_name_const_iterator_type;
    typedef marker_shape_type::const_iterator marker_shape_const_iterator_type;

    typedef std::vector<node_type> parameter_rectangle_type;
    typedef std::map<std::string, parameter_rectangle_type > parameter_name_type;
    typedef std::map<std::string, parameter_name_type > parameter_shape_type;
    typedef parameter_name_type::const_iterator parameter_name_const_iterator_type;
    typedef parameter_shape_type::const_iterator parameter_shape_const_iterator_type;


    // gestion des lignes : shape,name,value,meshSize
    typedef boost::tuple<std::string,std::string,uint16_type,double > ligne_type;
    typedef std::list< ligne_type > ligne_type_type;
    typedef std::list< ligne_type_type > ligne_name_type;
    typedef ligne_type_type::const_iterator ligne_type_const_iterator_type;
    typedef ligne_name_type::const_iterator ligne_name_const_iterator_type;

    // gestion des surfaces : shape,name,(numGlobSurface,valueOfLineloop),meshSize
    typedef boost::tuple<std::string,std::string,std::pair<int,int>,double > surface_type;
    typedef std::list< surface_type > surface_type_type;
    typedef std::list< surface_type_type > surface_name_type;
    typedef surface_type_type::const_iterator surface_type_const_iterator_type;
    typedef surface_name_type::const_iterator surface_name_const_iterator_type;

    // gestion des volumes : shape,name,(numGlobVolume,value),meshSize
    typedef boost::tuple<std::string,std::string, std::pair<int,int>,double > volume_type;
    typedef std::list< volume_type > volume_type_type;
    typedef std::list< volume_type_type > volume_name_type;
    typedef volume_type_type::const_iterator volume_type_const_iterator_type;
    typedef volume_name_type::const_iterator volume_name_const_iterator_type;

    // gestion des surfaceLoop : shape,name, numLoopLoc->list<value>
    typedef boost::tuple<std::string,std::string, std::map< int, std::list<int> > > surfaceloop_type;
    typedef std::list< surfaceloop_type > surfaceloop_type_type;
    typedef std::list< surfaceloop_type_type > surfaceloop_name_type;
    typedef surfaceloop_type_type::const_iterator surfaceloop_type_const_iterator_type;
    typedef surfaceloop_name_type::const_iterator surfaceloop_name_const_iterator_type;


    GeoGMSHTool( uint16_type __dim, std::string __shape="NO_SHAPE", std::string __name="NO_NAME", double __meshSize=0.1 );

    GeoGMSHTool( uint16_type __dim,  std::string const & geoUserStr, double __meshSize=0.1, std::string __shape="NO_SHAPE", std::string __name="NO_NAME" );

    GeoGMSHTool( GeoGMSHTool const & m );

    virtual ~GeoGMSHTool() {}

    void zeroCpt();

    void operator=( GeoGMSHTool const & m );

    GeoGMSHTool operator+( const GeoGMSHTool & m );
    GeoGMSHTool operator-( const GeoGMSHTool & m );

    GeoGMSHTool opFusion( const GeoGMSHTool & m,int __typeop );

    void init( int orderGeo,
               std::string gmshFormatVersion,
               double hmin=0,double hmax=1e22,
               int refine=0,
               bool optimize3dNetgen=true,
               GMSH_PARTITIONER partitioner=GMSH_PARTITIONER_CHACO,
               int partitions=1,
               bool partition_file=false );

    /*
     *
     */
    void initData( std::string __shape,
                   std::string __name,
                   double __meshSize,
                   std::vector<GeoTool::Node> & __param,
                   uint16_type dim,
                   uint16_type __nbligne,
                   uint16_type __nbsurface,
                   uint16_type __nbvolume );

    /*
     *Utile pour la fct geoStr()
     *Pas de maj pour cptSurface et cptVolume car traitement different
     */
    void updateData( GeoGMSHTool const & m );

    /*
     * Update the output stringstream wich generate the gmsh code
     */
    void updateOstr( std::string __str )
    {
        *M_ostr << __str;
    }

    /*
     * Generate the gmsh code
     */
    void geoStr();

    /*
     * Clean
     */
    void cleanOstr()
    {
        M_ostr.reset( new std::ostringstream() );
    }

    /**
     * Print information
     */
    void showMe() const;


    BOOST_PARAMETER_MEMBER_FUNCTION(
        (void),
        setMarker,
        tag,
        (required
         ( type, (std::string))
         ( name, (std::string)) )
        (optional
         (markerAll, (bool), false)
         (marker1, (bool), false)
         (marker2, (bool), false)
         (marker3, (bool), false)
         (marker4, (bool), false)
         (marker5, (bool), false)
         (marker6, (bool), false)
         (marker7, (bool), false)
         (marker8, (bool), false)
         (marker9, (bool), false)
         (marker10, (bool), false)
         (marker11, (bool), false)
         (marker12, (bool), false)
         ))
        {
            std::vector<bool> mymarkers(12,markerAll);
            if (!markerAll) {
                mymarkers[0]=marker1;
                mymarkers[1]=marker2;
                mymarkers[2]=marker3;
                mymarkers[3]=marker4;
                mymarkers[4]=marker5;
                mymarkers[5]=marker6;
                mymarkers[6]=marker7;
                mymarkers[7]=marker8;
                mymarkers[8]=marker9;
                mymarkers[9]=marker10;
                mymarkers[10]=marker11;
                mymarkers[11]=marker12;
            }
            this->setMarkerImpl(type,name,mymarkers);
        }
    virtual void setMarkerImpl( std::string type, std::string name, std::vector<bool> const& markers )
    {
        CHECK( false ) << "not implemented in GeoGMSHTool class, but in shape class";
    }


    BOOST_PARAMETER_MEMBER_FUNCTION(
        ( typename Feel::detail::mesh<Args>::ptrtype ), // return type
        createMesh, // function name
        tag,
        ( required
          ( mesh, * )
          ( name, ( std::string ) )
        ) //required
        ( optional
          ( format,         *, option(_name="gmsh.format").template as<int>() )
          ( straighten,     *( boost::is_integral<mpl::_> ), 1 )
          ( refine,          *( boost::is_integral<mpl::_> ), 0 )
          ( partitions,   *( boost::is_integral<mpl::_> ), Environment::worldComm().size() )
          ( partition_file,   *( boost::is_integral<mpl::_> ), 0 )
          ( partitioner,   *( boost::is_integral<mpl::_> ), GMSH_PARTITIONER_CHACO )
          ( worldcomm,      *, Environment::worldComm() )
          ( hmin,     ( double ), 0 )
          ( hmax,     ( double ), 1e22 )
          ( optimize3d_netgen, *( boost::is_integral<mpl::_> ), true )
        ) //optional
    )
    {
        typedef typename Feel::detail::mesh<Args>::type _mesh_type;
        typedef typename Feel::detail::mesh<Args>::ptrtype _mesh_ptrtype;

        _mesh_ptrtype _mesh( mesh );
        _mesh->setWorldComm( worldcomm );

        if ( worldcomm.isActive() )
        {

            this->cleanOstr();
            this->zeroCpt();
            Gmsh gmsh( _mesh_type::nDim, _mesh_type::nOrder, worldcomm );
            gmsh.setRecombine( _mesh_type::shape_type::is_hypercube );
            gmsh.setRefinementLevels( refine );
            gmsh.setFileFormat( (GMSH_FORMAT)format );
            gmsh.setNumberOfPartitions( partitions );
            gmsh.setPartitioner( partitioner );
            gmsh.setMshFileByPartition( partition_file );
            this->init( _mesh_type::nOrder,gmsh.version(),
                        hmin,hmax,refine,
                        optimize3d_netgen,
                        partitioner,partitions,partition_file );

            std::string geostring;

            if ( M_geoIsDefineByUser )
            {
                geostring= M_ostrDefineByUser->str();
            }

            else
            {
                this->geoStr();
                geostring = M_ostr->str();
            }


            std::string fname;
            bool gen;
            boost::tie( fname, gen ) = gmsh.generate( name,
                                                      geostring,
                                                      false,false,false );

            ImporterGmsh<_mesh_type> import( fname, FEELPP_GMSH_FORMAT_VERSION, worldcomm );
            _mesh->accept( import );
            _mesh->components().set ( MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
            _mesh->updateForUse();

            if ( straighten && _mesh_type::nOrder > 1 )
                return straightenMesh( _mesh, worldcomm.subWorldComm(), false, false );

        } // if (worldcomm.isActive())

        return _mesh;
    }


#if 0 // old function
    template<typename mesh_type>
    boost::shared_ptr<mesh_type>
    createMesh( std::string name, int straighten = 1, WorldComm const& worldcomm=Environment::worldComm() )
    {
        boost::shared_ptr<mesh_type> mesh( new mesh_type );
        mesh->setWorldComm( worldcomm );

        if ( worldcomm.isActive() )
        {
            this->cleanOstr();
            this->zeroCpt();

            Gmsh gmsh( mesh_type::nDim,mesh_type::nOrder, worldcomm );
            gmsh.setOrder( mesh_type::nOrder );
            gmsh.setRecombine( mesh_type::shape_type::is_hypercube );

            this->init( mesh_type::nOrder,gmsh.version() );

            std::string geostring;

            if ( M_geoIsDefineByUser )
            {
                geostring= M_ostrDefineByUser->str();
            }

            else
            {
                this->geoStr();
                geostring = M_ostr->str();
            }

            std::string fname = gmsh.generate( name,
                                               geostring,false,false,false );

            ImporterGmsh<mesh_type> import( fname, FEELPP_GMSH_FORMAT_VERSION, worldcomm );
            mesh->accept( import );
            mesh->components().set ( MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
            mesh->updateForUse();

            if ( straighten && mesh_type::nOrder > 1 )
                return straightenMesh( mesh );
        } // if (worldcomm.isActive())

        return mesh;
    }
#endif
    /*_________________________________________________*
     *_________________________________________________*
     * Accessor                                        *
     *_________________________________________________*
     *_________________________________________________*/

    uint16_type dim() const
    {
        return  M_dim;
    }
    std::string name() const
    {
        return M_name;
    }

    uint16_type cptPt() const
    {

        return M_cptPt;
    }
    uint16_type cptLine() const
    {
        return M_cptLine;
    }
    uint16_type cptLineLoop() const
    {
        return M_cptLineLoop;
    }
    uint16_type cptSurface() const
    {
        return M_cptSurface;
    }
    uint16_type cptTableau() const
    {
        return M_cptTableau;   //voir les extrudes par exemple
    }
    uint16_type cptSurfaceLoop() const
    {
        return M_cptSurfaceLoop;
    }
    uint16_type cptVolume() const
    {
        return M_cptVolume;
    }

    /*_________________________________________________*
     * Parameter
     *_________________________________________________*/

    parameter_shape_const_iterator_type
    paramShapeBegin() const
    {
        return M_paramShape->begin();
    }

    parameter_shape_const_iterator_type paramShapeEnd() const
    {
        return M_paramShape->end();
    }

    parameter_name_const_iterator_type
    paramNameBegin( std::string __shape ) const
    {
        return M_paramShape->find( __shape )->second.begin();
    }

    parameter_name_const_iterator_type
    paramNameEnd( std::string __shape ) const
    {
        return M_paramShape->find( __shape )->second.end();
    }

    parameter_rectangle_type
    getParameter( std::string __shape, std::string __name ) const
    {
        return M_paramShape->find( __shape )->second.find( __name )->second;
    }

    /*_________________________________________________*
     * Marker
     *_________________________________________________*/
    /*
    marker_shape_const_iterator_type
    markShapeBegin() const
    {
        return M_markShape->begin();
    }

    marker_shape_const_iterator_type
    markShapeEnd() const
    {
        return M_markShape->end();
    }*/


    marker_type_const_iterator_type
    markerTypeBegin( /*std::string __shape*/ ) const
    {
        //return M_markShape->find(__shape)->second.begin();
        return M_markShape->begin();
    }

    marker_type_const_iterator_type
    markerTypeEnd( /*std::string __shape*/ ) const
    {
        //return M_markShape->find(__shape)->second.end();
        return M_markShape->end();
    }
    /*
    marker_type_type
    markerType(std::string __shape) const
    {
        //return M_markShape->find(__shape)->second;
        }*/

    marker_markerName_const_iterator_type
    markerMarkerNameBegin( /*std::string __shape,*/ std::string __type ) const
    {
        //return M_markShape->find(__shape)->second.find(__type)->second.begin();
        return M_markShape->find( __type )->second.begin();
    }

    marker_markerName_const_iterator_type
    markerMarkerNameEnd( /*std::string __shape,*/ std::string __type ) const
    {
        //return M_markShape->find(__shape)->second.find(__type)->second.end();
        return M_markShape->find( __type )->second.end();
    }

    marker_markerName_type
    markerMarkerName( /*std::string __shape,*/ std::string __type ) const
    {
        //return M_markShape->find(__shape)->second.find(__type)->second;
        return M_markShape->find( __type )->second;
    }

    //std::list<marker_base_type>::const_iterator
    std::vector<marker_base_type>::const_iterator
    markerListIndiceBegin( /*std::string __shape,*/ std::string __type ,std::string __markerName ) const
    {
        return M_markShape->find( __type )->second.find( __markerName )->second.begin();
    }

    //std::list<marker_base_type>::const_iterator
    std::vector<marker_base_type>::const_iterator
    markerListIndiceEnd( /*std::string __shape,*/ std::string __type ,std::string __markerName ) const
    {
        //return M_markShape->find(__shape)->second.find(__type)->second.find(__markerName)->second.end();
        return M_markShape->find( __type )->second.find( __markerName )->second.end();
    }


    //std::list<marker_base_type>
    std::vector<marker_base_type>
    getMarkerName( /*std::string __shape,*/ std::string __type ,std::string __markerName ) const
    {
        //return M_markShape->find(__shape)->second.find(__type)->second.find(__markerName)->second;
        return M_markShape->find( __type )->second.find( __markerName )->second;
    }

    /*_________________________________________________*
     *_________________________________________________*
     * Members                                         *
     *_________________________________________________*
     *_________________________________________________*/

    uint16_type M_dim;
    std::string M_name;
    // memory
    uint16_type M_cptPt;
    uint16_type M_cptLine;
    uint16_type M_cptLineLoop;
    uint16_type M_cptSurface;
    uint16_type M_cptTableau;
    uint16_type M_cptSurfaceLoop;
    uint16_type M_cptVolume;

    // gestion des surface : shape,name,value
    // value is the marker associated to the planeSurface (init to 0 and to use when call geoStr())
    //std::list< std::list< boost::tuple<std::string,std::string, uint16_type > > > M_surfaceList;
    boost::shared_ptr<ligne_name_type> M_ligneList;
    boost::shared_ptr<surface_name_type> M_surfaceList;
    boost::shared_ptr<volume_name_type> M_volumeList;
    boost::shared_ptr<surfaceloop_name_type> M_surfaceLoopList;

    boost::shared_ptr<std::ostringstream> M_ostrExtrude;
    boost::shared_ptr<std::ostringstream> M_ostrSurfaceLoop;


    // data containers
    //boost::shared_ptr<map_shape_names_type> M_map_Shape;
    boost::shared_ptr<parameter_shape_type> M_paramShape;
    //boost::shared_ptr<marker_shape_type> M_markShape;
    boost::shared_ptr<marker_type_type> M_markShape;

    // output string
    boost::shared_ptr<std::ostringstream> M_ostr;

    boost::shared_ptr<std::ostringstream> M_ostrDefineByUser;
    bool M_geoIsDefineByUser;


    // usefull for writing geofile
    void addPoint(/*std::string name, */detail::GeoToolPoint const& pt)
    {
        //M_buildDataPoint[name].push_back(pt);
        M_buildDataPoint[pt.globalId()] =  pt;
    }
    void addLine(detail::GeoToolLine const& line)
    {
        M_buildDataLine[line.globalId()] = line;
    }
    void addLineLoop(detail::GeoToolLineLoop const& lineLoop)
    {
        M_buildDataLineLoop[lineLoop.globalId()] = lineLoop;
    }
    void addSurface(detail::GeoToolSurface const& surf)
    {
        M_buildDataSurface[surf.globalId()] = surf;
    }
    void addSurfaceLoop(detail::GeoToolSurfaceLoop const& surfLoop)
    {
        M_buildDataSurfaceLoop[surfLoop.globalId()] = surfLoop;
    }
    void addVolume(detail::GeoToolVolume const& vol)
    {
        M_buildDataVolume[vol.globalId()] = vol;
    }
    std::map<size_type, detail::GeoToolPoint> M_buildDataPoint;
    std::map<size_type, detail::GeoToolLine> M_buildDataLine;
    std::map<size_type, detail::GeoToolLineLoop> M_buildDataLineLoop;
    std::map<size_type, detail::GeoToolSurface> M_buildDataSurface;
    std::map<size_type, detail::GeoToolSurfaceLoop> M_buildDataSurfaceLoop;
    std::map<size_type, detail::GeoToolVolume> M_buildDataVolume;


    /*int
    lineIdFromNameAndLocalId( std::string name, int localId ) const
    {
        for ( auto const& line : M_buildDataLine )
            if ( line.second.name() == name && line.second.localId() == localId )
                return line.first;
        return 0;
     }*/

    int
    surfaceIdFromName( std::string name ) const
    {
        for ( auto const& surf : M_buildDataSurface )
            if ( surf.second.name() == name )
                return surf.first;
        return 0;
    }

    GeoGMSHTool&
    fusion( GeoGMSHTool const& gt1, int marker1, GeoGMSHTool const& gt2, int marker2, bool keepInterface=true)
    {
        if ( keepInterface )
            M_fusionMarkersLineWithInterface.push_back( FusionMarkers( gt1,marker1,gt2,marker2 ) );
        else
            M_fusionMarkersLineWithoutInterface.push_back( FusionMarkers( gt1,marker1,gt2,marker2 ) );
        return *this;
    }


    std::vector<FusionMarkers> const& fusionMarkersLineWithInterface() { return M_fusionMarkersLineWithInterface; }
    std::vector<FusionMarkers> const& fusionMarkersLineWithoutInterface() { return M_fusionMarkersLineWithoutInterface; }
    std::vector<FusionMarkers> M_fusionMarkersLineWithInterface, M_fusionMarkersLineWithoutInterface;

    std::vector<detail::ApplyFusionMarkersLine> M_buildApplyFusionMarker;

    bool hasSameOrientation( detail::GeoToolLine const& l1, detail::GeoToolLine const& l2 ) const
    {
        CHECK( l1.lineType() == l2.lineType() ) << "error not same type";
        if ( l1.lineType() == "line" )
        {
            CHECK( l1.listPt().size() == 2 ) << "error";
            auto thenode1p1 = M_buildDataPoint.find(l1.listPt().front())->second.node();
            auto thenode2p1 = M_buildDataPoint.find(l2.listPt().front())->second.node();

            if ( ( std::abs(thenode1p1[0] - thenode2p1[0]) < 1e-9 ) &&
                 ( std::abs(thenode1p1[1] - thenode2p1[1]) < 1e-9 ) &&
                 ( std::abs(thenode1p1[2] - thenode2p1[2]) < 1e-9 ) )
                return true;
            else
                return false;
        }
        else CHECK( false ) << "TODO";

    }


    bool lineLoopIsClosed( detail::GeoToolLineLoop const& lineloop ) const;
    bool lineLoopHasConnection( detail::GeoToolLineLoop const& lineloop1, detail::GeoToolLineLoop const& lineloop2 ) const;
    void lineLoopApplyConnection( detail::GeoToolLineLoop & lineloop1, detail::GeoToolLineLoop const& lineloop2 );
    std::set<int> lineLoopPointIdsNotConnected( detail::GeoToolLineLoop const& lineloop ) const;

    // name->(list of local id )
    //std::map<std::string,std::set<size_type> > M_fusionMarkerDoublonPoint,M_fusionMarkerDoublonLine,M_fusionMarkerDoublonLineLoop;


    std::pair<bool,std::string>
    findPhysicalMarker(std::string markerType, std::string __name, int __numLoc )
    {
        if ( this->M_markShape->find( markerType/*"line"*/) != this->M_markShape->end() )
        {
            for ( auto markName : this->markerMarkerName(markerType/*"line"*/) )
            {
                for ( auto lineMarked : markName.second )
                {
                    //std::string __shape = boost::get<2>( *__dg );
                    //std::string __name = boost::get<3>( *__dg );
                    if ( __name == lineMarked.get<1>() && __numLoc == lineMarked.get<2>() )
                    {
                        //myLine.setPhysicalMarker( markName.first );
                        return std::make_pair( true,markName.first );
                    }
                }
            }
        }
        return std::make_pair( false,std::string("") );
    }


};

/*_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*
 * Function on the namespace                       *
 *_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*/


void run( data_geo_ptrtype __dg );

template <uint16_type Numero>
node_type
param( data_geo_ptrtype __dg );


void
writePoint( uint16_type __numLoc, data_geo_ptrtype __dg ,double __x1,double __x2=0, double __x3=0 );

void
writeLine( uint16_type __numLoc, data_geo_ptrtype __dg ,uint16_type __n1, uint16_type __n2 );

void
writeCircle( uint16_type __numLoc, data_geo_ptrtype __dg ,uint16_type __n1, uint16_type __n2, uint16_type __n3 );

void
writeEllipse( uint16_type __numLoc, data_geo_ptrtype __dg ,uint16_type __n1, uint16_type __n2, uint16_type __n3, uint16_type __n4 );

void
writeSpline( uint16_type __numLoc, data_geo_ptrtype __dg ,Loop __loop );

void
writeBSpline( uint16_type __numLoc, data_geo_ptrtype __dg ,Loop __loop );

void
writeLineLoop( uint16_type __numLoc, data_geo_ptrtype __dg , Loop /*const*/ __loop );

void
writePlaneSurface( uint16_type __numLoc, data_geo_ptrtype __dg , uint16_type __ind );

void
writeRuledSurface( uint16_type __numLoc, data_geo_ptrtype __dg , uint16_type __ind );

void
writeExtrudeSurface( uint16_type __numLoc,data_geo_ptrtype __dg , uint16_type __ind,Loop /*const*/ __loop );

void
writePtInSurface( data_geo_ptrtype __dg , uint16_type __indPt,uint16_type __indSurf );

void
writeSurfaceLoop( uint16_type __numLoc, data_geo_ptrtype __dg , Loop /*const*/ __loop );

void
writeVolume( uint16_type __numLoc, data_geo_ptrtype __dg , uint16_type __ind );

boost::tuple<Node,Node,Node>
computeBasisOrthogonal( node_type dir,node_type centre );


/*_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*
 * PREPROCESSOR METHODS                            *
 *_________________________________________________*
 *_________________________________________________*
 *_________________________________________________*/




/*                                                 */
/**/
#define GEOTOOL_FOR_COMP(r, state)                                      \
        BOOST_PP_LESS( BOOST_PP_TUPLE_ELEM(2, 0, state),           \
                            BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, 1, state)) \
                            )                                           \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_FOR_INCR(r, state)                          \
        (                                                   \
         BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, 0, state)),	\
         BOOST_PP_TUPLE_ELEM(2, 1, state) )                 \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_SHAPE_FOR_PARAM_SIGNATURE(r, state)                     \
    Node BOOST_PP_CAT( __param, BOOST_PP_TUPLE_ELEM(2,0,state) ) = Node(0,0,0) BOOST_PP_COMMA() \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_SHAPE_PARAM_SIGNATURE(state)                           \
        BOOST_PP_FOR( (0, BOOST_PP_SUB(GEOTOOL_SHAPE_NBPARAM(BOOST_PP_TUPLE_ELEM(2,0,state)),1) ), \
                      GEOTOOL_FOR_COMP,                                 \
                      GEOTOOL_FOR_INCR,                                 \
                      GEOTOOL_SHAPE_FOR_PARAM_SIGNATURE)                \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_SHAPE_CLASS(r,state)                                    \
        class GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state)) : public GeoGMSHTool \
        {                                                               \
        public :                                                        \
                                                                        \
            typedef GeoGMSHTool::node_type node_type;                   \
            typedef GeoTool::Node Node;                                 \
                                                                        \
                                                                        \
            GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state))(double __meshSize, \
                                                                     std::string __name, \
                                                                     GEOTOOL_SHAPE_PARAM_SIGNATURE(state) \
                                                                     uint16_type type = 0 ); /*Ne sert a rien, juste a cause de la virgule au dessus)*/ \
                                                                        \
                                                                        \
                                                                        \
            GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state))(const GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state)) & m) \
                :                                                       \
                GeoGMSHTool(m),                                         \
                M_param(m.M_param)                                      \
                {}                                                      \
                                                                        \
                                                                        \
            void setMarkerImpl( std::string type, std::string name, std::vector<bool> const& markers ); \
                                                                        \
                                                                        \
            /*std::string M_name;*/                                     \
                                                                        \
            static const std::string shape() { return GEOTOOL_SHAPE_NAME_STR(BOOST_PP_TUPLE_ELEM(2,0, state));} \
            /*const std::string name() const {return M_name;}*/         \
                                                                        \
                                                                        \
            std::vector<GeoTool::Node> M_param;                        \
                                                                        \
        };                                                              \
        /**/
/*_________________________________________________*/
/*_________________________________________________*/
/*                                                 */
/**/



//creation des classes representants les objets geotool
BOOST_PP_FOR( ( 0, BOOST_PP_SUB( BOOST_PP_ARRAY_SIZE( GEOTOOL_SHAPE ),1 ) ),
              GEOTOOL_FOR_COMP,
              GEOTOOL_FOR_INCR,
              GEOTOOL_SHAPE_CLASS )



template<typename mesh_type>
boost::shared_ptr<mesh_type>
createMeshFromGeoFile( std::string geofile,std::string name,double meshSize,int straighten = 1,
                       int partitions=1, WorldComm worldcomm=Environment::worldComm(),
                       int partition_file = 0, GMSH_PARTITIONER partitioner = GMSH_PARTITIONER_CHACO )
{

    boost::shared_ptr<mesh_type> mesh( new mesh_type );
    mesh->setWorldComm( worldcomm );

    if ( !worldcomm.isActive() ) return mesh;

    Gmsh gmsh( mesh_type::nDim,mesh_type::nOrder,worldcomm );
    gmsh.setCharacteristicLength( meshSize );
    gmsh.setNumberOfPartitions( partitions );
    gmsh.setPartitioner( partitioner );
    gmsh.setMshFileByPartition( partition_file );
    gmsh.setRecombine( mesh_type::shape_type::is_hypercube );

    std::ostringstream ostr;

    // preambule :
    ostr << "Mesh.MshFileVersion = " << gmsh.version() << ";\n"
         << "Mesh.CharacteristicLengthExtendFromBoundary=1;\n"
         << "Mesh.CharacteristicLengthFromPoints=1;\n"
         << "Mesh.ElementOrder=" << gmsh.order() << ";\n"
         << "Mesh.SecondOrderIncomplete = 0;\n"
         << "Mesh.Algorithm = 6;\n" // 2D mesh algorithm (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad)
         << "Mesh.Algorithm3D = 4;\n" // 3D mesh algorithm (1=Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D)
         << "Mesh.OptimizeNetgen=1;\n"
         << "// partitioning data\n"
         << "Mesh.Partitioner=" << partitioner << ";\n"
         << "Mesh.NbPartitions=" << partitions << ";\n"
         << "Mesh.MshFilePartitioned=" << partition_file << ";\n";


    std::string contenu;
    std::ifstream ifstr( geofile.c_str(), std::ios::in );

    if ( ifstr )
    {
        // each line of the stream is appended in contenu
        while ( getline( ifstr, contenu ) )
            ostr << contenu<<"\n";

        ifstr.close();
    }


    std::string fname;
    bool generated;
    boost::tie( fname, generated ) = gmsh.generate( name,
                                                    ostr.str(),false,false,true );

    ImporterGmsh<mesh_type> import( fname, FEELPP_GMSH_FORMAT_VERSION, worldcomm );

    mesh->accept( import );
    mesh->components().set ( MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
    mesh->updateForUse();

    if ( straighten && mesh_type::nOrder > 1 )
        return straightenMesh( mesh, worldcomm.subWorldComm() );

    return mesh;
}



}//GeoTool

} //Feel
#endif /* FEELPP_GEOTOOL_HPP */

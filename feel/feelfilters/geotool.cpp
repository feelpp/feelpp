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
   \file geotool.cpp
   \author Vincent Chabannes <vincent.chabannes@imag.fr>
   \date 2011-03-03
 */



#include <feel/feelfilters/geotool.hpp>
#include <feel/feelfilters/geotoolshape.cpp>

#if defined( FEELPP_HAS_GMSH_H )
#include <GmshConfig.h>
#endif

namespace Feel
{

namespace GeoTool
{

FusionMarkers::FusionMarkers(GeoGMSHTool const& gt1, int marker1,GeoGMSHTool const& gt2, int marker2)
    :
    M_nameGT1( gt1.name() ), M_nameGT2( gt2.name() ),
    M_marker1( marker1 ),M_marker2( marker2 )
{}


bool
GeoGMSHTool::lineLoopIsClosed( detail::GeoToolLineLoop const& lineloop ) const
{
    if ( lineloop.listLine().size() == 0 )
        return false;

    // init with first line in lineloop
    int lineIdStart = lineloop.listLine().front();
    int lineIdStartPos = (lineIdStart>0)? lineIdStart: -lineIdStart;
    auto const& thelineStart = M_buildDataLine.find(lineIdStartPos)->second;
    int pointIdStart = (lineIdStart>0)? thelineStart.firstPointIdConnection() : thelineStart.secondPointIdConnection();
    int currentPointId = (lineIdStart>0)? thelineStart.secondPointIdConnection() : thelineStart.firstPointIdConnection();

    for (  int thelineId : lineloop.listLine() )
    {
        // ignore the first line
        if ( thelineId == lineIdStart ) continue;

        int thelineIdPos = ( thelineId>0 )? thelineId : -thelineId;
        auto const& theline = M_buildDataLine.find(thelineIdPos)->second;
        if ( thelineId > 0 )
        {
            if ( currentPointId != theline.firstPointIdConnection() )
                return false;
            currentPointId = theline.secondPointIdConnection();
        }
        else
        {
            if ( currentPointId != theline.secondPointIdConnection() )
                return false;
            currentPointId = theline.firstPointIdConnection();
        }
    }

    // check last point is the first point
    if ( currentPointId != pointIdStart )
        return false;

    // the lineloop is well closed
    return true;
}

std::set<int>
GeoGMSHTool::lineLoopPointIdsNotConnected( detail::GeoToolLineLoop const& lineloop ) const
{
    std::set<int> pointIdNotConnected;
    for ( int thelineIdA : lineloop.listLine() )
    {
        int thelineIdPosA = ( thelineIdA>0 )? thelineIdA : -thelineIdA;
        auto const& thelineA = M_buildDataLine.find(thelineIdPosA)->second;
        int firstPointIdA = (thelineIdA>0)? thelineA.firstPointIdConnection() : thelineA.secondPointIdConnection();
        int secondPointIdA = (thelineIdA>0)? thelineA.secondPointIdConnection() : thelineA.firstPointIdConnection();
        bool findInternalPointFirst = false, findInternalPointSecond = false;
        for ( int thelineIdB : lineloop.listLine() )
        {
            if ( thelineIdA == thelineIdB ) continue;
            int thelineIdPosB = ( thelineIdB>0 )? thelineIdB : -thelineIdB;
            auto const& thelineB = M_buildDataLine.find(thelineIdPosB)->second;
            int firstPointIdB = (thelineIdB>0)? thelineB.firstPointIdConnection() : thelineB.secondPointIdConnection();
            int secondPointIdB = (thelineIdB>0)? thelineB.secondPointIdConnection() : thelineB.firstPointIdConnection();

            if ( firstPointIdA == secondPointIdB )
                findInternalPointFirst=true;
            if ( secondPointIdA == firstPointIdB  )
                findInternalPointSecond=true;

            if ( findInternalPointFirst && findInternalPointSecond )
                break;
        }
        if ( !findInternalPointFirst )
            pointIdNotConnected.insert( firstPointIdA );
        if ( !findInternalPointSecond )
            pointIdNotConnected.insert( secondPointIdA );
    }
    return pointIdNotConnected;
}
bool
GeoGMSHTool::lineLoopHasConnection( detail::GeoToolLineLoop const& lineloop1, detail::GeoToolLineLoop const& lineloop2 ) const
{
    // list of point in lineloop1 not connected
    std::set<int> pointIdNotConnected1 = this->lineLoopPointIdsNotConnected(lineloop1);
    std::set<int> pointIdNotConnected2 = this->lineLoopPointIdsNotConnected(lineloop2);
    for( int ptId1 : pointIdNotConnected1 )
    {
        CHECK( M_buildDataPoint.find(ptId1) != M_buildDataPoint.end() ) << " not find point " << ptId1;
        auto const& thePt1 = M_buildDataPoint.find(ptId1)->second;

        //bool findRelatedPoint=false;
        for( int ptId2 : pointIdNotConnected2 )
        {
            CHECK( M_buildDataPoint.find(ptId2) != M_buildDataPoint.end() ) << " not find point " << ptId2;
            auto const& thePt2 = M_buildDataPoint.find(ptId2)->second;
            if ( thePt1.representSameNode( thePt2 ) )
                return true;
        }
    }

    return false;
}

void
GeoGMSHTool::lineLoopApplyConnection( detail::GeoToolLineLoop & lineloop1, detail::GeoToolLineLoop const& lineloop2 )
{
    std::list<int> myNewLinesId;
    int lineIdStart = lineloop1.listLine().front();
    int lineIdStartPos = (lineIdStart>0)? lineIdStart: -lineIdStart;
    auto const& thelineStart = M_buildDataLine.find(lineIdStartPos)->second;
    int pointIdFront = (lineIdStart>0)? thelineStart.firstPointIdConnection() : thelineStart.secondPointIdConnection();
    int pointIdBack = (lineIdStart>0)? thelineStart.secondPointIdConnection() : thelineStart.firstPointIdConnection();

    myNewLinesId.push_back( lineIdStart );

    std::list<int> lineIdsConcatenate;
    for (  int thelineId : lineloop1.listLine() )
        lineIdsConcatenate.push_back(thelineId);
    for (  int thelineId : lineloop2.listLine() )
        lineIdsConcatenate.push_back(thelineId);


    bool hasFinish=false;
    while ( !hasFinish )
    {
        int nLineAdded=0;
        for ( int thelineId : lineIdsConcatenate )
        {
            if ( std::find(myNewLinesId.begin(),myNewLinesId.end(), thelineId ) != myNewLinesId.end() )
                continue;

            int thelineIdPos = ( thelineId>0 )? thelineId : -thelineId;
            auto const& theline = M_buildDataLine.find(thelineIdPos)->second;

            CHECK( M_buildDataPoint.find(pointIdBack) != M_buildDataPoint.end() ) << " not find point " << pointIdBack;
            auto const& thePtBack = M_buildDataPoint.find(pointIdBack)->second;
            CHECK( M_buildDataPoint.find(pointIdFront) != M_buildDataPoint.end() ) << " not find point " << pointIdFront;
            auto const& thePtFront = M_buildDataPoint.find(pointIdFront)->second;
            auto const& thePtFirstPtConnection = M_buildDataPoint.find(theline.firstPointIdConnection())->second;
            auto const& thePtSecondPtConnection = M_buildDataPoint.find(theline.secondPointIdConnection())->second;

            if ( pointIdBack == theline.firstPointIdConnection() )
            {
                myNewLinesId.push_back( thelineIdPos/*thelineId*/ ); // positive sign
                pointIdBack = theline.secondPointIdConnection();
                ++nLineAdded;
            }
            else if ( pointIdBack == theline.secondPointIdConnection() )
            {
                myNewLinesId.push_back( -thelineIdPos/*-thelineId*/ ); // negative sign
                pointIdBack = theline.firstPointIdConnection();
                ++nLineAdded;
            }
            else if ( thePtBack.representSameNode( thePtFirstPtConnection ) )
            {
                myNewLinesId.push_back( thelineIdPos/*thelineId*/ ); // positive sign
                pointIdBack = theline.secondPointIdConnection();
                // replace point in Line
                M_buildDataLine[thelineIdPos].replacePointId( thePtFirstPtConnection.globalId(),thePtBack.globalId() );
                ++nLineAdded;
            }
            else if ( thePtBack.representSameNode( thePtSecondPtConnection ) )
            {
                myNewLinesId.push_back( -thelineIdPos/*-thelineId*/ ); // negative sign
                pointIdBack = theline.firstPointIdConnection();
                // replace point in Line
                M_buildDataLine[thelineIdPos].replacePointId( thePtSecondPtConnection.globalId(),thePtBack.globalId() );
                ++nLineAdded;
            }

            // if line has been added, go to next line
            if ( nLineAdded > 0 )
                continue;

            if ( pointIdFront == theline.secondPointIdConnection() )
            {
                myNewLinesId.push_front( thelineIdPos/*thelineId*/ ); // positive sign
                pointIdFront = theline.firstPointIdConnection();
                ++nLineAdded;
            }
            else if ( pointIdFront == theline.firstPointIdConnection() )
            {
                myNewLinesId.push_front( -thelineIdPos/*-thelineId*/ ); // negative sign
                pointIdFront = theline.secondPointIdConnection();
                ++nLineAdded;
            }
            else if ( thePtFront.representSameNode( thePtSecondPtConnection ) )
            {
                myNewLinesId.push_back( thelineIdPos/*thelineId*/ ); // positive sign
                pointIdFront = theline.firstPointIdConnection();
                // replace point in Line
                M_buildDataLine[thelineIdPos].replacePointId( thePtSecondPtConnection.globalId(),thePtFront.globalId() );
                ++nLineAdded;
            }
            else if ( thePtFront.representSameNode( thePtFirstPtConnection ) )
            {
                myNewLinesId.push_back( -thelineIdPos/*-thelineId*/ ); // negative sign
                pointIdFront = theline.secondPointIdConnection();
                // replace point in Line
                M_buildDataLine[thelineIdPos].replacePointId( thePtFirstPtConnection.globalId(),thePtFront.globalId() );
                ++nLineAdded;
            }

        } // for ( int thelineId : lineIdsConcatenate )

        if ( nLineAdded == 0 )
            hasFinish=true;

    } //  while ( !hasFinish )


    // if front point and back point have same node but not same id, fix numbering
    auto const& thePtBack = M_buildDataPoint.find(pointIdBack)->second;
    auto const& thePtFront = M_buildDataPoint.find(pointIdFront)->second;
    if ( pointIdFront != pointIdBack && thePtFront.representSameNode( thePtBack ) )
    {
        for ( int thelineId : myNewLinesId )
        {
            int thelineIdPos = ( thelineId>0 )? thelineId : -thelineId;
            M_buildDataLine[thelineIdPos].replacePointId( thePtFront.globalId(),thePtBack.globalId() );
        }
    }
#if 0
    std::cout << "RES : ";
    for ( auto res : myNewLinesId )
        std::cout << res << " ";
    std::cout << "\n";
#endif

    lineloop1.setLines( myNewLinesId );

}




GeoGMSHTool::GeoGMSHTool( uint16_type __dim,
                          std::string __shape,
                          std::string __name,
                          double __meshSize )
    :
    M_dim( __dim ),
    M_name(__name),
    M_cptPt( 1 ),
    M_cptLine( 1 ),
    M_cptLineLoop( 1 ),
    M_cptSurface( 1 ),
    M_cptTableau( 1 ),
    M_cptSurfaceLoop( 1 ),
    M_cptVolume( 1 ),
    M_ligneList( new ligne_name_type() ),
    M_surfaceList( new surface_name_type() ),
    M_volumeList( new volume_name_type() ),
    M_surfaceLoopList( new surfaceloop_name_type() ),
    M_ostrExtrude( new std::ostringstream() ),
    M_ostrSurfaceLoop( new std::ostringstream() ),
    M_paramShape( new parameter_shape_type() ),
    M_markShape( new marker_type_type() ),
    M_ostr( new std::ostringstream() ),
    M_ostrDefineByUser( new std::ostringstream() ),
    M_geoIsDefineByUser( false )
{
}

GeoGMSHTool::GeoGMSHTool( uint16_type __dim,
                          std::string const & geoUserStr,
                          double __meshSize,
                          std::string __shape,
                          std::string __name )
    :
    M_dim( __dim ),
    M_name(__name),
    M_cptPt( 1 ),
    M_cptLine( 1 ),
    M_cptLineLoop( 1 ),
    M_cptSurface( 1 ),
    M_cptTableau( 1 ),
    M_cptSurfaceLoop( 1 ),
    M_cptVolume( 1 ),
    M_ligneList( new ligne_name_type() ),
    M_surfaceList( new surface_name_type() ),
        M_volumeList( new volume_name_type() ),
        M_surfaceLoopList( new surfaceloop_name_type() ),
        M_ostrExtrude( new std::ostringstream() ),
        M_ostrSurfaceLoop( new std::ostringstream() ),
        M_paramShape( new parameter_shape_type() ),
        M_markShape( new marker_type_type() ),
        M_ostr( new std::ostringstream() ),
        M_ostrDefineByUser( new std::ostringstream() ),
        M_geoIsDefineByUser( true )
{
    *M_ostrDefineByUser << geoUserStr;
}

GeoGMSHTool::GeoGMSHTool( GeoGMSHTool const & m )
    :
    M_dim( m. M_dim ),
    M_name( m.M_name ),
    M_cptPt( m.M_cptPt ),
    M_cptLine( m.M_cptLine ),
    M_cptLineLoop( m.M_cptLineLoop ),
    M_cptSurface( m.M_cptSurface ),
    M_cptTableau( m.M_cptTableau ),
    M_cptSurfaceLoop( m.M_cptSurfaceLoop ),
    M_cptVolume( m.M_cptVolume ),
    M_ligneList( new ligne_name_type( *( m.M_ligneList ) ) ),
    M_surfaceList( new surface_name_type( *( m.M_surfaceList ) ) ),
    M_volumeList( new volume_name_type( *( m.M_volumeList ) ) ),
    M_surfaceLoopList( new surfaceloop_name_type( *m.M_surfaceLoopList ) ),
    M_ostrExtrude( new std::ostringstream() ),
    M_ostrSurfaceLoop( new std::ostringstream() ),
    M_paramShape( new parameter_shape_type( *( m.M_paramShape ) ) ),
    M_markShape( new marker_type_type( *( m.M_markShape ) ) ),
    M_ostr( new std::ostringstream() ),
    M_ostrDefineByUser( new std::ostringstream() ),
    M_geoIsDefineByUser( m.M_geoIsDefineByUser ),
    M_buildDataPoint( m.M_buildDataPoint),
    M_buildDataLine( m.M_buildDataLine ),
    M_buildDataLineLoop( m.M_buildDataLineLoop ),
                           M_buildDataSurface( m.M_buildDataSurface ),
                           M_buildDataSurfaceLoop( m.M_buildDataSurfaceLoop ),
                           M_buildDataVolume( m.M_buildDataVolume ),
                           M_fusionMarkersLineWithInterface( m.M_fusionMarkersLineWithInterface ),
                           M_fusionMarkersLineWithoutInterface( m.M_fusionMarkersLineWithoutInterface )
{
    updateOstr( ( m.M_ostr )->str() );
    *M_ostrExtrude << ( m.M_ostrExtrude )->str();
    *M_ostrSurfaceLoop << ( m.M_ostrSurfaceLoop )->str();

    if ( M_geoIsDefineByUser ) *M_ostrDefineByUser  << ( m.M_ostrDefineByUser )->str();
}

void
GeoGMSHTool::zeroCpt()
{
    M_cptPt=1;
    M_cptLine=1;
    M_cptLineLoop=1;
    M_cptSurface=1;
    M_cptTableau=1;
    M_cptSurfaceLoop=1;
    M_cptVolume=1;


    ligne_name_type::iterator itLigne = this->M_ligneList->begin();
    ligne_name_type::iterator itLigne_end = this->M_ligneList->end();

    for ( ; itLigne != itLigne_end; ++itLigne )
    {
        ligne_type_type::iterator itLigne2 = itLigne->begin();
        ligne_type_type::iterator itLigne2_end = itLigne->end();

        for ( ; itLigne2 != itLigne2_end; ++itLigne2 )
        {
            boost::get<2>( *itLigne2 )=0;
        }
    }


    auto itSurf = this->M_surfaceList->begin();
    auto itSurf_end = this->M_surfaceList->end();

    for ( ; itSurf != itSurf_end; ++itSurf )
    {
        auto itSurf2 = itSurf->begin();
        auto itSurf2_end = itSurf->end();

        for ( ; itSurf2 != itSurf2_end; ++itSurf2 )
        {
            itSurf2->get<2>() = std::make_pair( 0,0 ); //.clear();
            //boost::get<2>(*itSurf2)=0;
        }
    }

    M_ostrExtrude.reset( new std::ostringstream() );
    M_ostrSurfaceLoop.reset( new std::ostringstream() );

    auto itVol = this->M_volumeList->begin();
    auto itVol_end = this->M_volumeList->end();

    for ( ; itVol != itVol_end; ++itVol )
    {
        auto itVol2 = itVol->begin();
        auto itVol2_end = itVol->end();

        for ( ; itVol2 != itVol2_end; ++itVol2 )
        {
            itVol2->get<2>() = std::make_pair( 0,0 ); //.clear();
            //boost::get<2>(*itVol2)=0;
        }
    }

    auto surfaceLoop_it = this->M_surfaceLoopList->begin();
    auto surfaceLoop_en = this->M_surfaceLoopList->end();

    for ( ; surfaceLoop_it!=surfaceLoop_en ; ++surfaceLoop_it )
    {
        auto surfaceLoop2_it =surfaceLoop_it->begin();
        auto surfaceLoop2_en =surfaceLoop_it->end();

        for ( ; surfaceLoop2_it!=surfaceLoop2_en ; ++surfaceLoop2_it )
            surfaceLoop2_it->get<2>().clear();
    }

} // end zeroCpt

void
GeoGMSHTool::operator=( GeoGMSHTool const & m )
{
    M_dim = m.M_dim;
    M_name = m.M_name;
    M_cptPt = m.M_cptPt;
    M_cptLine = m.M_cptLine;
    M_cptLineLoop = m.M_cptLineLoop;
    M_cptSurface = m.M_cptSurface;
    M_cptTableau = m.M_cptTableau;
    M_cptSurfaceLoop = m.M_cptSurfaceLoop;
    M_cptVolume = m.M_cptVolume;

    M_ligneList.reset( new ligne_name_type( *( m.M_ligneList ) ) );
    M_surfaceList.reset( new surface_name_type( *( m.M_surfaceList ) ) );
    M_volumeList.reset( new volume_name_type( *( m.M_volumeList ) ) );
    M_surfaceLoopList.reset( new surfaceloop_name_type( *m.M_surfaceLoopList ) );
    M_ostrExtrude.reset( new std::ostringstream() );
    *M_ostrExtrude << ( m.M_ostrExtrude )->str();
    M_ostrSurfaceLoop.reset( new std::ostringstream() );
    *M_ostrSurfaceLoop << ( m.M_ostrSurfaceLoop )->str();

    M_paramShape.reset( new parameter_shape_type( *( m.M_paramShape ) ) );

    M_markShape.reset( new marker_type_type( *( m.M_markShape ) ) );
    M_ostr.reset( new std::ostringstream() );
    updateOstr( ( m.M_ostr )->str() );

    //new
    M_buildDataPoint = m.M_buildDataPoint;
    M_buildDataLine = m.M_buildDataLine;
    M_buildDataLineLoop = m.M_buildDataLineLoop;
    M_buildDataSurface = m.M_buildDataSurface;
    M_buildDataSurfaceLoop = m.M_buildDataSurfaceLoop;
    M_buildDataVolume = m.M_buildDataVolume;
    M_fusionMarkersLineWithInterface = m.M_fusionMarkersLineWithInterface;
    M_fusionMarkersLineWithoutInterface = m.M_fusionMarkersLineWithoutInterface;

}

void
GeoGMSHTool::initData( std::string __shape,
                       std::string __name,
                       double __meshSize,
                       std::vector<GeoTool::Node> & __param,
                       uint16_type dim,
                       uint16_type __nbligne,
                       uint16_type __nbsurface,
                       uint16_type __nbvolume )
{
    boost::tuple<std::string,double> __id = boost::make_tuple( __name, __meshSize );

    ( *( M_paramShape ) )[__shape][__name].resize( __param.size() );

    for ( uint16_type n=0; n<__param.size(); ++n )
    {
        ( *( M_paramShape ) )[__shape][__name][n] = __param[n].getNode();
    }



    if ( dim>=1 )
    {
        //Attention 0 par defaut pour dire que ce n'est pas initialiser
        for ( uint16_type n=0; n<__nbligne; ++n )
        {
            //std::list< boost::tuple<std::string,std::string, uint16_type  >	>__listTemp;
            ligne_type_type __listTemp;
            __listTemp.push_back( boost::make_tuple( __shape,__name,0,__meshSize ) );
            M_ligneList->push_back( __listTemp );
        }
    }

    if ( dim>=2 )
    {
        for ( uint16_type n=0; n<__nbsurface; ++n )
        {
            //Attention 0 par defaut pour dire que ce n'est pas initialiser
            std::pair<int,int> listEmpty = std::make_pair( 0,0 );
            surface_type_type __listTemp;
            __listTemp.push_back( boost::make_tuple( __shape,__name,listEmpty,__meshSize ) );
            M_surfaceList->push_back( __listTemp );
        }
    }

    if ( dim==3 )
    {
        //Attention 0 par defaut pour dire que ce n'est pas initialiser
        for ( uint16_type n=0; n<__nbvolume; ++n )
        {
            std::pair<int,int> listEmpty = std::make_pair( 0,0 );
            surface_type_type __listTemp;
            __listTemp.push_back( boost::make_tuple( __shape,__name,listEmpty,__meshSize ) );
            M_volumeList->push_back( __listTemp );
        }

        std::map<int,std::list<int> > listEmpty;
        listEmpty.clear();
        surfaceloop_type_type __listTemp;
        __listTemp.clear();
        __listTemp.push_back( boost::make_tuple( __shape,__name,listEmpty ) );
        M_surfaceLoopList->push_back( __listTemp );
    }

}

void
GeoGMSHTool::updateData( GeoGMSHTool const & m )
{
    M_cptPt = m.M_cptPt;
    M_cptLine = m.M_cptLine;
    M_cptLineLoop = m.M_cptLineLoop;
    M_cptTableau = m.M_cptTableau;
    M_cptSurfaceLoop = m.M_cptSurfaceLoop;

    M_cptSurface = m.M_cptSurface;
    M_cptVolume = m.M_cptVolume;

    M_paramShape = m.M_paramShape;
    M_markShape = m.M_markShape;

    M_ligneList.reset( new ligne_name_type( *( m.M_ligneList ) ) );
    M_surfaceList.reset( new surface_name_type( *( m.M_surfaceList ) ) );
    M_volumeList.reset( new volume_name_type( *( m.M_volumeList ) ) );
    M_surfaceLoopList.reset( new surfaceloop_name_type( *m.M_surfaceLoopList ) );

    M_ostrExtrude.reset( new std::ostringstream() );
    *M_ostrExtrude << ( m.M_ostrExtrude )->str();

    M_ostrSurfaceLoop.reset( new std::ostringstream() );
    *M_ostrSurfaceLoop << ( m.M_ostrSurfaceLoop )->str();

    // new
    M_buildDataPoint = m.M_buildDataPoint;
    M_buildDataLine = m.M_buildDataLine;
    M_buildDataLineLoop = m.M_buildDataLineLoop;
    M_buildDataSurface = m.M_buildDataSurface;
    M_buildDataSurfaceLoop = m.M_buildDataSurfaceLoop;
    M_buildDataVolume = m.M_buildDataVolume;
    M_fusionMarkersLineWithInterface = m.M_fusionMarkersLineWithInterface;
    M_fusionMarkersLineWithoutInterface = m.M_fusionMarkersLineWithoutInterface;

}

void
GeoGMSHTool::showMe() const
{
    std::cout << "cptPt = " << M_cptPt << "\n"
              << "cptLine = " << M_cptLine << "\n"
              << "cptLineLoop = " << M_cptLineLoop << "\n"
              << "cptTableau = " << M_cptTableau << "\n"
              << "cptSurfaceLoop = " << M_cptSurfaceLoop << "\n"
              << "cptSurface = " << M_cptSurface << "\n"
              << "cptVolume = " << M_cptVolume << "\n";

    if ( this->dim()==2 )
    {
        surface_name_const_iterator_type itSurfff = M_surfaceList->begin();
        surface_name_const_iterator_type itSurfff_end = M_surfaceList->end();
        for ( ; itSurfff != itSurfff_end; ++itSurfff )
        {
            std::cout << "[";
            surface_type_const_iterator_type itSurfff2 = itSurfff->begin();
            surface_type_const_iterator_type itSurfff2_end = itSurfff->end();
            for ( ; itSurfff2 != itSurfff2_end ; ++itSurfff2 )
            {
                std::string Qshape = boost::get<0>( *itSurfff2 );
                std::string Qname = boost::get<1>( *itSurfff2 );
                double QmeshSize=boost::get<3>( *itSurfff2 );
                std::cout << "(" << Qshape << "," << Qname << "," << QmeshSize << ")";
                //listPPP.push_back(boost::make_tuple(Qshape,Qname,QmeshSize));
                //setPPP.insert( boost::make_tuple( Qshape,Qname,QmeshSize ) );
            }
            std::cout << "]";
        }
        std::cout << "\n";
    }

    std::cout << "Physical Markers\n";

    for ( auto const& markShape : *M_markShape )
        for ( auto const& markName : markShape.second )
        {
            std::cout << "[" << markShape.first << "][" << markName.first << "] : ";
            for ( auto const& themark : markName.second )
                std::cout << " (" << themark.get<0>() << "," << themark.get<1>() << "," << themark.get<2>() << ") ";
            std::cout << "\n";
        }

}


GeoGMSHTool
GeoGMSHTool::operator-( const GeoGMSHTool & m )
{
    return this->opFusion( m,2 );
}

GeoGMSHTool
GeoGMSHTool::operator+( const GeoGMSHTool & m )
{
    return this->opFusion( m,1 );
}



GeoGMSHTool
GeoGMSHTool::opFusion( const GeoGMSHTool & m,int __typeOp )
{

    GeoGMSHTool __geoTool( this->dim() );

    if ( __typeOp==1 && this->dim()==1 )
    {
        //__geoTool.M_ligneList.reset(new surface_name_type(*(this->M_ligneList)));
        __geoTool.M_ligneList.reset( new ligne_name_type( *( this->M_ligneList ) ) );
        ligne_name_const_iterator_type itLine = m.M_ligneList->begin();
        ligne_name_const_iterator_type itLine_end = m.M_ligneList->end();

        for ( ; itLine != itLine_end; ++itLine )
        {
            ligne_type_const_iterator_type itLine2 = itLine->begin();
            ligne_type_const_iterator_type itLine2_end = itLine->end();
            ligne_type_type __listTemp;

            for ( ; itLine2 != itLine2_end; ++itLine2 )
            {
                __listTemp.push_back( *itLine2 );
            }

            __geoTool.M_ligneList->push_back( __listTemp );
        }
    }

    //Add Plane Surface for operator + : (((rect,u1,_)))+(((circ,u2,_))) -> (((rect,u1,_)),((circ,u2,_)))
    if ( ( __typeOp==1 && this->dim()>=2 ) || ( __typeOp==2 && this->dim()==3 ) )
    {
        __geoTool.M_surfaceList.reset( new surface_name_type( *( this->M_surfaceList ) ) );
        surface_name_const_iterator_type itSurf = m.M_surfaceList->begin();
        surface_name_const_iterator_type itSurf_end = m.M_surfaceList->end();

        for ( ; itSurf != itSurf_end; ++itSurf )
        {
            surface_type_const_iterator_type itSurf2 = itSurf->begin();
            surface_type_const_iterator_type itSurf2_end = itSurf->end();
            surface_type_type __listTemp;

            for ( ; itSurf2 != itSurf2_end; ++itSurf2 )
            {
                __listTemp.push_back( *itSurf2 );
            }

            __geoTool.M_surfaceList->push_back( __listTemp );
        }
    }

    // Add Plane Surface for operator - : (((rect,u1,_)))-(((circ,u2,_))) -> (((rect,u1,_),(circ,u2,_)))
    else if ( __typeOp==2 && this->dim()==2 )
    {
        __geoTool.M_surfaceList.reset( new surface_name_type( *( this->M_surfaceList ) ) );
        surface_name_const_iterator_type itSurf = m.M_surfaceList->begin();
        surface_name_const_iterator_type itSurf_end = m.M_surfaceList->end();

        for ( ; itSurf != itSurf_end; ++itSurf )
        {
            surface_type_const_iterator_type itSurf2 = itSurf->begin();
            surface_type_const_iterator_type itSurf2_end = itSurf->end();
            surface_type_type __listTemp;

            for ( ; itSurf2 != itSurf2_end; ++itSurf2 )
            {
                __geoTool.M_surfaceList->begin()->push_back( *itSurf2 );
            }
        }
    }


    // Add surfaceLoop
    if (this->dim()==3)
        {
            __geoTool.M_surfaceLoopList.reset( new surfaceloop_name_type( *( this->M_surfaceLoopList ) ) );
            surfaceloop_name_const_iterator_type itSurfLoop = m.M_surfaceLoopList->begin();
            surfaceloop_name_const_iterator_type itSurfLoop_end = m.M_surfaceLoopList->end();
            for ( ; itSurfLoop != itSurfLoop_end; ++itSurfLoop )
                {
                    surfaceloop_type_const_iterator_type itSurfLoop2 = itSurfLoop->begin();
                    surfaceloop_type_const_iterator_type itSurfLoop2_end = itSurfLoop->end();

                    for ( ; itSurfLoop2 != itSurfLoop2_end; ++itSurfLoop2 )
                        {
                            __geoTool.M_surfaceLoopList->begin()->push_back( *itSurfLoop2 );
                        }
                }
        }

    //Add Volume for operator + : (((rect,u1,_)))+(((circ,u2,_))) -> (((rect,u1,_)),((circ,u2,_)))
    if ( __typeOp==1 && this->dim()==3 )
    {
        __geoTool.M_volumeList.reset( new volume_name_type( *( this->M_volumeList ) ) );
        volume_name_const_iterator_type itVol = m.M_volumeList->begin();
        volume_name_const_iterator_type itVol_end = m.M_volumeList->end();

        for ( ; itVol != itVol_end; ++itVol )
        {
            volume_type_const_iterator_type itVol2 = itVol->begin();
            volume_type_const_iterator_type itVol2_end = itVol->end();
            volume_type_type __listTemp;

            for ( ; itVol2 != itVol2_end; ++itVol2 )
            {
                __listTemp.push_back( *itVol2 );
            }

            __geoTool.M_volumeList->push_back( __listTemp );
        }
    }

    // Add Volume for operator - : (((rect,u1,_)))-(((circ,u2,_))) -> (((rect,u1,_),(circ,u2,_)))
    else if ( __typeOp==2 && this->dim()==3 )
    {

        __geoTool.M_volumeList.reset( new volume_name_type( *( this->M_volumeList ) ) );
        volume_name_const_iterator_type itVol = m.M_volumeList->begin();
        volume_name_const_iterator_type itVol_end = m.M_volumeList->end();

        for ( ; itVol != itVol_end; ++itVol )
        {
            volume_type_const_iterator_type itVol2 = itVol->begin();
            volume_type_const_iterator_type itVol2_end = itVol->end();
            volume_type_type __listTemp;

            for ( ; itVol2 != itVol2_end; ++itVol2 )
            {
                __geoTool.M_volumeList->begin()->push_back( *itVol2 );
            }
        }
    }

    //get data from this (easy)
    __geoTool.M_paramShape.reset( new parameter_shape_type( *this->M_paramShape ) );
    __geoTool.M_markShape.reset( new marker_type_type ( *this->M_markShape ) );


    if ( this->dim()==1 )
    {
        //get data from (more hard because no duplication)
        ligne_name_const_iterator_type itShape = m.M_ligneList->begin();
        ligne_name_const_iterator_type itShape_end = m.M_ligneList->end();

        for ( ; itShape != itShape_end ; ++itShape )
        {
            auto itName = itShape->begin();
            auto itName_end = itShape->end();

            for ( ; itName != itName_end ; ++itName )
            {
                ( *( __geoTool. M_paramShape ) )[boost::get<0>( *itName )][ boost::get<1>( *itName )]
                    = m.getParameter( boost::get<0>( *itName ),boost::get<1>( *itName ) );
            }
        }
    }

    else if ( this->dim()==2 )
    {
        //get data from (more hard because no duplication)
        surface_name_const_iterator_type itShape = m.M_surfaceList->begin();
        surface_name_const_iterator_type itShape_end = m.M_surfaceList->end();

        for ( ; itShape != itShape_end ; ++itShape )
        {
            auto itName = itShape->begin();
            auto itName_end = itShape->end();

            for ( ; itName != itName_end ; ++itName )
            {
                ( *( __geoTool. M_paramShape ) )[boost::get<0>( *itName )][ boost::get<1>( *itName )]
                    = m.getParameter( boost::get<0>( *itName ),boost::get<1>( *itName ) );
            }
        }
    }

    else if ( this->dim()==3 )
    {
        //get data from (more hard because no duplication)
        volume_name_const_iterator_type itShape = m.M_volumeList->begin();
        volume_name_const_iterator_type itShape_end = m.M_volumeList->end();

        for ( ; itShape != itShape_end ; ++itShape )
        {
            auto itName = itShape->begin();
            auto itName_end = itShape->end();

            for ( ; itName != itName_end ; ++itName )
            {
                ( *( __geoTool. M_paramShape ) )[boost::get<0>( *itName )][ boost::get<1>( *itName )]
                    = m.getParameter( boost::get<0>( *itName ),boost::get<1>( *itName ) );
            }
        }
    }

    //update marker
    marker_type_const_iterator_type itMarkType = m.markerTypeBegin();
    marker_type_const_iterator_type itMarkType_end = m.markerTypeEnd();

    while ( itMarkType!=itMarkType_end )
    {
        if ( ( __typeOp==1 ) ||
                ( itMarkType->first=="point" ) ||
                ( itMarkType->first=="line" && this->dim()>=2 ) ||
                ( __typeOp==2 && itMarkType->first=="surface" && this->dim()==3 ) )
        {
            marker_markerName_const_iterator_type itMarkName = m.markerMarkerNameBegin( itMarkType->first );
            marker_markerName_const_iterator_type itMarkName_end = m.markerMarkerNameEnd( itMarkType->first );

            while ( itMarkName!=itMarkName_end )
            {
                //Est-ce UTILE?????
                if ( !m.getMarkerName( itMarkType->first, itMarkName->first ).empty() )
                {
                    std::vector/*list*/<marker_base_type>::const_iterator itLRef=m.markerListIndiceBegin( itMarkType->first,
                            itMarkName->first );
                    std::vector/*list*/<marker_base_type>::const_iterator itLRef_end=m.markerListIndiceEnd( itMarkType->first,
                            itMarkName->first );

                    while ( itLRef != itLRef_end )
                    {
                        ( *( __geoTool.M_markShape ) )[itMarkType->first][itMarkName->first].push_back( *itLRef );
                        ++itLRef;
                    }
                }

                ++itMarkName;
            }
        }

        ++itMarkType;
    }

    return __geoTool;
}


void//std::string
GeoGMSHTool::init( int orderGeo, std::string gmshFormatVersion,
                   double hmin,double hmax,int refine,
                   bool optimize3dNetgen,
                   GMSH_PARTITIONER partitioner, int partitions, bool partition_file )
{
    //fait dans gmsh.cpp
    *M_ostr << "Mesh.MshFileVersion = "<< gmshFormatVersion <<";\n"
             << "Mesh.CharacteristicLengthExtendFromBoundary=1;\n"
             << "Mesh.CharacteristicLengthFromPoints=1;\n"
             << "Mesh.ElementOrder="<< orderGeo <<";\n"
             << "Mesh.SecondOrderIncomplete = 0;\n"
             << "Mesh.Algorithm = 6;\n";
#if defined(HAVE_TETGEN)
    *M_ostr << "Mesh.Algorithm3D = 1;\n"; // 3D mesh algorithm (1=Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D)
#else
    *M_ostr << "Mesh.Algorithm3D = 4;\n";
#endif
    *M_ostr << "Mesh.RecombinationAlgorithm=0;\n"; // Mesh recombination algorithm (0=standard, 1=blossom)

    if (optimize3dNetgen)
        *M_ostr << "Mesh.OptimizeNetgen=1;\n";
    else
        *M_ostr << "Mesh.OptimizeNetgen=0;\n";

    *M_ostr << "Mesh.CharacteristicLengthMin=" << hmin << ";\n"
             << "Mesh.CharacteristicLengthMax=" << hmax << ";\n"
             << "// refine level " << refine << "\n" // just to rewrite file if refine change
             << "// partitioning data\n"
             << "Mesh.Partitioner=" << partitioner << ";\n"
             << "Mesh.NbPartitions=" << partitions << ";\n"
             << "Mesh.MshFilePartitioned=" << partition_file << ";\n";


    //<< "Mesh.CharacteristicLengthFromCurvature=1;\n";
}


void
GeoGMSHTool::geoStr()
{
    //this->showMe();
#if 0 //
    //data memory ( type->shape->name )
    std::vector<std::map<std::string,std::map<std::string, std::map<uint16_type,uint16_type> > > > __dataMemGlob( 6 );
    std::vector<std::map<std::string,std::map<std::string, std::map<uint16_type,bool> > > > __dataMemGlobSurf1( 2 );
    std::vector<std::map<std::string,std::map<std::string, std::map<uint16_type,std::string> > > > __dataMemGlobSurf2( 2 );
    std::vector<std::map<std::string,std::map<std::string, std::map<uint16_type,bool> > > > __dataMemGlobIsRuled( 1 );
    std::vector<std::map<std::string,std::map<std::string, std::map<uint16_type,std::list<uint16_type> > > > > __dataMemGlobPtsInSurf( 1 );
    // type -> name -> num surfLoop -> list de surfLoop
    std::map<std::string,std::map<std::string, std::map<int,std::list<int> > > > __dataMemGlobSurfaceLoop;
    __dataMemGlobSurfaceLoop.clear();
#endif //
    //construction list ordonne d'objet a construire
    std::list<boost::tuple<std::string,std::string,double> > listPPP;

    std::set<boost::tuple<std::string,std::string,double> > setPPP;

    if ( this->dim()==1 )
    {
        ligne_name_const_iterator_type itLigne = M_ligneList->begin();
        ligne_name_const_iterator_type itLigne_end = M_ligneList->end();

        for ( ; itLigne != itLigne_end; ++itLigne )
        {
            ligne_type_const_iterator_type itLigne2 = itLigne->begin();
            ligne_type_const_iterator_type itLigne2_end = itLigne->end();

            for ( ; itLigne2 != itLigne2_end ; ++itLigne2 )
            {
                std::string Qshape = boost::get<0>( *itLigne2 );
                std::string Qname = boost::get<1>( *itLigne2 );
                double QmeshSize=boost::get<3>( *itLigne2 );
                //listPPP.push_back(boost::make_tuple(Qshape,Qname,QmeshSize));
                setPPP.insert( boost::make_tuple( Qshape,Qname,QmeshSize ) );
            }
        }
    }


    if ( this->dim()==2 )
    {
        surface_name_const_iterator_type itSurfff = M_surfaceList->begin();
        surface_name_const_iterator_type itSurfff_end = M_surfaceList->end();

        for ( ; itSurfff != itSurfff_end; ++itSurfff )
        {
            surface_type_const_iterator_type itSurfff2 = itSurfff->begin();
            surface_type_const_iterator_type itSurfff2_end = itSurfff->end();

            for ( ; itSurfff2 != itSurfff2_end ; ++itSurfff2 )
            {
                std::string Qshape = boost::get<0>( *itSurfff2 );
                std::string Qname = boost::get<1>( *itSurfff2 );
                double QmeshSize=boost::get<3>( *itSurfff2 );
                //listPPP.push_back(boost::make_tuple(Qshape,Qname,QmeshSize));
                setPPP.insert( boost::make_tuple( Qshape,Qname,QmeshSize ) );
            }
        }
    }

    else if ( this->dim()==3 )
    {

        volume_name_const_iterator_type itSurfff = M_volumeList->begin();
        volume_name_const_iterator_type itSurfff_end = M_volumeList->end();

        for ( ; itSurfff != itSurfff_end; ++itSurfff )
        {
            volume_type_const_iterator_type itSurfff2 = itSurfff->begin();
            volume_type_const_iterator_type itSurfff2_end = itSurfff->end();

            for ( ; itSurfff2 != itSurfff2_end ; ++itSurfff2 )
            {
                std::string Qshape = boost::get<0>( *itSurfff2 );
                std::string Qname = boost::get<1>( *itSurfff2 );
                double QmeshSize=boost::get<3>( *itSurfff2 );
                //listPPP.push_back(boost::make_tuple(Qshape,Qname,QmeshSize));
                setPPP.insert( boost::make_tuple( Qshape,Qname,QmeshSize ) );
            }
        }


    }

#if 0
    auto itList = listPPP.begin();
    auto itList_end = listPPP.end();
#else
    auto itList = setPPP.begin();
    auto itList_end = setPPP.end();
#endif

    for ( ; itList!=itList_end; ++itList )
    {
        std::string Qshape = boost::get<0>( *itList );
        std::string Qname = boost::get<1>( *itList );
        //std::cout << "\n Qshape="<<Qshape <<" Qname="<<Qname<<std::endl;
        //*M_ostr << "h=" << boost::get<2>( *itList ) << ";\n";


        //local data memory
        vec_map_data_ptrtype __dataMem( new vec_map_data_type( 8 ) );
        vec_map_data_surf1_ptrtype __dataMemSurf1( new vec_map_data_surf1_type( 2 ) );
        vec_map_data_surf2_ptrtype __dataMemSurf2( new vec_map_data_surf2_type( 2 ) );
        vec_map_data_surf1_ptrtype __dataMemIsRuled( new vec_map_data_surf1_type( 1 ) );
        vec_map_data_ptsinsurf_ptrtype __dataMemPtsInSurf( new vec_map_data_ptsinsurf_type( 1 ) );

        map_surfaceLoop_type __dataMemLocSurfaceLoop;
        __dataMemLocSurfaceLoop.clear();

        GeoGMSHTool_ptrtype __geoTool( new GeoGMSHTool( this->dim() ) );
        __geoTool->updateData( *this );
#if 0 //
        __geoTool->cleanOstr();
#endif
        GeoTool::data_geo_ptrtype __data_geoTool( new GeoTool::data_geo_type( boost::make_tuple( __geoTool,
                __dataMem,
                Qshape,//itShape->first,
                Qname,//boost::get<0>(*itName),
                __dataMemSurf1,
                __dataMemSurf2,
                __dataMemIsRuled,
                __dataMemPtsInSurf,
                __dataMemLocSurfaceLoop,
                itList->get<2>() ) ) );

        // generate the code for the geometry
        run( __data_geoTool );
#if 0
        std::cout << "Qshape " << Qshape << " Qname" << Qname << "\n";
        for ( auto const& mapPt : ( *( boost::get<1>( *__data_geoTool ) ) )[0] )
        {
            std::cout << "(" << mapPt.first << "," << mapPt.second << ");";
        }
        std::cout << "\n\n";
#endif
#if 0 //
        __dataMemGlob[0][Qshape][Qname] = ( *( boost::get<1>( *__data_geoTool ) ) )[0]; //pts
        __dataMemGlob[1][Qshape][Qname] = ( *( boost::get<1>( *__data_geoTool ) ) )[1]; //lines
        __dataMemGlob[2][Qshape][Qname] = ( *( boost::get<1>( *__data_geoTool ) ) )[2]; //lineLoop
        __dataMemGlob[3][Qshape][Qname] = ( *( boost::get<1>( *__data_geoTool ) ) )[3]; //Surface
        __dataMemGlob[4][Qshape][Qname] = ( *( boost::get<1>( *__data_geoTool ) ) )[4]; //SurfaceLoop
        __dataMemGlob[5][Qshape][Qname] = ( *( boost::get<1>( *__data_geoTool ) ) )[5]; //Volume
        __dataMemGlobSurf1[0][Qshape][Qname] = ( *( boost::get<4>( *__data_geoTool ) ) )[0]; // bool : surface is tab gmsh
        __dataMemGlobSurf2[0][Qshape][Qname] = ( *( boost::get<5>( *__data_geoTool ) ) )[0]; // string : surface name tab :
        __dataMemGlobSurf1[1][Qshape][Qname] = ( *( boost::get<4>( *__data_geoTool ) ) )[1]; // bool : volume is tab gmsh
        __dataMemGlobSurf2[1][Qshape][Qname] = ( *( boost::get<5>( *__data_geoTool ) ) )[1]; // string : volume name tab :
        __dataMemGlobIsRuled[0][Qshape][Qname] = ( *( boost::get<6>( *__data_geoTool ) ) )[0]; // bool : surface is ruled
        __dataMemGlobPtsInSurf[0][Qshape][Qname] = ( *( boost::get<7>( *__data_geoTool ) ) )[0]; // list of uint16_type : pts in surface

        __dataMemGlobSurfaceLoop[Qshape][Qname] = __data_geoTool->get<8>();
#endif //
        // get infos
        this->updateData( *boost::get<0>( *__data_geoTool ) );
#if 0 //
        this->updateOstr( boost::get<0>( *__data_geoTool )->M_ostr->str() );
#endif

    }

    //this->showMe();

    //--------------------------------------------------------------------------//
    //--------------------------------------------------------------------------//
    // apply fusion of line/surface and operator +/-
    //--------------------------------------------------------------------------//
    //--------------------------------------------------------------------------//

    std::set<int> ptIdErased,lineIdErased,lineloopIdErased,surfaceIdErased,volumeIdErased;

    if ( this->dim() == 2 )
    {
        //---------------------------------------------------------------------//
        //---------------------------------------------------------------------//
        // CASE 1 : with interface
        M_buildApplyFusionMarker.resize( this->fusionMarkersLineWithInterface().size() );
        for ( auto const& mylinelooppair : this->M_buildDataLineLoop )
        {
            auto const& mylineloop = mylinelooppair.second;
            for ( int mylineGlobalId : mylineloop.listLine() )
            {
                int mylineGlobalIdPos = ( mylineGlobalId > 0 )? mylineGlobalId : -mylineGlobalId;
                auto const& myline = M_buildDataLine[mylineGlobalIdPos];
                int lineLocalId = myline.localId();

                //M_buildApplyFusionMarker
                //auto itfm = std::find_if( this->fusionMarkersLineWithInterface().begin(), this->fusionMarkersLineWithInterface().end(),
                //                          [lineLocalId](FusionMarkers const& fm ) { return fm.marker1() == lineLocalId; } );
                bool find1=false,find2=false;
                int fusId1=0,fusId2=0;
                for ( int k=0; k<this->fusionMarkersLineWithInterface().size() && ( !find1 || !find2 )  ; ++k )
                {
                    if ( this->fusionMarkersLineWithInterface()[k].nameGT1() == mylineloop.name() &&
                         this->fusionMarkersLineWithInterface()[k].marker1() == lineLocalId )
                    {
                        fusId1=k;find1=true;
                        //std::cout << "find marker1 " << lineLocalId << "\n";
                        M_buildApplyFusionMarker[k].setGlobalId1(mylineGlobalId, mylineloop.globalId() /*myline.globalId()*/);
                    }
                    if ( this->fusionMarkersLineWithInterface()[k].nameGT2() == mylineloop.name() &&
                         this->fusionMarkersLineWithInterface()[k].marker2() == lineLocalId )
                    {
                        fusId2=k;find2=true;
                        //std::cout << "find marker2 " << lineLocalId << "\n";
                        M_buildApplyFusionMarker[k].setGlobalId2( mylineGlobalId, mylineloop.globalId() /*myline.globalId()*/);
                    }
                }
            }
        }

        //this->updateFusionMarkersLineWithInterface();
        for ( int k=0; k<this->fusionMarkersLineWithInterface().size(); ++k )
        {
            int mylineGlobalId1 = M_buildApplyFusionMarker[k].globalId1Line();
            int mylineGlobalId1Pos = ( mylineGlobalId1 > 0 )? mylineGlobalId1 : -mylineGlobalId1;
            auto const& myline1 = M_buildDataLine.find(mylineGlobalId1Pos)->second;

            int mylineGlobalId2 = M_buildApplyFusionMarker[k].globalId2Line();
            int mylineGlobalId2Pos = ( mylineGlobalId2 > 0 )? mylineGlobalId2 : -mylineGlobalId2;
            auto const& myline2 = M_buildDataLine.find(mylineGlobalId2Pos)->second;
            //std::cout << "identify same line " << mylineGlobalId1Pos << " vs " << mylineGlobalId2Pos << "\n";
            int mylineloopGlobalId2 = M_buildApplyFusionMarker[k].globalId2LineLoop();

            std::map<int,int> mapPointToReplace;
            if ( this->hasSameOrientation( myline1,myline2 ) )
            {
                // add here a map for lineloop if 
                M_buildDataLineLoop[mylineloopGlobalId2].replaceLineId( mylineGlobalId2, mylineGlobalId1 );

                auto it = myline1.listPt().begin();
                auto en = myline1.listPt().end();
                auto it2 = myline2.listPt().begin();
                for ( ; it != en ;++it,++it2 )
                {
                    mapPointToReplace[*it2] = *it;
                }
                //std::cout << "has same orientation\n";
            }
            else
            {
                M_buildDataLineLoop[mylineloopGlobalId2].replaceLineId( mylineGlobalId2, -mylineGlobalId1 );
                //std::cout << "has NOT same orientation\n";

                auto it = myline1.listPt().begin();
                auto en = myline1.listPt().end();
                auto it2 = myline2.listPt().end();--it2;
                for ( ; it != en ;++it,--it2 )
                {
                    mapPointToReplace[*it2] = *it;
                }
            }
            //M_buildDataLine.erase( mylineGlobalId2Pos );
            lineIdErased.insert(mylineGlobalId2Pos);
            // replace point id in line ( TODO reduce search)
            for ( auto& mylinepair : M_buildDataLine )
                for ( auto const& ptToRep : mapPointToReplace )
                    mylinepair.second.replacePointId( ptToRep.first, ptToRep.second );

            for ( auto const& ptToRep : mapPointToReplace )
                ptIdErased.insert( ptToRep.first );
        }

        //---------------------------------------------------------------------//
        //---------------------------------------------------------------------//
        // CASE 2 : without interface

        // define new surface to build from initial surface (with localId of line interface)
        // newNameSurf -> ( name base surf -> list of local line Id )
        std::map<std::string,std::map<std::string,std::set<int> > > mapNewSurface;
        int cptNewSurfaceAdded=0;
        for ( int k=0; k<this->fusionMarkersLineWithoutInterface().size(); ++k )
        {
            std::string name1 = this->fusionMarkersLineWithoutInterface()[k].nameGT1();
            std::string name2 = this->fusionMarkersLineWithoutInterface()[k].nameGT2();
            int markerId1 = this->fusionMarkersLineWithoutInterface()[k].marker1();
            int markerId2 = this->fusionMarkersLineWithoutInterface()[k].marker2();
            // search if one of surface has already register
            std::string findName1,findName2;
            for ( auto const& newSurf : mapNewSurface )
            {
                if ( findName1.empty() && newSurf.second.find( name1 ) != newSurf.second.end() )
                    findName1 = newSurf.first;
                if ( findName2.empty() && newSurf.second.find( name2 ) != newSurf.second.end() )
                    findName2 = newSurf.first;
                if ( !findName1.empty() && !findName2.empty() ) break;
            }
            if ( findName1.empty() && findName2.empty() )
            {
                std::string newNameRegister = ( boost::format("ConcatenateSurface%1%")%cptNewSurfaceAdded ).str();
                //mapNewSurface[newNameRegister].insert(name1);
                mapNewSurface[newNameRegister][name1].insert(markerId1);
                mapNewSurface[newNameRegister][name2].insert(markerId2);
                ++cptNewSurfaceAdded;
            }
            else if ( !findName1.empty() && findName2.empty() )
            {
                //mapNewSurface[findName1].insert(name2);
                mapNewSurface[findName1][name1].insert(markerId1);
                mapNewSurface[findName1][name2].insert(markerId2);
            }
            else if ( findName1.empty() && !findName2.empty() )
            {
                //mapNewSurface[findName2].insert(name1);
                mapNewSurface[findName2][name1].insert(markerId1);
                mapNewSurface[findName2][name2].insert(markerId2);
            }
            else
            {
                CHECK( false ) << "TODO";
            }
        } // for ( int k=0; ... )

#if 0
        for ( auto const& newSurf : mapNewSurface )
        {
            std::cout << "[" << newSurf.first << "] -> ";
            for ( auto surfName : newSurf.second )
            {
                std::cout << "[" << surfName.first << " : ";
                for ( auto surfLocId : surfName.second )
                    std::cout << " " << surfLocId;
                std::cout << "]";
            }
            std::cout << "\n";
        }
#endif

        // create new surfaces
        for ( auto const& newSurf : mapNewSurface )
        {
            // create new surface
            int idNewSurf = this->M_cptSurface;
            ++this->M_cptSurface; //update counter
            std::string surfType = "plane";//surfaceRegisterFront.surfaceType(); // plane, ruled
            detail::GeoToolSurface myFusionSurface( newSurf.first,invalid_size_type_value,idNewSurf,surfType );

            std::vector<detail::GeoToolLineLoop> newLineLoopMemory;
            // first pass : compute new line loop (not closed) and identify lineloops and surfaces to erase
            // add all lineloops and take one surface marker
            for ( auto surfName : newSurf.second )
            {
                int idSurf = this->surfaceIdFromName( surfName.first );
                CHECK( idSurf > 0 ) << "surface not find";
                auto const& mysurf = M_buildDataSurface.find( idSurf )->second;
                if ( myFusionSurface.physicalMarker().empty() )
                    myFusionSurface.setPhysicalMarker( mysurf.physicalMarker() );

                //myFusionSurface.addLineLoop( mysurf.listLineLoop() );

                std::set<int> thelineIdErasedWithFusion;
                for ( int localIdLineFusion : surfName.second )
                {
                    // search global line id define in fusion with surface name and local id
                    // must be apply with initial object
                    for ( int thelineloopId : mysurf.listLineLoop() )
                    {
                        CHECK( M_buildDataLineLoop.find(thelineloopId) != M_buildDataLineLoop.end() ) << "lineloop " << thelineloopId << "not exist";
                        auto const& thelineloop = M_buildDataLineLoop.find(thelineloopId)->second;
                        for ( int thelineId : thelineloop.listLine() )
                        {
                            int thelineIdPos = ( thelineId > 0 )? thelineId : -thelineId;
                            CHECK( M_buildDataLine.find(thelineIdPos) != M_buildDataLine.end() ) << "line " << thelineIdPos << "not exist";
                            auto const& theline = M_buildDataLine.find(thelineIdPos)->second;

                            if ( theline.name() == mysurf.name() && // igore line already modified by a fusion with interface
                                 theline.localId() == localIdLineFusion )
                            {
                                //std::cout << "newSurf.first " << newSurf.first << " surfName.first " << surfName.first << " theline.globalId() "<< theline.globalId() << "\n";
                                thelineIdErasedWithFusion.insert( theline.globalId() );
                                lineIdErased.insert( theline.globalId() );
                            }
                        }
                    }
                }
#if 0
                std::cout << "thelineIdErasedWithFusion : ";
                for ( auto const& hola : thelineIdErasedWithFusion )
                    std::cout << hola << " ";
                std::cout << "\n";
#endif
                // build new lineloops without interface line ( lineloop are not closed and must be fix after)
                for ( int thelineloopId : mysurf.listLineLoop() )
                {
                    CHECK( M_buildDataLineLoop.find(thelineloopId) != M_buildDataLineLoop.end() ) << "lineloop " << thelineloopId << "not exist";
                    auto const& thelineloop = M_buildDataLineLoop.find(thelineloopId)->second;
                    std::list<int> myNewLinesId;
                    for ( int thelineId : thelineloop.listLine() )
                    {
                        int thelineIdPos = ( thelineId > 0 )? thelineId : -thelineId;
                        if ( thelineIdErasedWithFusion.find( thelineIdPos ) == thelineIdErasedWithFusion.end() )
                            myNewLinesId.push_back( thelineId );
                    }

                    if ( myNewLinesId.size() > 0 )
                    {
                        int idNewLineLoop = this->M_cptLineLoop;
                        ++this->M_cptLineLoop; //update counter
                        detail::GeoToolLineLoop myNewLineLoop(myFusionSurface.name(), invalid_size_type_value,idNewLineLoop );
                        myNewLineLoop.setLines( myNewLinesId );
                        newLineLoopMemory.push_back( myNewLineLoop );
                        //this->addLineLoop(myNewLineLoop);
                        //myFusionSurface.addLineLoop( mysurf.listLineLoop() );
                    }
                    lineloopIdErased.insert(thelineloop.globalId());
                }
                //thelineIdErasedWithFusion

                surfaceIdErased.insert( mysurf.globalId() );
            } // for ( auto surfName : )

            // remove interfaces lines in line loop
#if 0
            std::cout << " newLineLoopMemory : \n";
            for ( auto const& idMem : newLineLoopMemory )
                idMem.showMe();
            //std::cout << idMem.globalId() << " ";
            //std::cout << "\n";
#endif

            // reconnect line loops and reconnect points
            CHECK( newLineLoopMemory.size() > 0 ) << "Not Good";
            std::vector<bool> isDone( newLineLoopMemory.size(), false );
            //for ( auto const& thelineloop : newLineLoopMemory )
            for ( int k=0;k<newLineLoopMemory.size();++k )
            {
                if ( !isDone[k] )
                {
                    int idNewLineLoop = this->M_cptLineLoop;
                    ++this->M_cptLineLoop; //update counter
                    detail::GeoToolLineLoop myNewLineLoop(myFusionSurface.name(), invalid_size_type_value,idNewLineLoop );
                    myNewLineLoop.setLines( newLineLoopMemory[k].listLine() );
                    isDone[k]=true;

                    bool lineloopClosed = false;
                    while ( !lineloopClosed/*this->lineLoopIsClosed(myNewLineLoop)*/ )
                    {
                        for ( int k2=0;k2<newLineLoopMemory.size() /*&& !lineloopClosed*//*!this->lineLoopIsClosed(myNewLineLoop)*/;++k2 )
                        {
                            if( !isDone[k2] )
                            {
                                if ( this->lineLoopHasConnection( myNewLineLoop, newLineLoopMemory[k2] ) )
                                {
                                    //std::cout << "HasConnection\n";
                                    this->lineLoopApplyConnection( myNewLineLoop, newLineLoopMemory[k2] );
                                    isDone[k2]=true;

                                    if ( this->lineLoopIsClosed(myNewLineLoop) )
                                    {
                                        lineloopClosed=true;
                                        break;
                                    }
                                }
                            }
                            //myNewLineLoop.showMe();
                        }
                    }
                    CHECK( this->lineLoopIsClosed(myNewLineLoop) ) << "lineloop must be closed here";

                    //std::cout << "finish lineloop\n";
                    //myNewLineLoop.showMe();

                    this->addLineLoop( myNewLineLoop );
                    myFusionSurface.addLineLoop( myNewLineLoop.globalId()  );
                }
            }

            this->addSurface(myFusionSurface);
            //myFusionSurface.showMe();

        } // for ( auto const& newSurf : mapNewSurface )


        //---------------------------------------------------
        // update M_surfaceList after fusion without interface

        //std::map< std::pair<std::string,int>,int > surfaceListModified;
        std::map< std::pair<std::string,int>,std::pair<std::string,int> > surfaceListModified;
        std::map< std::pair<std::string,int>,std::pair<std::string,int> > surfaceListErased;
        if ( mapNewSurface.size() > 0 )
        {
            for ( auto const& newSurf : mapNewSurface )
            {
                std::string nameNewSurf = newSurf.first;
                int idNewSurf = this->surfaceIdFromName( nameNewSurf );
                std::string nameInitialUsedForInsertNewSurf;
                for ( auto surfName : newSurf.second )
                {
                    std::string nameInitialSurf = surfName.first;
                    int idInitialSurf = this->surfaceIdFromName( nameInitialSurf );

                    // use first surface for replacing in surface list
                    if ( nameInitialUsedForInsertNewSurf.empty() )
                    {
                        nameInitialUsedForInsertNewSurf = nameInitialSurf;
                        surfaceListModified[std::make_pair(nameInitialSurf,idInitialSurf)]=std::make_pair(nameNewSurf,idNewSurf);
                    }
                    else
                        surfaceListErased[std::make_pair(nameInitialSurf,idInitialSurf)]=std::make_pair(nameNewSurf,idNewSurf);

                    //surfaceListErased.insert( std::make_pair(nameInitialSurf,idInitialSurf) );
                }
            }
        }
        // 2 cases : front and others(diff)

        std::map<int,std::set<int> > surfaceListMovedSurfFromFrontFusion;
        for ( auto const& itSurf : *this->M_surfaceList )
        {
            CHECK( itSurf.size() > 0 ) << "no surface";
            if ( itSurf.size() == 1 ) continue;

            surface_type_const_iterator_type itSurf2front = itSurf.begin();
            if (itSurf2front->get<2>().second==0) continue; // surface useless
            CHECK( this->M_buildDataSurface.find( itSurf2front->get<2>().first ) != this->M_buildDataSurface.end() ) << "error";

            std::string nameSurf = itSurf2front->get<1>();
            int idSurf = itSurf2front->get<2>().first;
            if ( surfaceListErased.find( std::make_pair( nameSurf,idSurf ) ) != surfaceListErased.end() )
            {
                int idNewSurf = surfaceListErased.find( std::make_pair( nameSurf,idSurf ) )->second.second;
                int thecpt=0;
                //++itSurf2front;
                for ( auto const& itSurf2 : itSurf )
                {
                    if ( thecpt > 0 )
                    {
                        int idSurfDiff = itSurf2.get<2>().first;
                        surfaceListMovedSurfFromFrontFusion[idNewSurf].insert( idSurfDiff );
                    }
                    ++thecpt;
                }
            }
        }

        // move maybe diff surface
        auto itSurf = this->M_surfaceList->begin();
        auto enSurf = this->M_surfaceList->end();
        for ( ; itSurf != enSurf ; ++itSurf )
        {
            //typedef surface_type_type::iterator surface_type_iterator_type;
            //surface_type_iterator_type itSurf2front = itSurf->begin();
            surface_type_const_iterator_type itSurf2front = itSurf->begin();
            if (itSurf2front->get<2>().second==0) continue; // surface useless

            std::string nameSurf = itSurf2front->get<1>();
            int idSurf = itSurf2front->get<2>().first;
            if ( surfaceListModified.find( std::make_pair( nameSurf,idSurf ) ) != surfaceListModified.end() )
            {
                auto findSurfToModify = surfaceListModified.find( std::make_pair( nameSurf,idSurf ) );
                int idNewSurf = findSurfToModify->second.second;
                if ( surfaceListMovedSurfFromFrontFusion.find(idNewSurf) !=surfaceListMovedSurfFromFrontFusion.end() )
                {
                    for ( int theidSurf : surfaceListMovedSurfFromFrontFusion.find(idNewSurf)->second )
                    {
                        auto const& mysurf = this->M_buildDataSurface.find( theidSurf )->second;
                        CHECK( mysurf.listLineLoop().size() == 1 ) << "diff can be used only with surface with one lineloop";
                        itSurf->push_back( boost::make_tuple( "blabla",mysurf.name(),std::make_pair( mysurf.globalId(),mysurf.listLineLoop().front() ),0. ) );
                    }
                }


            }
        }

        // modified and delete fusion surface
        /*auto*/ itSurf = this->M_surfaceList->begin();
        //auto enSurf = this->M_surfaceList->end();
        for ( ; itSurf != enSurf ; /*++itSurf*/ )
        {
            //surface_type_const_iterator_type itSurf2front = itSurf->begin();
            // get first surface (which is not a diff surface)
            typedef surface_type_type::iterator surface_type_iterator_type;
            surface_type_iterator_type itSurf2front = itSurf->begin();
            if ( itSurf2front->get<2>().second==0 ) continue; // surface useless

            std::string nameSurf = itSurf2front->get<1>();
            int idSurf = itSurf2front->get<2>().first;

            bool hasErasedAllSurf=false;

            auto findSurfToModify = surfaceListModified.find( std::make_pair( nameSurf,idSurf ) );
            if ( findSurfToModify != surfaceListModified.end() )
            {
                int idNewSurf = findSurfToModify->second.second;
                CHECK( this->M_buildDataSurface.find( idNewSurf ) != this->M_buildDataSurface.end() ) << "error";
                auto const& mysurf = this->M_buildDataSurface.find( idNewSurf )->second;
                boost::get<0>(*itSurf2front) = "blabla"; // shape (usefull??)
                boost::get<1>(*itSurf2front) = mysurf.name();
                CHECK( mysurf.listLineLoop().size() == 1 ) << "diff can be used only with surface with one lineloop";
                boost::get<2>(*itSurf2front) = std::make_pair( mysurf.globalId(),mysurf.listLineLoop().front() ); // surf id, lineloop id
                boost::get<3>(*itSurf2front) = 0.; // mesh size (usefull??)
            }
            else if ( surfaceListErased.find( std::make_pair( nameSurf,idSurf ) ) != surfaceListErased.end() )
            {
                // delete front surface with all diff surface
                itSurf = this->M_surfaceList->erase( itSurf );
                hasErasedAllSurf=true;
            }


            // treat diff surface
            if ( !hasErasedAllSurf && itSurf->size() > 1 )
            {
                auto itSurfDiff = ++itSurf2front;
                auto enSurfDiff = itSurf->end();
                while ( itSurfDiff != enSurfDiff )
                {
                    std::string nameSurfDiff = itSurfDiff->get<1>();
                    int idSurfDiff = itSurfDiff->get<2>().first;
                    if ( surfaceListErased.find( std::make_pair( nameSurfDiff,idSurfDiff ) ) != surfaceListErased.end() )
                    {
                        // delete diff surface
                        itSurfDiff = itSurf->erase( itSurfDiff );
                    }
                    else
                    {
                        auto findSurfDiffToModify = surfaceListModified.find( std::make_pair( nameSurfDiff,idSurfDiff ) );
                        if ( findSurfDiffToModify != surfaceListModified.end() )
                        {

                            // modify diff surface
                            int idNewSurfDiff = findSurfDiffToModify->second.second;
                            CHECK( this->M_buildDataSurface.find( idNewSurfDiff ) != this->M_buildDataSurface.end() ) << "error";
                            auto const& mysurf = this->M_buildDataSurface.find( idNewSurfDiff )->second;
                            boost::get<0>(*itSurfDiff) = "blablaDiff"; // shape (usefull??)
                            boost::get<1>(*itSurfDiff) = mysurf.name();
                            CHECK( mysurf.listLineLoop().size() == 1 ) << "diff can be used only with surface with one lineloop";
                            boost::get<2>(*itSurfDiff) = std::make_pair( mysurf.globalId(),mysurf.listLineLoop().front() ); // surf id, lineloop id
                            boost::get<3>(*itSurfDiff) = 0.; // mesh size (usefull??)
                        }
                        ++itSurfDiff;
                    }
                }
            }
            if( !hasErasedAllSurf ) ++itSurf;
        } // for ( ; itSurf != enSurf ; ++itSurf )

        //this->showMe();

    } // if ( dim == 2 )

    // create diff surface
    for ( auto const& itSurf : *this->M_surfaceList )
    {
        CHECK( itSurf.size() > 0 ) << "no surface";
        if ( itSurf.size() == 1 ) continue;

        surface_type_const_iterator_type itSurf2front = itSurf.begin();
        if (itSurf2front->get<2>().second==0) continue; // surface useless
        CHECK( this->M_buildDataSurface.find( itSurf2front->get<2>().first ) != this->M_buildDataSurface.end() ) << "error";

        //std::string theshape = itSurf2->get<0>();
        //std::string thename = itSurf2->get<1>();

        int idNewSurf = this->M_cptSurface;
        ++this->M_cptSurface; //update counter

        auto const& surfaceRegisterFront = this->M_buildDataSurface.find( itSurf2front->get<2>().first )->second;

        std::string newSurfName;
        int thecpt=0;
        for ( auto const& itSurf2 : itSurf )
        {
            if ( thecpt > 0 ) newSurfName += "-";
            newSurfName += itSurf2.get<1>();
            ++thecpt;
        }
        std::string surfType = surfaceRegisterFront.surfaceType(); // plane, ruled
        detail::GeoToolSurface mySurf( newSurfName,invalid_size_type_value,idNewSurf,surfType );
        mySurf.setPhysicalMarker( surfaceRegisterFront.physicalMarker() );

        for ( auto const& itSurf2 : itSurf )
        {
            CHECK( this->M_buildDataSurface.find( itSurf2.get<2>().first ) != this->M_buildDataSurface.end() ) << "error";
            for ( auto llId : this->M_buildDataSurface.find( itSurf2.get<2>().first )->second.listLineLoop() )
                mySurf.addLineLoop( llId );

            // search if the surface is present in surface list and there is only this one
            // in this this case, not erase
            bool doEraseSurf = true;
            for ( auto const& itSurfb : *this->M_surfaceList )
                if ( itSurfb.size() == 1 && itSurfb.front().get<2>().first == itSurf2.get<2>().first)
                    doEraseSurf = false;
            if ( doEraseSurf )
                surfaceIdErased.insert( itSurf2.get<2>().first );

        }
        this->addSurface(mySurf);
    }


    // create diff volume
    for ( auto const& itVol : *this->M_volumeList )
    {
        CHECK( itVol.size() > 0 ) << "no surface";
        if ( itVol.size() == 1 ) continue;

        auto itVol2front = itVol.front();
        if (itVol2front.get<2>().second==0) continue; // surface useless
        CHECK( this->M_buildDataVolume.find( itVol2front.get<2>().first ) != this->M_buildDataVolume.end() ) << "error";

        auto const& volumeRegisterFront = this->M_buildDataVolume.find( itVol2front.get<2>().first )->second;

        int idNewVol = this->M_cptVolume;
        ++this->M_cptVolume; //update counter
        detail::GeoToolVolume myVol( invalid_size_type_value,idNewVol );
        myVol.setPhysicalMarker( volumeRegisterFront.physicalMarker() );

        for ( auto const& itVol2 : itVol )
        {
            CHECK( this->M_buildDataVolume.find( itVol2.get<2>().first ) != this->M_buildDataVolume.end() ) << "error";
            for ( auto slId : this->M_buildDataVolume.find( itVol2.get<2>().first )->second.listSurfaceLoop() )
                myVol.addSurfaceLoop( slId );

            // search if the volume is present in volume list and there is only this one
            // in this this case, not erase
            bool doEraseVol = true;
            for ( auto const& itVolb : *this->M_volumeList )
                if ( itVolb.size() == 1 && itVolb.front().get<2>().first == itVol2.get<2>().first)
                    doEraseVol = false;
            if ( doEraseVol )
                volumeIdErased.insert( itVol2.get<2>().first );

        }
        this->addVolume(myVol);
    }




    for ( auto pid : ptIdErased )
        M_buildDataPoint.erase( pid );

    for ( auto lid : lineIdErased )
        M_buildDataLine.erase( lid );

    for ( auto llid : lineloopIdErased )
        M_buildDataLineLoop.erase( llid );

    for ( auto sid : surfaceIdErased )
        M_buildDataSurface.erase( sid );

    for ( auto vid : volumeIdErased )
        M_buildDataVolume.erase( vid );


    //--------------------------------------------------------------------------//
    //--------------------------------------------------------------------------//
    // generate code for all entities in geo file
    //--------------------------------------------------------------------------//
    //--------------------------------------------------------------------------//

    // save markers in these containers
    std::map<std::string,std::set<int> > markerPoints,markerLines,markerSurf,markerVol;

    for ( auto const& mypt : this->M_buildDataPoint )
    {
        //mypt.second.showMe();
        CHECK( mypt.first == mypt.second.globalId() ) << "error";
        std::ostringstream __ostr;
        __ostr << "Point(" << mypt.second.globalId()
               << ") = {"
               << std::scientific << std::setprecision( 16 ) << mypt.second.node()[0] << ","
               << std::scientific << std::setprecision( 16 ) << mypt.second.node()[1] << ","
               << std::scientific << std::setprecision( 16 ) << mypt.second.node()[2] <<","
               << mypt.second.hSize() <<"};\n";
        this->updateOstr( __ostr.str() );

        // save point markers
        if ( !mypt.second.physicalMarker().empty() )
            markerPoints[mypt.second.physicalMarker()].insert( mypt.second.globalId() );
    }
    for ( auto const& myline : this->M_buildDataLine )
    {
        //myline.second.showMe();
        std::ostringstream __ostr;
        if ( myline.second.lineType() == "line" )
            __ostr << "Line";
        else if ( myline.second.lineType() == "circle" )
            __ostr << "Circle";
        else if ( myline.second.lineType() == "ellipse" )
            __ostr << "Ellipse";
        else if ( myline.second.lineType() == "spline" )
            __ostr << "Spline";
        else if ( myline.second.lineType() == "bspline" )
            __ostr << "BSpline";
        else CHECK( false ) << "error";

        __ostr << "(" << myline.second.globalId() << ") = { ";
        int thecptPt=0;
        for ( size_type ptId : myline.second.listPt() )
        {
            if ( thecptPt > 0 )
                __ostr << ",";
            __ostr << ptId;
            ++thecptPt;
        }
        __ostr << "};\n";
        this->updateOstr( __ostr.str() );

        // save line markers
        if ( !myline.second.physicalMarker().empty() )
            markerLines[myline.second.physicalMarker()].insert( myline.second.globalId() );
    }
    for ( auto const& mylineloop : this->M_buildDataLineLoop )
    {
        //mylineloop.second.showMe();
        std::ostringstream __ostr;
        __ostr << "Line Loop(" << mylineloop.second.globalId()
               << ") = {" ;
        int thecptLine=0;
        for ( int lineId : mylineloop.second.listLine() )
        {
            if ( thecptLine > 0 )
                __ostr << ",";
            __ostr << lineId;
            ++thecptLine;
        }
        __ostr << "};\n";
        this->updateOstr( __ostr.str() );
    }
    for ( auto const& mysurf : this->M_buildDataSurface )
    {
        //mysurf.second.showMe();
        std::ostringstream __ostr;
        if ( mysurf.second.surfaceType() == "plane" )
            __ostr << "Plane Surface(";
        else if ( mysurf.second.surfaceType() == "ruled" )
            __ostr << "Ruled Surface(";
        __ostr << mysurf.first << ") = {" ;

        int thecpt=0;
        for ( auto const& lineloopId : mysurf.second.listLineLoop() )
        {
            if ( thecpt > 0 ) __ostr << ",";
            __ostr << lineloopId;
            ++thecpt;
        }
        __ostr << "};\n";

        for ( int ptId : mysurf.second.ptsInSurface() )
            __ostr << "Point{" << ptId << "} In Surface{" << mysurf.first << "};\n";

        this->updateOstr( __ostr.str() );

        // save surface markers
        if ( !mysurf.second.physicalMarker().empty() )
            markerSurf[mysurf.second.physicalMarker()].insert( mysurf.second.globalId() );
    }
    for ( auto const& mysurfloop : this->M_buildDataSurfaceLoop )
    {
        //mysurfloop.second.showMe();
        std::ostringstream __ostr;
        __ostr << "Surface Loop(" <<  mysurfloop.first
               << ") = {" ;
        int thecpt=0;
        for ( int surfId : mysurfloop.second.listSurface() )
        {
            if ( thecpt > 0 )
                __ostr << ",";
            __ostr << surfId;
            ++thecpt;
        }
        __ostr << "};\n";
        this->updateOstr( __ostr.str() );
    }
    for ( auto const& myvol : this->M_buildDataVolume )
    {
        //myvol.second.showMe();
        std::ostringstream __ostr;
        __ostr << "Volume(" <<  myvol.first
               << ") = {" ;
        int thecpt=0;
        for ( int surfLoopId : myvol.second.listSurfaceLoop() )
        {
            if ( thecpt > 0 )
                __ostr << ",";
            __ostr << surfLoopId;
            ++thecpt;
        }
        __ostr << "};\n";
        this->updateOstr( __ostr.str() );

        // save volume markers
        if ( !myvol.second.physicalMarker().empty() )
            markerVol[myvol.second.physicalMarker()].insert( myvol.second.globalId() );
    }

    //--------------------------------------------------------------------------//
    //--------------------------------------------------------------------------//
    // generate code for physical markers in geo file
    //--------------------------------------------------------------------------//
    //--------------------------------------------------------------------------//

    for ( auto const& mymarkPoint : markerPoints )
    {
        *M_ostr << "Physical Point(\"" << mymarkPoint.first << "\") = {";
        int thecpt=0;
        for ( int mypointId : mymarkPoint.second )
        {
            if ( thecpt > 0 ) *M_ostr << ",";
            *M_ostr << mypointId;
            ++thecpt;
        }
        *M_ostr << "};\n";
    }

    for ( auto const& mymarkLine : markerLines )
    {
        *M_ostr << "Physical Line(\"" << mymarkLine.first << "\") = {";
        int thecpt=0;
        for ( int mylineId : mymarkLine.second )
        {
            if ( thecpt > 0 ) *M_ostr << ",";
            *M_ostr << mylineId;
            ++thecpt;
        }
        *M_ostr << "};\n";
    }

    for ( auto const& mymarkSurf : markerSurf )
    {
        *M_ostr << "Physical Surface(\"" << mymarkSurf.first << "\") = {";
        int thecpt=0;
        for ( int mysurfId : mymarkSurf.second )
        {
            if ( thecpt > 0 ) *M_ostr << ",";
            *M_ostr << mysurfId;
            ++thecpt;
        }
        *M_ostr << "};\n";
    }

    for ( auto const& myvol : markerVol )
    {
        *M_ostr << "Physical Volume(\"" << myvol.first << "\") = {";
        int thecpt=0;
        for ( int myvolId : myvol.second )
        {
            if ( thecpt > 0 ) *M_ostr << ",";
            *M_ostr << myvolId;
            ++thecpt;
        }
        *M_ostr << "};\n";
    }


} // geoStr()





#define GEOTOOL_GENERATE_RUN(r,state)                                   \
        if( boost::get<2>(*__dg) ==  GEOTOOL_SHAPE_NAME_STR(BOOST_PP_TUPLE_ELEM(2,0,state)) ) \
            {                                                           \
                BOOST_PP_CAT(run,GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state))) (__dg); \
            }                                                           \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/


// A refaire avec les boost pp
void run( data_geo_ptrtype __dg )
{
    BOOST_PP_FOR( ( 0, BOOST_PP_SUB( BOOST_PP_ARRAY_SIZE( GEOTOOL_SHAPE ),1 ) ),
                  GEOTOOL_INSTANTIATES_FOR_COMP,
                  GEOTOOL_INSTANTIATES_FOR_INCR,
                  GEOTOOL_GENERATE_RUN );
}

/*_________________________________________________*
 *_________________________________________________*
 * Function on the namespace                       *
 *_________________________________________________*
 *_________________________________________________*/

template <uint16_type Numero>
node_type
param( data_geo_ptrtype __dg )
{
    std::string __shape = boost::get<2>( *__dg );
    std::string __name = boost::get<3>( *__dg );
    //node_type __node = boost::get<Numero>(boost::get<0>(*__dg)->M_paramShape->find(__shape)->second[__name]);
    node_type __node =( boost::get<0>( *__dg )->M_paramShape->find( __shape )->second[__name] )[Numero];

    __node.resize( 3 );

    if ( __node.size()<2 )
        __node( 1 )=0.0;

    if ( __node.size()<3 )
        __node( 2 )=0.0;

    return __node;
}

/*_________________________________________________*
 *_________________________________________________*
 * Function on the namespace                       *
 *_________________________________________________*
 *_________________________________________________*/


void writePoint( uint16_type __numLoc, data_geo_ptrtype __dg ,double __x1,double __x2, double __x3 )
{
    ( *( boost::get<1>( *__dg ) ) )[0][__numLoc] = boost::get<0>( *__dg )->cptPt(); //            __mapPt[0][__numLoc] = boost::get<0>(*__dg)->cptPt();

    //auto name = __dg->get<3>();
    detail::GeoToolPoint myPt(__x1,__x2,__x3, __dg->get<9>(), __numLoc, __dg->get<0>()->cptPt() );

    auto mymark = __dg->get<0>()->findPhysicalMarker( "point", __dg->get<3>()/*name*/, __numLoc );
    if ( mymark.first )
        myPt.setPhysicalMarker( mymark.second );

    __dg->get<0>()->addPoint(/*name,*/myPt);

#if 0
    std::ostringstream __ostr;
    __ostr << "Point(" << boost::get<0>( *__dg )->cptPt()
           << ") = {"
           << std::scientific << std::setprecision( 16 ) << __x1 << ","
           << std::scientific << std::setprecision( 16 ) << __x2 << ","
           << std::scientific << std::setprecision( 16 ) << __x3 <<", h};\n";
    boost::get<0>( *__dg )->updateOstr( __ostr.str() );
#endif
    ++( boost::get<0>( *__dg )->M_cptPt );
}

/*_________________________________________________*/

void
writeLine( uint16_type __numLoc, data_geo_ptrtype __dg ,uint16_type __n1, uint16_type __n2 )
{
    ( *( boost::get<1>( *__dg ) ) )[1][__numLoc] = boost::get<0>( *__dg )->cptLine();

    std::string name = __dg->get<3>();
    detail::GeoToolLine myLine(name,__numLoc, __dg->get<0>()->cptLine(), "line" );
    myLine.setPoints( { (*(__dg->get<1>()))[0][__n1], (*(__dg->get<1>()))[0][__n2] } );

    auto mymark = __dg->get<0>()->findPhysicalMarker( "line", __dg->get<3>()/*name*/, __numLoc );
    if ( mymark.first )
        myLine.setPhysicalMarker( mymark.second );
#if 0
    if ( __dg->get<0>()->M_markShape->find( "line") != __dg->get<0>()->M_markShape->end() )
    {
        for ( auto markName : __dg->get<0>()->markerMarkerName("line") )
        {
            for ( auto lineMarked : markName.second )
            {
                //std::string __shape = boost::get<2>( *__dg );
                std::string __name = boost::get<3>( *__dg );
                if ( __name == lineMarked.get<1>() && __numLoc == lineMarked.get<2>() )
                {
                    myLine.setPhysicalMarker( markName.first );
                }
            }
        }
    }
#endif
    __dg->get<0>()->addLine(myLine);

#if 0
    std::ostringstream __ostr;
    __ostr << "Line(" << boost::get<0>( *__dg )->cptLine()
           << ") = {"
           << ( *( boost::get<1>( *__dg ) ) )[0][__n1] << ","
           << ( *( boost::get<1>( *__dg ) ) )[0][__n2] << "};\n";
    boost::get<0>( *__dg )->updateOstr( __ostr.str() );
#endif
    ++boost::get<0>( *__dg )->M_cptLine;
}

/*_________________________________________________*/

void
writeCircle( uint16_type __numLoc, data_geo_ptrtype __dg ,uint16_type __n1, uint16_type __n2, uint16_type __n3 )
{
    ( *( boost::get<1>( *__dg ) ) )[1][__numLoc] = boost::get<0>( *__dg )->cptLine();

    std::string name = __dg->get<3>();
    detail::GeoToolLine myLine(name,__numLoc, __dg->get<0>()->cptLine(), "circle" );
    myLine.setPoints( { (*(__dg->get<1>()))[0][__n1], (*(__dg->get<1>()))[0][__n2], (*(__dg->get<1>()))[0][__n3] } );

    auto mymark = __dg->get<0>()->findPhysicalMarker( "line", __dg->get<3>()/*name*/, __numLoc );
    if ( mymark.first )
        myLine.setPhysicalMarker( mymark.second );

    __dg->get<0>()->addLine(myLine);
#if 0
    std::ostringstream __ostr;
    __ostr << "Circle(" << boost::get<0>( *__dg )->cptLine() //cptCircle
           << ") = {"
           << ( *( boost::get<1>( *__dg ) ) )[0][__n1] << ","
           << ( *( boost::get<1>( *__dg ) ) )[0][__n2] << ","
           << ( *( boost::get<1>( *__dg ) ) )[0][__n3] << "};\n";

    boost::get<0>( *__dg )->updateOstr( __ostr.str() );
#endif
    ++boost::get<0>( *__dg )->M_cptLine;
}

/*_________________________________________________*/

void
writeEllipse( uint16_type __numLoc, data_geo_ptrtype __dg ,uint16_type __n1, uint16_type __n2, uint16_type __n3, uint16_type __n4 )
{
    ( *( boost::get<1>( *__dg ) ) )[1][__numLoc] = boost::get<0>( *__dg )->cptLine();

    std::string name = __dg->get<3>();
    detail::GeoToolLine myLine(name,__numLoc, __dg->get<0>()->cptLine(), "ellipse" );
    myLine.setPoints( { (*(__dg->get<1>()))[0][__n1], (*(__dg->get<1>()))[0][__n2], (*(__dg->get<1>()))[0][__n3], (*(__dg->get<1>()))[0][__n4] } );

    auto mymark = __dg->get<0>()->findPhysicalMarker( "line", __dg->get<3>()/*name*/, __numLoc );
    if ( mymark.first )
        myLine.setPhysicalMarker( mymark.second );

    __dg->get<0>()->addLine(myLine);
#if 0
    std::ostringstream __ostr;
    __ostr << "Ellipse(" << boost::get<0>( *__dg )->cptLine() //cptCircle
           << ") = {"
           << ( *( boost::get<1>( *__dg ) ) )[0][__n1] << ","
           << ( *( boost::get<1>( *__dg ) ) )[0][__n2] << ","
           << ( *( boost::get<1>( *__dg ) ) )[0][__n3] << ","
           << ( *( boost::get<1>( *__dg ) ) )[0][__n4] << "};\n";
    boost::get<0>( *__dg )->updateOstr( __ostr.str() );
#endif

    ++boost::get<0>( *__dg )->M_cptLine;
}

/*_________________________________________________*/

void
writeSpline( uint16_type __numLoc, data_geo_ptrtype __dg ,Loop __loop )
{
    ( *( boost::get<1>( *__dg ) ) )[1][__numLoc] = boost::get<0>( *__dg )->cptLine();

    std::string name = __dg->get<3>();
    detail::GeoToolLine myLine(name,__numLoc, __dg->get<0>()->cptLine(), "spline" );
    std::list<size_type> myptId;
    for ( auto it=__loop.begin(), en=__loop.end() ; it!=en ; ++it )
        myptId.push_back( *it );
    myLine.setPoints( myptId );

    auto mymark = __dg->get<0>()->findPhysicalMarker( "line", __dg->get<3>()/*name*/, __numLoc );
    if ( mymark.first )
        myLine.setPhysicalMarker( mymark.second );

    __dg->get<0>()->addLine(myLine);

#if 0
    std::ostringstream __ostr;
    __ostr << "Spline(" << boost::get<0>( *__dg )->cptLine() //cptCircle
           << ") = {";

    std::list<int>::const_iterator it= __loop.begin();
    std::list<int>::const_iterator it_end= --__loop.end();

    while ( it!=it_end )
    {
        __ostr << ( *( boost::get<1>( *__dg ) ) )[0][*it] <<"," ;
        ++it;
    }

    __ostr << ( *( boost::get<1>( *__dg ) ) )[0][*it] << "};\n";

    boost::get<0>( *__dg )->updateOstr( __ostr.str() );
#endif
    ++boost::get<0>( *__dg )->M_cptLine;
}

/*_________________________________________________*/

void
writeBSpline( uint16_type __numLoc, data_geo_ptrtype __dg ,Loop __loop )
{
    ( *( boost::get<1>( *__dg ) ) )[1][__numLoc] = boost::get<0>( *__dg )->cptLine();

    std::string name = __dg->get<3>();
    detail::GeoToolLine myLine(name,__numLoc, __dg->get<0>()->cptLine(), "bspline" );
    std::list<size_type> myptId;
    for ( auto it=__loop.begin(), en=__loop.end() ; it!=en ; ++it )
        myptId.push_back( *it );
    myLine.setPoints( myptId );

    auto mymark = __dg->get<0>()->findPhysicalMarker( "line", __dg->get<3>()/*name*/, __numLoc );
    if ( mymark.first )
        myLine.setPhysicalMarker( mymark.second );

    __dg->get<0>()->addLine(myLine);
#if 0
    std::ostringstream __ostr;
    __ostr << "BSpline(" << boost::get<0>( *__dg )->cptLine() //cptCircle
           << ") = {";

    std::list<int>::const_iterator it= __loop.begin();
    std::list<int>::const_iterator it_end= --__loop.end();

    while ( it!=it_end )
    {
        __ostr << ( *( boost::get<1>( *__dg ) ) )[0][*it] <<"," ;
        ++it;
    }

    __ostr << ( *( boost::get<1>( *__dg ) ) )[0][*it] << "};\n";

    boost::get<0>( *__dg )->updateOstr( __ostr.str() );
#endif
    ++boost::get<0>( *__dg )->M_cptLine;
}

/*_________________________________________________*/

void
writeLineLoop( uint16_type __numLoc, data_geo_ptrtype __dg , Loop /*const*/ __loop )
{
    ( *( boost::get<1>( *__dg ) ) )[2][__numLoc] = boost::get<0>( *__dg )->cptLineLoop();


    std::string name = __dg->get<3>();
    detail::GeoToolLineLoop myLineLoop(name,__numLoc, __dg->get<0>()->cptLineLoop() );
    std::list<int> myLinesId;
    for ( auto itl= __loop.begin(),enl= __loop.end();itl!=enl;++itl )
    {
        if ( *itl>0 )
            myLinesId.push_back( (*__dg->get<1>() )[1][*itl] );
        else
            myLinesId.push_back( -(*__dg->get<1>() )[1][-*itl] );
    }
    myLineLoop.setLines( myLinesId );
    __dg->get<0>()->addLineLoop(myLineLoop);

#if 0
    std::ostringstream __ostr;
    __ostr << "Line Loop(" << boost::get<0>( *__dg )->cptLineLoop()
           << ") = {" ;
    std::list<int>::const_iterator it= __loop.begin();
    std::list<int>::const_iterator it_end= --__loop.end();

    while ( it!=it_end )
    {
        if ( *it>0 )
            __ostr << ( *( boost::get<1>( *__dg ) ) )[1][*it] <<"," ;

        else
            __ostr << "-" << ( *( boost::get<1>( *__dg ) ) )[1][-*it] <<"," ;

        ++it;
    }

    if ( *it>0 )
        __ostr << ( *( boost::get<1>( *__dg ) ) )[1][*it] << "};\n";

    else
        __ostr << "-" << ( *( boost::get<1>( *__dg ) ) )[1][-*it] << "};\n";

    boost::get<0>( *__dg )->updateOstr( __ostr.str() );
#endif
    ++boost::get<0>( *__dg )->M_cptLineLoop;
}

/*_________________________________________________*/

void
writePtInSurface( data_geo_ptrtype __dg , uint16_type __indLocPt,uint16_type __indLocSurf )
{

    auto indPtGlob = ( *( boost::get<1>( *__dg ) ) )[0][__indLocPt];
    // A SUPP la ligne dessous
    ( *( boost::get<7>( *__dg ) ) )[0][( *( boost::get<1>( *__dg ) ) )[3][__indLocSurf]].push_back( indPtGlob );

    auto indSurfGlob = (*__dg->get<1>() )[3][__indLocSurf];
    CHECK(__dg->get<0>()->M_buildDataSurface.find(indSurfGlob) != __dg->get<0>()->M_buildDataSurface.end() ) << "not find surface";
    __dg->get<0>()->M_buildDataSurface[indSurfGlob].addPtInSurface( indPtGlob );
}


//ici on n'ecrit pas, on memorise cause des operations de difference
//l'ecriture est realise dans geoStr()
void
writePlaneSurface( uint16_type __numLoc, data_geo_ptrtype __dg , uint16_type __ind )
{
    ( *( boost::get<1>( *__dg ) ) )[3][__numLoc] = boost::get<0>( *__dg )->cptSurface(); //num local to global
    ( *( boost::get<4>( *__dg ) ) )[0][( *( boost::get<1>( *__dg ) ) )[3][__numLoc]] = false; //is tab gmsh
    ( *( boost::get<5>( *__dg ) ) )[0][( *( boost::get<1>( *__dg ) ) )[3][__numLoc]] = ""; //name of tab
    ( *( boost::get<6>( *__dg ) ) )[0][( *( boost::get<1>( *__dg ) ) )[3][__numLoc]] = false; //isRuled

    bool __find=false;
    //Memorize in surfaceList
    GeoGMSHTool::surface_name_type::iterator itSurf = __dg->get<0>()->M_surfaceList->begin();
    GeoGMSHTool::surface_name_type::iterator itSurf_end = __dg->get<0>()->M_surfaceList->end();

    //uint cptSurf=0, refLineLoop=0;
    uint16_type cptSurf = __dg->get<0>()->cptSurface();
    uint16_type refLineLoop= ( *( boost::get<1>( *__dg ) ) )[2][__ind];


    std::string name = __dg->get<3>();
    detail::GeoToolSurface mySurf( name,__numLoc,cptSurf,"plane" );
    mySurf.setLineLoop( refLineLoop );

    auto mymark = __dg->get<0>()->findPhysicalMarker( "surface", __dg->get<3>()/*name*/, __numLoc );
    if ( mymark.first )
        mySurf.setPhysicalMarker( mymark.second );

    __dg->get<0>()->addSurface(mySurf);


#if 1
    for ( ; itSurf !=itSurf_end; ++itSurf )
    {
        auto itSurf2 = itSurf->begin();
        if (  boost::get<0>( *itSurf2 ) == boost::get<2>( *__dg ) &&  // same shape
              boost::get<1>( *itSurf2 ) == boost::get<3>( *__dg ) ) // same name
        {
            if ( itSurf2->get<2>().second == 0 && !__find )
            {
                itSurf2->get<2>().first = cptSurf;
                itSurf2->get<2>().second = refLineLoop;
                __find=true;
            }
        }
    }

    itSurf = __dg->get<0>()->M_surfaceList->begin();
    bool __find2=false;
#endif
    for ( ; itSurf !=itSurf_end; ++itSurf )
    {
        if (itSurf->size()<2) continue;
        GeoGMSHTool::surface_type_type::iterator itSurf2 = itSurf->begin();++itSurf2;
        GeoGMSHTool::surface_type_type::iterator itSurf2_end = itSurf->end();
        while ( itSurf2 !=itSurf2_end )
        {
            if (  boost::get<0>( *itSurf2 ) == boost::get<2>( *__dg ) &&  // same shape
                  boost::get<1>( *itSurf2 ) == boost::get<3>( *__dg ) ) // same name
                {
                    //on cherche la 1ere surface non init
                    if ( itSurf2->get<2>().second == 0 && !__find2 )
                    {
                        itSurf2->get<2>().first = cptSurf;//__dg->get<0>()->cptSurface();
                        itSurf2->get<2>().second = refLineLoop;// ( *( boost::get<1>( *__dg ) ) )[2][__ind]; // get the lineloop
                        __find2=true;
                    }
                }
            ++itSurf2;
        } //
    } // for ( ; itSurf !=itSurf_end; ++itSurf )

    ++boost::get<0>( *__dg )->M_cptSurface; //update counter
}

//ici on n'ecrit pas, on memorise cause des operations de difference
//l'ecriture est realise dans geoStr()
void
writeRuledSurface( uint16_type __numLoc, data_geo_ptrtype __dg , uint16_type __ind )
{
    ( *( boost::get<1>( *__dg ) ) )[3][__numLoc] = boost::get<0>( *__dg )->cptSurface(); //num local to global
    ( *( boost::get<4>( *__dg ) ) )[0][( *( boost::get<1>( *__dg ) ) )[3][__numLoc]] = false; //is tab gmsh
    ( *( boost::get<5>( *__dg ) ) )[0][( *( boost::get<1>( *__dg ) ) )[3][__numLoc]] = ""; //name of tab
    ( *( boost::get<6>( *__dg ) ) )[0][( *( boost::get<1>( *__dg ) ) )[3][__numLoc]] = true; //isRuled

    uint16_type cptSurf = __dg->get<0>()->cptSurface();
    uint16_type refLineLoop= ( *( boost::get<1>( *__dg ) ) )[2][__ind];

    std::string name = __dg->get<3>();
    detail::GeoToolSurface mySurf( name,__numLoc,cptSurf,"ruled" );
    mySurf.setLineLoop( refLineLoop );

    auto mymark = __dg->get<0>()->findPhysicalMarker( "surface", __dg->get<3>()/*name*/, __numLoc );
    if ( mymark.first )
        mySurf.setPhysicalMarker( mymark.second );

    __dg->get<0>()->addSurface(mySurf);



    bool __find=false;
    //Memorize in surfaceList
    GeoGMSHTool::surface_name_type::iterator itSurf = boost::get<0>( *__dg )->M_surfaceList->begin();
    GeoGMSHTool::surface_name_type::iterator itSurf_end = boost::get<0>( *__dg )->M_surfaceList->end();

#if 0
    for ( ; itSurf !=itSurf_end; ++itSurf )
    {
        GeoGMSHTool::surface_type_type::iterator itSurf2 = itSurf->begin();
        GeoGMSHTool::surface_type_type::iterator itSurf2_end = itSurf->end();

        while ( itSurf2 !=itSurf2_end )
        {
            if ( boost::get<0>( *itSurf2 ) == boost::get<2>( *__dg ) )
            {
                if ( boost::get<1>( *itSurf2 ) == boost::get<3>( *__dg ) )
                {
                    //on cherche la 1ere surface non init
                    if ( itSurf2->get<2>().second == 0 && !__find )
                    {
                        itSurf2->get<2>().first = __dg->get<0>()->cptSurface();
                        itSurf2->get<2>().second = ( *( boost::get<1>( *__dg ) ) )[2][__ind];
                        __find=true;
                    }
                }
            }

            ++itSurf2;
        }
    }
#elif 0
    uint16_type cptSurf=0, refLineLoop=0;
    for ( ; itSurf !=itSurf_end; ++itSurf )
    {
        GeoGMSHTool::surface_type_type::iterator itSurf2 = itSurf->begin();
        GeoGMSHTool::surface_type_type::iterator itSurf2_end = itSurf->end();
        while ( itSurf2 !=itSurf2_end )
        {
            if (  boost::get<0>( *itSurf2 ) == boost::get<2>( *__dg ) &&  // same shape
                  boost::get<1>( *itSurf2 ) == boost::get<3>( *__dg ) ) // same name
                {
                    if ( itSurf2->get<2>().second == 0 )
                    {
                        //on cherche la 1ere surface non init
                        if (!__find)
                        {
                            cptSurf = __dg->get<0>()->cptSurface();
                            refLineLoop= ( *( boost::get<1>( *__dg ) ) )[2][__ind];
                            itSurf2->get<2>().first = cptSurf;//__dg->get<0>()->cptSurface();
                            itSurf2->get<2>().second = refLineLoop;// ( *( boost::get<1>( *__dg ) ) )[2][__ind]; // get the lineloop
                            __find=true;
                        }
                        else if (__dg->get<0>()->dim()==2)
                        {
                            itSurf2->get<2>().first = cptSurf;
                            itSurf2->get<2>().second = refLineLoop;
                        }
                    }
                }

            ++itSurf2;
        } //
    } // for ( ; itSurf !=itSurf_end; ++itSurf )
#else

    for ( ; itSurf !=itSurf_end; ++itSurf )
    {
        auto itSurf2 = itSurf->begin();
        if (  boost::get<0>( *itSurf2 ) == boost::get<2>( *__dg ) &&  // same shape
              boost::get<1>( *itSurf2 ) == boost::get<3>( *__dg ) ) // same name
        {
            if ( itSurf2->get<2>().second == 0 && !__find )
            {
                itSurf2->get<2>().first = cptSurf;
                itSurf2->get<2>().second = refLineLoop;
                __find=true;
            }
        }
    }

    itSurf = __dg->get<0>()->M_surfaceList->begin();
    bool __find2=false;

    for ( ; itSurf !=itSurf_end; ++itSurf )
    {
        if (itSurf->size()<2) continue;
        GeoGMSHTool::surface_type_type::iterator itSurf2 = itSurf->begin();++itSurf2;
        GeoGMSHTool::surface_type_type::iterator itSurf2_end = itSurf->end();
        while ( itSurf2 !=itSurf2_end )
        {
            if (  boost::get<0>( *itSurf2 ) == boost::get<2>( *__dg ) &&  // same shape
                  boost::get<1>( *itSurf2 ) == boost::get<3>( *__dg ) ) // same name
                {
                    //on cherche la 1ere surface non init
                    if ( itSurf2->get<2>().second == 0 && !__find2 )
                    {
                        itSurf2->get<2>().first = cptSurf;//__dg->get<0>()->cptSurface();
                        itSurf2->get<2>().second = refLineLoop;// ( *( boost::get<1>( *__dg ) ) )[2][__ind]; // get the lineloop
                        __find2=true;
                    }
                }
            ++itSurf2;
        } //
    } // for ( ; itSurf !=itSurf_end; ++itSurf )
#endif

    ++boost::get<0>( *__dg )->M_cptSurface;
}

void
writeExtrudeSurface( uint16_type __numLoc,data_geo_ptrtype __dg , uint16_type __ind,Loop /*const*/ __loop )
{
    Loop __loopDef = Loop()>>0;

    for ( uint16_type i=2; i< __loop.size()+2; ++i )
        __loopDef>>i;

    std::list<int>::const_iterator itDef = __loopDef.begin();
    std::ostringstream __nametab;
    __nametab << "tableau"<<boost::get<0>( *__dg )->M_cptTableau;

    //<< "out23[] = Extrude {-delta_gap,0,0} { Surface{60}; } ;\n"
    std::ostringstream __ostr;
    __ostr << __nametab.str()
           << "[]= Extrude { 4 ,0,0} { Surface{"
           << ( *( boost::get<1>( *__dg ) ) )[3][__ind]
           <<"}; } ;\n";

    std::list<int>::const_iterator it= __loop.begin();
    std::list<int>::const_iterator it_end= __loop.end();

    while ( it!=it_end )
    {
        ( *( boost::get<1>( *__dg ) ) )[3][*it] = *itDef; //boost::get<0>(*__dg)->cptTableau();
        ( *( boost::get<4>( *__dg ) ) )[0][( *( boost::get<1>( *__dg ) ) )[3][*it]] = true;
        ( *( boost::get<5>( *__dg ) ) )[0][( *( boost::get<1>( *__dg ) ) )[3][*it]] = __nametab.str();
        ( *( boost::get<6>( *__dg ) ) )[0][( *( boost::get<1>( *__dg ) ) )[3][*it]] = false;
        ++it;
        ++itDef;
    }

    //volume
    //gmsh semble creer les volumes dans l'ordre 1,2,3,...
    //boost::get<0>(*__dg)->cptVolume() == tableau[1]
    ( *( boost::get<1>( *__dg ) ) )[5][__numLoc] = boost::get<0>( *__dg )->cptVolume();
    //(*(boost::get<1>(*__dg)))[3][*it] = 1;
    ( *( boost::get<4>( *__dg ) ) )[1][__numLoc] = true;
    ( *( boost::get<5>( *__dg ) ) )[1][__numLoc] = __nametab.str();

    *( boost::get<0>( *__dg )->M_ostrExtrude ) << __ostr.str();

    ++boost::get<0>( *__dg )->M_cptTableau;
    ++boost::get<0>( *__dg )->M_cptVolume;

}


void
writeSurfaceLoop( uint16_type __numLoc, data_geo_ptrtype __dg , Loop /*const*/ __loop )
{
    ( *( boost::get<1>( *__dg ) ) )[4][__numLoc] = __dg->get<0>()->cptSurfaceLoop(); // num local to global


    std::string name = __dg->get<3>();
    detail::GeoToolSurfaceLoop mySurfaceLoop(name,__numLoc, __dg->get<0>()->cptSurfaceLoop() );
    std::list<int> mySurfaceId;
    for ( auto itl= __loop.begin(),enl= __loop.end();itl!=enl;++itl )
    {
        if ( *itl>0 )
            mySurfaceId.push_back( (*__dg->get<1>() )[3][*itl] );
        else
            mySurfaceId.push_back( -(*__dg->get<1>() )[3][-*itl] );
    }
    mySurfaceLoop.setSurfaces( mySurfaceId );
    __dg->get<0>()->addSurfaceLoop(mySurfaceLoop);



    auto surfaceLoop_it = __dg->get<0>()->M_surfaceLoopList->begin();
    auto surfaceLoop_en =  __dg->get<0>()->M_surfaceLoopList->end();

    for ( ; surfaceLoop_it!=surfaceLoop_en ; ++surfaceLoop_it )
    {
        auto surfaceLoop2_it =surfaceLoop_it->begin();
        auto surfaceLoop2_en =surfaceLoop_it->end();

        for ( ; surfaceLoop2_it!=surfaceLoop2_en ; ++surfaceLoop2_it )
        {
            if ( ( surfaceLoop2_it->get<0>() == __dg->get<2>() ) &&
                    ( surfaceLoop2_it->get<1>() == __dg->get<3>() ) ) // search shape and name
            {
                surfaceLoop2_it->get<2>()[__dg->get<0>()->cptSurfaceLoop()].clear();
                auto loop_it= __loop.begin();
                auto loop_en= __loop.end();

                for ( ; loop_it!=loop_en ; ++loop_it )
                    surfaceLoop2_it->get<2>()[__dg->get<0>()->cptSurfaceLoop()].push_back( *loop_it );
            }

        }

    }

    ++boost::get<0>( *__dg )->M_cptSurfaceLoop;

}

//ici on n'ecrit pas, on memorise cause des operations de difference
//l'ecriture est realise dans geoStr()
void
writeVolume( uint16_type __numLoc, data_geo_ptrtype __dg , uint16_type __ind )
{
    ( *( boost::get<1>( *__dg ) ) )[5][__numLoc] = boost::get<0>( *__dg )->cptVolume(); //num local to global



    uint16_type cptVol = __dg->get<0>()->cptVolume();
    uint16_type refSurfaceLoop= ( *( boost::get<1>( *__dg ) ) )[4][__ind];

    detail::GeoToolVolume myVol(__numLoc,cptVol );
    myVol.setSurfaceLoop( refSurfaceLoop );

    auto mymark = __dg->get<0>()->findPhysicalMarker( "volume", __dg->get<3>()/*name*/, __numLoc );
    if ( mymark.first )
        myVol.setPhysicalMarker( mymark.second );

    __dg->get<0>()->addVolume(myVol);



    bool __find=false;
    //Memorize in volumeList
    GeoGMSHTool::volume_name_type::iterator itVol = boost::get<0>( *__dg )->M_volumeList->begin();
    GeoGMSHTool::volume_name_type::iterator itVol_end = boost::get<0>( *__dg )->M_volumeList->end();

#if 0
    for ( ; itSurf !=itSurf_end; ++itSurf )
    {
        GeoGMSHTool::volume_type_type::iterator itSurf2 = itSurf->begin();
        GeoGMSHTool::volume_type_type::iterator itSurf2_end = itSurf->end();

        while ( itSurf2 !=itSurf2_end )
        {
            if ( boost::get<0>( *itSurf2 ) == boost::get<2>( *__dg ) ) // same shape
            {
                if ( boost::get<1>( *itSurf2 ) == boost::get<3>( *__dg ) ) // same name
                {
                    //on cherche la 1ere surface non init
                    if ( itSurf2->get<2>().second == 0 && !__find )
                    {
                        itSurf2->get<2>().first = __dg->get<0>()->cptVolume();
                        itSurf2->get<2>().second = ( *( boost::get<1>( *__dg ) ) )[4][__ind];

                        //boost::get<2>(*itSurf2) = (*(boost::get<1>(*__dg)))[4][__ind];
                        __find=true;
                    }
                }
            }

            ++itSurf2;
        }
    }
#else

    for ( ; itVol !=itVol_end; ++itVol )
    {
        auto itVol2 = itVol->begin();
        if (  boost::get<0>( *itVol2 ) == boost::get<2>( *__dg ) &&  // same shape
              boost::get<1>( *itVol2 ) == boost::get<3>( *__dg ) ) // same name
        {
            if ( itVol2->get<2>().second == 0 && !__find )
            {
                itVol2->get<2>().first = cptVol;
                itVol2->get<2>().second = refSurfaceLoop;
                __find=true;
            }
        }
    }

    itVol = __dg->get<0>()->M_volumeList->begin();
    bool __find2=false;

    for ( ; itVol !=itVol_end; ++itVol )
    {
        if (itVol->size()<2) continue;
        auto itVol2 = itVol->begin();++itVol2;
        auto itVol2_end = itVol->end();
        while ( itVol2 !=itVol2_end )
        {
            if (  boost::get<0>( *itVol2 ) == boost::get<2>( *__dg ) &&  // same shape
                  boost::get<1>( *itVol2 ) == boost::get<3>( *__dg ) ) // same name
                {
                    //on cherche la 1ere surface non init
                    if ( itVol2->get<2>().second == 0 && !__find2 )
                    {
                        itVol2->get<2>().first = cptVol;
                        itVol2->get<2>().second = refSurfaceLoop;
                        __find2=true;
                    }
                }
            ++itVol2;
        }
    } // for ( ; itVol !=itVol_end; ++itVol )


#endif

    ++boost::get<0>( *__dg )->M_cptVolume; //update counter
}

/*_________________________________________________*/

boost::tuple<Node,Node,Node>
computeBasisOrthogonal( node_type dir,node_type centre )
{

    double norm_dir=std::sqrt( dir( 0 )*dir( 0 )+dir( 1 )*dir( 1 )+dir( 2 )*dir( 2 ) );
    dir( 0 )=dir( 0 )/norm_dir;
    dir( 1 )=dir( 1 )/norm_dir;
    dir( 2 )=dir( 2 )/norm_dir;

    // plane coefficient equation
    double a=dir( 0 );
    double b=dir( 1 );
    double c=dir( 2 );
    //double d=-dir( 0 )*centre( 0 )-dir( 1 )*centre( 1 )-dir( 2 )*centre( 2 ); //-N scalaire OA
    double d=0;// consider origin (0,0,0)
    double rayon=1;

    // found a point in plane other that the center
    Node ptBis( 0,0,0 );
    Node ptPerturb( dir( 0 )+1,dir( 1 ),dir( 2 ) );
    double normPtPerturb = math::sqrt(math::pow(ptPerturb(0),2)+math::pow(ptPerturb(1),2)+math::pow(ptPerturb(2),2) );
    double vdotn = ptPerturb(0)*dir(0)+ptPerturb(1)*dir(1)+ptPerturb(2)*dir(2);

    if ( std::abs(std::abs(vdotn/normPtPerturb)-1) < 1e-9 )
    {
        // if colinear else change pt perturbation
        ptPerturb(0)=dir(0);ptPerturb(1)=dir(1)+1;ptPerturb(2)=dir(2);
        normPtPerturb = math::sqrt(math::pow(ptPerturb(0),2)+math::pow(ptPerturb(1),2)+math::pow(ptPerturb(2),2) );
        vdotn = ptPerturb(0)*dir(0)+ptPerturb(1)*dir(1)+ptPerturb(2)*dir(2);
    }
    CHECK( std::abs(std::abs(vdotn/normPtPerturb)-1) > 1e-9 ) << "colinear to normal : vdotn=" << vdotn
                                                              << " normal=" << dir(0) << "," << dir(1) << "," << dir(2)
                                                              << " v=" << ptPerturb(0) << "," << ptPerturb(1) << "," << ptPerturb(2)
                                                              << "\n";
    ptBis(0)=ptPerturb(0)-vdotn*dir(0);
    ptBis(1)=ptPerturb(1)-vdotn*dir(1);
    ptBis(2)=ptPerturb(2)-vdotn*dir(2);

    CHECK( std::abs(a*ptBis(0)+b*ptBis(1)+c*ptBis(2)+d) < 1e-9 ) << "point is not on plane\n";

    // first base vector in plane
    Node u( ptBis( 0 )/*-centre( 0 )*/,ptBis( 1 )/*-centre( 1 )*/,ptBis( 2 )/*-centre( 2 )*/ );

    //a=u2v3-u3v2 ; b=u3v1-u1v3 ; c=u1v2-u2v1
    // second base vector in plane (orthogonal basis (u,v,dir) avec v=u ProdVect dir )
    Node v( u( 1 )*dir( 2 ) - u( 2 )*dir( 1 ),
            u( 2 )*dir( 0 ) - u( 0 )*dir( 2 ), //u( 0 )*dir( 2 ) - u( 2 )*dir( 0 ),
            u( 0 )*dir( 1 ) - u( 1 )*dir( 0 ) );

    double norm_u = std::sqrt( u( 0 )*u( 0 )+u( 1 )*u( 1 )+u( 2 )*u( 2 ) );
    double norm_v = std::sqrt( v( 0 )*v( 0 )+v( 1 )*v( 1 )+v( 2 )*v( 2 ) );
    u( 0 )/=norm_u;
    u( 1 )/=norm_u;
    u( 2 )/=norm_u;
    v( 0 )/=norm_v;
    v( 1 )/=norm_v;
    v( 2 )/=norm_v;

    Node D( dir( 0 ),dir( 1 ),dir( 2 ) );

    /*std::cout << "dir : " << D(0) << " " << D(1) << " " << D(2) << "\n"
              << "u : " << u(0) << " " << u(1) << " " << u(2) << "\n"
              << "v : " << v(0) << " " << v(1) << " " << v(2) << "\n";*/

    return boost::make_tuple( D,u,v );

}




































#define GEOTOOL_FOR_COMP2(r, state)                                     \
        BOOST_PP_NOT_EQUAL( BOOST_PP_TUPLE_ELEM(4, 0, state),           \
                            BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(4, 1, state)) \
                            )                                           \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_FOR_INCR2(r, state)                         \
        (                                                   \
         BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(4, 0, state)),    \
         BOOST_PP_TUPLE_ELEM(4, 1, state),                  \
         BOOST_PP_TUPLE_ELEM(4, 2, state),                  \
         BOOST_PP_TUPLE_ELEM(4, 3, state) )                 \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_FOR_MARKER_POINT_MACRO2(r, state)                       \
        __listMarker.push_back(boost::make_tuple(this->shape(),this->name(), \
                                                 GEOTOOL_MARKER_POINT_MARKVALUE( BOOST_PP_CAT(GEOTOOL_MARKER_POINT_, \
                                                                                             GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(4,3,state) )), \
                                                                                BOOST_PP_TUPLE_ELEM(4, 2, state), \
                                                                                BOOST_PP_TUPLE_ELEM(4, 0, state) ) \
                                                 )                      \
                               );                                       \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_FOR_MARKER_LINE_MACRO2(r, state)                        \
        __listMarker.push_back(boost::make_tuple(this->shape(),this->name(), \
                                                 GEOTOOL_MARKER_LINE_MARKVALUE( BOOST_PP_CAT(GEOTOOL_MARKER_LINE_, \
                                                                                             GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(4,3,state) )), \
                                                                                BOOST_PP_TUPLE_ELEM(4, 2, state), \
                                                                                BOOST_PP_TUPLE_ELEM(4, 0, state) ) \
                                                 )                      \
                               );                                       \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#if 1
#define GEOTOOL_FOR_MARKER_SURFACE_MACRO2(r, state)                     \
        __listMarker.push_back(boost::make_tuple(this->shape(),this->name(), \
                                                 GEOTOOL_MARKER_SURFACE_MARKVALUE( BOOST_PP_CAT(GEOTOOL_MARKER_SURFACE_, \
                                                                                                GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(4,3,state) )), \
                                                                                   BOOST_PP_TUPLE_ELEM(4, 2, state), \
                                                                                   BOOST_PP_TUPLE_ELEM(4, 0, state) ) \
                                                 )                      \
                               );                                       \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#else
#define GEOTOOL_FOR_MARKER_SURFACE_MACRO2(r, state)                     \
        __listMarker.push_back(boost::make_tuple(this->shape(),this->name(), \
                                                 GEOTOOL_MARKER_SURFACE_MARKVALUE( BOOST_PP_CAT(GEOTOOL_MARKER_SURFACE_, \
                                                                                                BOOST_PP_IF(BOOST_PP_GREATER(GEOTOOL_SHAPE_NBSURFACE(BOOST_PP_TUPLE_ELEM(4,3,state)) ,0), \
                                                                                                GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(4,3,state) ), \
                                                                                                DEFAULT)), \
                                                                                   BOOST_PP_TUPLE_ELEM(4, 2, state), \
                                                                                   BOOST_PP_TUPLE_ELEM(4, 0, state) ) \
                                                 )                      \
                               );                                       \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/

#endif


#if 0
#define GEOTOOL_FOR_MARKER_VOLUME_MACRO2(r, state)                      \
        __listMarker.push_back(boost::make_tuple(this->shape(),this->name(), \
                                                 GEOTOOL_MARKER_VOLUME_MARKVALUE( BOOST_PP_CAT(GEOTOOL_MARKER_VOLUME_, \
                                                                                               BOOST_PP_IF(BOOST_PP_GREATER(GEOTOOL_SHAPE_NBVOLUME(BOOST_PP_TUPLE_ELEM(4,3,state)) ,0), \
                                                                                                           GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(4,3,state)), \
                                                                                                           DEFAULT)), \
                                                                                  BOOST_PP_TUPLE_ELEM(4, 2, state), \
                                                                                  BOOST_PP_TUPLE_ELEM(4, 0, state) ) \
                                                 )                      \
                               );                                       \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#else
#define GEOTOOL_FOR_MARKER_VOLUME_MACRO2(r, state)                      \
        __listMarker.push_back(boost::make_tuple(this->shape(),this->name(), \
                                                 GEOTOOL_MARKER_VOLUME_MARKVALUE( BOOST_PP_CAT(GEOTOOL_MARKER_VOLUME_, \
                                                                                               GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(4,3,state) )), \
                                                                                  BOOST_PP_TUPLE_ELEM(4, 2, state), \
                                                                                  BOOST_PP_TUPLE_ELEM(4, 0, state) ) \
                                                 )                      \
                               );                                       \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#endif


#define GEOTOOL_FOR_MARKER_POINT_MACRO(r, state)                         \
        if (BOOST_PP_CAT(marker,                                        \
                         BOOST_PP_ADD(BOOST_PP_TUPLE_ELEM(3, 0, state),	\
                                      1 ) ) )                           \
            {                                                           \
                BOOST_PP_FOR( (0,                                       \
                               BOOST_PP_SUB(GEOTOOL_MARKER_POINT_NBMARK(BOOST_PP_CAT(GEOTOOL_MARKER_POINT_, \
                                                                                    GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(3,2,state))), \
                                                                       BOOST_PP_TUPLE_ELEM(3, 0, state) ),1), \
                               BOOST_PP_TUPLE_ELEM(3, 0, state),		\
                               BOOST_PP_TUPLE_ELEM(3, 2, state)			\
                               ),                                       \
                              GEOTOOL_FOR_COMP2, GEOTOOL_FOR_INCR2, GEOTOOL_FOR_MARKER_POINT_MACRO2) \
                    }                                                   \

/**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_FOR_MARKER_LINE_MACRO(r, state)                         \
        if (BOOST_PP_CAT(marker,                                        \
                         BOOST_PP_ADD(BOOST_PP_TUPLE_ELEM(3, 0, state),	\
                                      1 ) ) )                           \
            {                                                           \
                BOOST_PP_FOR( (0,                                       \
                               BOOST_PP_SUB(GEOTOOL_MARKER_LINE_NBMARK(BOOST_PP_CAT(GEOTOOL_MARKER_LINE_, \
                                                                                    GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(3,2,state))), \
                                                                       BOOST_PP_TUPLE_ELEM(3, 0, state) ),1), \
                               BOOST_PP_TUPLE_ELEM(3, 0, state),		\
                               BOOST_PP_TUPLE_ELEM(3, 2, state)			\
                               ),                                       \
                              GEOTOOL_FOR_COMP2, GEOTOOL_FOR_INCR2, GEOTOOL_FOR_MARKER_LINE_MACRO2) \
                    }                                                   \

/**/
/*_________________________________________________*/
/*                                                 */
/**/
#if 1
#define GEOTOOL_FOR_MARKER_SURFACE_MACRO(r, state)                      \
        if (BOOST_PP_CAT(marker,                                        \
                         BOOST_PP_ADD(BOOST_PP_TUPLE_ELEM(3, 0, state),	\
                                      1 ) ) )                           \
            {                                                           \
                BOOST_PP_FOR( (0,                                       \
                               BOOST_PP_SUB(GEOTOOL_MARKER_SURFACE_NBMARK(BOOST_PP_CAT(GEOTOOL_MARKER_SURFACE_, \
                                                                                       GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(3,2,state)) ), \
                                                                          BOOST_PP_TUPLE_ELEM(3, 0, state) ),1), \
                               BOOST_PP_TUPLE_ELEM(3, 0, state),		\
                               BOOST_PP_TUPLE_ELEM(3, 2, state)			\
                               ),                                       \
                              GEOTOOL_FOR_COMP2, GEOTOOL_FOR_INCR2, GEOTOOL_FOR_MARKER_SURFACE_MACRO2) \
                    }                                                   \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#else
#define GEOTOOL_FOR_MARKER_SURFACE_MACRO(r, state)                      \
        if (BOOST_PP_CAT(marker,                                        \
                         BOOST_PP_ADD(BOOST_PP_TUPLE_ELEM(3, 0, state),	\
                                      1 ) ) )                           \
            {                                                           \
                BOOST_PP_FOR( (0,                                       \
                               BOOST_PP_SUB(GEOTOOL_MARKER_SURFACE_NBMARK(BOOST_PP_CAT(GEOTOOL_MARKER_SURFACE_, \
                                                                                       BOOST_PP_IF(BOOST_PP_GREATER(GEOTOOL_SHAPE_NBSURFACE(BOOST_PP_TUPLE_ELEM(3,2,state)),0), \
                                                                                                   GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(3,2,state)), \
                                                                                                   DEFAULT )), \
                                                                          BOOST_PP_TUPLE_ELEM(3, 0, state) ),1) , \
                               BOOST_PP_TUPLE_ELEM(3, 0, state),		\
                               BOOST_PP_TUPLE_ELEM(3, 2, state)			\
                               ),                                       \
                              GEOTOOL_FOR_COMP2, GEOTOOL_FOR_INCR2, GEOTOOL_FOR_MARKER_SURFACE_MACRO2) \
                   }                                                  \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/

#endif

#if 0
#define GEOTOOL_FOR_MARKER_VOLUME_MACRO(r, state)                       \
        if (BOOST_PP_CAT(marker,                                        \
                         BOOST_PP_ADD(BOOST_PP_TUPLE_ELEM(3, 0, state),	\
                                      1 ) ) )                           \
            {                                                           \
                BOOST_PP_FOR( (0,                                       \
                               BOOST_PP_SUB(GEOTOOL_MARKER_VOLUME_NBMARK(BOOST_PP_CAT(GEOTOOL_MARKER_VOLUME_, \
                                                                                      BOOST_PP_IF(BOOST_PP_GREATER(GEOTOOL_SHAPE_NBVOLUME(BOOST_PP_TUPLE_ELEM(3,2,state)) ,0), \
                                                                                                  GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(3,2,state)),\
                                                                                                  DEFAULT )), \
                                                                          BOOST_PP_TUPLE_ELEM(3, 0, state) ),1), \
                               BOOST_PP_TUPLE_ELEM(3, 0, state),		\
                               BOOST_PP_TUPLE_ELEM(3, 2, state)			\
                               ),                                       \
                              GEOTOOL_FOR_COMP2, GEOTOOL_FOR_INCR2, GEOTOOL_FOR_MARKER_VOLUME_MACRO2) \
                    }                                                   \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#else
#define GEOTOOL_FOR_MARKER_VOLUME_MACRO(r, state)                       \
        if (BOOST_PP_CAT(marker,                                        \
                         BOOST_PP_ADD(BOOST_PP_TUPLE_ELEM(3, 0, state),	\
                                      1 ) ) )                           \
            {                                                           \
                BOOST_PP_FOR( (0,                                       \
                               BOOST_PP_SUB(GEOTOOL_MARKER_VOLUME_NBMARK(BOOST_PP_CAT(GEOTOOL_MARKER_VOLUME_, \
                                                                                      GEOTOOL_SHAPE_NAME_MACRO(BOOST_PP_TUPLE_ELEM(3,2,state)) ), \
                                                                         BOOST_PP_TUPLE_ELEM(3, 0, state) ),1), \
                               BOOST_PP_TUPLE_ELEM(3, 0, state),		\
                               BOOST_PP_TUPLE_ELEM(3, 2, state)			\
                               ),                                       \
                              GEOTOOL_FOR_COMP2, GEOTOOL_FOR_INCR2, GEOTOOL_FOR_MARKER_VOLUME_MACRO2) \
                    }                                                   \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/

#endif


#define GEOTOOL_FOR_COMP1(r, state)                                     \
        BOOST_PP_NOT_EQUAL( BOOST_PP_TUPLE_ELEM(3, 0, state),           \
                            BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(3, 1, state)) \
                            )                                           \
        /**/
/*_________________________________________________*/

#define GEOTOOL_FOR_INCR1(r, state)                         \
        (                                                   \
         BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(3, 0, state)),	\
         BOOST_PP_TUPLE_ELEM(3, 1, state),                  \
         BOOST_PP_TUPLE_ELEM(3, 2, state) )                 \
        /**/
/*_________________________________________________*/
/*_________________________________________________*/
#define GEOTOOL_SHAPE_PARAM(r, state)                                   \
    M_param[BOOST_PP_TUPLE_ELEM(2,0,state)] = BOOST_PP_CAT( __param,    \
                                                            BOOST_PP_TUPLE_ELEM(2,0,state) ); \
    /**/


/*_________________________________________________*/
/*_________________________________________________*/


#define GEOTOOL_SHAPE_FOR_PARAM_SIGNATURE_NO_DEFAULT_ARG(r, state)      \
    Node BOOST_PP_CAT( __param, BOOST_PP_TUPLE_ELEM(2,0,state) ) BOOST_PP_COMMA() \
        /**/
/*_________________________________________________*/
/*                                                 */
/**/
#define GEOTOOL_SHAPE_PARAM_SIGNATURE_NO_DEFAULT_ARG(state)             \
    BOOST_PP_FOR( (0, BOOST_PP_SUB(GEOTOOL_SHAPE_NBPARAM(BOOST_PP_TUPLE_ELEM(2,0,state)),1) ), \
                  GEOTOOL_FOR_COMP,                                     \
                  GEOTOOL_FOR_INCR,                                     \
                  GEOTOOL_SHAPE_FOR_PARAM_SIGNATURE_NO_DEFAULT_ARG)     \
                                                                          /**/



#define GEOTOOL_SHAPE_CLASS_IMPL(r,state)                               \
    GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state))::          \
    GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state))(double __meshSize, \
                                                             std::string __name, \
                                                             GEOTOOL_SHAPE_PARAM_SIGNATURE_NO_DEFAULT_ARG(state) \
                                                             uint16_type type ) /*Ne sert a rien, juste a cause de la virgule au dessus)*/ \
    :                                                                   \
    GeoGMSHTool( GEOTOOL_SHAPE_DIM(BOOST_PP_TUPLE_ELEM(2,0,state)),shape(), __name, __meshSize) \
    /*M_name(__name)*/                                                  \
    {                                                                   \
        M_param.resize( GEOTOOL_SHAPE_NBPARAM(BOOST_PP_TUPLE_ELEM(2,0,state))); \
        BOOST_PP_FOR( (0, BOOST_PP_SUB(GEOTOOL_SHAPE_NBPARAM(BOOST_PP_TUPLE_ELEM(2,0,state)),1) ), \
                      GEOTOOL_FOR_COMP,                                 \
                      GEOTOOL_FOR_INCR,                                 \
                      GEOTOOL_SHAPE_PARAM);                             \
                                                                        \
        this->initData(shape(),                                         \
                       __name,                                          \
                       __meshSize,                                      \
                       M_param,                                         \
                       GEOTOOL_SHAPE_DIM(BOOST_PP_TUPLE_ELEM(2,0,state)), \
                       1,                                               \
                       GEOTOOL_SHAPE_NBSURFACE(BOOST_PP_TUPLE_ELEM(2,0,state)), \
                       GEOTOOL_SHAPE_NBVOLUME(BOOST_PP_TUPLE_ELEM(2,0,state))); \
    }                                                                   \
                                                                        \
    void                                                                \
    GEOTOOL_SHAPE_NAME_CLASS(BOOST_PP_TUPLE_ELEM(2,0,state))            \
        ::setMarkerImpl( std::string type, std::string name, std::vector<bool> const& markers ) \
    {                                                                   \
        bool marker1=markers[0];                                        \
        bool marker2=markers[1];                                        \
        bool marker3=markers[2];                                        \
        bool marker4=markers[3];                                        \
        bool marker5=markers[4];                                        \
        bool marker6=markers[5];                                        \
        bool marker7=markers[6];                                        \
        bool marker8=markers[7];                                        \
        bool marker9=markers[8];                                        \
        bool marker10=markers[9];                                       \
        bool marker11=markers[10];                                      \
        bool marker12=markers[11];                                      \
                                                                        \
        std::vector<marker_base_type> __listMarker = (*(M_markShape))[type][name]; \
                                                                        \
    if (type=="point")                                                  \
        {                                                               \
            BOOST_PP_FOR( (0, BOOST_PP_SUB(BOOST_PP_ARRAY_SIZE( BOOST_PP_CAT(GEOTOOL_MARKER_POINT_, \
                                                                             GEOTOOL_SHAPE_NAME_MACRO( BOOST_PP_TUPLE_ELEM(2,0,state)))), \
                                           1), BOOST_PP_TUPLE_ELEM(2,0,state)), \
                          GEOTOOL_FOR_COMP1,                            \
                      GEOTOOL_FOR_INCR1,                                \
                      GEOTOOL_FOR_MARKER_POINT_MACRO)                   \
            }                                                           \
    else if (type=="line")                                              \
    {                                                                   \
        BOOST_PP_FOR( (0, BOOST_PP_SUB(BOOST_PP_ARRAY_SIZE( BOOST_PP_CAT(GEOTOOL_MARKER_LINE_, \
                                                                         GEOTOOL_SHAPE_NAME_MACRO( BOOST_PP_TUPLE_ELEM(2,0,state)))), \
                                       1), BOOST_PP_TUPLE_ELEM(2,0,state)), \
                      GEOTOOL_FOR_COMP1,                                \
                      GEOTOOL_FOR_INCR1,                                \
                      GEOTOOL_FOR_MARKER_LINE_MACRO)                    \
            }                                                           \
    else if (type=="surface")                                           \
    {                                                                   \
        BOOST_PP_IF(BOOST_PP_NOT_EQUAL(GEOTOOL_SHAPE_NBSURFACE(BOOST_PP_TUPLE_ELEM(2,0,state)),0), \
                    BOOST_PP_FOR( (0, BOOST_PP_SUB(BOOST_PP_ARRAY_SIZE( BOOST_PP_CAT(GEOTOOL_MARKER_SURFACE_, \
                                                                                     GEOTOOL_SHAPE_NAME_MACRO( BOOST_PP_TUPLE_ELEM(2,0,state)))), \
                                                   1), BOOST_PP_TUPLE_ELEM(2,0,state)), \
                                  GEOTOOL_FOR_COMP1,                    \
                                  GEOTOOL_FOR_INCR1,                    \
                                  GEOTOOL_FOR_MARKER_SURFACE_MACRO),    \
                    )                                                   \
            }                                                           \
    else if (type=="volume")                                            \
    {                                                                   \
        BOOST_PP_IF(BOOST_PP_NOT_EQUAL(GEOTOOL_SHAPE_NBVOLUME(BOOST_PP_TUPLE_ELEM(2,0,state)),0), \
                    BOOST_PP_FOR( (0, BOOST_PP_SUB(BOOST_PP_ARRAY_SIZE( BOOST_PP_CAT(GEOTOOL_MARKER_VOLUME_, \
                                                                                     GEOTOOL_SHAPE_NAME_MACRO( BOOST_PP_TUPLE_ELEM(2,0,state)))), \
                                                   1), BOOST_PP_TUPLE_ELEM(2,0,state)), \
                                  GEOTOOL_FOR_COMP1,                    \
                                  GEOTOOL_FOR_INCR1,                    \
                                  GEOTOOL_FOR_MARKER_VOLUME_MACRO),     \
                    )                                                   \
            }                                                           \
                                                                        \
    (*(M_markShape))[type][name] = __listMarker;                        \
                                                                        \
    }                                                                   \

    /**/

//creation des classes representants les objets geotool
BOOST_PP_FOR( ( 0, BOOST_PP_SUB( BOOST_PP_ARRAY_SIZE( GEOTOOL_SHAPE ),1 ) ),
              GEOTOOL_FOR_COMP,
              GEOTOOL_FOR_INCR,
              GEOTOOL_SHAPE_CLASS_IMPL )



















} //namespace GeoTool

} //namespace Feel

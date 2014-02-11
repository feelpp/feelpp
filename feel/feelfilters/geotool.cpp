
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
#if defined( FEELPP_HAS_GMSH_H )
#include <GmshConfig.h>
#endif

namespace Feel
{

namespace GeoTool
{

GeoGMSHTool::GeoGMSHTool( uint16_type __dim,
                          std::string __shape,
                          std::string __name,
                          double __meshSize )
    :
    M_dim( __dim ),
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
    M_geoIsDefineByUser( m.M_geoIsDefineByUser )
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
                    std::list<marker_base_type>::const_iterator itLRef=m.markerListIndiceBegin( itMarkType->first,
                            itMarkName->first );
                    std::list<marker_base_type>::const_iterator itLRef_end=m.markerListIndiceEnd( itMarkType->first,
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

    //data memory ( type->shape->name )
    std::vector<std::map<std::string,std::map<std::string, std::map<uint16_type,uint16_type> > > > __dataMemGlob( 6 );
    std::vector<std::map<std::string,std::map<std::string, std::map<uint16_type,bool> > > > __dataMemGlobSurf1( 2 );
    std::vector<std::map<std::string,std::map<std::string, std::map<uint16_type,std::string> > > > __dataMemGlobSurf2( 2 );
    std::vector<std::map<std::string,std::map<std::string, std::map<uint16_type,bool> > > > __dataMemGlobIsRuled( 1 );
    std::vector<std::map<std::string,std::map<std::string, std::map<uint16_type,std::list<uint16_type> > > > > __dataMemGlobPtsInSurf( 1 );
    // type -> name -> num surfLoop -> list de surfLoop
    std::map<std::string,std::map<std::string, std::map<int,std::list<int> > > > __dataMemGlobSurfaceLoop;
    __dataMemGlobSurfaceLoop.clear();

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

        *M_ostr << "h=" << boost::get<2>( *itList ) << ";\n";
        //data memory
        vec_map_data_ptrtype __dataMem( new vec_map_data_type( 8 ) );
        vec_map_data_surf1_ptrtype __dataMemSurf1( new vec_map_data_surf1_type( 2 ) );
        vec_map_data_surf2_ptrtype __dataMemSurf2( new vec_map_data_surf2_type( 2 ) );
        vec_map_data_surf1_ptrtype __dataMemIsRuled( new vec_map_data_surf1_type( 1 ) );
        vec_map_data_ptsinsurf_ptrtype __dataMemPtsInSurf( new vec_map_data_ptsinsurf_type( 1 ) );

        map_surfaceLoop_type __dataMemLocSurfaceLoop;
        __dataMemLocSurfaceLoop.clear();

        GeoGMSHTool_ptrtype __geoTool( new GeoGMSHTool( this->dim() ) );
        __geoTool->updateData( *this );
        __geoTool->cleanOstr();

        GeoTool::data_geo_ptrtype __data_geoTool( new GeoTool::data_geo_type( boost::make_tuple( __geoTool,
                __dataMem,
                Qshape,//itShape->first,
                Qname,//boost::get<0>(*itName),
                __dataMemSurf1,
                __dataMemSurf2,
                __dataMemIsRuled,
                __dataMemPtsInSurf,
                __dataMemLocSurfaceLoop
                                                                                               ) ) );

        // generate the code for the geometry
        run( __data_geoTool );


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

        // get infos
        this->updateData( *boost::get<0>( *__data_geoTool ) );
        this->updateOstr( boost::get<0>( *__data_geoTool )->M_ostr->str() );

    }


    //std::map<int,int> mapSurfaceRenumbering;


    //Write the planes surfaces
    //Fait ici a cause des opertateurs (+,-)
    surface_name_const_iterator_type itSurf = this->M_surfaceList->begin();
    surface_name_const_iterator_type itSurf_end = this->M_surfaceList->end();
    std::ostringstream __surface_str;
    //counter of surface
    uint16_type __surfnumber=1;//referencier une liste des surf dans les writePlaneSurface
    // build map between local id surf with name -> global id surf
    std::map<std::string,std::map<std::string, std::map<uint16_type,uint16_type> > > __dataSurfacePost;
    // counter of local surf
    std::map<std::string,std::map<std::string, uint16_type> > __dataSurfacePostCpt;

    for ( ; itSurf != itSurf_end; ++itSurf )
    {
        surface_type_const_iterator_type itSurf2 = itSurf->begin();

        if (itSurf2->get<2>().second==0) continue; // surface useless


        // utile pour les physical marker
        //if ((__dataSurfacePostCpt[boost::get<0>(*itSurf2)][boost::get<1>(*itSurf2)]).empty())
        if ( __dataSurfacePostCpt.find( boost::get<0>( *itSurf2 ) )==__dataSurfacePostCpt.end() )
            __dataSurfacePostCpt[boost::get<0>( *itSurf2 )][boost::get<1>( *itSurf2 )]=1;

        else if ( __dataSurfacePostCpt[boost::get<0>( *itSurf2 )].find( boost::get<1>( *itSurf2 ) ) == __dataSurfacePostCpt[boost::get<0>( *itSurf2 )].end() )
            __dataSurfacePostCpt[boost::get<0>( *itSurf2 )][boost::get<1>( *itSurf2 )]=1;

        else
            ++__dataSurfacePostCpt[itSurf2->get<0>()][itSurf2->get<1>()];

        __dataSurfacePost[itSurf2->get<0>()][itSurf2->get<1>()]
        [__dataSurfacePostCpt[itSurf2->get<0>()][itSurf2->get<1>()]]=__surfnumber;

        // si la surface est issue d'un extrude => stoker dans un tab gmsh
        if ( true/* ! __dataMemGlobSurf1[0][itSurf2->get<0>()][itSurf2->get<1>()][__surfnumber]*/ )
        {
            //Attention : On fait a cause des op - : sinon les markers surfaces sont incorrectes(l'idee est de marquer la 1ere sous-surface)
            //__dataMemGlob[3][boost::get<0>(*itSurf2)][boost::get<1>(*itSurf2)][__surfnumber]=__surfnumber;
            //std::cout << "\n taille " << __dataMemGlobIsRuled[0][itSurf2->get<0>()][itSurf2->get<1>()].size() << " surnumber "<< __surfnumber << std::endl;

            if ( ! __dataMemGlobIsRuled[0][itSurf2->get<0>()][itSurf2->get<1>()][itSurf2->get<2>().first/*__surfnumber*/] )
                __surface_str << "Plane Surface(" << __surfnumber << ") = {" ;

            else
                __surface_str << "Ruled Surface(" << __surfnumber << ") = {" ;

            //std::cout << "name " << itSurf2->get<1>() << " map realsurf " << itSurf2->get<2>().first << " with locsurf  " <<  __dataSurfacePostCpt[itSurf2->get<0>()][itSurf2->get<1>()] << std::endl;
            //mapSurfaceRenumbering[itSurf2->get<2>().first] = __dataSurfacePostCpt[itSurf2->get<0>()][itSurf2->get<1>()];//__surfnumber;

            surface_type_const_iterator_type itSurf2_end = --itSurf->end();
            for ( ; itSurf2 != itSurf2_end; ++itSurf2 )
            {
                //std::cout << "\n num Glob SURF " << itSurf2->get<2>().first << " __surfnumber " << __surfnumber << std::endl;
                __surface_str << itSurf2->get<2>().second << ",";
            }

            // ATTTENTION : FAIT ICI CAR 1 SEUL!!!!!!!!!!!!!!!!!!
            //mapSurfaceRenumbering[itSurf2->get<2>().first] = __surfnumber;
            //std::cout << "\n num Glob SURF " << itSurf2->get<2>().first << " "<< itSurf2->get<2>().second << " __surfnumber " << __surfnumber << std::endl;

            __surface_str << itSurf2->get<2>().second;

            __surface_str << "};\n";

            //maybe add more pts in surface
            itSurf2 = itSurf->begin();
            auto ptInSurf_it = __dataMemGlobPtsInSurf[0][boost::get<0>( *itSurf2 )][boost::get<1>( *itSurf2 )][itSurf2->get<2>().first].begin();
            auto ptInSurf_en = __dataMemGlobPtsInSurf[0][boost::get<0>( *itSurf2 )][boost::get<1>( *itSurf2 )][itSurf2->get<2>().first].end();
            for ( ; ptInSurf_it != ptInSurf_en ; ++ptInSurf_it )
            {
                __surface_str << "Point{" << *ptInSurf_it << "} In Surface{" << __surfnumber << "};\n";
            }

        }

        ++__surfnumber;

    }

    this->updateOstr( __surface_str.str() );

    //---------------------------------------------------------------------------------//
    //Write the extrude surfaces
    this->updateOstr( M_ostrExtrude->str() );

    //---------------------------------------------------------------------------------//
    //Write the surfaces loops

    int __nSurfaceLoop=1;
    std::ostringstream __ostrSurfaceLoop;
    auto surfaceLoop_it = this->M_surfaceLoopList->begin();
    auto surfaceLoop_en = this->M_surfaceLoopList->end();

    for ( ; surfaceLoop_it!=surfaceLoop_en ; ++surfaceLoop_it )
    {
        auto surfaceLoop2_it =surfaceLoop_it->begin();
        auto surfaceLoop2_en =surfaceLoop_it->end();

        for ( ; surfaceLoop2_it!=surfaceLoop2_en ; ++surfaceLoop2_it )
        {
            /*std::cout << "\n SurfaceLoop shape : " << surfaceLoop2_it->get<0>()
                      << " name " << surfaceLoop2_it->get<1>()
                      << " size " << surfaceLoop2_it->get<2>().size()
                      << std::endl;*/

            auto surfaceLoop3_it =surfaceLoop2_it->get<2>().begin();
            auto surfaceLoop3_en =surfaceLoop2_it->get<2>().end();

            for ( ; surfaceLoop3_it!=surfaceLoop3_en ; ++surfaceLoop3_it )
            {
                //numLoc = surfaceLoop3_it->first
                __ostrSurfaceLoop << "Surface Loop(" << __nSurfaceLoop
                                  << ") = {" ;
                auto surfaceLoop4_it =surfaceLoop3_it->second.begin();
                auto surfaceLoop4_en =--surfaceLoop3_it->second.end();

                for ( ; surfaceLoop4_it!=surfaceLoop4_en ; ++surfaceLoop4_it )
                {
                    //std::cout << "\n HOLA " << *surfaceLoop4_it << " map " << mapSurfaceRenumbering[*surfaceLoop4_it] <<std::endl;
                    __ostrSurfaceLoop << __dataSurfacePost[surfaceLoop2_it->get<0>()][surfaceLoop2_it->get<1>()][ /*mapSurfaceRenumbering[*/ *surfaceLoop4_it/*]*/ ] << ",";
                }

                __ostrSurfaceLoop << __dataSurfacePost[surfaceLoop2_it->get<0>()][surfaceLoop2_it->get<1>()][ /*mapSurfaceRenumbering[*/ *surfaceLoop4_it/*]*/ ] <<"};\n";
                ++__nSurfaceLoop;

            } // surfaceLoop
        } // name
    } // shape

    this->updateOstr( __ostrSurfaceLoop.str() );

    //---------------------------------------------------------------------------------//
    //Write the volumes
    //Fait ici a cause des opertateurs (+,-)
    volume_name_const_iterator_type itVol = this->M_volumeList->begin();
    volume_name_const_iterator_type itVol_end = this->M_volumeList->end();
    std::ostringstream __volume_str;
    //counter of volume
    uint16_type __volnumber=1;
    std::map<std::string,std::map<std::string, std::map<uint16_type,uint16_type> > > __dataVolumePost;
    std::map<std::string,std::map<std::string, uint16_type> > __dataVolumePostCpt;

    for ( ; itVol != itVol_end; ++itVol )
    {
        volume_type_const_iterator_type itVol2 = itVol->begin();
        volume_type_const_iterator_type itVol2_end = --itVol->end();

        // utile pour les physical marker
        //if ((__dataSurfacePostCpt[boost::get<0>(*itSurf2)][boost::get<1>(*itSurf2)]).empty())
        if ( __dataVolumePostCpt.find( boost::get<0>( *itVol2 ) )==__dataVolumePostCpt.end() )
            __dataVolumePostCpt[boost::get<0>( *itVol2 )][boost::get<1>( *itVol2 )]=1;

        else if ( __dataVolumePostCpt[boost::get<0>( *itVol2 )].find( boost::get<1>( *itVol2 ) ) == __dataVolumePostCpt[boost::get<0>( *itVol2 )].end() )
            __dataVolumePostCpt[boost::get<0>( *itVol2 )][boost::get<1>( *itVol2 )]=1;

        else
            ++__dataVolumePostCpt[boost::get<0>( *itVol2 )][boost::get<1>( *itVol2 )];

        __dataVolumePost[boost::get<0>( *itVol2 )]
        [boost::get<1>( *itVol2 )]
        [__dataVolumePostCpt[boost::get<0>( *itVol2 )][boost::get<1>( *itVol2 )]]=__volnumber;


        // si le volume est issue d'un extrude => stocker dans un tab gmsh => pas d'affichage
        if ( true /*! __dataMemGlobSurf1[1][boost::get<0>( *itVol2 )][boost::get<1>( *itVol2 )][__volnumber]*/ )
        {
            __volume_str << "Volume(" << __volnumber << ") = {" ;

            for ( ; itVol2 != itVol2_end; ++itVol2 )
            {
                __volume_str << itVol2->get<2>().second << ",";
            }

            __volume_str << itVol2->get<2>().second;
            __volume_str << "};\n";
        }

        ++__volnumber;

    }


    this->updateOstr( __volume_str.str() );


    // generate the code for the marker
    marker_type_const_iterator_type itMarkType= ( *( M_markShape ) ).begin();
    marker_type_const_iterator_type itMarkType_end=( *( M_markShape ) ).end();

    while ( itMarkType!=itMarkType_end )
    {
        marker_markerName_const_iterator_type itMarkName = ( *( M_markShape ) )[itMarkType->first].begin();
        marker_markerName_const_iterator_type itMarkName_end=( *( M_markShape ) )[itMarkType->first].end();

        while ( itMarkName!=itMarkName_end )
        {
            if ( itMarkType->first=="point" )
            {
                //on cree un nouvelle list dont on enleve les doublons(aux cas où!)
                std::list<marker_base_type> newListMark;
                std::list<marker_base_type>::const_iterator itMark = itMarkName->second.begin();
                std::list<marker_base_type>::const_iterator itMark_end = itMarkName->second.end();

                for ( ; itMark!=itMark_end; ++itMark )
                {
                    auto value = __dataMemGlob[0][boost::get<0>( *itMark )][boost::get<1>( *itMark )][boost::get<2>( *itMark )];
                    std::list<marker_base_type>::const_iterator itnewMark = newListMark.begin();
                    std::list<marker_base_type>::const_iterator itnewMark_end = newListMark.end();
                    bool find=false;

                    while ( itnewMark != itnewMark_end && !find )
                    {
                        auto value2 = __dataMemGlob[0][boost::get<0>( *itnewMark )][boost::get<1>( *itnewMark )][boost::get<2>( *itnewMark )];

                        if ( value == value2 ) find=true;

                        ++itnewMark;
                    }

                    if ( !find ) newListMark.push_back( *itMark );
                }


                *M_ostr << "Physical Point(\"" << itMarkName->first << "\") = {";

                ///*std::list<marker_base_type>::const_iterator*/ itMark = itMarkName->second.begin();
                ///*std::list<marker_base_type>::const_iterator*/ itMark_end = --itMarkName->second.end();
                auto itMarkTTT=newListMark.begin();
                auto itMarkTTT_end=--newListMark.end();

                while ( itMarkTTT!=itMarkTTT_end )
                {
                    *M_ostr << __dataMemGlob[0][boost::get<0>( *itMarkTTT )][boost::get<1>( *itMarkTTT )][boost::get<2>( *itMarkTTT )]
                             <<",";
                    ++itMarkTTT;
                }

                *M_ostr << __dataMemGlob[0][boost::get<0>( *itMarkTTT )][boost::get<1>( *itMarkTTT )][boost::get<2>( *itMarkTTT )] << "};\n";
            }

            else if ( itMarkType->first=="line" )
            {
                //on cree un nouvelle list dont on enleve les doublons(aux cas où!)
                std::list<marker_base_type> newListMark;
                std::list<marker_base_type>::const_iterator itMark = itMarkName->second.begin();
                std::list<marker_base_type>::const_iterator itMark_end = itMarkName->second.end();

                for ( ; itMark!=itMark_end; ++itMark )
                {
                    auto value = __dataMemGlob[1][boost::get<0>( *itMark )][boost::get<1>( *itMark )][boost::get<2>( *itMark )];
                    std::list<marker_base_type>::const_iterator itnewMark = newListMark.begin();
                    std::list<marker_base_type>::const_iterator itnewMark_end = newListMark.end();
                    bool find=false;

                    while ( itnewMark != itnewMark_end && !find )
                    {
                        auto value2 = __dataMemGlob[1][boost::get<0>( *itnewMark )][boost::get<1>( *itnewMark )][boost::get<2>( *itnewMark )];

                        if ( value == value2 ) find=true;

                        ++itnewMark;
                    }

                    if ( !find ) newListMark.push_back( *itMark );
                }


                *M_ostr << "Physical Line(\"" << itMarkName->first << "\") = {";

                ///*std::list<marker_base_type>::const_iterator*/ itMark = itMarkName->second.begin();
                ///*std::list<marker_base_type>::const_iterator*/ itMark_end = --itMarkName->second.end();
                auto itMarkTTT=newListMark.begin();
                auto itMarkTTT_end=--newListMark.end();

                while ( itMarkTTT!=itMarkTTT_end )
                {
                    *M_ostr << __dataMemGlob[1][boost::get<0>( *itMarkTTT )][boost::get<1>( *itMarkTTT )][boost::get<2>( *itMarkTTT )]
                             <<",";
                    ++itMarkTTT;
                }

                *M_ostr << __dataMemGlob[1][boost::get<0>( *itMarkTTT )][boost::get<1>( *itMarkTTT )][boost::get<2>( *itMarkTTT )] << "};\n";
            }

            else if ( itMarkType->first=="surface" )
            {

                //on cree un nouvelle list dont on enleve les doublons(aux cas où!)
                std::list<marker_base_type> newListMark;
                std::list<marker_base_type>::const_iterator itMark = itMarkName->second.begin();
                std::list<marker_base_type>::const_iterator itMark_end = itMarkName->second.end();

                for ( ; itMark!=itMark_end; ++itMark )
                {
                    auto value = __dataSurfacePost[boost::get<0>( *itMark )][boost::get<1>( *itMark )][boost::get<2>( *itMark )];
                    std::list<marker_base_type>::const_iterator itnewMark = newListMark.begin();
                    std::list<marker_base_type>::const_iterator itnewMark_end = newListMark.end();
                    bool find=false;

                    while ( itnewMark != itnewMark_end && !find )
                    {
                        auto value2 = __dataSurfacePost[boost::get<0>( *itnewMark )][boost::get<1>( *itnewMark )][boost::get<2>( *itnewMark )];

                        if ( value == value2 ) find=true;

                        ++itnewMark;
                    }

                    if ( !find ) newListMark.push_back( *itMark );
                }

                *M_ostr << "Physical Surface(\"" << itMarkName->first << "\") = {";

                //std::list<marker_base_type>::const_iterator itMark = itMarkName->second.begin();
                //std::list<marker_base_type>::const_iterator itMark_end = --itMarkName->second.end();

                auto itMarkTTT=newListMark.begin();
                auto itMarkTTT_end=--newListMark.end();

                while ( itMarkTTT!=itMarkTTT_end )
                {
                    if ( ! __dataMemGlobSurf1[0][boost::get<0>( *itMarkTTT )][boost::get<1>( *itMarkTTT )][boost::get<2>( *itMarkTTT )] )
                        *M_ostr << /*__dataMemGlob[3]*/__dataSurfacePost[boost::get<0>( *itMarkTTT )][boost::get<1>( *itMarkTTT )][boost::get<2>( *itMarkTTT )]
                                 <<",";

                    else
                        *M_ostr << __dataMemGlobSurf2[0][boost::get<0>( *itMarkTTT )][boost::get<1>( *itMarkTTT )][boost::get<2>( *itMarkTTT )]
                                 <<"["
                                 << /*__dataMemGlob[3]*/__dataSurfacePost[boost::get<0>( *itMarkTTT )][boost::get<1>( *itMarkTTT )][boost::get<2>( *itMarkTTT )]
                                 <<"],";

                    ++itMarkTTT;
                }

                if ( ! __dataMemGlobSurf1[0][boost::get<0>( *itMarkTTT )][boost::get<1>( *itMarkTTT )][boost::get<2>( *itMarkTTT )] )
                    *M_ostr << /*__dataMemGlob[3]*/__dataSurfacePost[boost::get<0>( *itMarkTTT )][boost::get<1>( *itMarkTTT )][boost::get<2>( *itMarkTTT )]
                             << "};\n";

                else
                    *M_ostr << __dataMemGlobSurf2[0][boost::get<0>( *itMarkTTT )][boost::get<1>( *itMarkTTT )][boost::get<2>( *itMarkTTT )]
                             <<"["
                             << /*__dataMemGlob[3]*/__dataSurfacePost[boost::get<0>( *itMarkTTT )][boost::get<1>( *itMarkTTT )][boost::get<2>( *itMarkTTT )]
                             <<"] };\n";
            }

            else if ( itMarkType->first=="volume" )
            {
#if 0
                *M_ostr << "Physical Volume(\"" << itMarkName->first << "\") = {";

                std::list<marker_base_type>::const_iterator itMark = itMarkName->second.begin();
                std::list<marker_base_type>::const_iterator itMark_end = --itMarkName->second.end();

                while ( itMark!=itMark_end )
                {
                    if ( ! __dataMemGlobSurf1[1][boost::get<0>( *itMark )][boost::get<1>( *itMark )][boost::get<2>( *itMark )] )
                        *M_ostr << __dataMemGlob[3][boost::get<0>( *itMark )][boost::get<1>( *itMark )][boost::get<2>( *itMark )]<<",";

                    else
                        *M_ostr << __dataMemGlobSurf2[1][boost::get<0>( *itMark )][boost::get<1>( *itMark )][boost::get<2>( *itMark )]
                                 <<"[1],";

                    ++itMark;
                }

                if ( ! __dataMemGlobSurf1[1][boost::get<0>( *itMark )][boost::get<1>( *itMark )][boost::get<2>( *itMark )] )
                    *M_ostr << __dataMemGlob[3][boost::get<0>( *itMark )][boost::get<1>( *itMark )][boost::get<2>( *itMark )] << "};\n";

                else
                    *M_ostr << __dataMemGlobSurf2[1][boost::get<0>( *itMark )][boost::get<1>( *itMark )][boost::get<2>( *itMark )]
                             <<"[1]};\n";

#else
                //on cree un nouvelle list dont on enleve les doublons(aux cas où!)
                std::list<marker_base_type> newListMark;
                std::list<marker_base_type>::const_iterator itMark = itMarkName->second.begin();
                std::list<marker_base_type>::const_iterator itMark_end = itMarkName->second.end();

                for ( ; itMark!=itMark_end; ++itMark )
                {
                    auto value = __dataVolumePost[boost::get<0>( *itMark )][boost::get<1>( *itMark )][boost::get<2>( *itMark )];
                    std::list<marker_base_type>::const_iterator itnewMark = newListMark.begin();
                    std::list<marker_base_type>::const_iterator itnewMark_end = newListMark.end();
                    bool find=false;

                    while ( itnewMark != itnewMark_end && !find )
                    {
                        auto value2 = __dataVolumePost[boost::get<0>( *itnewMark )][boost::get<1>( *itnewMark )][boost::get<2>( *itnewMark )];

                        if ( value == value2 ) find=true;

                        ++itnewMark;
                    }

                    if ( !find ) newListMark.push_back( *itMark );
                }

                *M_ostr << "Physical Volume(\"" << itMarkName->first << "\") = {";

                auto itMarkTTT=newListMark.begin();
                auto itMarkTTT_end=--newListMark.end();

                while ( itMarkTTT!=itMarkTTT_end )
                {
                    if ( ! __dataMemGlobSurf1[1][boost::get<0>( *itMarkTTT )][boost::get<1>( *itMarkTTT )][boost::get<2>( *itMarkTTT )] )
                        *M_ostr << __dataVolumePost[boost::get<0>( *itMarkTTT )][boost::get<1>( *itMarkTTT )][boost::get<2>( *itMarkTTT )]
                                 <<",";

                    else
                        *M_ostr << __dataMemGlobSurf2[1][boost::get<0>( *itMarkTTT )][boost::get<1>( *itMarkTTT )][boost::get<2>( *itMarkTTT )]
                                 <<"["
                                 << __dataVolumePost[boost::get<0>( *itMarkTTT )][boost::get<1>( *itMarkTTT )][boost::get<2>( *itMarkTTT )]
                                 <<"],";

                    ++itMarkTTT;
                }

                if ( ! __dataMemGlobSurf1[1][boost::get<0>( *itMarkTTT )][boost::get<1>( *itMarkTTT )][boost::get<2>( *itMarkTTT )] )
                    *M_ostr << __dataVolumePost[boost::get<0>( *itMarkTTT )][boost::get<1>( *itMarkTTT )][boost::get<2>( *itMarkTTT )]
                             << "};\n";

                else
                    *M_ostr << __dataMemGlobSurf2[1][boost::get<0>( *itMarkTTT )][boost::get<1>( *itMarkTTT )][boost::get<2>( *itMarkTTT )]
                             <<"["
                             << __dataVolumePost[boost::get<0>( *itMarkTTT )][boost::get<1>( *itMarkTTT )][boost::get<2>( *itMarkTTT )]
                             <<"] };\n";

#endif
            }

            ++itMarkName;
        }

        ++itMarkType;
    }

    //++itShape;
    //}



    //std::cout << "\n HOLA "<< M_ostr->str()<<std::endl;
    //return M_ostr->str();

}





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
    std::ostringstream __ostr;
    __ostr << "Point(" << boost::get<0>( *__dg )->cptPt()
           << ") = {"
           << std::scientific << std::setprecision( 16 ) << __x1 << ","
           << std::scientific << std::setprecision( 16 ) << __x2 << ","
           << std::scientific << std::setprecision( 16 ) << __x3 <<", h};\n";
    boost::get<0>( *__dg )->updateOstr( __ostr.str() );
    ++( boost::get<0>( *__dg )->M_cptPt );
}

/*_________________________________________________*/

void
writeLine( uint16_type __numLoc, data_geo_ptrtype __dg ,uint16_type __n1, uint16_type __n2 )
{
    ( *( boost::get<1>( *__dg ) ) )[1][__numLoc] = boost::get<0>( *__dg )->cptLine();
    std::ostringstream __ostr;
    __ostr << "Line(" << boost::get<0>( *__dg )->cptLine()
           << ") = {"
           << ( *( boost::get<1>( *__dg ) ) )[0][__n1] << ","
           << ( *( boost::get<1>( *__dg ) ) )[0][__n2] << "};\n";

    boost::get<0>( *__dg )->updateOstr( __ostr.str() );

    ++boost::get<0>( *__dg )->M_cptLine;
}

/*_________________________________________________*/

void
writeCircle( uint16_type __numLoc, data_geo_ptrtype __dg ,uint16_type __n1, uint16_type __n2, uint16_type __n3 )
{
    ( *( boost::get<1>( *__dg ) ) )[1][__numLoc] = boost::get<0>( *__dg )->cptLine();

    std::ostringstream __ostr;
    __ostr << "Circle(" << boost::get<0>( *__dg )->cptLine() //cptCircle
           << ") = {"
           << ( *( boost::get<1>( *__dg ) ) )[0][__n1] << ","
           << ( *( boost::get<1>( *__dg ) ) )[0][__n2] << ","
           << ( *( boost::get<1>( *__dg ) ) )[0][__n3] << "};\n";

    boost::get<0>( *__dg )->updateOstr( __ostr.str() );

    ++boost::get<0>( *__dg )->M_cptLine;
}

/*_________________________________________________*/

void
writeEllipse( uint16_type __numLoc, data_geo_ptrtype __dg ,uint16_type __n1, uint16_type __n2, uint16_type __n3, uint16_type __n4 )
{
    ( *( boost::get<1>( *__dg ) ) )[1][__numLoc] = boost::get<0>( *__dg )->cptLine();

    std::ostringstream __ostr;
    __ostr << "Ellipse(" << boost::get<0>( *__dg )->cptLine() //cptCircle
           << ") = {"
           << ( *( boost::get<1>( *__dg ) ) )[0][__n1] << ","
           << ( *( boost::get<1>( *__dg ) ) )[0][__n2] << ","
           << ( *( boost::get<1>( *__dg ) ) )[0][__n3] << ","
           << ( *( boost::get<1>( *__dg ) ) )[0][__n4] << "};\n";

    boost::get<0>( *__dg )->updateOstr( __ostr.str() );

    ++boost::get<0>( *__dg )->M_cptLine;
}

/*_________________________________________________*/

void
writeSpline( uint16_type __numLoc, data_geo_ptrtype __dg ,Loop __loop )
{
    ( *( boost::get<1>( *__dg ) ) )[1][__numLoc] = boost::get<0>( *__dg )->cptLine();

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

    ++boost::get<0>( *__dg )->M_cptLine;
}

/*_________________________________________________*/

void
writeBSpline( uint16_type __numLoc, data_geo_ptrtype __dg ,Loop __loop )
{
    ( *( boost::get<1>( *__dg ) ) )[1][__numLoc] = boost::get<0>( *__dg )->cptLine();

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

    ++boost::get<0>( *__dg )->M_cptLine;
}

/*_________________________________________________*/

void
writeLineLoop( uint16_type __numLoc, data_geo_ptrtype __dg , Loop /*const*/ __loop )
{
    ( *( boost::get<1>( *__dg ) ) )[2][__numLoc] = boost::get<0>( *__dg )->cptLineLoop();

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

    ++boost::get<0>( *__dg )->M_cptLineLoop;
}

/*_________________________________________________*/

void
writePtInSurface( data_geo_ptrtype __dg , uint16_type __indLocPt,uint16_type __indLocSurf )
{

    auto indPtGlob = ( *( boost::get<1>( *__dg ) ) )[0][__indLocPt];
    ( *( boost::get<7>( *__dg ) ) )[0][( *( boost::get<1>( *__dg ) ) )[3][__indLocSurf]].push_back( indPtGlob );
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
    uint16_type cptSurf = __dg->get<0>()->cptSurface();
    uint16_type refLineLoop= ( *( boost::get<1>( *__dg ) ) )[2][__ind];

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
    uint16_type cptVol = __dg->get<0>()->cptVolume();
    uint16_type refSurfaceLoop= ( *( boost::get<1>( *__dg ) ) )[4][__ind];

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

    //coefficient du plan qui a pour normal dir
    double a=dir( 0 );
    double b=dir( 1 );
    double c=dir( 2 );
    double d=-dir( 0 )*centre( 0 )-dir( 1 )*centre( 1 )-dir( 2 )*centre( 2 ); //-N scalaire OA

    double rayon=1;

    //un point du plan autre que centre
    Node ptBis( 0,0,0 );

    if ( a==0 )
    {
        ptBis( 0 )=centre( 0 )+rayon;

        if ( !( std::abs( c )<1e-8 ) )
        {
            ptBis( 1 )=centre( 1 );
            ptBis( 2 )=( -b*ptBis( 1 )-d )/c;
        }

        else if ( !( std::abs( b )<1e-8 ) )
        {
            ptBis( 2 )=centre( 2 );
            ptBis( 1 )=( -c*ptBis( 2 )-d )/b;
        }
    }

    else if ( b==0 )
    {
        ptBis( 1 )=centre( 1 )+rayon;

        if ( ! (std::abs( c )<1e-8 ) )
        {
            ptBis( 0 )=centre( 0 );
            ptBis( 2 )=( -a*ptBis( 0 )-d )/c;
        }

        else if ( ! (std::abs( a )<1e-8 ) )
        {
            ptBis( 2 )=centre( 2 );
            ptBis( 0 )=( -c*ptBis( 2 )-d )/a;
        }
    }

    else if ( c==0 )
    {
        double xtemp=centre( 0 );
        double ztemp=centre( 1 )+rayon;
        double ytemp=( -a*xtemp-c*ztemp-d )/b;
        ptBis = Node( xtemp,ytemp,ztemp );
    }

    else
        ptBis = Node( centre( 0 )+rayon,centre( 1 ),centre( 2 ) );

    //un veteur du plan
    Node u( ptBis( 0 )-centre( 0 ),ptBis( 1 )-centre( 1 ),ptBis( 2 )-centre( 2 ) );

    //a=u2v3-u3v2 ; b=u3v1-u1v3 ; c=u1v2-u2v1
    //deuxieme vecteur qui forme une base orthogonal (u,v,dir) avec v=u ProdVect dir
    Node v( u( 1 )*dir( 2 ) - u( 2 )*dir( 1 ),
            u( 0 )*dir( 2 ) - u( 2 )*dir( 0 ),
            u( 0 )*dir( 1 ) - u( 1 )*dir( 0 ) );

    double norm_u = std::sqrt( u( 0 )*u( 0 )+u( 1 )*u( 1 )+u( 2 )*u( 2 ) );
    double norm_v = std::sqrt( v( 0 )*v( 0 )+v( 1 )*u( 1 )+v( 2 )*v( 2 ) );
    u( 0 )/=norm_u;
    u( 1 )/=norm_u;
    u( 2 )/=norm_u;
    v( 0 )/=norm_v;
    v( 1 )/=norm_v;
    v( 2 )/=norm_v;

    Node D( dir( 0 ),dir( 1 ),dir( 2 ) );

    return boost::make_tuple( D,u,v );

}

/*_________________________________________________*/



/*_________________________________________________*
 *_________________________________________________*
 * Function user                                   *
 *_________________________________________________*
 *_________________________________________________*/


void
runLine( data_geo_ptrtype dg )
{
    node_type PtA = param<0>( dg );
    node_type PtB = param<1>( dg );

    writePoint( 1, dg , PtA( 0 ), PtA( 1 ) );
    writePoint( 2, dg , PtB( 0 ), PtB( 1 ) );
    writeLine( 1, dg , 1 , 2 );
}


void
runTriangle( data_geo_ptrtype dg )
{
    node_type PtA = param<0>( dg );
    node_type PtB = param<1>( dg );
    node_type PtC = param<2>( dg );

    writePoint( 1, dg , PtA( 0 ), PtA( 1 ) );
    writePoint( 2, dg , PtB( 0 ), PtB( 1 ) );
    writePoint( 3, dg , PtC( 0 ), PtC( 1 ) );

    writeLine( 1, dg , 1 , 2 );
    writeLine( 2, dg , 2 , 3 );
    writeLine( 3, dg , 3 , 1 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3 );

    writePlaneSurface( 1, dg, 1 );
}



void
runRectangle( data_geo_ptrtype dg )
{

    node_type PtA = param<0>( dg );
    node_type PtB = param<1>( dg );

    writePoint( 1, dg , PtA( 0 ), PtA( 1 ) );
    writePoint( 2, dg , PtB( 0 ), PtA( 1 ) );
    writePoint( 3, dg , PtB( 0 ), PtB( 1 ) );
    writePoint( 4, dg , PtA( 0 ), PtB( 1 ) );

    writeLine( 1, dg , 1 , 2 );
    writeLine( 2, dg , 2 , 3 );
    writeLine( 3, dg , 3 , 4 );
    writeLine( 4, dg , 4 , 1 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4 );

    writePlaneSurface( 1, dg, 1 );

}

void
runQuadrangle( data_geo_ptrtype dg )
{

    node_type PtA = param<0>( dg );
    node_type PtB = param<1>( dg );
    node_type PtC = param<2>( dg );
    node_type PtD = param<3>( dg );

    writePoint( 1, dg , PtA( 0 ), PtA( 1 ) );
    writePoint( 2, dg , PtB( 0 ), PtB( 1 ) );
    writePoint( 3, dg , PtC( 0 ), PtC( 1 ) );
    writePoint( 4, dg , PtD( 0 ), PtD( 1 ) );

    writeLine( 1, dg , 1 , 2 );
    writeLine( 2, dg , 2 , 3 );
    writeLine( 3, dg , 3 , 4 );
    writeLine( 4, dg , 4 , 1 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4 );

    writePlaneSurface( 1, dg, 1 );

}

void
runPentagon( data_geo_ptrtype dg )
{

    node_type PtA = param<0>( dg );
    node_type PtB = param<1>( dg );
    node_type PtC = param<2>( dg );
    node_type PtD = param<3>( dg );
    node_type PtE = param<4>( dg );

    writePoint( 1, dg , PtA( 0 ), PtA( 1 ) );
    writePoint( 2, dg , PtB( 0 ), PtB( 1 ) );
    writePoint( 3, dg , PtC( 0 ), PtC( 1 ) );
    writePoint( 4, dg , PtD( 0 ), PtD( 1 ) );
    writePoint( 5, dg , PtE( 0 ), PtE( 1 ) );

    writeLine( 1, dg , 1 , 2 );
    writeLine( 2, dg , 2 , 3 );
    writeLine( 3, dg , 3 , 4 );
    writeLine( 4, dg , 4 , 5 );
    writeLine( 5, dg , 5 , 1 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4>>5 );

    writePlaneSurface( 1, dg, 1 );

}


void
runHexagon( data_geo_ptrtype dg )
{

    node_type PtA = param<0>( dg );
    node_type PtB = param<1>( dg );
    node_type PtC = param<2>( dg );
    node_type PtD = param<3>( dg );
    node_type PtE = param<4>( dg );
    node_type PtF = param<5>( dg );

    writePoint( 1, dg , PtA( 0 ), PtA( 1 ) );
    writePoint( 2, dg , PtB( 0 ), PtB( 1 ) );
    writePoint( 3, dg , PtC( 0 ), PtC( 1 ) );
    writePoint( 4, dg , PtD( 0 ), PtD( 1 ) );
    writePoint( 5, dg , PtE( 0 ), PtE( 1 ) );
    writePoint( 6, dg , PtF( 0 ), PtF( 1 ) );

    writeLine( 1, dg , 1 , 2 );
    writeLine( 2, dg , 2 , 3 );
    writeLine( 3, dg , 3 , 4 );
    writeLine( 4, dg , 4 , 5 );
    writeLine( 5, dg , 5 , 6 );
    writeLine( 6, dg , 6 , 1 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4>>5>>6 );

    writePlaneSurface( 1, dg, 1 );

}


void
runCircle( data_geo_ptrtype dg )
{
    node_type PtA = param<0>( dg );
    node_type PtB = param<1>( dg );

    writePoint( 1, dg , PtB( 0 ), PtB( 1 ) );
    writePoint( 2, dg , PtA( 0 ), PtA( 1 ) );
    writePoint( 3, dg , 2*PtA( 0 )-PtB( 0 ), 2*PtA( 1 )-PtB( 1 ) );

    writeCircle( 1, dg, 1, 2, 3 );
    writeCircle( 2, dg, 3, 2, 1 );

    writeLineLoop( 1, dg, Loop()>>1>>2 );

    writePlaneSurface( 1, dg, 1 );
    writePtInSurface( dg,2,1 );
}


void
runEllipse( data_geo_ptrtype dg )
{
    node_type PtC = param<0>( dg );
    node_type PtMinor = param<1>( dg );
    node_type PtMajor = param<2>( dg );

    writePoint( 1, dg , PtC( 0 ), PtC( 1 ) );
    writePoint( 2, dg , PtMinor( 0 ), PtMinor( 1 ) );
    writePoint( 3, dg , 2*PtC( 0 )-PtMinor( 0 ), 2*PtC( 1 )-PtMinor( 1 ) );
    writePoint( 4, dg , PtMajor( 0 ), PtMajor( 1 ) );
    writePoint( 5, dg , 2*PtC( 0 )-PtMajor( 0 ), 2*PtC( 1 )-PtMajor( 1 ) );

    writeEllipse( 1, dg, 2, 1, 4, 4 );
    writeEllipse( 2, dg, 4, 1, 3, 3 );
    writeEllipse( 3, dg, 3, 1, 5, 5 );
    writeEllipse( 4, dg, 5, 1, 2, 2 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4 );

    writePlaneSurface( 1, dg, 1 );
    writePtInSurface( dg,2,1 );
}



void
runPie( data_geo_ptrtype dg )
{
    node_type PtA = param<0>( dg );
    node_type PtB = param<1>( dg );
    node_type PtC = param<2>( dg );

    writePoint( 1, dg , PtB( 0 ), PtB( 1 ) );
    writePoint( 2, dg , PtA( 0 ), PtA( 1 ) );
    writePoint( 3, dg , PtC( 0 ), PtC( 1 ) );

    writeCircle( 1, dg, 1, 2, 3 );
    writeLine( 2, dg , 3 , 2 );
    writeLine( 3, dg , 2 , 1 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3 );

    writePlaneSurface( 1, dg, 1 );
    writePtInSurface( dg,2,1 );

}

void
runSpecial_1a( data_geo_ptrtype dg )
{
    double yh=1.0;
    Node a1( 0.0, yh );
    Node a2( 3.0, yh );
    Node a3( 6.0, yh+0.7 );
    Node a4( 7.3, yh-0.5 );
    Node a5( 8.5, yh );
    Node a6( 11.0, yh );
    double ep=0.3;
    //_______________________________________________//
    writePoint( 1, dg , a1( 0 ), a1( 1 ) );
    writePoint( 2, dg , a2( 0 ), a2( 1 ) );
    writePoint( 3, dg , a3( 0 ), a3( 1 ) );
    writePoint( 4, dg , a4( 0 ), a4( 1 ) );
    writePoint( 5, dg , a5( 0 ), a5( 1 ) );
    writePoint( 6, dg , a6( 0 ), a6( 1 ) );

    writeSpline( 1, dg, Loop()>>1>>2>>3>>4>>5>>6 );

    writePoint( 7, dg , a1( 0 ), a1( 1 )+ep );
    writePoint( 8, dg , a2( 0 ), a2( 1 )+ep );
    writePoint( 9, dg , a3( 0 ), a3( 1 )+ep );
    writePoint( 10, dg , a4( 0 ), a4( 1 )+ep );
    writePoint( 11, dg , a5( 0 ), a5( 1 )+ep );
    writePoint( 12, dg , a6( 0 ), a6( 1 )+ep );

    writeSpline( 2, dg, Loop()>>7>>8>>9>>10>>11>>12 );

    writeLine( 3, dg, 1,7 );
    writeLine( 4, dg, 6,12 );

    writeLineLoop( 1, dg, Loop()>>1>>4>>-2>>-3 );
    writePlaneSurface( 1, dg, 1 );

    writePoint( 13, dg , a1( 0 ), -a1( 1 ) );
    writePoint( 14, dg , a2( 0 ), -a2( 1 ) );
    writePoint( 15, dg , a3( 0 ), -a3( 1 ) );
    writePoint( 16, dg , a4( 0 ), -a4( 1 ) );
    writePoint( 17, dg , a5( 0 ), -a5( 1 ) );
    writePoint( 18, dg , a6( 0 ), -a6( 1 ) );

    writeSpline( 5, dg, Loop()>>13>>14>>15>>16>>17>>18 );

    writePoint( 19, dg , a1( 0 ), -a1( 1 )-ep );
    writePoint( 20, dg , a2( 0 ), -a2( 1 )-ep );
    writePoint( 21, dg , a3( 0 ), -a3( 1 )-ep );
    writePoint( 22, dg , a4( 0 ), -a4( 1 )-ep );
    writePoint( 23, dg , a5( 0 ), -a5( 1 )-ep );
    writePoint( 24, dg , a6( 0 ), -a6( 1 )-ep );

    writeSpline( 6, dg, Loop()>>19>>20>>21>>22>>23>>24 );

    writeLine( 7, dg, 13,19 );
    writeLine( 8, dg, 18,24 );

    writeLineLoop( 2, dg, Loop()>>5>>8>>-6>>-7 );
    writePlaneSurface( 2, dg, 2 );

}

void
runSpecial_1b( data_geo_ptrtype dg )
{
    double yh=1.0;
    Node a1( 0.0, yh );
    Node a2( 3.0, yh );
    Node a3( 6.0, yh+0.7 );
    Node a4( 7.3, yh-0.5 );
    Node a5( 8.5, yh );
    Node a6( /*11.0*/16.0, yh );
    //_______________________________________________//
    writePoint( 1, dg , a1( 0 ), a1( 1 ) );
    writePoint( 2, dg , a2( 0 ), a2( 1 ) );
    writePoint( 3, dg , a3( 0 ), a3( 1 ) );
    writePoint( 4, dg , a4( 0 ), a4( 1 ) );
    writePoint( 5, dg , a5( 0 ), a5( 1 ) );
    writePoint( 6, dg , a6( 0 ), a6( 1 ) );

    writeSpline( 1, dg, Loop()>>1>>2>>3>>4>>5>>6 );

    writePoint( 7, dg , a1( 0 ), -a1( 1 ) );
    writePoint( 8, dg , a2( 0 ), -a2( 1 ) );
    writePoint( 9, dg , a3( 0 ), -a3( 1 ) );
    writePoint( 10, dg , a4( 0 ), -a4( 1 ) );
    writePoint( 11, dg , a5( 0 ), -a5( 1 ) );
    writePoint( 12, dg , a6( 0 ), -a6( 1 ) );

    writeSpline( 2, dg, Loop()>>7>>8>>9>>10>>11>>12 );

    writeLine( 3, dg, 1,7 );
    writeLine( 4, dg, 6,12 );

    writeLineLoop( 1, dg, Loop()>>1>>4>>-2>>-3 );

    writePlaneSurface( 1, dg, 1 );

}

void
runPeanut( data_geo_ptrtype dg )
{
    node_type PtCenter = param<0>( dg );
    node_type majorRadiusParam = param<1>( dg ); //2
    node_type minorRadiusParam = param<2>( dg ); //1
    node_type penautRadiusParam = param<3>( dg ); //0.1
    double majorRadius = majorRadiusParam( 0 );
    double minorRadius = minorRadiusParam( 0 );
    double penautRadius = penautRadiusParam( 0 );
    writePoint( 1, dg , PtCenter( 0 )               , PtCenter( 1 )+penautRadius, 0. );
    writePoint( 2, dg , PtCenter( 0 )               , PtCenter( 1 )-penautRadius, 0. );
    writePoint( 3, dg , PtCenter( 0 )-majorRadius   , PtCenter( 1 )             , 0. );
    writePoint( 4, dg , PtCenter( 0 )+majorRadius   , PtCenter( 1 )             , 0. );
    writePoint( 5, dg , PtCenter( 0 )-majorRadius/4., PtCenter( 1 )+minorRadius , 0. );
    writePoint( 6, dg , PtCenter( 0 )-majorRadius/4., PtCenter( 1 )-minorRadius , 0. );
    writePoint( 7, dg , PtCenter( 0 )+majorRadius/4., PtCenter( 1 )+minorRadius , 0. );
    writePoint( 8, dg , PtCenter( 0 )+majorRadius/4., PtCenter( 1 )-minorRadius , 0. );
    writeBSpline( 1, dg, Loop()>>1>>5>>3>>6>>2>>8>>4>>7>>1 );

    writeLineLoop( 1, dg, Loop()>>1 );
    writePlaneSurface( 1, dg, 1 );
}

void
runHexahedron( data_geo_ptrtype dg )
{

    node_type Pt1 = param<0>( dg );
    node_type Pt2 = param<1>( dg );
    node_type Pt3 = param<2>( dg );
    node_type Pt4 = param<3>( dg );
    node_type Pt5 = param<4>( dg );
    node_type Pt6 = param<5>( dg );
    node_type Pt7 = param<6>( dg );
    node_type Pt8 = param<7>( dg );

    writePoint( 1, dg , Pt1( 0 ), Pt1( 1 ), Pt1( 2 ) );
    writePoint( 2, dg , Pt2( 0 ), Pt2( 1 ), Pt2( 2 ) );
    writePoint( 3, dg , Pt3( 0 ), Pt3( 1 ), Pt3( 2 ) );
    writePoint( 4, dg , Pt4( 0 ), Pt4( 1 ), Pt4( 2 ) );
    writePoint( 5, dg , Pt5( 0 ), Pt5( 1 ), Pt5( 2 ) );
    writePoint( 6, dg , Pt6( 0 ), Pt6( 1 ), Pt6( 2 ) );
    writePoint( 7, dg , Pt7( 0 ), Pt7( 1 ), Pt7( 2 ) );
    writePoint( 8, dg , Pt8( 0 ), Pt8( 1 ), Pt8( 2 ) );

    writeLine( 1, dg , 1 , 2 );
    writeLine( 2, dg , 2 , 3 );
    writeLine( 3, dg , 3 , 4 );
    writeLine( 4, dg , 4 , 1 );
    writeLine( 5, dg , 5 , 6 );
    writeLine( 6, dg , 6 , 7 );
    writeLine( 7, dg , 7 , 8 );
    writeLine( 8, dg , 8 , 5 );
    writeLine( 9, dg , 1 , 5 );
    writeLine( 10, dg , 2 , 6 );
    writeLine( 11, dg , 3 , 7 );
    writeLine( 12, dg , 4 , 8 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4 );
    writePlaneSurface( 1, dg, 1 );
    writeLineLoop( 2, dg, Loop()>>-5>>-8>>-7>>-6 );
    writePlaneSurface( 2, dg, 2 );
    writeLineLoop( 3, dg, Loop()>>-1>>9>>5>>-10 );
    writePlaneSurface( 3, dg, 3 );
    writeLineLoop( 4, dg, Loop()>>10>>6>>-11>>-2 );
    writePlaneSurface( 4, dg, 4 );
    writeLineLoop( 5, dg, Loop()>>11>>7>>-12>>-3 );
    writePlaneSurface( 5, dg, 5 );
    writeLineLoop( 6, dg, Loop()>>8>>-9>>-4>>12 );
    writePlaneSurface( 6, dg, 6 );

    writeSurfaceLoop( 1, dg, Loop()>>1>>2>>3>>4>>5>>6 );
    writeVolume( 1, dg, 1 );
}

void
runTetrahedron( data_geo_ptrtype dg )
{
    node_type Pt1 = param<0>( dg );
    node_type Pt2 = param<1>( dg );
    node_type Pt3 = param<2>( dg );
    node_type Pt4 = param<3>( dg );

    writePoint( 1, dg , Pt1( 0 ), Pt1( 1 ), Pt1( 2 ) );
    writePoint( 2, dg , Pt2( 0 ), Pt2( 1 ), Pt2( 2 ) );
    writePoint( 3, dg , Pt3( 0 ), Pt3( 1 ), Pt3( 2 ) );
    writePoint( 4, dg , Pt4( 0 ), Pt4( 1 ), Pt4( 2 ) );

    writeLine( 1, dg , 1 , 2 );
    writeLine( 2, dg , 2 , 3 );
    writeLine( 3, dg , 3 , 1 );
    writeLine( 4, dg , 1 , 4 );
    writeLine( 5, dg , 2 , 4 );
    writeLine( 6, dg , 3 , 4 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3 );
    writePlaneSurface( 1, dg, 1 );

    //writeLineLoop( 2, dg, Loop()>>5>>-4>>1 );
    writeLineLoop( 2, dg, Loop()>>-1>>4>>-5 );
    writePlaneSurface( 2, dg, 2 );
    //writeLineLoop( 3, dg, Loop()>>2>>6>>-5 );
    writeLineLoop( 3, dg, Loop()>>5>>-6>>-2 );
    writePlaneSurface( 3, dg, 3 );
    //writeLineLoop( 4, dg, Loop()>>3>>4>>-6 );
    writeLineLoop( 4, dg, Loop()>>6>>-4>>-3 );
    writePlaneSurface( 4, dg, 4 );

    writeSurfaceLoop( 1, dg, Loop()>>1>>2>>3>>4 );

    writeVolume( 1, dg, 1 );
}


void
runCube( data_geo_ptrtype dg )
{

    node_type Pt1 = param<0>( dg );
    node_type Pt7 = param<1>( dg );
    double lenx = Pt7( 0 )-Pt1( 0 );
    double leny = Pt7( 1 )-Pt1( 1 );
    double lenz = Pt7( 2 )-Pt1( 2 );
    node_type Pt2 = Pt1;
    Pt2( 0 )+=lenx;
    node_type Pt3 = Pt1;
    Pt3( 0 )+=lenx;
    Pt3( 1 )+=leny;
    node_type Pt4 = Pt1;
    Pt4( 1 )+=leny;
    node_type Pt5 = Pt1;
    Pt5( 2 )+=lenz;
    node_type Pt6 = Pt1;
    Pt6( 0 )+=lenx;
    Pt6( 2 )+=lenz;
    // Pt7 was give above by the user
    node_type Pt8 = Pt1;
    Pt8( 1 )+=leny;
    Pt8( 2 )+=lenz;

    writePoint( 1, dg , Pt1( 0 ), Pt1( 1 ), Pt1( 2 ) );
    writePoint( 2, dg , Pt2( 0 ), Pt2( 1 ), Pt2( 2 ) );
    writePoint( 3, dg , Pt3( 0 ), Pt3( 1 ), Pt3( 2 ) );
    writePoint( 4, dg , Pt4( 0 ), Pt4( 1 ), Pt4( 2 ) );
    writePoint( 5, dg , Pt5( 0 ), Pt5( 1 ), Pt5( 2 ) );
    writePoint( 6, dg , Pt6( 0 ), Pt6( 1 ), Pt6( 2 ) );
    writePoint( 7, dg , Pt7( 0 ), Pt7( 1 ), Pt7( 2 ) );
    writePoint( 8, dg , Pt8( 0 ), Pt8( 1 ), Pt8( 2 ) );

    writeLine( 1, dg , 1 , 2 );
    writeLine( 2, dg , 2 , 3 );
    writeLine( 3, dg , 3 , 4 );
    writeLine( 4, dg , 4 , 1 );
    writeLine( 5, dg , 5 , 6 );
    writeLine( 6, dg , 6 , 7 );
    writeLine( 7, dg , 7 , 8 );
    writeLine( 8, dg , 8 , 5 );
    writeLine( 9, dg , 1 , 5 );
    writeLine( 10, dg , 2 , 6 );
    writeLine( 11, dg , 3 , 7 );
    writeLine( 12, dg , 4 , 8 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4 );
    writePlaneSurface( 1, dg, 1 );
    writeLineLoop( 2, dg, Loop()>>-5>>-8>>-7>>-6 );
    writePlaneSurface( 2, dg, 2 );
    writeLineLoop( 3, dg, Loop()>>-1>>9>>5>>-10 );
    writePlaneSurface( 3, dg, 3 );
    writeLineLoop( 4, dg, Loop()>>10>>6>>-11>>-2 );
    writePlaneSurface( 4, dg, 4 );
    writeLineLoop( 5, dg, Loop()>>11>>7>>-12>>-3 );
    writePlaneSurface( 5, dg, 5 );
    writeLineLoop( 6, dg, Loop()>>8>>-9>>-4>>12 );
    writePlaneSurface( 6, dg, 6 );

    writeSurfaceLoop( 1, dg, Loop()>>1>>2>>3>>4>>5>>6 );

    writeVolume( 1, dg, 1 );
}

void
runCylindre( data_geo_ptrtype dg )
{
#if 0
    node_type centre = param<0>( dg );
    double par = param<1>( dg )( 0 );

    Node c1( -1,-par,0 );
    Node c2( -1, 0,-par );
    Node c3( -1, par,0 );
    Node c4( -1, 0, par );
    writePoint( 1, dg , centre( 0 ), centre( 1 ), centre( 2 ) );
    writePoint( 2, dg , c1( 0 ), c1( 1 ), c1( 2 ) );
    writePoint( 3, dg , c2( 0 ), c2( 1 ), c2( 2 ) );
    writePoint( 4, dg , c3( 0 ), c3( 1 ), c3( 2 ) );
    writePoint( 5, dg , c4( 0 ), c4( 1 ), c4( 2 ) );

    writeCircle( 1, dg, 2, 1, 3 );
    writeCircle( 2, dg, 3, 1, 4 );
    writeCircle( 3, dg, 4, 1, 5 );
    writeCircle( 4, dg, 5, 1, 2 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4 );
    writePlaneSurface( 1, dg, 1 );

    //volume 1, extrude surface 1, creation surface 2,3,4,5,6
    writeExtrudeSurface( 1,dg,1,Loop()>>2>>3>>4>>5>>6 );
#endif

    node_type centre = param<0>( dg );
    node_type direction = param<1>( dg );
    double rayon = param<2>( dg )( 0 );
    double longueur=param<3>( dg )( 0 );

    auto basis = computeBasisOrthogonal( direction,centre );
    Node dir( boost::get<0>( basis ) );
    Node u( boost::get<1>( basis ) );
    Node v( boost::get<2>( basis ) );

    Node geoX1( centre( 0 )+u( 0 )*rayon,
                centre( 1 )+u( 1 )*rayon,
                centre( 2 )+u( 2 )*rayon  );

    Node geoX2( centre( 0 )+v( 0 )*rayon,
                centre( 1 )+v( 1 )*rayon,
                centre( 2 )+v( 2 )*rayon  );

    Node geoX3( centre( 0 )-u( 0 )*rayon,
                centre( 1 )-u( 1 )*rayon,
                centre( 2 )-u( 2 )*rayon  );

    Node geoX4( centre( 0 )-v( 0 )*rayon,
                centre( 1 )-v( 1 )*rayon,
                centre( 2 )-v( 2 )*rayon  );

    Node centre2( centre( 0 )+longueur*dir( 0 ),
                  centre( 1 )+longueur*dir( 1 ),
                  centre( 2 )+longueur*dir( 2 ) );

    Node geoX5( geoX1( 0 )+longueur*dir( 0 ),
                geoX1( 1 )+longueur*dir( 1 ),
                geoX1( 2 )+longueur*dir( 2 ) );

    Node geoX6( geoX2( 0 )+longueur*dir( 0 ),
                geoX2( 1 )+longueur*dir( 1 ),
                geoX2( 2 )+longueur*dir( 2 ) );

    Node geoX7( geoX3( 0 )+longueur*dir( 0 ),
                geoX3( 1 )+longueur*dir( 1 ),
                geoX3( 2 )+longueur*dir( 2 ) );

    Node geoX8( geoX4( 0 )+longueur*dir( 0 ),
                geoX4( 1 )+longueur*dir( 1 ),
                geoX4( 2 )+longueur*dir( 2 ) );

    //--------------------------------------------------------------------------//

    writePoint( 1, dg, centre( 0 ), centre( 1 ) ,centre( 2 ) );
    writePoint( 2, dg, geoX1( 0 ), geoX1( 1 ) ,geoX1( 2 ) );
    writePoint( 3, dg, geoX2( 0 ), geoX2( 1 ) ,geoX2( 2 ) );
    writePoint( 4, dg, geoX3( 0 ), geoX3( 1 ) ,geoX3( 2 ) );
    writePoint( 5, dg, geoX4( 0 ), geoX4( 1 ) ,geoX4( 2 ) );

    writeCircle( 1, dg, 2,1,3 );
    writeCircle( 2, dg, 3,1,4 );
    writeCircle( 3, dg, 4,1,5 );
    writeCircle( 4, dg, 5,1,2 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4 );

    writePoint( 6, dg, centre2( 0 ), centre2( 1 ) ,centre2( 2 ) );
    writePoint( 7, dg, geoX5( 0 ), geoX5( 1 ) ,geoX5( 2 ) );
    writePoint( 8, dg, geoX6( 0 ), geoX6( 1 ) ,geoX6( 2 ) );
    writePoint( 9, dg, geoX7( 0 ), geoX7( 1 ) ,geoX7( 2 ) );
    writePoint( 10, dg, geoX8( 0 ), geoX8( 1 ) ,geoX8( 2 ) );

    writeCircle( 5, dg, 7,6,8 );
    writeCircle( 6, dg,8,6,9 );
    writeCircle( 7, dg,9,6,10 );
    writeCircle( 8 ,dg,10,6,7 );

    writeLineLoop( 2, dg, Loop()>>-5>>-8>>-7>>-6 );

    writeLine( 9, dg, 4, 9 );
    writeLine( 10, dg, 5, 10 );
    writeLine( 11, dg, 2, 7 );
    writeLine( 12, dg, 8, 3 );
    writeLineLoop( 13, dg, Loop()>>-9>>-2>>-12>>6 );
    writeLineLoop( 15, dg, Loop()>>9>>7>>-10>>-3 );
    writeLineLoop( 17, dg, Loop()>>10>>8>>-11>>-4 );
    writeLineLoop( 19, dg, Loop()>>11>>5>>12>>-1 );

    writePlaneSurface( 1, dg, 1 );
    writePlaneSurface( 2, dg, 2 );
    writePtInSurface(dg,1,1 );
    writePtInSurface(dg,6,2 );

    writeRuledSurface( 3, dg, 13 );
    writeRuledSurface( 4, dg, 15 );
    writeRuledSurface( 5, dg, 17 );
    writeRuledSurface( 6, dg, 19 );

    writeSurfaceLoop( 1, dg, Loop()>>5>>4>>3>>2>>6>>1 );
    writeVolume( 1, dg, 1 );

}

void
runSphere( data_geo_ptrtype dg )
{

    node_type centre = param<0>( dg );
    double R = param<1>( dg )( 0 );

    double x_center = centre( 0 );
    double y_center = centre( 1 );
    double z_center = centre( 2 );

    writePoint( 1, dg, x_center, y_center, z_center );
    writePoint( 2, dg, x_center - R, y_center, z_center );
    writePoint( 4, dg, x_center, y_center - R, z_center );
    writePoint( 5, dg, x_center + R, y_center, z_center );
    writePoint( 8, dg, x_center, y_center, z_center - R );
    writePoint( 11, dg, x_center, y_center + R, z_center );
    writePoint( 14, dg, x_center, y_center, z_center + R );
    writeCircle( 1, dg, 2, 1, 4 );
    writeCircle( 2, dg, 4, 1, 5 );
    writeCircle( 3, dg, 2, 1, 8 );
    writeCircle( 4, dg, 4, 1, 8 );
    writeCircle( 6, dg, 2, 1, 11 );
    writeCircle( 7, dg, 8, 1, 11 );
    writeCircle( 9, dg, 2, 1, 14 );
    writeCircle( 10, dg, 11, 1, 14 );
    writeCircle( 13, dg, 14, 1, 4 );
    writeCircle( 15, dg, 8, 1, 5 );
    writeCircle( 18, dg, 11, 1, 5 );
    writeCircle( 21, dg, 14, 1, 5 );

#if 0
    writeLineLoop( 1, dg, Loop()>>1>>4>>-3 );
    writeRuledSurface( 1, dg,1 );
    writeLineLoop( 2, dg, Loop()>>3>>7>>-6 );
    writeRuledSurface( 2, dg,2 );
    writeLineLoop( 3, dg, Loop()>>6>>10>>-9 );
    writeRuledSurface( 3,dg,3 );
    writeLineLoop( 4, dg, Loop()>>9>>13>>-1 );
    writeRuledSurface( 4, dg, 4 );
    writeLineLoop( 5, dg, Loop()>>-15>>-4>>2 );
    writeRuledSurface( 5, dg, 5 );
    writeLineLoop( 6, dg, Loop()>>-18>>-7>>15 );
    writeRuledSurface( 6, dg, 6 );
    writeLineLoop( 7, dg, Loop()>>-21>>-10>>18 );
    writeRuledSurface( 7,dg,7 );
    writeLineLoop( 8, dg, Loop()>>-2>>-13>>21 );
    writeRuledSurface( 8, dg,8 );
    writeSurfaceLoop( 1, dg, Loop()>>1>>4>>3>>2>>6>>7>>8>>5 );

#else
    writeLineLoop( 1, dg, Loop()>>1>>4>>-3 );
    writeRuledSurface( 1, dg,1 );
    writeLineLoop( 2, dg, Loop()>>7>>-6>>3 );//
    writeRuledSurface( 2, dg,2 );
    writeLineLoop( 3, dg, Loop()>>6>>10>>-9 );//
    writeRuledSurface( 3,dg,3 );
    writeLineLoop( 4, dg, Loop()>>13>>-1>>9 );//
    writeRuledSurface( 4, dg, 4 );
    writeLineLoop( 5, dg, Loop()>>-4>>2>>-15 );//
    writeRuledSurface( 5, dg, 5 );
    writeLineLoop( 6, dg, Loop()>>-18>>-7>>15 );//
    writeRuledSurface( 6, dg, 6 );
    writeLineLoop( 7, dg, Loop()>>-10>>18>>-21 );//
    writeRuledSurface( 7,dg,7 );
    writeLineLoop( 8, dg, Loop()>>-2>>-13>>21 );
    writeRuledSurface( 8, dg,8 );
    writeSurfaceLoop( 1, dg, Loop()>>3>>2>>6>>7>>8>>5>>1>>4 );
#endif

    writeVolume( 1, dg, 1 );
}



void
runTube( data_geo_ptrtype dg )
{

    node_type centre = param<0>( dg );
    node_type direction = param<1>( dg );
    double rayon = param<2>( dg )( 0 );
    double longueur=param<3>( dg )( 0 );
    double epaisseur=param<4>( dg )( 0 );

    auto basis = computeBasisOrthogonal( direction,centre );
    Node dir( boost::get<0>( basis ) );
    Node u( boost::get<1>( basis ) );
    Node v( boost::get<2>( basis ) );

    Node geoX1( centre( 0 )+u( 0 )*rayon,
                centre( 1 )+u( 1 )*rayon,
                centre( 2 )+u( 2 )*rayon  );

    Node geoX2( centre( 0 )+v( 0 )*rayon,
                centre( 1 )+v( 1 )*rayon,
                centre( 2 )+v( 2 )*rayon  );

    Node geoX3( centre( 0 )-u( 0 )*rayon,
                centre( 1 )-u( 1 )*rayon,
                centre( 2 )-u( 2 )*rayon  );

    Node geoX4( centre( 0 )-v( 0 )*rayon,
                centre( 1 )-v( 1 )*rayon,
                centre( 2 )-v( 2 )*rayon  );

    Node centre2( centre( 0 )+longueur*dir( 0 ),
                  centre( 1 )+longueur*dir( 1 ),
                  centre( 2 )+longueur*dir( 2 ) );

    Node geoX5( geoX1( 0 )+longueur*dir( 0 ),
                geoX1( 1 )+longueur*dir( 1 ),
                geoX1( 2 )+longueur*dir( 2 ) );

    Node geoX6( geoX2( 0 )+longueur*dir( 0 ),
                geoX2( 1 )+longueur*dir( 1 ),
                geoX2( 2 )+longueur*dir( 2 ) );

    Node geoX7( geoX3( 0 )+longueur*dir( 0 ),
                geoX3( 1 )+longueur*dir( 1 ),
                geoX3( 2 )+longueur*dir( 2 ) );

    Node geoX8( geoX4( 0 )+longueur*dir( 0 ),
                geoX4( 1 )+longueur*dir( 1 ),
                geoX4( 2 )+longueur*dir( 2 ) );


    Node geoX1B( centre( 0 )+u( 0 )*( rayon+epaisseur ),
                 centre( 1 )+u( 1 )*( rayon+epaisseur ),
                 centre( 2 )+u( 2 )*( rayon+epaisseur )  );

    Node geoX2B( centre( 0 )+v( 0 )*( rayon+epaisseur ),
                 centre( 1 )+v( 1 )*( rayon+epaisseur ),
                 centre( 2 )+v( 2 )*( rayon+epaisseur ) );

    Node geoX3B( centre( 0 )-u( 0 )*( rayon+epaisseur ),
                 centre( 1 )-u( 1 )*( rayon+epaisseur ),
                 centre( 2 )-u( 2 )*( rayon+epaisseur ) );

    Node geoX4B( centre( 0 )-v( 0 )*( rayon+epaisseur ),
                 centre( 1 )-v( 1 )*( rayon+epaisseur ),
                 centre( 2 )-v( 2 )*( rayon+epaisseur ) );

    Node geoX5B( geoX1B( 0 )+longueur*dir( 0 ),
                 geoX1B( 1 )+longueur*dir( 1 ),
                 geoX1B( 2 )+longueur*dir( 2 ) );

    Node geoX6B( geoX2B( 0 )+longueur*dir( 0 ),
                 geoX2B( 1 )+longueur*dir( 1 ),
                 geoX2B( 2 )+longueur*dir( 2 ) );

    Node geoX7B( geoX3B( 0 )+longueur*dir( 0 ),
                 geoX3B( 1 )+longueur*dir( 1 ),
                 geoX3B( 2 )+longueur*dir( 2 ) );

    Node geoX8B( geoX4B( 0 )+longueur*dir( 0 ),
                 geoX4B( 1 )+longueur*dir( 1 ),
                 geoX4B( 2 )+longueur*dir( 2 ) );

    //--------------------------------------------------------------------------//

    writePoint( 1, dg, centre( 0 ), centre( 1 ) ,centre( 2 ) );
    writePoint( 2, dg, geoX1( 0 ), geoX1( 1 ) ,geoX1( 2 ) );
    writePoint( 3, dg, geoX2( 0 ), geoX2( 1 ) ,geoX2( 2 ) );
    writePoint( 4, dg, geoX3( 0 ), geoX3( 1 ) ,geoX3( 2 ) );
    writePoint( 5, dg, geoX4( 0 ), geoX4( 1 ) ,geoX4( 2 ) );

    writeCircle( 1, dg, 2,1,3 );
    writeCircle( 2, dg, 3,1,4 );
    writeCircle( 3, dg, 4,1,5 );
    writeCircle( 4, dg, 5,1,2 );

    writeLineLoop( 1, dg, Loop()>>1>>2>>3>>4 );

    writePoint( 6, dg, centre2( 0 ), centre2( 1 ) ,centre2( 2 ) );
    writePoint( 7, dg, geoX5( 0 ), geoX5( 1 ) ,geoX5( 2 ) );
    writePoint( 8, dg, geoX6( 0 ), geoX6( 1 ) ,geoX6( 2 ) );
    writePoint( 9, dg, geoX7( 0 ), geoX7( 1 ) ,geoX7( 2 ) );
    writePoint( 10, dg, geoX8( 0 ), geoX8( 1 ) ,geoX8( 2 ) );

    writeCircle( 5, dg, 7,6,8 );
    writeCircle( 6, dg,8,6,9 );
    writeCircle( 7, dg,9,6,10 );
    writeCircle( 8 ,dg,10,6,7 );

    writeLineLoop( 2, dg,Loop()>>5>>6>>7>>8 );

    writeLine( 9, dg, 4, 9 );
    writeLine( 10, dg, 5, 10 );
    writeLine( 11, dg, 2, 7 );
    writeLine( 12, dg, 8, 3 );

    writeLineLoop( 3, dg, Loop()>>9>>-6>>12>>2 );
    writeLineLoop( 4, dg, Loop()>>9>>7>>-10>>-3 );
    writeLineLoop( 5, dg, Loop()>>10>>8>>-11>>-4 );
    writeLineLoop( 6, dg, Loop()>>11>>5>>12>>-1 );

    //writePlaneSurface(1, dg, 1);
    //writePlaneSurface(2, dg, 2);

    writeRuledSurface( 3, dg, 3 );
    writeRuledSurface( 4, dg, 4 );
    writeRuledSurface( 5, dg, 5 );
    writeRuledSurface( 6, dg, 6 );


    writePoint( 11, dg, geoX1B( 0 ), geoX1B( 1 ) ,geoX1B( 2 ) );
    writePoint( 12, dg, geoX2B( 0 ), geoX2B( 1 ) ,geoX2B( 2 ) );
    writePoint( 13, dg, geoX3B( 0 ), geoX3B( 1 ) ,geoX3B( 2 ) );
    writePoint( 14, dg, geoX4B( 0 ), geoX4B( 1 ) ,geoX4B( 2 ) );

    writeCircle( 13, dg, 11,1,12 );
    writeCircle( 14, dg, 12,1,13 );
    writeCircle( 15, dg, 13,1,14 );
    writeCircle( 16, dg, 14,1,11 );

    writePoint( 15, dg, geoX5B( 0 ), geoX5B( 1 ) ,geoX5B( 2 ) );
    writePoint( 16, dg, geoX6B( 0 ), geoX6B( 1 ) ,geoX6B( 2 ) );
    writePoint( 17, dg, geoX7B( 0 ), geoX7B( 1 ) ,geoX7B( 2 ) );
    writePoint( 18, dg, geoX8B( 0 ), geoX8B( 1 ) ,geoX8B( 2 ) );

    writeCircle( 17, dg, 15,6,16 );
    writeCircle( 18, dg, 16,6,17 );
    writeCircle( 19, dg, 17,6,18 );
    writeCircle( 20 ,dg, 18,6,15 );

    writeLine( 21, dg, 13, 17 );
    writeLine( 22, dg, 14, 18 );
    writeLine( 23, dg, 11, 15 );
    writeLine( 24, dg, 16, 12 );

    writeLineLoop( 7, dg, Loop()>>21>>-18>>24>>14 );
    writeLineLoop( 8, dg, Loop()>>21>>19>>-22>>-15 );
    writeLineLoop( 9, dg, Loop()>>22>>20>>-23>>-16 );
    writeLineLoop( 10, dg, Loop()>>23>>17>>24>>-13 );

    writeRuledSurface( 7, dg, 7 );
    writeRuledSurface( 8, dg, 8 );
    writeRuledSurface( 9, dg, 9 );
    writeRuledSurface( 10, dg, 10 );

    writeLine( 25, dg, 2, 11 );
    writeLine( 26, dg, 3, 12 );
    writeLine( 27, dg, 4, 13 );
    writeLine( 28, dg, 5, 14 );
    writeLine( 29, dg, 7, 15 );
    writeLine( 30, dg, 8, 16 );
    writeLine( 31, dg, 9, 17 );
    writeLine( 32, dg, 10, 18 );


    writeLineLoop( 11, dg, Loop()>>19>>-32>>-7>>31 );
    writeLineLoop( 12, dg, Loop()>>20>>-29>>-8>>32 );
    writeLineLoop( 13, dg, Loop()>>6>>31>>-18>>-30 );
    writeLineLoop( 14, dg, Loop()>>5>>30>>-17>>-29 );
    writeLineLoop( 15, dg, Loop()>>15>>-28>>-3>>27 );
    writeLineLoop( 16, dg, Loop()>>27>>-14>>-26>>2 );
    writeLineLoop( 17, dg, Loop()>>26>>-13>>-25>>1 );
    writeLineLoop( 18, dg, Loop()>>25>>-16>>-28>>4 );

    // internal surface
    writeLineLoop( 19 ,dg, Loop()>>23>>-29>>-11>>25 );
    writeLineLoop( 20 ,dg, Loop()>>10>>32>>-22>>-28 );
    writeLineLoop( 21 ,dg, Loop()>>9>>31>>-21>>-27 );
    writeLineLoop( 22 ,dg, Loop()>>24>>-26>>-12>>30 );

    // internal surface
    writePlaneSurface( 19,dg,19 );
    writeRuledSurface( 20,dg,20 );
    writePlaneSurface( 21,dg,21 );
    writePlaneSurface( 22,dg,22 );

    // Inlet or outlet?
    writeRuledSurface( 11,dg,11 );
    writeRuledSurface( 12,dg,12 );
    writeRuledSurface( 13,dg,13 );
    writeRuledSurface( 14,dg,14 );

    // Inlet or outlet?
    writeRuledSurface( 15,dg,15 );
    writeRuledSurface( 16,dg,16 );
    writeRuledSurface( 17,dg,17 );
    writeRuledSurface( 18,dg,18 );



    //writeSurfaceLoop(1, dg, Loop()>>3>>4>>5>>6>>7>>8>>9>>10>>11>>12>>13>>14>>15>>16>>17>>18);
    //writeSurfaceLoop(1, dg, Loop()>>3>>4>>5>>6>>7>>8>>9>>10>>11>>12>>13>>14>>15>>16>>17>>18>>19>>20>>21>>22);
    //writeVolume(1, dg, 1);
#if 0
    writeSurfaceLoop( 1, dg, Loop()>>6>>9>>2>>13>>18>>19 );
    writeVolume( 1,dg,1 );
    writeSurfaceLoop( 2, dg, Loop()>>18>>7>>10>>3>>16>>17 );
    writeVolume( 2,dg,2 );
    writeSurfaceLoop( 37, dg, Loop()>>20>>12>>4>>15>>8>>17 );
    writeVolume( 3,dg,3 );
    writeSurfaceLoop( 4, dg, Loop()>>11>>1>>14>>5>>19>>20 );
    writeVolume( 4,dg,4 );
#else
    writeSurfaceLoop( 1, dg, Loop()>>4>>16>>8>>19>>12>>9 );
    writeVolume( 1,dg,1 );
    writeSurfaceLoop( 2, dg, Loop()>>5>>15>>1>>18>>11>>12 );
    writeVolume( 2,dg,2 );
    writeSurfaceLoop( 3, dg, Loop()>>6>>13>>2>>17>>11>>10 );
    writeVolume( 3,dg,3 );
    writeSurfaceLoop( 4, dg, Loop()>>10>>7>>14>>3>>20>>9 );
    writeVolume( 4,dg,4 );
#endif


} // runTube


void
runSpecial3D_1( data_geo_ptrtype dg )
{
    double lgstruct=0.35101;
    double xL = 0.6-lgstruct;
    double yMin = -0.12,yMax=0.12;
    writePoint(1,dg,  xL, yMin, 0.19 );
    writePoint(2,dg, 0.6, yMin, 0.19 );
    writePoint(3,dg, 0.6, yMin, 0.21 );
    writePoint(4,dg,  xL, yMin, 0.21 );
    writePoint(5,dg,  xL, yMax, 0.19 );
    writePoint(6,dg, 0.6, yMax, 0.19 );
    writePoint(7,dg, 0.6, yMax, 0.21 );
    writePoint(8,dg,  xL, yMax, 0.21 );
    // point sup
    writePoint(9,dg,  0.2     , yMin, 0.2 );// center
    writePoint(10,dg, 0.2-0.05, yMin, 0.2 ); // on circle
    writePoint(11,dg, 0.2     , yMax, 0.2 );// center
    writePoint(12,dg, 0.2-0.05, yMax, 0.2); // on circle

    writeLine(1,dg, 1,2);
    writeLine(2,dg, 2,3);
    writeLine(3,dg, 3,4);
    writeLine(4,dg, 4,1);
    writeLine(5,dg, 5,6);
    writeLine(6,dg, 6,7);
    writeLine(7,dg, 7,8);
    writeLine(8,dg, 8,5);
    writeLine(9,dg, 1,5);
    writeLine(10,dg,2,6);
    writeLine(11,dg,3,7);
    writeLine(12,dg,4,8);
    // line sup
    writeCircle(13,dg, 1,9,10 );
    writeCircle(14,dg, 10,9,4 );
    writeCircle(15,dg, 5,11,12 );
    writeCircle(16,dg, 12,11,8 );
    // line on cylinder
    writeLine(17,dg, 12, 10);

    writeLineLoop(1,dg, Loop()>>1>>2>>3>>4);
    writeLineLoop(2,dg, Loop()>>-8>>-7>>-6>>-5);
    writeLineLoop(3,dg, Loop()>>9>>5>>-10>>-1);
    writeLineLoop(4,dg, Loop()>>10>>6>>-11>>-2);
    writeLineLoop(5,dg, Loop()>>11>>7>>-12>>-3);
    //writeLineLoop(6,dg, Loop()>>9>>-8>>-12>>4);
    writePlaneSurface(1,dg,1);
    writePlaneSurface(2,dg,2);
    writePlaneSurface(3,dg,3);
    writePlaneSurface(4,dg,4);
    writePlaneSurface(5,dg,5);
    //writePlaneSurface(6,dg,6);

    writeLineLoop(6,dg,  Loop()>>17>>14>>12>>-16 );
    writeLineLoop(7,dg,  Loop()>>-15>>-9>>13>>-17);
    writeLineLoop(8,dg,  Loop()>>-13>>-4>>-14);
    writeLineLoop(9,dg, Loop()>>15>>16>>8);
    writeRuledSurface(6,dg,6);
    writeRuledSurface(7,dg,7);
    writePlaneSurface(8,dg,8);
    writePlaneSurface(9,dg,9);
    writePtInSurface(dg,9,8 );
    writePtInSurface(dg,11,9 );

#if 0
    writeSurfaceLoop( 1, dg, Loop()>>7>>10>>8>>9>>6 );
    writeVolume( 1,dg,1 );
    writeSurfaceLoop( 2, dg, Loop()>>5>>4>>3>>1>>2>>6 );
    writeVolume( 2,dg,2 );
#else
    writeSurfaceLoop( 1, dg, Loop()>>6>>9>>7>>8>>5>>4>>3>>1>>2 );
    writeVolume( 1,dg,1 );
#endif

} // runSpecial3D_1


} //namespace GeoTool

} //namespace Feel

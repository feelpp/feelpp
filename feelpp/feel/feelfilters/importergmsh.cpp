/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 28 Dec 2016

 Copyright (C) 2016 Feel++ Consortium

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
#include <feel/feelfilters/importergmsh.hpp>

namespace Feel
{
#if defined(FEELPP_HAS_GMSH_H)
#if !defined( FEELPP_HAS_GMSH_API )
// EDF MED format
const auto med = GmshReaderFactory::instance().emplace( ".med",
                                                        []( std::string const& fname,std::string const& gname )
                                                        {
                                                            GModel *m = new GModel( gname );
                                                            int status = 1;
#ifdef FEELPP_HAS_GMSH_HAS_MED
                                                            status = m->readMED(fname);
#else
                                                            throw std::logic_error("Gmsh MED support is not available. Cannot load MED file");
#endif
                                                            return status;
                                                        });
// Medit inria format
const auto mesh = GmshReaderFactory::instance().emplace( ".mesh",
                                                         []( std::string const& fname,std::string const& gname )
                                                         {
                                                             GModel *m = new GModel( gname );
                                                             int status = m->readMESH(fname);
                                                             return status;
                                                         });
// Nastran Bulk Data File format
const auto reader_bdf = GmshReaderFactory::instance().emplace( ".bdf",
                                                               []( std::string const& fname,std::string const& gname )
                                                               {
                                                                   GModel *m = new GModel( gname );
                                                                   int status = m->readBDF(fname);
                                                                   return status;
                                                               });
// Plot3D files
const auto reader_p3d = GmshReaderFactory::instance().emplace( ".p3d",
                                                               []( std::string const& fname,std::string const& gname )
                                                               {
                                                                   GModel *m = new GModel( gname );
                                                                   int status = m->readP3D(fname);
                                                                   return status;
                                                               });
// CFD General Notation System files
const auto reader_cgns = GmshReaderFactory::instance().emplace( ".cgns",
                                                                []( std::string const& fname,std::string const& gname )
                                                                {
                                                                    GModel *m = new GModel( gname );
                                                                    int status = 1;
#ifdef FEELPP_HAS_GMSH_HAS_CGNS
                                                                    status = m->readCGNS(fname);
#else
                                                                    throw std::logic_error("Gmsh CGNS(HDF) support is not available. Cannot load CGNS file");
#endif
                                                                    return status;
                                                                });
#endif

#endif // FEELPP_HAS_GMSH_H
namespace detail
{
bool isOnProcessor( std::vector<rank_type> ghosts, rank_type partition, rank_type worldcommrank, rank_type worldcommsize)
{
    // maybe proc id does not start at 0
    std::for_each( ghosts.begin(), ghosts.end(), [&worldcommsize]( rank_type& g ) { g = (g % worldcommsize); } );
    
    if ( worldcommsize == 1 )
        return true;

    if ( worldcommrank == partition )
        return true;
    
    // is the element a ghost cell
    // look into ghosts if 'partition' is present
    auto it = std::find( ghosts.begin(), ghosts.end(), worldcommrank );
    if ( it != ghosts.end() )
        return true;
    return false;
}

} // detail


// From Gmsh - Common/StringUtils.h
void SwapBytes(char *array, int size, int n)
{
    char *x = new char[size];
    for(int i = 0; i < n; i++) {
        char *a = &array[i * size];
        memcpy(x, a, size);
        for(int c = 0; c < size; c++)
            a[size - 1 - c] = x[c];
    }
    delete [] x;
}

} // Feel




#if defined( FEELPP_HAS_GMSH_API )
int getInfoMSH(const int typeMSH, std::string & elementName)
{
    int dim, order, numNodes;
    std::vector<double> parametricCoord;
#if GMSH_VERSION_LESS_THAN( 4,5,0 )
    gmsh::model::mesh::getElementProperties( typeMSH,elementName,dim,order,numNodes,parametricCoord );
#else
    int numPrimaryNodes;
    gmsh::model::mesh::getElementProperties( typeMSH,elementName,dim,order,numNodes,parametricCoord,numPrimaryNodes );
#endif
    return numNodes;
}
#else
int getInfoMSH(const int typeMSH, const char **const name)
{
    CHECK(typeMSH != MSH_POLYG_ && typeMSH != MSH_POLYH_ && typeMSH != MSH_POLYG_B)
        << "GMSH Element type " << typeMSH << " not supported by Feel++\n";
    CHECK( typeMSH < MSH_NUM_TYPE ) << "Invalid GMSH element type " << typeMSH << "\n";

#if defined( FEELPP_HAS_GMSH_H)
    return MElement::getInfoMSH(typeMSH, name);
#else
    // from Gmsh - Geo/MElement.h
    switch(typeMSH){
    case MSH_PNT     : if(name) *name = "Point";            return 1;
    case MSH_LIN_1   : if(name) *name = "Line 1";           return 1;
    case MSH_LIN_2   : if(name) *name = "Line 2";           return 2;
    case MSH_LIN_3   : if(name) *name = "Line 3";           return 2 + 1;
    case MSH_LIN_4   : if(name) *name = "Line 4";           return 2 + 2;
    case MSH_LIN_5   : if(name) *name = "Line 5";           return 2 + 3;
    case MSH_LIN_6   : if(name) *name = "Line 6";           return 2 + 4;
    case MSH_LIN_7   : if(name) *name = "Line 7";           return 2 + 5;
    case MSH_LIN_8   : if(name) *name = "Line 8";           return 2 + 6;
    case MSH_LIN_9   : if(name) *name = "Line 9";           return 2 + 7;
    case MSH_LIN_10  : if(name) *name = "Line 10";          return 2 + 8;
    case MSH_LIN_11  : if(name) *name = "Line 11";          return 2 + 9;
    case MSH_LIN_B   : if(name) *name = "Line Border";      return 2;
    case MSH_LIN_C   : if(name) *name = "Line Child";       return 2;
    case MSH_TRI_1   : if(name) *name = "Triangle 1";       return 1;
    case MSH_TRI_3   : if(name) *name = "Triangle 3";       return 3;
    case MSH_TRI_6   : if(name) *name = "Triangle 6";       return 3 + 3;
    case MSH_TRI_9   : if(name) *name = "Triangle 9";       return 3 + 6;
    case MSH_TRI_10  : if(name) *name = "Triangle 10";      return 3 + 6 + 1;
    case MSH_TRI_12  : if(name) *name = "Triangle 12";      return 3 + 9;
    case MSH_TRI_15  : if(name) *name = "Triangle 15";      return 3 + 9 + 3;
    case MSH_TRI_15I : if(name) *name = "Triangle 15I";     return 3 + 12;
    case MSH_TRI_21  : if(name) *name = "Triangle 21";      return 3 + 12 + 6;
    case MSH_TRI_28  : if(name) *name = "Triangle 28";      return 3 + 15 + 10;
    case MSH_TRI_36  : if(name) *name = "Triangle 36";      return 3 + 18 + 15;
    case MSH_TRI_45  : if(name) *name = "Triangle 45";      return 3 + 21 + 21;
    case MSH_TRI_55  : if(name) *name = "Triangle 55";      return 3 + 24 + 28;
    case MSH_TRI_66  : if(name) *name = "Triangle 66";      return 3 + 27 + 36;
    case MSH_TRI_18  : if(name) *name = "Triangle 18";      return 3 + 15;
    case MSH_TRI_21I : if(name) *name = "Triangle 21I";     return 3 + 18;
    case MSH_TRI_24  : if(name) *name = "Triangle 24";      return 3 + 21;
    case MSH_TRI_27  : if(name) *name = "Triangle 27";      return 3 + 24;
    case MSH_TRI_30  : if(name) *name = "Triangle 30";      return 3 + 27;
    case MSH_TRI_B   : if(name) *name = "Triangle Border";  return 3;
    case MSH_QUA_1   : if(name) *name = "Quadrilateral 1";  return 1;
    case MSH_QUA_4   : if(name) *name = "Quadrilateral 4";  return 4;
    case MSH_QUA_8   : if(name) *name = "Quadrilateral 8";  return 4 + 4;
    case MSH_QUA_9   : if(name) *name = "Quadrilateral 9";  return 9;
    case MSH_QUA_16  : if(name) *name = "Quadrilateral 16"; return 16;
    case MSH_QUA_25  : if(name) *name = "Quadrilateral 25"; return 25;
    case MSH_QUA_36  : if(name) *name = "Quadrilateral 36"; return 36;
    case MSH_QUA_49  : if(name) *name = "Quadrilateral 49"; return 49;
    case MSH_QUA_64  : if(name) *name = "Quadrilateral 64"; return 64;
    case MSH_QUA_81  : if(name) *name = "Quadrilateral 81"; return 81;
    case MSH_QUA_100 : if(name) *name = "Quadrilateral 100";return 100;
    case MSH_QUA_121 : if(name) *name = "Quadrilateral 121";return 121;
    case MSH_QUA_12  : if(name) *name = "Quadrilateral 12"; return 12;
    case MSH_QUA_16I : if(name) *name = "Quadrilateral 16I";return 16;
    case MSH_QUA_20  : if(name) *name = "Quadrilateral 20"; return 20;
    case MSH_QUA_24  : if(name) *name = "Quadrilateral 24"; return 24;
    case MSH_QUA_28  : if(name) *name = "Quadrilateral 28"; return 28;
    case MSH_QUA_32  : if(name) *name = "Quadrilateral 32"; return 32;
    case MSH_QUA_36I : if(name) *name = "Quadrilateral 36I";return 36;
    case MSH_QUA_40  : if(name) *name = "Quadrilateral 40"; return 40;
    case MSH_POLYG_  : if(name) *name = "Polygon";          return 0;
    case MSH_POLYG_B : if(name) *name = "Polygon Border";   return 0;
    case MSH_TET_1   : if(name) *name = "Tetrahedron 1";    return 1;
    case MSH_TET_4   : if(name) *name = "Tetrahedron 4";    return 4;
    case MSH_TET_10  : if(name) *name = "Tetrahedron 10";   return 4 + 6;
    case MSH_TET_20  : if(name) *name = "Tetrahedron 20";   return 4 + 12 + 4;
    case MSH_TET_35  : if(name) *name = "Tetrahedron 35";   return 4 + 18 + 12 + 1;
    case MSH_TET_56  : if(name) *name = "Tetrahedron 56";   return 4 + 24 + 24 + 4;
    case MSH_TET_84  : if(name) *name = "Tetrahedron 84";   return (7*8*9)/6;
    case MSH_TET_120 : if(name) *name = "Tetrahedron 120";  return (8*9*10)/6;
    case MSH_TET_165 : if(name) *name = "Tetrahedron 165";  return (9*10*11)/6;
    case MSH_TET_220 : if(name) *name = "Tetrahedron 220";  return (10*11*12)/6;
    case MSH_TET_286 : if(name) *name = "Tetrahedron 286";  return (11*12*13)/6;
    case MSH_TET_16  : if(name) *name = "Tetrahedron 16";   return 4 + 6*2;
    case MSH_TET_22  : if(name) *name = "Tetrahedron 22";   return 4 + 6*3;
    case MSH_TET_28  : if(name) *name = "Tetrahedron 28";   return 4 + 6*4;
    case MSH_TET_34  : if(name) *name = "Tetrahedron 34";   return 4 + 6*5;
    case MSH_TET_40  : if(name) *name = "Tetrahedron 40";   return 4 + 6*6;
    case MSH_TET_46  : if(name) *name = "Tetrahedron 46";   return 4 + 6*7;
    case MSH_TET_52  : if(name) *name = "Tetrahedron 52";   return 4 + 6*8;
    case MSH_TET_58  : if(name) *name = "Tetrahedron 58";   return 4 + 6*9;
    case MSH_HEX_1   : if(name) *name = "Hexahedron 1";     return 1;
    case MSH_HEX_8   : if(name) *name = "Hexahedron 8";     return 8;
    case MSH_HEX_20  : if(name) *name = "Hexahedron 20";    return 8 + 12;
    case MSH_HEX_27  : if(name) *name = "Hexahedron 27";    return 8 + 12 + 6 + 1;
    case MSH_HEX_64  : if(name) *name = "Hexahedron 64";    return 64;
    case MSH_HEX_125 : if(name) *name = "Hexahedron 125";   return 125;
    case MSH_HEX_216 : if(name) *name = "Hexahedron 216";   return 216;
    case MSH_HEX_343 : if(name) *name = "Hexahedron 343";   return 343;
    case MSH_HEX_512 : if(name) *name = "Hexahedron 512";   return 512;
    case MSH_HEX_729 : if(name) *name = "Hexahedron 729";   return 729;
    case MSH_HEX_1000: if(name) *name = "Hexahedron 1000";  return 1000;
    case MSH_HEX_32  : if(name) *name = "Hexahedron 32";    return 8 + 12*2;
    case MSH_HEX_44  : if(name) *name = "Hexahedron 44";    return 8 + 12*3;
    case MSH_HEX_56  : if(name) *name = "Hexahedron 56";    return 8 + 12*4;
    case MSH_HEX_68  : if(name) *name = "Hexahedron 68";    return 8 + 12*5;
    case MSH_HEX_80  : if(name) *name = "Hexahedron 80";    return 8 + 12*6;
    case MSH_HEX_92  : if(name) *name = "Hexahedron 92";    return 8 + 12*7;
    case MSH_HEX_104 : if(name) *name = "Hexahedron 104";   return 8 + 12*8;
    case MSH_PRI_1   : if(name) *name = "Prism 1";          return 1;
    case MSH_PRI_6   : if(name) *name = "Prism 6";          return 6;
    case MSH_PRI_15  : if(name) *name = "Prism 15";         return 6 + 9;
    case MSH_PRI_18  : if(name) *name = "Prism 18";         return 6 + 9 + 3;
    case MSH_PRI_40  : if(name) *name = "Prism 40";         return 6 + 18 + 12+2 + 2*1;
    case MSH_PRI_75  : if(name) *name = "Prism 75";         return 6 + 27 + 27+6 + 3*3;
    case MSH_PRI_126 : if(name) *name = "Prism 126";        return 6 + 36 + 48+12 + 4*6;
    case MSH_PRI_196 : if(name) *name = "Prism 196";        return 6 + 45 + 75+20 + 5*10;
    case MSH_PRI_288 : if(name) *name = "Prism 288";        return 6 + 54 + 108+30 + 6*15;
    case MSH_PRI_405 : if(name) *name = "Prism 405";        return 6 + 63 + 147+42 + 7*21;
    case MSH_PRI_550 : if(name) *name = "Prism 550";        return 6 + 72 + 192+56 + 8*28;
    case MSH_PRI_24  : if(name) *name = "Prism 24";         return 6 + 9*2;
    case MSH_PRI_33  : if(name) *name = "Prism 33";         return 6 + 9*3;
    case MSH_PRI_42  : if(name) *name = "Prism 42";         return 6 + 9*4;
    case MSH_PRI_51  : if(name) *name = "Prism 51";         return 6 + 9*5;
    case MSH_PRI_60  : if(name) *name = "Prism 60";         return 6 + 9*6;
    case MSH_PRI_69  : if(name) *name = "Prism 69";         return 6 + 9*7;
    case MSH_PRI_78  : if(name) *name = "Prism 78";         return 6 + 9*8;
    case MSH_PYR_1   : if(name) *name = "Pyramid 1";        return 1;
    case MSH_PYR_5   : if(name) *name = "Pyramid 5";        return 5;
    case MSH_PYR_13  : if(name) *name = "Pyramid 13";       return 5 + 8;
    case MSH_PYR_14  : if(name) *name = "Pyramid 14";       return 5 + 8 + 1;
    case MSH_PYR_30  : if(name) *name = "Pyramid 30";       return 5 + 8*2 + 4*1  + 1*4  + 1;
    case MSH_PYR_55  : if(name) *name = "Pyramid 55";       return 5 + 8*3 + 4*3  + 1*9  + 5;
    case MSH_PYR_91  : if(name) *name = "Pyramid 91";       return 5 + 8*4 + 4*6  + 1*16 + 14;
    case MSH_PYR_140 : if(name) *name = "Pyramid 140";      return 5 + 8*5 + 4*10 + 1*25 + 30;
    case MSH_PYR_204 : if(name) *name = "Pyramid 204";      return 5 + 8*6 + 4*15 + 1*36 + 55;
    case MSH_PYR_285 : if(name) *name = "Pyramid 285";      return 5 + 8*7 + 4*21 + 1*49 + 91;
    case MSH_PYR_385 : if(name) *name = "Pyramid 385";      return 5 + 8*8 + 4*28 + 1*64 + 140;
    case MSH_PYR_21  : if(name) *name = "Pyramid 21";       return 5 + 8*2;
    case MSH_PYR_29  : if(name) *name = "Pyramid 29";       return 5 + 8*3;
    case MSH_PYR_37  : if(name) *name = "Pyramid 37";       return 5 + 8*4;
    case MSH_PYR_45  : if(name) *name = "Pyramid 45";       return 5 + 8*5;
    case MSH_PYR_53  : if(name) *name = "Pyramid 53";       return 5 + 8*6;
    case MSH_PYR_61  : if(name) *name = "Pyramid 61";       return 5 + 8*7;
    case MSH_PYR_69  : if(name) *name = "Pyramid 69";       return 5 + 8*8;
    case MSH_TRIH_4 : if(name) *name = "Trihedron 4";       return 4;
    case MSH_POLYH_  : if(name) *name = "Polyhedron";       return 0;
    case MSH_PNT_SUB : if(name) *name = "Point Xfem";       return 1;
    case MSH_LIN_SUB : if(name) *name = "Line Xfem";        return 2;
    case MSH_TRI_SUB : if(name) *name = "Triangle Xfem";    return 3;
    case MSH_TET_SUB : if(name) *name = "Tetrahedron Xfem"; return 4;
    default:
        if(name) *name = "Unknown";
        return 0;
    }

#endif // FEELPP_HAS_GMSH_H
}

#endif // !FEELPP_HAS_GMSH_API

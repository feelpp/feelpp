/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-06-30

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
   \file gambit.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-06-30
 */
#include <iostream>
#include <fstream>
#include <iomanip>


#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>


#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

#include <feel/feelcore/feel.hpp>

#include <feel/feeldiscr/mesh.hpp>

#include <feel/feelfilters/importergambit.hpp>


namespace Feel
{
namespace gambit
{
std::istream & eatline( std::istream & s )    // eat a whole line from std::istream
{
    while ( s.get() != '\n' && s . good() )
    {}

    return s ;
}

std::istream & eat_comments( std::istream & s )    //eat lines starting with '!%#;$'
{
    char c = 'a';
    s.get( c ) ;

    while ( c == '!' ||
            c == '%' ||
            c == '#' ||
            c == ';' ||
            c == '$' )
    {
        s >> eatline ;
        s.get( c ) ;
    }

    return s.putback( c ) ;
}

std::istream & next_good_line( std::istream & s, std::string & line )
{
    s >> eat_comments;
    getline( s, line );
    return s;
}

const int quad4::face[][2] = { {0, 1},
    {1, 2},
    {2, 3},
    {3, 0}
};

const int quad8::face[][3] = { {0, 1, 2},
    {2, 3, 4},
    {4, 5, 6},
    {6, 7, 0}
};
const int tria3::face[][2] = { {0, 1},
    {1, 2},
    {2, 0}
};
const int tria6::face[][3] = { {0, 1, 2},
    {2, 3, 4},
    {4, 5, 0}
};

const int tetra4::face[][3] = { {1, 0, 2},
    {0, 1, 3},
    {1, 2, 3},
    {2, 0, 3}
};
const int tetra10::face[][6] =  { {2, 1, 0, 3, 5, 4},
    {0, 1, 2, 7, 9, 6},
    {2, 4, 5, 8, 9, 7},
    {5, 3, 0, 6, 9, 8}
};
bool
read( std::string const& filename,
      nodes_type& nodes,
      nodes_boundary_type& boundary,
      elements_type& elements )
{
    Debug( 8012 ) << "Reading Gambit mesh file   (" << filename << ")"
                  << ":" << "\n";


    std::ifstream file( filename.c_str() );

    if ( file.fail() )
    {
        Debug( 8012 ) << "Reading mesh file " << filename
                      << " impossible" << "\n";
        throw std::logic_error( std::string( "cannot open gambit mesh file: " ) + filename );
    }

    std::string sdummy;
    std::string line;

    int         numnp;
    int         nelem;
    int         ngrps;
    int         nbsets;
    int         ndfcd;
    int         ndfvl;

    std::map<std::string, int> groups;
    int thegroupid = 1;

    while ( next_good_line( file, line ).good() )
    {
        if ( line.find( "GAMBIT NEUTRAL FILE" ) != std::string::npos )
        {
            // Reading header file
            file >> sdummy;
            Debug( 8012 ) << "Mesh name        : " << sdummy << "\n";
            file >> sdummy >> sdummy >> sdummy >> sdummy;
            Debug( 8012 ) << "Version          : " << sdummy << "\n";
            file >> sdummy >> sdummy >> sdummy >> sdummy >> sdummy;
            file >> sdummy >> sdummy >> sdummy >> sdummy >> sdummy;
            file >> numnp >> nelem >> ngrps >> nbsets >> ndfcd >> ndfvl;
            Debug( 8012 ) << "  Number of nodes         : " << numnp << "\n";
            Debug( 8012 ) << "  Number of elements      : " << nelem << "\n";
            Debug( 8012 ) << "  Number of elem groups   : " << ngrps << "\n";
            Debug( 8012 ) << "  Number of BC sets       : " << nbsets << "\n";
            Debug( 8012 ) << "  Number of directions    : " << ndfcd << "\n";
            Debug( 8012 ) << "  Number of velocity comp : " << ndfvl<< "\n";

            if ( ndfcd == 3 )
                nodes.resize( 3*numnp );

            else
                nodes.resize( 2*numnp );

            boundary.resize( numnp );

            elements.resize( nelem );
            file >> sdummy;
        }

        if ( line.find( "NODAL" ) != std::string::npos )
        {
            //            Debug( 8012 ) << sdummy << "\n";
            Debug( 8012 ) << "Reading nodes coordinates ... " << "\n";

            int k = 0;


            for ( int inode = 0; inode < numnp; ++inode )
            {
                file >> k;
                k--;

                if ( ndfcd == 3 ) // 3D
                    file >> nodes[3*k+0] >> nodes[3*k+1] >> nodes[3*k+2];

                else // 2D
                    file >> nodes[2*k+0] >> nodes[2*k+1];

                boundary[k].get<0>() = false;
                boundary[k].get<1>() = 0;
            }

            Debug( 8012 ) << "done." << "\n";
        }

        if ( line.find( "ELEMENTS/CELLS" ) != std::string::npos )
        {
            /**
               NE      Global element number (not required to be sequential or
               continuous)
               NTYPE   Element geometry type: 1 = Edge
               2 = Quadrilateral
               3 = Triangle
               4 = Brick
               5 = Wedge (Prism)
               6 = Tetrahedron
               7 = Pyramid
               NDP   Number of nodes that define the element
               NODE  List of nodes that define the element
            */
            Debug( 8012 ) << "Reading connectivity      ... " << "\n";

            for ( int ielem = 0; ielem < nelem; ++ielem )
            {
                int k, ntype, ndp;

                // NE
                file >> k;
                k--;

                // NTYPE
                file >> ntype;
                elements[k].get<0>() = ntype;
                elements[k].get<1>() = 0;

                // NDP
                file >> ndp;
                elements[k].get<2>().resize( ndp );

                // NODES
                for ( int i = 0; i < ndp; ++i )
                {
                    int inode;
                    file >> inode;
                    --inode;
                    boost::get<2>( elements[k] )[i] = inode;
                }

                elements[k].get<3>() = boost::make_tuple( -1, -1 );
            }

            Debug( 8012 ) << "Reading connectivity done." << "\n";
        }

        if ( line.find( "ELEMENT GROUP" ) != std::string::npos )
        {
            Debug( 8012 ) << "Reading ELEMENT GROUP...." << "\n";

            /*
              NGP    Element group number
              NELGP  Number of elements in group
              MTYP   Material type (NOTE: Interpretation of this flag is solver-
              dependent.)
              0 = Undefined
              1 = Conjugate
              2 = Fluid
              3 = Porous
              4 = Solid
              5 = Deformable
              NFLAGS Number of solver-dependent flags
            */
            std::string dummy;
            int ngroup, nelgp, mtyp, nflags ;
            std::string thegroup;
            file >> dummy >> ngroup
                 >> dummy >> nelgp
                 >> dummy >> mtyp
                 >> dummy >> nflags;
            file >> thegroup >> nflags;
            //groups[thegroup] = thegroupid++;
            Debug( 8012 ) << "group : " << ngroup
                          << " with id " << thegroupid
                          << " and  " << nelgp << " elements\n";

            for ( int i = 0; i < nelgp; ++i )
            {
                int nel;
                file >> nel;
                --nel;
                elements[nel].get<1>() = thegroupid;
            }

            Debug( 8012 ) << "Reading ELEMENT GROUP done" << "\n";
            thegroupid++;
        }

        if ( line.find( "BOUNDARY CONDITIONS" ) != std::string::npos )
        {
            /*
              NAME    Name of boundary-condition set
              ITYPE   Data type (0 = node; 1 = element/cell)
              NENTRY  Number of data records in boundary-condition set
              NVALUES Number of values for each data record
              IBCODE1 (Optional) Boundary condition code 1
              IBCODE2 (Optional) Boundary condition code 2
              IBCODE3 (Optional) Boundary condition code 3
              IBCODE4 (Optional) Boundary condition code 4
              IBCODE5 (Optional) Boundary condition code 5
            */
            std::string frname;
            file >> frname;
            Debug( 8012 ) << "Reading " << frname << " BC(with id=" << thegroupid << ") ... " << "\n";

            int itype;

            int nentry;
            int nvalues;
            int ibcode;

            file >> itype >> nentry >> nvalues >> ibcode;

            for ( int kfr = 0; kfr < nentry; ++kfr )
            {
                int nfr;
                file >> nfr;
                nfr--;


                if ( itype == 0 )
                {
                    /*
                      NODE                    Node number
                      (VALUES(I),I=1,NVALUES) Nodal values
                    */
                    boundary[kfr].get<0>() = true;
                    boundary[kfr].get<1>() = thegroupid;
                }

                else if ( itype == 1 )
                {
                    /*
                      ELEM                    Element/cell number
                      ELEMENT TYPE            Element type
                      FACE                    Face number
                      (VALUES(I),I=1,NVALUES) Element/cell values
                    */
                    int type, face;
                    file >> type >> face;
                    --face;


                    switch ( type )
                    {
                    case Element::TRIANGLE:
                        elements[nfr].get<3>() = boost::make_tuple( face, thegroupid );

                        for ( int n = 0; n < 2; ++n )
                        {
                            Debug( 8013 )  << "node " << tria3::face[face][n] << " in element " << nfr
                                           << " has  condition " << thegroupid << "\n";
                            boundary[boost::get<2>( elements[nfr] )[tria3::face[face][n]]].get<0>() = true;
                            boundary[boost::get<2>( elements[nfr] )[tria3::face[face][n]]].get<1>() = thegroupid;
                        }

                        break;

                    case Element::TETRAHEDRON:
                        for ( int n = 0; n < 3; ++n )
                        {
                            int g2l[]= { 3, 2, 1, 0 };
                            elements[nfr].get<3>() = boost::make_tuple( g2l[face], thegroupid );

                            int local_index = Feel::details::tetra<1>::f2p( g2l[face], n );
                            Debug( 8013 )  << "node " << local_index << " in element " << nfr
                                           << " has  condition " << thegroupid << "\n";
                            boundary[boost::get<2>( elements[nfr] )[local_index]].get<0>() = true;
                            boundary[boost::get<2>( elements[nfr] )[local_index]].get<1>() = thegroupid;
                        }

                        break;

                    }

                    FEELPP_ASSERT( elements[nfr].get<0>() == type )( type )( elements[nfr].get<0>() ).error( "[gambit] inconsistent element type for boundary conditions specifications" );
                }

            }

            Debug( 8012 )  << "Reading " << frname << " BC(with id=" << thegroupid << ")   done." << "\n";
            ++thegroupid;

        }

    }


    Debug( 8012 ) << "completed." << "\n";

    return true;
}
} // gambit
} // Feel

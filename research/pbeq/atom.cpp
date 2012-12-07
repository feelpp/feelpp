/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

This file is part of the Feel library

Author(s): Simone Deparis <simone.deparis@epfl.ch>
Date: 2007-07-10

Copyright (C) 2007 EPFL

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
   \file atom.cpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2007-07-10
*/

#include "atom.hpp"


namespace Feel
{
Atom::Atom( )
    :
    M_center( dim )
{
    std::fill( M_center.begin(), M_center.end(), value_type( 0 ) );
}

Atom::Atom( Atom const& atom )
    :
    M_center  ( atom.M_center ),
    M_id      ( atom.M_id ),
    M_name    ( atom.M_name ),
    M_peptide ( atom.M_peptide ),
    M_pepId   ( atom.M_pepId ),
    M_charge  ( atom.M_charge ),
    M_radius  ( atom.M_radius ),
    M_radius2  ( atom.M_radius2 )
{
}

Atom::~Atom( )
{
}


Atom&
Atom::operator=( Atom const& atom )
{
    if ( this != &atom )
    {
        M_center  = atom.M_center;
        M_id      = atom.M_id;
        M_name    = atom.M_name;
        M_peptide = atom.M_peptide;
        M_pepId   = atom.M_pepId;
        M_charge  = atom.M_charge;
        M_radius  = atom.M_radius;
        M_radius2  = atom.M_radius2;
    }

    return *this;
}



int
Atom::readPQRline( std::istream& PQRline, bool const neglectChrRad )
{
    PQRline >> M_id
            >> M_name
            >> M_peptide
            >> M_pepId;

    for ( int i( 0 ); i< dim; i++ )
        PQRline >> M_center[i];

    if ( neglectChrRad )
    {
        std::string tmp;
        std::getline( PQRline,tmp );

    }

    else
    {
        PQRline >> M_charge
                >> M_radius;

        M_radius2 = M_radius * M_radius;

        M_ID = -1;
        M_pName = M_peptide;
    }

    if ( PQRline.eof() )
        return 1;

    if ( PQRline.fail() )
    {
        showMe();
        return -1;
    }

    return 0;

}


/**
 * read one cha.crd  and one rad.crd lines which looks like
ID id  peptide name   x       y        z       pName  pID charge
1    1 GLY  N    -24.90000   6.97500   4.91400 SEG1  -2   -0.30000

ID id  peptide name   x       y        z       pName  pID radius
1    1 GLY  N    -24.90000   6.97500   4.91400 SEG1  -2   1.8500

 * retuns 0 if line read susscessfully, otherwise 1
 */

int
Atom::readCRDlines( std::istream& chaline, std::istream& radline )
{
    M_charge = readCRDlineCenter( chaline );

    M_radius = readCRDlineCenter( radline );

    M_radius2 = M_radius * M_radius;

    if ( chaline.eof() || radline.eof() )
        return 1;

    if ( chaline.fail() || radline.fail() )
    {
        showMe();
        return -1;
    }

    return 0;

}


value_type
Atom::readCRDlineCenter( std::istream& line )
{
    value_type last;
    line >> M_ID
         >> M_id
         >> M_peptide
         >> M_name;

    for ( int i( 0 ); i< dim; i++ )
        line >> M_center[i];

    line >> M_pName
         >> M_pepId
         >> last;

    return last;
}




void Atom::showMe() const
{
    std::ostringstream ostr;
    ostr << M_id << " "
         << M_name << " "
         << M_peptide << " "
         << M_pepId << " ";

    for ( int i( 0 ); i< dim; i++ )
        ostr << M_center[i] << " ";

    ostr << M_charge << " "
         << M_radius;

    LOG(INFO) << ostr.str() << "\n";

}

}

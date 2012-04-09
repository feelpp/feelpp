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
   \file atom.hpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2007-07-10
*/

#ifndef _ATOM_HPP
#define _ATOM_HPP

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/glas.hpp>

#include <vector>
#include <istream>
#include <string>

#include "pbeqapplication.hpp"

namespace Feel
{

/**
 * \class Atom
 * \brief Atom's definition and functions
 *
 * This class defines an atom through its position, charge and radius.
 * This is supposed to work with the Molecule class.
 * \ingroup ??
 * @author Simone Deparis
 * @see Im, Beglov, and Roux, "Continuum Solvation Model: ..."
 */

class Atom
{
public:
    /** @name Static values
     */
    //@{

    static const int dim = 3;

    //@}

    /** @name Typedefs
     */
    //@{

    typedef node<value_type>::type node_type;

    //@}


    /** @name Constructors, destructor
     */
    //@{


    /**
     * default atom constructor: atom at origin with zero radius and charges
     */
    Atom();

    Atom( Atom const& atom );

    ~Atom();

    /** @name Operator overloads
     */
    //@{

    Atom& operator=( Atom const& atom );

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * accessors
     */


    node_type               const& center()  const
    {
        return M_center;
    }
    uint16_type             const& id()      const
    {
        return M_id;
    }
    std::string             const& name()    const
    {
        return M_name;
    }
    std::string             const& peptide() const
    {
        return M_peptide;
    }
    std::string             const& pepId()   const
    {
        return M_pepId;
    }
    value_type              const& charge()  const
    {
        return M_charge;
    }
    value_type              const& radius()  const
    {
        return M_radius;
    }
    value_type              const& radius2()  const
    {
        return M_radius2;
    }

    /**
     * modifiers
     */

    void setCenter ( node_type               const& center )
    {
        M_center = center;
    }
    void setId     ( uint16_type             const& id     )
    {
        M_id     = id;
    }
    void setName   ( std::string             const& name   )
    {
        M_name   = name;
    }
    void setPeptide( std::string             const& peptide )
    {
        M_peptide=peptide;
    }
    void setPepId  ( std::string             const& pepId  )
    {
        M_pepId  = pepId;
    }
    void setCharge ( value_type              const& charge )
    {
        M_charge = charge;
    }
    void setRadius ( value_type              const& radius )
    {
        M_radius = radius;
        M_radius2 = M_radius*M_radius;
    }


    //@}


    /**
     * read one PQR line which looks like
    id name peptide pID      x       y       z     charge radius
    1  N   PRO     1     -15.656  10.553  -3.283 -0.0700 1.8500
     * retuns 0 if line read susscessfully, 1 if at end of file, otherwise -1
     * if neglectChrRad = true, charge and radius are read but neglected
     */
    int readPQRline( std::istream& PQRline, bool const neglectChrRad = false );

    /**
     * read one cha.crd  and one rad.crd lines which looks like
    ID id  peptide name   x       y        z       pName  pID charge
    1    1 GLY  N    -24.90000   6.97500   4.91400 SEG1  -2   -0.30000

    ID id  peptide name   x       y        z       pName  pID radius
    1    1 GLY  N    -24.90000   6.97500   4.91400 SEG1  -2   1.8500

     * retuns 0 if line read susscessfully, 1 if at end of file, otherwise -1
     */
    int readCRDlines( std::istream& chaline, std::istream& radline );

    void showMe() const;


private:

    value_type readCRDlineCenter( std::istream& line );

    node_type               M_center;
    int16_type              M_ID; //
    uint16_type             M_id;
    std::string             M_name;
    std::string             M_peptide;
    std::string             M_pName; //
    std::string             M_pepId;
    value_type              M_charge;
    value_type              M_radius;
    value_type              M_radius2;

}; // end class Atom


} // end namespace Feel


#endif

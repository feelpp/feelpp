/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2010-01-21

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
/**
   \file stokesbase.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-01-21
 */
#ifndef __StokesBase_H
#define __StokesBase_H 1

#include <feel/feelcore/factory.hpp>
#include <feel/feelcore/singleton.hpp>

namespace Feel
{
/**
 * \class StokesBase
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
    class StokesBase
    {
    public:


        /** @name Constants
         */
        //@{


        //@}

        /** @name Typedefs
         */
        //@{

        //! Stokes factory
        typedef Feel::Singleton< Feel::Factory< StokesBase, std::string > > factory_type;

        //@}

        /** @name Constructors, destructor
         */
        //@{

        //! default constructor
        StokesBase();
        //! copy constructor
        StokesBase( StokesBase const & );
        //! destructor
        virtual ~StokesBase();

        /**
         * initialize the polynomial visualisation class
         * \param vm variables map from the command line
         */
        virtual void init( po::variables_map const& vm );

        /**
         * create a new polynomial visualisation class
         * \param vm variables map from the command line
         */
        static StokesBase* New( po::variables_map const& vm );


        //@}

        /** @name Operator overloads
         */
        //@{

        //! copy operator
        StokesBase& operator=( StokesBase const & o);

        //@}

        /** @name Accessors
         */
        //@{

        //! return the variables map
        po::variables_map vm() const { return M_vm; }

        //! return the polynomial name
        std::string name() const { return M_name; }

        //@}

        /** @name  Mutators
         */
        //@{


        //@}

        /** @name  Methods
         */
        //@{


        //@}



    protected:

    private:

        //! variables map
        po::variables_map M_vm;

        //! stokes family name
        std::string M_name;

    };
}
#endif /* __StokesBase_H */

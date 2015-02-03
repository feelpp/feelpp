/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2004-11-09

  Copyright (C) 2004,2005 EPFL
  Copyright (C) 2007-2012 Universite Joseph Fourier (Grenoble I)

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
/*!
 * \file exporterhdf5.hpp
 * \brief HDF5 and XDMF exporter
 * \author Benjamin Vanthong <benjamin.vanthong@gmail.com>
 * \author Alexandre Ancel <alexandre.ancel@cemosis.fr>
 * \date 2014-08-28
 */
#ifndef __Exporterhdf5_H
#define __Exporterhdf5_H 1

#if defined(FEELPP_HAS_HDF5)

#include <feel/feelfilters/exporter.hpp>
#include <feel/feelcore/hdf5.hpp>


namespace Feel 
{
namespace fs = boost::filesystem;

template <typename MeshType, int N>
class Exporterhdf5
    : 
        public Exporter <MeshType, N>
{
    typedef Exporter<MeshType, N> super;
public: 
    typedef MeshType mesh_type;
    typedef typename mesh_type::value_type value_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;    
    typedef typename super::timeset_type timeset_type;
    typedef typename super::timeset_ptrtype timeset_ptrtype;
    typedef typename super::timeset_iterator timeset_iterator;
    typedef typename super::timeset_const_iterator timeset_const_iterator;

    Exporterhdf5( WorldComm const& worldComm = Environment::worldComm() );
    Exporterhdf5( std::string const& __p = "default", int freq = 1, WorldComm const& worldComm = Environment::worldComm() );
    Exporterhdf5( po::variables_map const& vm=Environment::vm(), std::string const& exp_prefix = "", WorldComm const& worldComm = Environment::worldComm() ) FEELPP_DEPRECATED;
    Exporterhdf5( std::string const& exp_prefix, WorldComm const& worldComm = Environment::worldComm() );

    Exporterhdf5( Exporterhdf5 const & __ex );

    ~Exporterhdf5(); 

    /** @name  Mutators
     */
    //@{

    Exporter<MeshType,N>* setOptions( std::string const& exp_prefix = "" )
        {
            super::setOptions( exp_prefix );

            return this;
        }
    Exporter<MeshType,N>* setOptions( po::variables_map const& vm, std::string const& exp_prefix = "" ) FEELPP_DEPRECATED
        {
            super::setOptions( exp_prefix );

            return this;
        }

        //@}


        void save() const;
        void visit( mesh_type* mesh);

    private :
        /*!
         * \brief Fonction used in almost all constructor to initialize the element's type
         */
        void init();

        /*!
         * \brief write .xmf and .h5 files for each process
         */
        void write() const;

        /*!
         * \brief write .h5 files for each process
         */
        void writeHDF5() const;

        /*!
         * \brief write .xmf file
         */
        void writeXDMF() const;

        /*!
         * \brief write mesh data
         */
        void saveMesh(mesh_ptrtype mesh, int stepIndex) const;

        /*!
         * \brief write informations of the mesh in .h5 file (unused for now)
         */
        void writeStats() const;

        /*!
         * \brief save solutions on nodes 
         * \param __step a time step
         * \param __var  iterator on solutions (begin)
         * \param en     iterator on solutions (end)
         */
        template<typename Iterator>
            void saveNodal( typename timeset_type::step_ptrtype __step, Iterator __var, Iterator en ) const;

        /*!
         * \brief save solutions on elements 
         * \param __step   a time step
         * \param __evar   iterator on solutions (begin)
         * \param __evaren iterator on solutions (end)
         */
        template<typename Iterator>
            void saveElement( typename timeset_type::step_ptrtype __step, Iterator __evar, Iterator __evaren ) const;

    private :
        mutable int tabCount;                   /*!< Number of tabs to print for Xdmf */
        mutable std::ostringstream M_fileName;        /*!< file name */
        mutable HDF5 M_HDF5;                   /*!< HDF5 IO */

        // Mesh geometry
        mutable std::string M_element_type;    /*!< element's type */

        mutable std::ofstream M_xmf;          /*!< Out stream to write the .xmf file */
        mutable std::ostringstream M_XDMFContent;               /*!< Content of Xdmf file */
};
} // Feel

#include <feel/feelfilters/exporterhdf5_impl.hpp>

#endif /* FEELL_HAS_HDF5 */
#endif /* __Exporterhdf5_H */


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
 * \author VANTHONG Benjamin <benjamin.vanthong@gmail.com>
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
    typedef Exporter<MeshType, N> super ;
    public: 
        typedef MeshType mesh_type ;
        typedef typename mesh_type::value_type value_type ;
        typedef boost::shared_ptr<mesh_type> mesh_ptrtype ;    
        typedef typename super::timeset_type timeset_type;
        typedef typename super::timeset_ptrtype timeset_ptrtype;
        typedef typename super::timeset_iterator timeset_iterator;
        typedef typename super::timeset_const_iterator timeset_const_iterator;

        Exporterhdf5 ( WorldComm const& worldComm = Environment::worldComm() ) ;
        Exporterhdf5 ( std::string const& __p = "default", int freq = 1, WorldComm const& worldComm = Environment::worldComm() ) ;
        Exporterhdf5( po::variables_map const& vm=Environment::vm(), std::string const& exp_prefix = "", WorldComm const& worldComm = Environment::worldComm() );

        Exporterhdf5 ( Exporterhdf5 const & __ex ) ;

        ~Exporterhdf5 () ; 

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


        void save () const ;
        void visit ( mesh_type* mesh) ;

    private :
        /*!
         * \brief Fonction used in almost all constructor to initialize the element's type
         */
        void init() ;

        /*!
         * \brief write .xmf and .h5 files for each process
         */
        void write ()  const ;

        /*!
         * \brief only one write .xmf file and one .h5 file for each time step  
         */
        void writeMerge ()  const ;

        /*!
         * \brief write points' coordonates   
         */
        void writePoints () const ;

        /*!
         * \brief write points' coordonates (merge version)
         */
        void writePointsMerge () const ;

        /*!
         * \brief write elements (an element is formed by several nodes) 
         */
        void writeElements () const ;

        /*!
         * \brief write elements (merge version) 
         */
        void writeElementsMerge () const ;

        /*!
         * \brief write informations of the mesh in .h5 file (unused for now)
         */
        void writeStats () const ;

        /*!
         * \brief save solutions on nodes 
         * \param __step a time step
         * \param __var  iterator on solutions (begin)
         * \param en     iterator on solutions (end)
         */
        template<typename Iterator>
            void saveNodal ( typename timeset_type::step_ptrtype __step, Iterator __var, Iterator en ) const ;

            /*!
             * \brief save solutions on nodes (merge version)
             * \param __step a time step
             * \param __var  iterator on solutions (begin)
             * \param en     iterator on solutions (end)
             */
        template<typename Iterator>
            void saveNodalMerge ( typename timeset_type::step_ptrtype __step, Iterator __var, Iterator en ) const ;

            /*!
             * \brief save solutions on elements  )
             * \param __step   a time step
             * \param __evar   iterator on solutions (begin)
             * \param __evaren iterator on solutions (end)
             */
        template<typename Iterator>
            void saveElement ( typename timeset_type::step_ptrtype __step, Iterator __evar, Iterator __evaren ) const ;

            /*!
             * \brief save solutions on elements (merge version) 
             * \param __step   a time step
             * \param __evar   iterator on solutions (begin)
             * \param __evaren iterator on solutions (end)
             */
        template<typename Iterator>
            void saveElementMerge ( typename timeset_type::step_ptrtype __step, Iterator __evar, Iterator __evaren ) const ;

        /*!
         * \brief open XDMF file (used only in version without merge)
         */
        void open_xdmf_xml () const ;

        /*!
         * \brief close properly XDMF file (used only in version without merge)
         */
        void close_xdmf_xml () const ;


        /*!
         * \brief a bubble sort of 2 arrays in the same time
         * \param ids    array of points' identifier
         * \param coords coordinates of points
         * \param n      size of those array
         */
        void bubbleSort (size_type * ids, value_type * coords, size_type n) const ;

    private :
        mutable std::string M_fileName ;        /*!< file name */
        mutable std::string M_fileNameStep ;    /*!< file name + time step */
        mutable HDF5 M_HDF5 ;                   /*!< HDF5 IO */
        mutable mesh_ptrtype M_meshOut ;        /*!< pointer on current mesh in a time step */

        // Mesh geometry
        mutable size_type M_elementNodes ;      /*!< number of nodes for one element */ 
        mutable size_type M_maxNumElements ;    /*!< number of elements for the current process */
        mutable size_type M_maxNumPoints ;      /*!< number of points for the current process */
        mutable size_type M_numParts ;          /*!< number of partitions */
        mutable std::string M_element_type ;    /*!< element's type */

        mutable size_type M_step = 0 ;          /*!< number of the current step */
        mutable std::ofstream M_xmf  ;          /*!< Out stream to write the .xmf file */

        mutable std::vector<size_type> M_uintBuffer ;           /*!< buffer of integer */
        mutable std::vector<value_type> M_realBuffer ;          /*!< buffer of double */
        mutable std::map<size_type, size_type> M_newPointId ;   /*!< new point identifier after sort */
        mutable std::ostringstream M_str ;                      /*!< buffer of string */
};
} // Feel

#include <feel/feelfilters/exporterhdf5_impl.hpp>

#endif /* FEELL_HAS_HDF5 */
#endif /* __Exporterhdf5_H */


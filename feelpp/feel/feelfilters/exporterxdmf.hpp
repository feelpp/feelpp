/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
#ifndef FEELPP_FILTERS_EXPORTERXDMF_HPP
#define FEELPP_FILTERS_EXPORTERXDMF_HPP 1

#if defined(FEELPP_HAS_HDF5)

#include <feel/feelfilters/exporter.hpp>
#include <feel/feelcore/hdf5.hpp>
#include <feel/feelfilters/detail/meshcontiguousnumberingmapping.hpp>


namespace Feel
{

template <typename MeshType, int N>
class ExporterXDMF
    :
        public Exporter <MeshType, N>
{
    typedef Exporter<MeshType, N> super;
public:
    typedef MeshType mesh_type;
    typedef typename mesh_type::value_type value_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    typedef typename super::timeset_type timeset_type;
    typedef typename super::timeset_ptrtype timeset_ptrtype;
    typedef typename super::timeset_iterator timeset_iterator;
    typedef typename super::timeset_const_iterator timeset_const_iterator;
protected :
    using steps_write_on_disk_type = typename super::steps_write_on_disk_type;
public :

    explicit ExporterXDMF( worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    ExporterXDMF( std::string const& __p = "default", int freq = 1, worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    ExporterXDMF( po::variables_map const& vm=Environment::vm(), std::string const& exp_prefix = "", worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() ) FEELPP_DEPRECATED;
    explicit ExporterXDMF( std::string const& exp_prefix, worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );

    ExporterXDMF( ExporterXDMF const & __ex );

    ~ExporterXDMF() override;

    /** @name  Mutators
     */
    //@{

    //@}


        void visit( mesh_type* mesh) override;

protected :

        //!  save the timeset
        void save( steps_write_on_disk_type const& stepsToWriteOnDisk ) const override;

    private :
        /*!
         * \brief Function used in almost all constructor to initialize the element's type
         */
        void init();

        /*!
         * \brief write .xmf and .h5 files for each process
         */
        void write() const;

        /*!
         * \brief write .h5 files for each process
         */
        void writeHDF5( steps_write_on_disk_type const& stepsToWriteOnDisk ) const;

        /*!
         * \brief write .xmf file
         */
        void writeXDMF( steps_write_on_disk_type const& stepsToWriteOnDisk ) const;

        /*!
         * \brief write mesh data
         */
        void saveMesh(mesh_ptrtype mesh, std::string const& geofilename ) const;

        /*!
         * \brief write information of the mesh in .h5 file (unused for now)
         */
    //void writeStats() const;
    void saveFields( typename timeset_type::step_ptrtype __step, std::string const& fieldsfilename ) const;
    template <bool IsNodal,typename MeshContiguousType, typename FieldDataType>
    void saveFields( std::string const& fieldsfilename, std::string const& groupName, MeshContiguousType const& mp, int part, FieldDataType const& fieldData,  std::vector<float> & realBuffer ) const;

        /*!
         * \brief save solutions on nodes or elements
         * \param __step a time step
         * \param __var  iterator on solutions (begin)
         * \param en     iterator on solutions (end)
         */
    //template<bool IsNodal,typename Iterator>
    //    void saveFields( typename timeset_type::step_ptrtype __step, Iterator __var, Iterator en ) const;

    private :
    //mutable int tabCount;                   /*!< Number of tabs to print for Xdmf */
    //mutable std::ostringstream M_fileName;        /*!< file name */
    mutable HDF5 M_hdf5Fields/*, M_hdf5Geo*/;                   /*!< HDF5 IO */

    // Mesh geometry
    mutable std::string M_element_type;    /*!< element's type */

    //    mutable std::ofstream M_xmf;          /*!< Out stream to write the .xmf file */
    //    mutable std::ostringstream M_XDMFContent;               /*!< Content of Xdmf file */
    mutable std::map<int,std::ostringstream> M_XDMFContent;
    mutable std::map<int,std::tuple<std::string,std::ostringstream>> M_XDMFGeoContent;

    mutable std::unordered_map<int, Feel::detail::MeshContiguousNumberingMapping<mesh_type,float>> M_cache_mp;
    mutable std::map<int,std::vector<size_type>> M_mapNodalArrayToDofId;
    mutable std::map<int,std::vector<size_type>> M_mapElementArrayToDofId;
};
} // Feel

#include <feel/feelfilters/exporterxdmf_impl.hpp>

#endif /* FEELL_HAS_HDF5 */
#endif

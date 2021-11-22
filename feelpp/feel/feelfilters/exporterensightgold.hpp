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
/**
   \file ExporterEnsightGold.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Alexandre Ancel <alexandre.ancel@cemosis.fr>
   \date 2006-11-26
 */
#ifndef FEELPP_FILTERS_EXPORTERENSIGHTGOLD_HPP
#define FEELPP_FILTERS_EXPORTERENSIGHTGOLD_HPP 1

#include <iostream>
#include <fstream>


#include <boost/lambda/lambda.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include <feel/feelmesh/filters.hpp>
#include <feel/feelfilters/detail/fileindex.hpp>
#include <feel/feelfilters/detail/meshcontiguousnumberingmapping.hpp>

namespace Feel
{
namespace fs = boost::filesystem;

/**
 * \class ExporterEnsightGold
 * \brief exporter to EnsightGold format
 *
 * \ingroup Exporter
 * @author Christophe Prud'homme
 * @author Alexandre Ancel
 */
template<typename MeshType, int N>
class ExporterEnsightGold
    :
public Exporter<MeshType, N>
{
    typedef Exporter<MeshType, N> super;
public:


    /** @name Typedefs
     */
    //@{

    typedef MeshType mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    using index_type = typename mesh_type::index_type;
    typedef typename super::timeset_type timeset_type;
    typedef typename super::timeset_ptrtype timeset_ptrtype;
    typedef typename super::timeset_iterator timeset_iterator;
    typedef typename super::timeset_const_iterator timeset_const_iterator;
protected :
    using mesh_contiguous_numbering_mapping_type = Feel::detail::MeshContiguousNumberingMapping<mesh_type,float>;
    using mesh_contiguous_numbering_mapping_ptrtype = std::shared_ptr<mesh_contiguous_numbering_mapping_type>;
    using step_ptrtype = typename super::step_ptrtype;
    using steps_write_on_disk_type = typename super::steps_write_on_disk_type;
public :
    //@}

    /** @name Constructors, destructor
     */
    //@{
    /**
     *
     The elements that are supported by the EnSight6 format are:

     \htmlonly
     <pre>
     1                 1------------------2        1----------2--------3
     point                   two node bar                three node bar


     7
     4-------------3          4-------------3
     3                 |             |          |             |
     3                        /\                |             |          |             |
     /\                      /  \               |             |        8 |             | 6
     /  \               6    /    \  5           |             |          |             |
     /    \                  /      \             |             |          |             |
     /      \                /        \            |             |          |             |
     /        \              /          \           |             |          |      5      |
     /          \            /    4       \          1-------------2          1-------------2
     1------------2           1------------2
     three node triangle       six node triangle       four node quadrangle     eight node quadrangle


     /\
     / |\
     /  |4\
     /   |  \
     /    |   \
     /     |    \
     1------|-----\
     \     |    3/
     \    |    /
     \  2|   /
     \  |  /
     \ | /
     \\2/

     four node tetrahedron
     </pre>
     \endhtmlonly

    */
    explicit ExporterEnsightGold( worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    ExporterEnsightGold( std::string const& __p = "default", int freq = 1, worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    explicit ExporterEnsightGold( std::string const& exp_prefix, worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );

    ExporterEnsightGold( ExporterEnsightGold const & __ex ) = default;

    ~ExporterEnsightGold() override;


    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the ensight element type
     */
    std::string const& elementType() const
    {
        return M_element_type;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    //@}

    /** @name  Methods
     */
    //@{

    void visit( mesh_type* mesh ) override;

    //@}



protected:

     //! save the timeset
    void save( steps_write_on_disk_type const& stepsToWriteOnDisk ) const override;

private:

    //! init the ensight exporter
    FEELPP_NO_EXPORT void init();

    //! write the '' file for ensight
    FEELPP_NO_EXPORT void writeSoSFile() const;

    //! write case file variables
    template<typename Iterator, typename TSt>
    void writeCaseFileVariables( Iterator it, Iterator end,
                                 std::string const& loc,
                                 TSt const& __ts,
                                 std::ostream& __out ) const;

    //! write the 'case' file for ensight
    FEELPP_NO_EXPORT void writeCaseFile() const;

    //! write the 'geo' file for ensight
    FEELPP_NO_EXPORT void writeGeoFiles( timeset_ptrtype __ts, mesh_ptrtype mesh, int timeIndex, bool isFirstStep ) const;
    FEELPP_NO_EXPORT void writeGeoMarkers( MPI_File fh, mesh_contiguous_numbering_mapping_type const& mp, bool writeHeaderBeginFile, bool writeBeginEndTimeSet, Feel::detail::FileIndex & index ) const;
    FEELPP_NO_EXPORT void writeGeoMarkedFaces(MPI_File fh, mesh_ptrtype mesh, std::pair<const std::string, std::vector<index_type> > & m) const;
    FEELPP_NO_EXPORT void writeGeoMarkedElements(MPI_File fh, mesh_contiguous_numbering_mapping_type const& mp, int part ) const;

    //! write the variables file for ensight
    FEELPP_NO_EXPORT void writeVariableFiles( timeset_ptrtype __ts, step_ptrtype step ) const;
    template<bool IsNodal,typename Iterator>
    FEELPP_NO_EXPORT void saveFields( timeset_ptrtype __ts, typename timeset_type::step_ptrtype __step, bool writeNewFile, std::string const& filenameStepIndex, bool isFirstStep, Iterator __var, Iterator en ) const;

private:
    mutable std::string M_filename;
    std::string M_element_type;
    std::string M_face_type;
    bool M_mergeTimeSteps;
    int M_packTimeSteps;

    // mapping allow to get ordering between Feel++ and Ensight format with curve element
    std::map<std::string,std::vector<uint16_type>> M_nodesOrderingInElementToEnsight;

    /* Number of digits used in timesteps */
    /* Set to 4 by default: range [0000; 9999] for timesteps */
    mutable int M_timeExponent;
    // file position for explicit pointers
    mutable MPI_Offset posInFile;
    mutable std::map<std::string, mesh_contiguous_numbering_mapping_ptrtype > M_cache_mp;
    mutable std::map<int,std::vector<size_type>> M_mapNodalArrayToDofId;
    mutable std::map<int,std::vector<size_type>> M_mapElementArrayToDofId;
};


} // Feel

//#if !defined( FEELPP_INSTANTIATION_MODE )
#include <feel/feelfilters/exporterensightgold_impl.hpp>
//#endif // FEELPP_INSTANTIATION_MODE

#endif /* __ExporterEnsightGold_H */

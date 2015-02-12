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
/**
   \file ExporterEnsightGold.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Alexandre Ancel <alexandre.ancel@cemosis.fr>
   \date 2006-11-26
 */
#ifndef __ExporterEnsightGold_H
#define __ExporterEnsightGold_H 1

#include <iostream>
#include <fstream>


#include <boost/lambda/lambda.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include <feel/feelmesh/filters.hpp>

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
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef typename super::timeset_type timeset_type;
    typedef typename super::timeset_ptrtype timeset_ptrtype;
    typedef typename super::timeset_iterator timeset_iterator;
    typedef typename super::timeset_const_iterator timeset_const_iterator;

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
    ExporterEnsightGold( WorldComm const& worldComm = Environment::worldComm() );
    ExporterEnsightGold( std::string const& __p = "default", int freq = 1, WorldComm const& worldComm = Environment::worldComm() );
    ExporterEnsightGold( po::variables_map const& vm=Environment::vm(), std::string const& exp_prefix = "", WorldComm const& worldComm = Environment::worldComm() ) FEELPP_DEPRECATED;
    ExporterEnsightGold( std::string const& exp_prefix, WorldComm const& worldComm = Environment::worldComm() );

    ExporterEnsightGold( ExporterEnsightGold const & __ex );

    ~ExporterEnsightGold();


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

    /**
     * \return the worldcomm passed in parameters of the exporter constructor (as the worldComm() can return a sequentialized one)
     */
    WorldComm const& worldCommBase() const
    {
        return M_worldCommBase;
    }


    //@}

    /** @name  Mutators
     */
    //@{

    Exporter<MeshType,N>* setOptions( po::variables_map const& vm, std::string const& exp_prefix = "" ) FEELPP_DEPRECATED
    {
        super::setOptions( exp_prefix );

        return this;
    }

    Exporter<MeshType,N>* setOptions( std::string const& exp_prefix = "" )
    {
        super::setOptions( exp_prefix );

        return this;
    }


    //@}

    /** @name  Methods
     */
    //@{

    /**
       save the timeset
    */
    void save() const;

    void visit( mesh_type* mesh );

    //@}



protected:

private:

    /**
     * init the ensight exporter
     */
    void init();

    /**
       write the '' file for ensight
    */
    void writeSoSFile() const;

    /**
       write the 'case' file for ensight
    */
    void writeCaseFile() const;

    /**
       updates the markers to be written by he exporters
    */
    void computeMarkersToWrite(mesh_ptrtype mesh) const;

    /**
       write the 'geo' file for ensight
    */
    void writeGeoFiles() const;
    void writeGeoMarkers(MPI_File fh, mesh_ptrtype mesh) const;
    void writeGeoHeader(MPI_File fh) const;
    void writeGeoMarkedFaces(MPI_File fh, mesh_ptrtype mesh, std::pair<const std::string, std::vector<size_type> > & m) const;
    void writeGeoMarkedElements(MPI_File fh, mesh_ptrtype mesh, size_type markerid) const;

    /**
       write the variables file for ensight
    */
    void writeVariableFiles() const;

    template<typename Iterator>
    void saveNodal( typename timeset_type::step_ptrtype __step, bool isFirstStep, Iterator __var, Iterator en ) const;

    template<typename Iterator>
    void saveElement( typename timeset_type::step_ptrtype __step, bool isFirstStep, Iterator __evar, Iterator __evaren ) const;

private:

    /* The purpose of this variable is to keep track */
    /* of the initial woldcomm that is given through the constructor */
    /* This is useful, when we don't want to merge the markers, */
    /* as the base worldComm is replaced with a sequential one */
    /* (and it also yields less code changes that having to use an alternate variable */
    /* containing the sequential worlcomm and leaving M_worldComm as the one in param) */
    WorldComm M_worldCommBase;

    mutable std::string M_filename;
    std::string M_element_type;
    std::string M_face_type;
    mutable std::set<int> M_markersToWrite;
    /* Number of digits used in timesteps */
    /* Set to 4 by default: range [0000; 9999] for timesteps */
    mutable int M_timeExponent;
};


} // Feel

//#if !defined( FEELPP_INSTANTIATION_MODE )
# include <feel/feelfilters/exporterensightgold_impl.hpp>
//#endif // FEELPP_INSTANTIATION_MODE

#endif /* __ExporterEnsightGold_H */

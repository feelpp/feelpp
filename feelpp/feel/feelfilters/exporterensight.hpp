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
   \file ExporterEnsight.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-11-26
 */
#ifndef FEELPP_FILTERS_EXPORTERENSIGHT_HPP
#define FEELPP_FILTERS_EXPORTERENSIGHT_HPP 1

#include <iostream>
#include <fstream>


#include <boost/lambda/lambda.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include <feel/feelmesh/filters.hpp>
#include <feel/feelfilters/detail/meshcontiguousnumberingmapping.hpp>

namespace Feel
{
namespace fs = boost::filesystem;

/**
 * \class ExporterEnsight
 * \brief exporter to Ensight format
 *
 * \ingroup Exporter
 * @author Christophe Prud'homme
 */
template<typename MeshType, int N>
class ExporterEnsight
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

    typedef typename super::timeset_type timeset_type;
    typedef typename super::timeset_ptrtype timeset_ptrtype;
    typedef typename super::timeset_iterator timeset_iterator;
    typedef typename super::timeset_const_iterator timeset_const_iterator;
protected :
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
    explicit ExporterEnsight( worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    ExporterEnsight( std::string const& __p = "default", int freq = 1, worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    ExporterEnsight( po::variables_map const& vm=Environment::vm(), std::string const& exp_prefix = "", worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() ) FEELPP_DEPRECATED;
    ExporterEnsight( std::string const& exp_prefix = "", worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() ) ;

    ExporterEnsight( ExporterEnsight const & __ex );

    ~ExporterEnsight() override;


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

    /**
       save the timeset
    */
    void save( steps_write_on_disk_type const& stepsToWriteOnDisk ) const override;

private:

    /**
     * init the ensight exporter
     */
    void init();

    /**
       write the '' file for ensight
    */
    void _F_writeSoSFile() const;

    /**
       write the 'case' file for ensight
    */
    void _F_writeCaseFile() const;

    /**
       write the 'geo' file for ensight
    */
    void _F_writeGeoFiles() const;

    /**
       write the variables file for ensight
    */
    void _F_writeVariableFiles() const;

    template<bool IsNodal,typename Iterator>
    void saveFields( typename timeset_type::step_ptrtype __step, Iterator __var, Iterator en ) const;


private:

    mutable std::string M_filename;
    std::string M_element_type;
    mutable std::unordered_map<int, Feel::detail::MeshPoints<float>> M_cache_mp;
    mutable std::optional<std::vector<size_type>> M_mapNodalArrayToDofId;
    mutable std::map<int,std::vector<size_type>> M_mapElementArrayToDofId;
};


} // Feel

//#if !defined( FEELPP_INSTANTIATION_MODE )
# include <feel/feelfilters/exporterensight_impl.hpp>
//#endif // FEELPP_INSTANTIATION_MODE

#endif /* __ExporterEnsight_H */

/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-03-30

  Copyright (C) 2005-2006 EPFL
  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file ExporterGmsh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-03-30
 */
#ifndef __ExporterGmsh_H
#define __ExporterGmsh_H 1

#include <iostream>
#include <fstream>


#include <boost/lambda/lambda.hpp>

#include <feel/feelcore/debug.hpp>

#include <feel/feelfilters/exporter.hpp>

namespace Feel
{
extern const char* FEELPP_GMSH_FORMAT_VERSION;
/**
 * \class ExporterGmsh
 * \brief Exporter to GMSH format
 *
 * \ingroup Exporter
 * @author Christophe Prud'homme
 */
template<typename MeshType, int N>
class ExporterGmsh
    :
public Exporter<MeshType,N>
{
public:


    /** @name Typedefs
     */
    //@{

    typedef MeshType mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef typename mesh_type::point_const_iterator point_const_iterator;

    typedef Exporter<MeshType,N> super;
    typedef typename mesh_type::value_type value_type;
    typedef typename super::timeset_type timeset_type;
    typedef typename super::timeset_ptrtype timeset_ptrtype;
    typedef typename super::timeset_iterator timeset_iterator;
    typedef typename super::timeset_const_iterator timeset_const_iterator;

    typedef typename timeset_type::step_type step_type;
    typedef typename timeset_type::step_ptrtype step_ptrtype;
    typedef typename timeset_type::step_const_iterator step_const_iterator;

    typedef typename matrix_node<value_type>::type matrix_node_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    ExporterGmsh( WorldComm const& worldComm = Environment::worldComm() );
    ExporterGmsh( std::string const& __p = "default", int freq = 1, WorldComm const& worldComm = Environment::worldComm() );
    ExporterGmsh( po::variables_map const& vm, std::string const& exp_prefix = "", WorldComm const& worldComm = Environment::worldComm() );

    ExporterGmsh( ExporterGmsh const & __ex );

    ~ExporterGmsh();

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //@}

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

    /** @name  Methods
     */
    //@{

    /**
       save the timeset
     */
    void save() const;

    /**
     * export mesh
     */
    void visit( mesh_type* mesh );

    /**
     * save the \p mesh to the file \p filename
     */
    void saveMesh( std::string const& filename, mesh_ptrtype mesh, bool parametric = false ) const;

    void gmshSaveAscii() const;

    void gmshSaveFormat( std::ostream& out, std::string const& version = FEELPP_GMSH_FORMAT_VERSION ) const;

    void gmshSavePhysicalNames( std::ostream& out, mesh_ptrtype mesh ) const;

    void gmshSaveNodesStart( std::ostream& out, mesh_ptrtype mesh, size_type nGlobPt, bool parametric = false ) const;
    void gmshSaveNodes( std::ostream& out, mesh_ptrtype mesh, bool parametric = false ) const;
    void gmshSaveNodesEnd( std::ostream& out, mesh_ptrtype mesh, bool parametric = false ) const;

    void gmshSaveElementsStart( std::ostream& out,size_type nGlobElt ) const;
    void gmshSaveElements( std::ostream& out, mesh_ptrtype __mesh, size_type indexEltStart ) const;
    void gmshSaveElementsEnd( std::ostream& out ) const;

    void gmshSaveNodeData( std::ostream& out, step_ptrtype __step ) const;

    void computeMinMax(step_ptrtype __step, std::map<std::string, std::vector<double> > & minMaxValues) const;
    void gmshSaveElementNodeData( std::ostream& out, step_ptrtype __step, size_type indexEltStart) const;

    template<typename ConvexType=typename mesh_type::shape_type>
    void gmshSaveOneElementAsMesh( std::string const& filename,
                                   typename mesh_type::element_type::super const& elt,
                                   PointSet<ConvexType,typename MeshType::value_type> const& ptset
                                   //=PointSet<ConvexType/*typename mesh_type::shape_type*/,typename MeshType::value_type>()
                                   ) const;

    template<typename ConvexRefType, typename ConvexPtSetType>
    void gmshSaveOneElementAsMesh( std::string const& filename,
                                   Reference<ConvexRefType,ConvexRefType::nDim,ConvexRefType::nOrder,ConvexRefType::nRealDim >  const& elt,
                                   PointSet<ConvexPtSetType,typename MeshType::value_type> const& ptset ) const;

    //@}

private:

    size_type numberOfGlobalPtAndIndex( mesh_ptrtype mesh ) const;

    boost::tuple<size_type,size_type> numberOfGlobalEltAndIndex( mesh_ptrtype mesh ) const;

    std::string M_element_type;

};

} // Feel

//#if !defined( FEELPP_INSTANTIATION_MODE )
# include <feel/feelfilters/exportergmsh_impl.hpp>
//#endif // FEELPP_INSTANTIATION_MODE

#endif /* __ExporterGmsh_H */

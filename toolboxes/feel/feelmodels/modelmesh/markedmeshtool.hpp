/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2011-07-05

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file markedmeshtool.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2011-07-05
 */

#ifndef FEELPP_MODELS_MARKEDMESHTOOL_H
#define FEELPP_MODELS_MARKEDMESHTOOL_H 1


#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>

namespace Feel
{
namespace FeelModels
{


template<typename IteratorRange,typename SpaceP0Type,typename VectorType>
void
updateP0EltMarkerFromFaceRange( IteratorRange const& range, std::shared_ptr<SpaceP0Type> const& XhP0,
                                std::shared_ptr<VectorType> & markEltVec )
{
    for ( auto itr = range.begin(), enr = range.end() ; itr!=enr ; ++itr )
    {
        auto const& face = boost::unwrap_ref( *itr );
        if ( face.isConnectedTo0() )
        {
            auto const& elt = face.element0();
            const size_type thedof = XhP0->dof()->localToGlobal(elt,0,0).index();
            markEltVec->add(thedof, 1.);
        }
        if ( face.isConnectedTo1() )
        {
            auto const& elt = face.element1();
            const size_type thedof = XhP0->dof()->localToGlobal(elt,0,0).index();
            markEltVec->add(thedof, 1.);
        }
    }
}



template< class MeshType >
class MarkedMeshTool
{
public :

    typedef MeshType mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    typedef FunctionSpace<mesh_type, bases<Lagrange<0,Scalar,Discontinuous> > > space_P0_type;
    typedef std::shared_ptr<space_P0_type> space_P0_ptrtype;
    typedef typename space_P0_type::element_type element_P0_type;
    typedef std::shared_ptr<element_P0_type> element_P0_ptrtype;


    MarkedMeshTool( mesh_ptrtype mesh )
        :
        M_mesh( mesh )
    {}

    void addFaceMarker( std::string mark ) { M_listMarkers.push_back( mark ); }
    void setFaceMarker( std::list<std::string> listmark ) { M_listMarkers = listmark; }
    void setFaceMarker( std::vector<std::string> vecmark )
    {
        M_listMarkers.clear();
        std::set<std::string> setmark( vecmark.begin(), vecmark.end() );
        for ( std::string const& mark : setmark )
            M_listMarkers.push_back( mark );
    }

    void
    buildSpaceP0()
    {
        M_XhP0 = space_P0_type::New(_mesh=M_mesh,
                                    _worldscomm=makeWorldsComm(1,M_mesh->worldCommPtr()),
                                    _extended_doftable=std::vector<bool>(1,true) );
    }

    void updateFaceMarker3FromInternalFaces()
    {
        M_mesh->updateMarker3WithRange(internalfaces(M_mesh),1);
    }

    void updateFaceMarker3FromFaceMarker()
    {
        for ( std::string mark : M_listMarkers )
            M_mesh->updateMarker3WithRange(markedfaces(M_mesh,mark),1);
    }

    void updateFaceMarker3FromEltConnectedToFaceMarker()
    {
        this->updateP0EltMarkerFromFaceMarker();
        M_mesh->updateMarker3( *M_markP0Elt );
        M_mesh->updateMarkersFromElements();
    }

    template < typename ExprT >
    void updateFaceMarker3FromExpr( vf::Expr<ExprT> const& expr, bool doUpdateForUse=true )
    {
        this->updateP0EltMarkerFromExpr(expr);
        if ( doUpdateForUse )
            this->updateForUseFaceMarker3();
    }
    void updateForUseFaceMarker3()
    {
        M_mesh->updateMarker3( *M_markP0Elt );
        M_mesh->updateMarkersFromElements();
    }

    void verbose()
    {
        size_type d = std::distance( marked3faces( M_mesh, 1 ).template get<1>(),
                                     marked3faces( M_mesh, 1 ).template get<2>() );
        size_type gd = d;

        if ( M_mesh->worldComm().localSize()>1 )
            mpi::all_reduce( M_mesh->worldComm().localComm(),
                             d,
                             gd,
                             std::plus<size_type>());

        if (M_mesh->worldComm().isMasterRank() )
            std::cout << "[ToolBoxMarkedMesh] : number of marked 3 faces: " << gd << std::endl;
    }

    void saveSubMeshFromMarked3Faces()
    {
        auto meshMark1 = createSubmesh( _mesh=M_mesh, _range=marked3faces(M_mesh,1),_context=EXTRACTION_KEEP_MESH_RELATION, _update=0, _only_on_boundary_faces=false);
        saveGMSHMesh(_mesh=meshMark1,_filename="submesh-marked3facesBy1.msh");
        auto meshMark0 = createSubmesh( _mesh=M_mesh, _range=marked3faces(M_mesh,0),_context=EXTRACTION_KEEP_MESH_RELATION, _update=0, _only_on_boundary_faces=false);
        saveGMSHMesh(_mesh=meshMark0,_filename="submesh-marked3facesBy0.msh");
    }

    void saveSubMeshFromMarked3Elements()
    {
        auto meshMark1 = createSubmesh( _mesh=M_mesh, _range=marked3elements(M_mesh,1),_context=EXTRACTION_KEEP_MESH_RELATION, _update=0, _only_on_boundary_faces=false);
        saveGMSHMesh(_mesh=meshMark1,_filename="submesh-marked3elementsBy1.msh");
        auto meshMark0 = createSubmesh( _mesh=M_mesh, _range=marked3elements(M_mesh,0),_context=EXTRACTION_KEEP_MESH_RELATION, _update=0, _only_on_boundary_faces=false);
        saveGMSHMesh(_mesh=meshMark0,_filename="submesh-marked3elementsBy0.msh");
    }


    void exportP0EltMarkerFromFaceMarker()
    {
        this->updateP0EltMarkerFromFaceMarker();
        auto myexporter = exporter( _mesh=M_mesh, _name="MyExportMarkP0Elt" );
        myexporter->step(0)->add( "markP0Elt", *M_markP0Elt );
        myexporter->save();
    }
private :

    void
    updateP0EltMarkerFromFaceMarker()
    {
        if (!M_XhP0)
            this->buildSpaceP0();
        if (!M_markP0Elt)
            M_markP0Elt = M_XhP0->elementPtr();

        auto markEltVec = backend()->newVector( M_XhP0 );

        for ( std::string mark : M_listMarkers )
            updateP0EltMarkerFromFaceRange( markedfaces(M_mesh,mark), M_XhP0, markEltVec );

        markEltVec->close();

        //*M_markP0Elt = *markEltVec;

        // because we add values, we need to push only the positives values
        for ( size_type k=0 ; k < M_XhP0->nLocalDof() ; ++k )
            if ( markEltVec->operator()(k) > 1e-8 ) M_markP0Elt->set( k,1. );
        //if ( M_markP0Elt->operator()(k) > 1e-8 ) M_markP0Elt->set( k,1. );
    }

    template < typename ExprT >
    void
    updateP0EltMarkerFromExpr( vf::Expr<ExprT> const& expr )
    {
        if (!M_XhP0)
            this->buildSpaceP0();
        if (!M_markP0Elt)
            M_markP0Elt = M_XhP0->elementPtr();

        auto markEltTemp = vf::project(_space=M_XhP0,
                                       _range=elements(M_mesh,EntityProcessType::ALL),
                                       _expr=expr );

        // use petsc vector for communicate values
        auto markEltVec = backend()->newVector( M_XhP0 );
        for ( size_type k=0 ; k < M_XhP0->nLocalDof() ; ++k )
            markEltVec->add(k,markEltTemp(k));
        //*markEltVec = markEltTemp;
        markEltVec->close();

        //*M_markP0Elt = *markEltVec;

        // because we add values, we need to push one the positives values
        for ( size_type k=0 ; k < M_XhP0->nLocalDof() ; ++k )
            if ( markEltVec->operator()(k) > 1e-8 ) M_markP0Elt->set( k,1. );
            //if ( M_markP0Elt->operator()(k) > 1e-8 ) M_markP0Elt->set( k,1. );
    }


private :
    mesh_ptrtype M_mesh;
    std::list<std::string> M_listMarkers;

    space_P0_ptrtype M_XhP0;
    element_P0_ptrtype M_markP0Elt;
};

} // namespace FeelModels
} // namespace Feel


#endif // FEELPP_MODELS_MARKEDMESHTOOL_H

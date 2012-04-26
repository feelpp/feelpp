/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   Date: 2007-07-21

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
   \file exportergmsh.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-07-21
*/
#ifndef __EXPORTERGMSH_CPP
#define __EXPORTERGMSH_CPP 1

#include <feel/feelcore/feel.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exportergmsh.hpp>
#include <feel/feelfilters/gmshenums.hpp>

namespace Feel
{
template<typename MeshType, int N>
ExporterGmsh<MeshType,N>::ExporterGmsh( std::string const& __p, int freq,
                                        WorldComm const& worldComm )
    :
    super( "gmsh", __p, freq, worldComm )
{

}
template<typename MeshType, int N>
ExporterGmsh<MeshType,N>::ExporterGmsh( po::variables_map const& vm, std::string const& exp_prefix,
                                        WorldComm const& worldComm )
    :
    super( vm, exp_prefix, worldComm )
{
}
template<typename MeshType, int N>
ExporterGmsh<MeshType,N>::ExporterGmsh( ExporterGmsh const & __ex )
    :
    super( __ex )
{}
template<typename MeshType, int N>
ExporterGmsh<MeshType,N>::~ExporterGmsh()
{}

template<typename MeshType, int N>
void
ExporterGmsh<MeshType,N>::save() const
{
    Debug( 8007 ) << "[ExporterGmsh] checking if frequency is ok\n";

    if ( this->cptOfSave() % this->freq()  )
    {
        this->saveTimeSet();
        return;
    }

    Debug( 8007 ) << "[ExporterGmsh] frequency is ok\n";

    Debug( 8007 ) << "[ExporterGmsh] save()...\n";

    gmshSaveAscii();

    Debug( 8007 ) << "[ExporterGmsh] saving done\n";

    this->saveTimeSet();
}

template<typename MeshType, int N>
void
ExporterGmsh<MeshType,N>::visit( mesh_type* )
{
}

template<typename MeshType, int N>
void
ExporterGmsh<MeshType,N>::gmshSaveAscii() const
{
    Debug( 8007 ) << "[gmshSaveascii] saving in gmsh ascii file format\n";

    timeset_const_iterator __ts_it = this->beginTimeSet();
    timeset_const_iterator __ts_en = this->endTimeSet();

    while ( __ts_it != __ts_en )
    {
        timeset_ptrtype __ts = *__ts_it;

        std::ostringstream __fname;

        __fname << __ts->name()  //<< this->prefix() //this->path()
                << "-" << M_comm.size() << "_" << M_comm.rank()
                << ".msh";

        /*std::string filename =  this->prefix()
          + __ts->name()
          + "-" + M_comm.size() + "_" + M_comm.rank()
          + ".msh";*/
        std::ofstream out;

        typename timeset_type::step_const_iterator __stepIt = __ts->beginStep();
        typename timeset_type::step_const_iterator __stepIt_end = __ts->endStep();
        __stepIt = boost::prior( __stepIt_end );

        if ( __stepIt != __stepIt_end )
        {
            step_ptrtype __step = *__stepIt;

            if ( __step->isInMemory() )
            {

                if ( __step->index()==1 )
                {
                    out.open( __fname.str().c_str(), std::ios::out );
                }

                else
                {
                    out.open( __fname.str().c_str(), std::ios::out | std::ios::app );
                }

                if ( out.fail() )
                {
                    Debug( 8007 ) << "cannot open " << __fname.str().c_str() << "\n";
                    exit( 0 );
                }

                Debug( 8007 ) << "[ExporterGmsh] saving model "
                              << __ts->name() << " at time step "
                              << __ts->index() << " in "
                              << __fname.str() << "\n";

                // save mesh only at first iteration
                if ( __stepIt == __ts->beginStep() )
                {
                    gmshSaveFormat( out );

                    gmshSavePhysicalNames( out, __step->mesh() );

                    auto nPointAndIndex = numberOfGlobalPtAndIndex( __step->mesh() );
                    auto nGlobPoint = nPointAndIndex.template get<0>();
                    auto indexPointStart = nPointAndIndex.template get<1>();
                    gmshSaveNodesStart( out, __step->mesh(), nGlobPoint );
                    gmshSaveNodes( out,__step->mesh(),indexPointStart );
                    gmshSaveNodesEnd( out, __step->mesh() );

                    auto nEltAndIndex = numberOfGlobalEltAndIndex( __step->mesh() );
                    auto nGlobElement = nEltAndIndex.template get<0>();
                    auto indexElementStart = nEltAndIndex.template get<1>();
                    gmshSaveElementsStart( out, nGlobElement );
                    gmshSaveElements( out, __step->mesh(), indexElementStart, indexPointStart );
                    gmshSaveElementsEnd( out );
                }

                //gmshSaveNodeData( out, __step);
                gmshSaveElementNodeData( out, __step );
            }
        }

        ++__ts_it;
    }
}


template<typename MeshType, int N>
void
ExporterGmsh<MeshType,N>::saveMesh( std::string const& filename, mesh_ptrtype mesh, bool parametric ) const
{

    if (  M_comm.rank() == 0 )
    {
        std::ofstream out( filename.c_str(), std::ios::out );

        if ( out.fail() )
        {
            Log() << "[ExporterGmsh::SaveMesh] cannot open " << filename << "\n";
            exit( 0 );
        }

        gmshSaveFormat( out );
        gmshSavePhysicalNames( out, mesh );
        out.close();
    }

    M_comm.barrier();

    //-----------------------------------------------------------------//

    auto nPointAndIndex = numberOfGlobalPtAndIndex( mesh );
    auto nGlobPoint = nPointAndIndex.template get<0>();
    auto indexPointStart = nPointAndIndex.template get<1>();

    if (  M_comm.rank() == 0 )
    {
        std::ofstream out( filename.c_str(), std::ios::app );
        gmshSaveNodesStart( out, mesh, nGlobPoint, parametric );
        out.close();
    }

    for ( int therank=0; therank<M_comm.size(); ++therank )
    {
        if ( therank == M_comm.rank() )
        {
            std::ofstream out( filename.c_str(), std::ios::app );
            gmshSaveNodes( out, mesh, indexPointStart, parametric );
            out.close();
        }

        M_comm.barrier();
    }

    if (  M_comm.rank() == 0 )
    {
        std::ofstream out( filename.c_str(), std::ios::app );
        gmshSaveNodesEnd( out, mesh, parametric );
        out.close();
    }

    //-----------------------------------------------------------------//

    auto nEltAndIndex = numberOfGlobalEltAndIndex( mesh );
    auto nGlobElement = nEltAndIndex.template get<0>();
    auto indexElementStart = nEltAndIndex.template get<1>();

    if (  M_comm.rank() == 0 )
    {
        std::ofstream out( filename.c_str(), std::ios::app );
        gmshSaveElementsStart( out, nGlobElement );
        out.close();
    }

    for ( int therank=0; therank<M_comm.size(); ++therank )
    {
        if ( therank == M_comm.rank() )
        {
            std::ofstream out( filename.c_str(), std::ios::app );
            gmshSaveElements( out, mesh, indexElementStart, indexPointStart );
            out.close();
        }

        M_comm.barrier();
    }

    if (  M_comm.rank() == 0 )
    {
        std::ofstream out( filename.c_str(), std::ios::app );
        gmshSaveElementsEnd( out );
        out.close();
    }
}

template<typename MeshType, int N>
void
ExporterGmsh<MeshType,N>::gmshSaveFormat( std::ostream& out, std::string const& version ) const
{
    out << "$MeshFormat\n"
        << version << " 0 " << sizeof( double ) << "\n"
        << "$EndMeshFormat\n";

}
template<typename MeshType, int N>
void
ExporterGmsh<MeshType,N>::gmshSavePhysicalNames( std::ostream& out, mesh_ptrtype mesh ) const
{
    // save Physical Names
    out << "$PhysicalNames\n";
    out << mesh->markerNames().size() << "\n";
    // physical-dimension physical-number "physical-name" ...
    BOOST_FOREACH( auto data, mesh->markerNames() )
    {
        out << data.second.template get<1>() << " "
            << data.second.template get<0>() << " "
            << "\"" << data.first << "\"" << "\n";
    }
    out << "$EndPhysicalNames\n";

}

template<typename MeshType, int N>
boost::tuple<size_type,size_type>
ExporterGmsh<MeshType,N>::numberOfGlobalPtAndIndex( mesh_ptrtype mesh ) const
{
    auto local_numberPoints = mesh->numPoints();
    auto global_numberPoints=local_numberPoints;

    mpi::all_reduce( M_comm,
                     local_numberPoints,
                     global_numberPoints,
                     std::plus<size_type>() );

    std::vector<size_type> all_localnumberPoint;

    mpi::all_gather( M_comm,
                     local_numberPoints,
                     all_localnumberPoint );

    size_type indexPtStart = 0;

    for ( int i=0; i<M_comm.rank(); ++i )
        indexPtStart+=all_localnumberPoint[i];

    return boost::make_tuple( global_numberPoints,indexPtStart );
}

template<typename MeshType, int N>
void
ExporterGmsh<MeshType,N>::gmshSaveNodesStart( std::ostream& out, mesh_ptrtype mesh, size_type nGlobPt, bool parametric ) const
{
    if ( parametric && mesh->isParametric() )
        out << "$ParametricNodes\n";

    else
        out << "$Nodes\n";

    out << nGlobPt << "\n";//number points
}

template<typename MeshType, int N>
void
ExporterGmsh<MeshType,N>::gmshSaveNodesEnd( std::ostream& out, mesh_ptrtype mesh, bool parametric ) const
{
    if ( parametric && mesh->isParametric() )
        out << "$EndParametricNodes\n";

    else
        out << "$EndNodes\n";
}

template<typename MeshType, int N>
void
ExporterGmsh<MeshType,N>::gmshSaveNodes( std::ostream& out, mesh_ptrtype mesh, size_type indexPtStart, bool parametric ) const
{

    point_const_iterator pt_it = mesh->beginPoint();
    point_const_iterator pt_en = mesh->endPoint();

    for ( ; pt_it!=pt_en ; ++pt_it )
    {
        out << pt_it->id()+1+indexPtStart
            << " "  << std::setw( 20 ) << std::setprecision( 16 ) << pt_it->node()[0];

        if ( mesh_type::nRealDim >= 2 )
            out << " "  << std::setw( 20 ) << std::setprecision( 16 ) << pt_it->node()[1];

        else
            out << " 0";

        if ( mesh_type::nRealDim >= 3 )
            out << " "  << std::setw( 20 ) << std::setprecision( 16 ) << pt_it->node()[2];

        else
            out << " 0";

        if ( parametric && mesh->isParametric() )
        {
            out << " " << pt_it->gDim() << " " << pt_it->gTag();

            if ( pt_it->gDim() == 1 )
                out << " " << std::setw( 20 ) << std::setprecision( 16 )    << pt_it->u();

            else if ( pt_it->gDim() == 2 )
                out << " " << std::setw( 20 ) << std::setprecision( 16 )    << pt_it->u()
                    << " " << std::setw( 20 ) << std::setprecision( 16 )    << pt_it->v();
        }

        out << "\n";
    }

}

template<typename MeshType, int N>
boost::tuple<size_type,size_type>
ExporterGmsh<MeshType,N>::numberOfGlobalEltAndIndex( mesh_ptrtype mesh ) const
{
    //auto allmarkedfaces = markedfaces( mesh );
    auto allmarkedfaces = boundaryfaces( mesh );
    size_type number_markedfaces= std::distance( allmarkedfaces.template get<1>(),allmarkedfaces.template get<2>() );

    auto eltOnProccess = elements( mesh );
    size_type number_elements= std::distance( eltOnProccess.template get<1>(),eltOnProccess.template get<2>() );

    auto local_numberElements = number_markedfaces+number_elements;
    auto global_numberElements=local_numberElements;

    mpi::all_reduce( M_comm,
                     local_numberElements,
                     global_numberElements,
                     std::plus<size_type>() );

    std::vector<size_type> all_localnumberElements;

    mpi::all_gather( M_comm,
                     local_numberElements,
                     all_localnumberElements );

    size_type indexEltStart = 0;

    for ( int i=0; i<M_comm.rank(); ++i )
        indexEltStart+=all_localnumberElements[i];

    return boost::make_tuple( global_numberElements,indexEltStart );
}

template<typename MeshType, int N>
void
ExporterGmsh<MeshType,N>::gmshSaveElementsStart( std::ostream& out, size_type nGlobElt ) const
{
    out << "$Elements\n"
        << nGlobElt << "\n";//number element
}

template<typename MeshType, int N>
void
ExporterGmsh<MeshType,N>::gmshSaveElementsEnd( std::ostream& out ) const
{
    out << "$EndElements\n";
}


template<typename MeshType, int N>
void
ExporterGmsh<MeshType,N>::gmshSaveElements( std::ostream& out, mesh_ptrtype mesh, size_type indexEltStart, size_type indexPtStart ) const
{

    //auto allmarkedfaces = markedfaces( mesh );
    auto allmarkedfaces = boundaryfaces( mesh );
    auto face_it = allmarkedfaces.template get<1>();
    auto face_end = allmarkedfaces.template get<2>();

    auto elem_number=indexEltStart+1;

    typedef typename MeshType::element_type element_type;
    GmshOrdering<element_type> ordering;

    typedef typename MeshType::face_type face_type;
    GmshOrdering<face_type> ordering_face;
    // save the faces

    for ( ; face_it != face_end; ++face_it )
    {
        // elm-number elm-type number-of-tags < tag > ... node-number-list
        out<< elem_number++ <<" ";
        out << ordering_face.type();
        // number-of-tags < tag >

        if ( FEELPP_GMSH_FORMAT_VERSION==std::string( "2.1" ) )
        {
            // out<<" 2 " << face_it->marker().value() << " " << face_it->marker2().value();
            out<<" 3 " << face_it->marker().value() << " " << face_it->marker2().value() << " " << face_it->processId()+1;
        }

        else if ( FEELPP_GMSH_FORMAT_VERSION==std::string( "2.2" ) )
        {
            uint16_type nbTag = 3 + face_it->numberOfPartitions();
            out << " " << nbTag
                << " " << face_it->marker().value()
                << " " << face_it->marker2().value()
                << " " << face_it->numberOfPartitions()
                << " " << face_it->processId()+1;

            for ( size_type i=0 ; i<face_it->numberOfNeighborPartitions(); ++i )
                out << " " << -( face_it->neighborPartitionIds()[i]+1 );
        }

        // node-number-list
        for ( uint16_type p=0; p<face_type::numPoints; ++p )
            out << " " << face_it->point( ordering_face.fromGmshId( p ) ).id()+1+indexPtStart;

        out<<"\n";
    } // faces


    auto eltOnProccess = elements( mesh );
    auto elt_it = eltOnProccess.template get<1>();
    auto elt_en = eltOnProccess.template get<2>();

    for ( ; elt_it != elt_en; ++elt_it )
    {
        out << elem_number++ <<" ";
        out << ordering.type();

        if ( FEELPP_GMSH_FORMAT_VERSION==std::string( "2.1" ) )
        {
            //out<<" 2 " << elt_it->marker().value() << " " << elt_it->marker2().value();
            out<<" 3 " << elt_it->marker().value() << " " << elt_it->marker2().value() << " " << elt_it->processId()+1;
        }

        else if ( FEELPP_GMSH_FORMAT_VERSION== std::string( "2.2" ) )
        {
            uint16_type nbTag = 3 + elt_it->numberOfPartitions();
            out << " " << nbTag
                << " " << elt_it->marker().value()
                << " " << elt_it->marker2().value()
                << " " << elt_it->numberOfPartitions()
                << " " << elt_it->processId()+1;

            for ( size_type i=0 ; i<elt_it->numberOfNeighborPartitions(); ++i )
                out << " " << -( elt_it->neighborPartitionIds()[i]+1 );
        }

        for ( uint16_type p=0; p<element_type::numPoints; ++p )
        {
            //std::cout << "index " << p << " -> " << ordering.fromGmshId(p) << " -> " << elt_it->point( ordering.fromGmshId(p) ).id()+1 << " : " << elt_it->point( ordering.fromGmshId(p) ).node() << "\n";
            out << " " << elt_it->point( ordering.fromGmshId( p ) ).id()+1+indexPtStart;
        }

        out<<"\n";
    } // elements

    //out << "$EndElements\n";

}



template<typename MeshType, int N>
void
ExporterGmsh<MeshType,N>::gmshSaveNodeData( std::ostream& out, step_ptrtype __step ) const
{
#if 0
    //!!!Not functionnal for curve element!!!

    typedef typename step_type::nodal_scalar_type nodal_scalar_type;
    typedef typename step_type::nodal_scalar_const_iterator nodal_scalar_const_iterator;

    //on parcourt le temp
    step_const_iterator __it = timeset->beginStep();
    step_const_iterator __end = timeset->endStep();

    for ( ; __it != __end ; ++__it )
    {
        nodal_scalar_const_iterator __var = ( *__it )->beginNodalScalar();
        nodal_scalar_const_iterator __varen = ( *__it )->endNodalScalar();

        out << "$NodeData\n";

        nodal_scalar_type const& __u = __var->second;
        //mesh_ptrtype mesh =__u.mesh();
        mesh_ptrtype mesh = ( *__it )->mesh();

        out << "1\n";//number of string tag
        out << "a scalar node\n";
        out << "1\n";//number of real tag
        out << "0.0\n";
        out << "3\n";//number of integer tags:
        out << "0\n";//the time step (0; time steps always start at 0)
        out << "1\n";//n-component (1 is scalar) field
        out << mesh->numPoints() << "\n";//number associated nodal values

        point_const_iterator pt_it = mesh->beginPoint();
        point_const_iterator pt_en = mesh->endPoint();

        for ( ; pt_it!=pt_en ; ++pt_it )
        {
            out << pt_it->id()+1
                <<" ";
            out << std::setw( 17 ) << std::setprecision( 16 ) <<__u( pt_it->id() );
            //__u(pt_it->node());
            out << "\n";
        }

        out << "$EndNodeData\n";
    }

#endif
}


template<typename MeshType, int N>
void
ExporterGmsh<MeshType,N>::gmshSaveElementNodeData( std::ostream& out,
        step_ptrtype __step ) const
{

    typedef typename mesh_type::element_const_iterator element_mesh_const_iterator;

    typedef typename step_type::nodal_scalar_type nodal_scalar_type;
    typedef typename step_type::nodal_scalar_const_iterator nodal_scalar_const_iterator;

    typedef typename step_type::nodal_vector_type nodal_vectorial_type;
    typedef typename step_type::nodal_vector_const_iterator nodal_vectorial_const_iterator;

    typedef typename step_type::element_scalar_type element_scalar_type;

    mesh_ptrtype mesh = __step->mesh();

    auto allmarkedfaces = boundaryfaces( mesh );
    int number_markedfaces= std::distance( allmarkedfaces.template get<1>(),allmarkedfaces.template get<2>() );

    uint16_type nLocalDof;
    size_type globaldof;
    typedef typename MeshType::element_type element_type;
    GmshOrdering<element_type> ordering;
    nodal_scalar_const_iterator __varScal = __step->beginNodalScalar();
    nodal_scalar_const_iterator __varScal_end = __step->endNodalScalar();

    for ( ; __varScal!=__varScal_end ; ++__varScal )
    {
        out << "\n$ElementNodeData\n";

        nodal_scalar_type const& __u = __varScal->second;

        //uint __nbCompFieldGMSH;
        //    if (nodal_scalar_type::functionspace_type::is_scalar)         { __nbCompFieldGMSH=1; }
        //else if (nodal_scalar_type::functionspace_type::is_vectorial) { __nbCompFieldGMSH=3; }
        //else if (nodal_scalar_type::functionspace_type::is_tensor2)   { __nbCompFieldGMSH=9; }

        out << "1\n";//number of string tag
        out << "\"" << __varScal->first <<"\"\n";//a scalar node\n";
        out << "1\n";//number of real tag
        out << __step->time() << "\n";//"0.0\n";//the time value (0.0)
        out << "3\n";//number of integer tags:
        out << __step->index()-1 << "\n";//"0\n";//the time step (0; time steps always start at 0)
        out << "1\n";//n-component (1 is scalar) field
        out << mesh->numElements() << "\n";//number associated nodal values

        element_mesh_const_iterator elt_it = mesh->beginElement();
        element_mesh_const_iterator elt_en = mesh->endElement();

        if ( !__u.areGlobalValuesUpdated() )
            __u.updateGlobalValues();

        for ( ; elt_it!=elt_en ; ++elt_it )
        {
            out << number_markedfaces + elt_it->id()+1;
            //nLocGeoPt = elt_it->nPoints();
            //nLocalDof = mesh->numLocalVertices();
            nLocalDof = nodal_scalar_type::functionspace_type::basis_type::nLocalDof;
            out << " " << nLocalDof;

            for ( uint16_type l = 0; l < nLocalDof; ++l )
            {
                uint16_type gmsh_l = ordering.fromGmshId( l );
                globaldof = boost::get<0>( __u.functionSpace()->dof()->localToGlobal( elt_it->id(), gmsh_l, 0 ) ); //l,c

                // verify that the dof points and mesh points coincide
#if !defined(NDEBUG)

                if ( ublas::norm_2( boost::get<0>( __u.functionSpace()->dof()->dofPoint( globaldof ) )-elt_it->point( ordering.fromGmshId( l ) ).node() ) > 1e-10 )
                {
                    std::cout << "------------------------------------------------------------\n";
                    std::cout << "invalid dof/mesh points\n";
                    std::cout << "dof global id:" << globaldof << " | local id:" << gmsh_l << "\n";
                    std::cout << "point global id:" <<  elt_it->point( ordering.fromGmshId( l ) ).id() << " | local id:" << gmsh_l << "\n";
                    std::cout << "node dof:  " << boost::get<0>( __u.functionSpace()->dof()->dofPoint( globaldof ) ) << "\n";
                    std::cout << "node element:  " << elt_it->point( ordering.fromGmshId( l ) ).node() << "\n";
                }

#endif // NDEBUG
                //out << " " << __u( globaldof);
                out << " " <<__u.container()( globaldof );
            }

            out << "\n";
        }

        out << "$EndElementNodeData\n";
    }

    nodal_vectorial_const_iterator __varVec = __step->beginNodalVector();
    nodal_vectorial_const_iterator __varVec_end = __step->endNodalVector();

    for ( ; __varVec!=__varVec_end ; ++__varVec )
    {
        out << "\n$ElementNodeData\n";

        nodal_vectorial_type const& __uVec = __varVec->second;

        uint16_type nComponents = __uVec.nComponents;

        out << "1\n";//number of string tag
        out << "\"" << __varVec->first <<"\"\n";//"a vectorial field\n";
        out << "1\n";//number of real tag
        out << __step->time() << "\n";//"0.0\n";//the time value (0.0)
        out << "3\n";//number of integer tags:
        out << __step->index() << "\n";//"0\n";//the time step (0; time steps always start at 0)
        out << "3\n";//n-component (3 is vectorial) field
        out << mesh->numElements() << "\n";//number associated nodal values

        element_mesh_const_iterator elt_it = mesh->beginElement();
        element_mesh_const_iterator elt_en = mesh->endElement();

        for ( ; elt_it!=elt_en ; ++elt_it )
        {
            out << number_markedfaces + elt_it->id()+1;
            //nLocalDof = mesh->numLocalVertices();
            nLocalDof = nodal_vectorial_type::functionspace_type::basis_type::nLocalDof;
            out << " " << nLocalDof;

            for ( uint16_type l = 0; l < nLocalDof; ++l )
            {
                uint16_type gmsh_l = ordering.fromGmshId( l );

                for ( uint16_type c = 0; c < 3; ++c )
                {
                    out << " ";

                    if ( c < nComponents )
                    {
                        globaldof = boost::get<0>( __uVec.functionSpace()->dof()->localToGlobal( elt_it->id(),
                                                   gmsh_l,
                                                   c ) );
                        //out << __uVec( globaldof);
                        out << __uVec.container()( globaldof );
                    }

                    else out << "0.0";
                }
            }

            out << "\n";
        }

        out << "$EndElementNodeData\n";
    }

    auto __ElmScal = __step->beginElementScalar();
    auto __ElmScal_end = __step->endElementScalar();

    for ( ; __ElmScal!=__ElmScal_end ; ++__ElmScal )
    {
        out << "\n$ElementData\n";

        element_scalar_type const& __u = __ElmScal->second;

        //uint __nbCompFieldGMSH;
        //    if (nodal_scalar_type::functionspace_type::is_scalar)         { __nbCompFieldGMSH=1; }
        //else if (nodal_scalar_type::functionspace_type::is_vectorial) { __nbCompFieldGMSH=3; }
        //else if (nodal_scalar_type::functionspace_type::is_tensor2)   { __nbCompFieldGMSH=9; }

        out << "1\n";//number of string tag
        out << "\"" << __ElmScal->first <<"\"\n";//a scalar node\n";
        out << "1\n";//number of real tag
        out << __step->time() << "\n";//"0.0\n";//the time value (0.0)
        out << "3\n";//number of integer tags:
        out << __step->index()-1 << "\n";//"0\n";//the time step (0; time steps always start at 0)
        out << "1\n";//n-component (1 is scalar) field
        out << mesh->numElements() << "\n";//number associated nodal values

        element_mesh_const_iterator elt_it = mesh->beginElement();
        element_mesh_const_iterator elt_en = mesh->endElement();

        if ( !__u.areGlobalValuesUpdated() )
            __u.updateGlobalValues();

        for ( ; elt_it!=elt_en ; ++elt_it )
        {
            globaldof = boost::get<0>( __u.functionSpace()->dof()->localToGlobal( elt_it->id(), 0, 0 ) ); //l,c

            out << number_markedfaces+elt_it->id()+1 << " " << /*__u( globaldof)*/__u.container()( globaldof ) << "\n";
        }

        out << "$EndElementData\n";
    }

}

#if 0
#if defined( FEELPP_INSTANTIATION_MODE )


//
// explicit instances
//

# define DIMS BOOST_PP_TUPLE_TO_LIST(3,(1,2,3))
# define ORDERS BOOST_PP_TUPLE_TO_LIST(5,(1,2,3,4,5))
# define ORDERS_FUN_GMSH BOOST_PP_TUPLE_TO_LIST(5,(1,2,3,4,5))

// exporter gmsh
# define FACTORY(LDIM,LORDER,ORDERFUN) template class ExporterGmsh<Mesh<Simplex<LDIM,LORDER,LDIM> >, ORDERFUN >;
# define FACTORY_OP(_, GDO) FACTORY GDO

BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_OP, 3, ( DIMS, ORDERS, ORDERS_FUN_GMSH ) )

#endif // FEELPP_INSTANTIATION_MODE
#endif
}
#endif // __EXPORTERGMSH_CPP

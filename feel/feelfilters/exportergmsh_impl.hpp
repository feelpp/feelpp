/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
ExporterGmsh<MeshType,N>::ExporterGmsh( WorldComm const& worldComm )
:
super( worldComm ),
M_element_type()
{
}
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
    DVLOG(2) << "[ExporterGmsh] checking if frequency is ok\n";

    if ( this->cptOfSave() % this->freq()  )
    {
        this->saveTimeSet();
        return;
    }

    DVLOG(2) << "[ExporterGmsh] frequency is ok\n";

    DVLOG(2) << "[ExporterGmsh] save()...\n";

    gmshSaveAscii();

    DVLOG(2) << "[ExporterGmsh] saving done\n";

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
    DVLOG(2) << "[gmshSaveascii] saving in gmsh ascii file format\n";

    timeset_const_iterator __ts_it = this->beginTimeSet();
    timeset_const_iterator __ts_en = this->endTimeSet();

    while ( __ts_it != __ts_en )
    {
        timeset_ptrtype __ts = *__ts_it;

        std::ostringstream __fname;

        /* If we want only one file, specify the same filename for each process */
        if(boption(_name="exporter.gmsh.merge") == true)
        {
            __fname << __ts->name()  //<< this->prefix() //this->path()
                    << "-" << this->worldComm().size()
                    << ".msh";
        }
        else
        {
            __fname << __ts->name()  //<< this->prefix() //this->path()
                    << "-" << this->worldComm().size() << "_" << this->worldComm().rank()
                    << ".msh";
        }

        /*std::string filename =  this->prefix()
          + __ts->name()
          + "-" + this->worldComm().size() + "_" + this->worldComm().rank()
          + ".msh";*/
        std::ofstream out;

        typename timeset_type::step_const_iterator __stepIt = __ts->beginStep();
        typename timeset_type::step_const_iterator __stepIt_end = __ts->endStep();
        __stepIt = boost::prior( __stepIt_end );

        if ( __stepIt != __stepIt_end )
        {
            step_ptrtype __step = *__stepIt;

            std::map<std::string, std::vector<double> > minMaxValues;

            if ( __step->isInMemory() )
            {
                auto nEltAndIndex = numberOfGlobalEltAndIndex( __step->mesh() );
                auto nGlobElement = nEltAndIndex.template get<0>();
                auto indexElementStart = nEltAndIndex.template get<1>();
                this->worldComm().barrier();

                /* Saving multiple files */
                if(boption(_name="exporter.gmsh.merge") == false)
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
                        DVLOG(2) << "cannot open " << __fname.str().c_str() << "\n";
                        exit( 0 );
                    }

                    DVLOG(2) << "[ExporterGmsh] saving model "
                                  << __ts->name() << " at time step "
                                  << __ts->index() << " in "
                                  << __fname.str() << "\n";

                    // save mesh only at first iteration
                    if ( __stepIt == __ts->beginStep() )
                    {
                        gmshSaveFormat( out );

                        gmshSavePhysicalNames( out, __step->mesh() );

                        auto pt_it = __step->mesh()->beginPointWithProcessId();
                        auto const pt_en = __step->mesh()->endPointWithProcessId();

                        gmshSaveNodesStart( out, __step->mesh(), std::distance(pt_it, pt_en) );
                        gmshSaveNodes( out,__step->mesh() );
                        gmshSaveNodesEnd( out, __step->mesh() );

                        auto eltOnProccess = elements( __step->mesh() );
                        auto elt_it = eltOnProccess.template get<1>();
                        auto elt_en = eltOnProccess.template get<2>();

                        auto allmarkedfaces = boundaryfaces( __step->mesh() );
                        auto face_it = allmarkedfaces.template get<1>();
                        auto face_end = allmarkedfaces.template get<2>();

                        gmshSaveElementsStart( out, std::distance(elt_it, elt_en) + std::distance(face_it, face_end) );
                        gmshSaveElements( out, __step->mesh(), indexElementStart );
                        gmshSaveElementsEnd( out );
                    }

                    //gmshSaveNodeData( out, __step);
                    gmshSaveElementNodeData( out, __step, indexElementStart);
                }
                /* saving data to one file */
                else
                {
                    DVLOG(2) << "[ExporterGmsh] saving model "
                                  << __ts->name() << " at time step "
                                  << __ts->index() << " in "
                                  << __fname.str() << "\n";

                    // save mesh only at first iteration
                    if ( __stepIt == __ts->beginStep() )
                    {
                        this->worldComm().barrier();
                        if(this->worldComm().isMasterRank())
                        {
                            out.open( __fname.str().c_str(), std::ios::out );
                            gmshSaveFormat( out );
                            gmshSavePhysicalNames( out, __step->mesh() );
                        }
                        
                        this->worldComm().barrier();
                        size_type nGlobPoint = numberOfGlobalPtAndIndex( __step->mesh() );
                        this->worldComm().barrier();

                        if(this->worldComm().isMasterRank())
                        {
                            gmshSaveNodesStart( out, __step->mesh(), nGlobPoint );
                            out.close();
                        }

                        for(int i = 0; i < this->worldComm().size(); i++)
                        {
                            if(i == this->worldComm().rank())
                            {
                                out.open( __fname.str().c_str(), std::ios::out | std::ios::app );
                                gmshSaveNodes( out,__step->mesh() );
                                out.close();
                            }
                            this->worldComm().barrier();
                        }

                        if(this->worldComm().isMasterRank())
                        {
                            out.open( __fname.str().c_str(), std::ios::out | std::ios::app );
                            gmshSaveNodesEnd( out, __step->mesh() );
                            gmshSaveElementsStart( out, nGlobElement );
                            out.close();
                        }
                        this->worldComm().barrier();

                        for(int i = 0; i < this->worldComm().size(); i++)
                        {
                            if(i == this->worldComm().rank())
                            {
                                out.open( __fname.str().c_str(), std::ios::out | std::ios::app );
                                gmshSaveElements( out, __step->mesh(), indexElementStart );
                                out.close();
                            }
                            this->worldComm().barrier();
                        }

                        if(this->worldComm().isMasterRank())
                        {
                            out.open( __fname.str().c_str(), std::ios::out | std::ios::app );
                            gmshSaveElementsEnd( out );
                            out.close();
                        }
                        this->worldComm().barrier();
                    }

                    for(int i = 0; i < this->worldComm().size(); i++)
                    {
                        if(i == this->worldComm().rank())
                        {
                            out.open( __fname.str().c_str(), std::ios::out | std::ios::app );
                            //gmshSaveNodeData( out, __step);
                            gmshSaveElementNodeData( out, __step, indexElementStart);
                            out.close();
                        }
                        this->worldComm().barrier();
                    }
                }

                // Correction of transfer function for scalar values through a geo script
                // only if we have more than 1 process
                if(this->worldComm().size() > 1)
                {
                    // The computation of the min/max, for readjusting the transfer function in Gmsh, 
                    // needs to be done with all processes available as it requires some communication
                    this->computeMinMax( __step, minMaxValues );

                    // the master rank then writes the corresponding geo file
                    if(this->worldComm().rank() == 0)
                    {
                        // open a geo file for output and correcting transfer function scale
                        std::ostringstream __geofname;
                        std::ostringstream __mshfname;

                        __geofname << __ts->name()  //<< this->prefix() //this->path()
                            << "-" << this->worldComm().size()
                            << ".geo";

                        std::ofstream geoout;
                        geoout.open(__geofname.str().c_str(), std::ios::out);

                        // merge the msh files, depending on the fact that we have 1 or several data files
                        if(boption(_name="exporter.gmsh.merge") == true)
                        {
                            __mshfname << __ts->name()  //<< this->prefix() //this->path()
                                << "-" << this->worldComm().size()
                                << ".msh";
                            geoout << "Merge \"" << __fname.str() << "\";" << std::endl; 
                        }
                        else
                        {
                            for(int i = 0; i < this->worldComm().size(); i++)
                            {
                                __mshfname.str("");
                                __mshfname << __ts->name()  //<< this->prefix() //this->path()
                                    << "-" << this->worldComm().size() << "_" << i
                                    << ".msh";
                                geoout << "Merge \"" << __mshfname.str() << "\";" << std::endl; 
                            }
                        }

                        /* if we have min/max values, we correct the scale of the view in Gmsh */
                        if(!minMaxValues.empty())
                        {
                            geoout << "nv = PostProcessing.NbViews-1;" << std::endl;
                            for(int i = 0; i < this->worldComm().size(); i++)
                            {
                                int j = 0;
                                for(std::map<std::string, std::vector<double> >::iterator it = minMaxValues.begin(); it!=minMaxValues.end() ; it++, j++)
                                {
                                    // make visible only first function values, declutter the interface
                                    geoout << "View[nv-" << ((this->worldComm().size() - 1 - i) * minMaxValues.size() + (minMaxValues.size() - 1 - j)) << "].";
                                    if(j == 0)
                                    {
                                        geoout << "Visible=1;" << std::endl;
                                    }
                                    else
                                    {
                                        geoout << "Visible=0;" << std::endl;
                                    }

                                    // if we have min-max values
                                    // we correct the range for the transfer function (for scalar values)
                                    if(it->second.size() == 2)
                                    {   
                                        geoout << "View[nv-" << ((this->worldComm().size() - 1 - i) * minMaxValues.size() + (minMaxValues.size() - 1 - j)) 
                                            << "].RangeType=2;" << std::endl;
                                        geoout << "View[nv-" << ((this->worldComm().size() - 1 - i) * minMaxValues.size() + (minMaxValues.size() - 1 - j))
                                            << "].CustomMin=" << it->second[0] << ";" << std::endl;
                                        geoout << "View[nv-" << ((this->worldComm().size() - 1 - i) * minMaxValues.size() + (minMaxValues.size() - 1 - j)) 
                                            << "].CustomMax=" << it->second[1] << ";" << std::endl;
                                    }
                                }
                            }
                        }

                        geoout.close();

                        /* If onelab is enabled, we register this filename to be loaded */
                        if(ioption(_name="onelab.enable" ) == 2)
                        {
                            Environment::olLoadInGmsh(__geofname.str());
                        }
                    }
                }
                /* if we have only one process, we register only the msh file for Onelab */
                else
                {
                    /* If onelab is enabled, we register the msh file to be loaded */
                    if(ioption(_name="onelab.enable" ) == 2)
                    {
                        Environment::olLoadInGmsh(__fname.str());
                    }
                }

            }
        }

        ++__ts_it;

    }
}

template<typename MeshType, int N>
void
ExporterGmsh<MeshType,N>::saveMesh( std::string const& filename, mesh_ptrtype mesh, bool parametric ) const
{

    if (  this->worldComm().rank() == 0 )
    {
        std::ofstream out( filename.c_str(), std::ios::out );

        if ( out.fail() )
        {
            LOG(INFO) << "[ExporterGmsh::SaveMesh] cannot open " << filename << "\n";
            exit( 0 );
        }

        gmshSaveFormat( out );
        gmshSavePhysicalNames( out, mesh );
        out.close();
    }

    this->worldComm().barrier();

    //-----------------------------------------------------------------//

    size_type nGlobPoint = numberOfGlobalPtAndIndex( mesh );

    if (  this->worldComm().rank() == 0 )
    {
        std::ofstream out( filename.c_str(), std::ios::app );
        gmshSaveNodesStart( out, mesh, nGlobPoint, parametric );
        out.close();
    }

    for ( int therank=0; therank<this->worldComm().size(); ++therank )
    {
        if ( therank == this->worldComm().rank() )
        {
            std::ofstream out( filename.c_str(), std::ios::app );
            gmshSaveNodes( out, mesh, parametric );
            out.close();
        }

        this->worldComm().barrier();
    }

    if (  this->worldComm().rank() == 0 )
    {
        std::ofstream out( filename.c_str(), std::ios::app );
        gmshSaveNodesEnd( out, mesh, parametric );
        out.close();
    }

    //-----------------------------------------------------------------//

    auto nEltAndIndex = numberOfGlobalEltAndIndex( mesh );
    auto nGlobElement = nEltAndIndex.template get<0>();
    auto indexElementStart = nEltAndIndex.template get<1>();

    if (  this->worldComm().rank() == 0 )
    {
        std::ofstream out( filename.c_str(), std::ios::app );
        gmshSaveElementsStart( out, nGlobElement );
        out.close();
    }

    for ( int therank=0; therank<this->worldComm().size(); ++therank )
    {
        if ( therank == this->worldComm().rank() )
        {
            std::ofstream out( filename.c_str(), std::ios::app );
            gmshSaveElements( out, mesh, indexElementStart );
            out.close();
        }

        this->worldComm().barrier();
    }

    if (  this->worldComm().rank() == 0 )
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
    if ( mesh->markerNames().size() == 0 ) return;

    // save Physical Names
    out << "$PhysicalNames\n";
    out << mesh->markerNames().size() << "\n";
    // physical-dimension physical-number "physical-name" ...
    BOOST_FOREACH( auto data, mesh->markerNames() )
    {
        out << data.second[1] << " "
            << data.second[0] << " "
            << "\"" << data.first << "\"" << "\n";
    }
    out << "$EndPhysicalNames\n";

}

template<typename MeshType, int N>
size_type
ExporterGmsh<MeshType,N>::numberOfGlobalPtAndIndex( mesh_ptrtype mesh ) const
{
    size_type nPointToWriteOnProcess=0;

    auto itPt = mesh->beginPointWithProcessId();
    auto const enPt = mesh->endPointWithProcessId();

    
    /* If we want only one file, specify the same filename for each process */
    if(boption(_name="exporter.gmsh.merge") == true)
    {
        for ( ; itPt!=enPt ; ++itPt )
        {
            if ( itPt->isLinkedToOtherPartitions() )
            {
                // add if the processId() is the min rank
                if (itPt->processId() < *std::min_element( itPt->neighborPartitionIds().begin(),itPt->neighborPartitionIds().end() ) )
                    ++nPointToWriteOnProcess;
            }
            else ++nPointToWriteOnProcess;
        }
    }
    // If we want several files, we need to copy the point that are shared with other partitions
    // even if the process id is not the minimal (otherwise we would be missing points)
    else
    {
        //nPointToWriteOnProcess = std::distance(itPt, enPt);
        for ( ; itPt!=enPt ; ++itPt )
        {
            ++nPointToWriteOnProcess;
        }
    }

    //size_type local_numberPoints = mesh->numPoints();
    size_type local_numberPoints = nPointToWriteOnProcess;
    size_type global_numberPoints=0;//local_numberPoints;

    mpi::all_reduce( this->worldComm(),
                     local_numberPoints,
                     global_numberPoints,
                     std::plus<size_type>() );

#if 0
    std::vector<size_type> all_localnumberPoint;

    mpi::all_gather( this->worldComm(),
                     local_numberPoints,
                     all_localnumberPoint );

    size_type indexPtStart = 0;

    for ( int i=0; i<this->worldComm().rank(); ++i )
        indexPtStart+=all_localnumberPoint[i];

    return boost::make_tuple( global_numberPoints,indexPtStart);
#endif
    return global_numberPoints;
}

template<typename MeshType, int N>
void
ExporterGmsh<MeshType,N>::gmshSaveNodesStart( std::ostream& out, mesh_ptrtype mesh, size_type nGlobPt, bool parametric ) const
{
    if ( parametric && mesh->isParametric() )
        out << "$ParametricNodes\n";

    else
        out << "$Nodes\n";

    // Save number of Nodes
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
ExporterGmsh<MeshType,N>::gmshSaveNodes( std::ostream& out, mesh_ptrtype mesh, bool parametric ) const
{
    auto pt_it = mesh->beginPointWithProcessId();
    auto const pt_en = mesh->endPointWithProcessId();

    for ( ; pt_it!=pt_en ; ++pt_it )
    {
        // if we want only one file, we discard nodes shared by different processes
        if(boption(_name="exporter.gmsh.merge") == true)
        {
            if ( pt_it->isLinkedToOtherPartitions() )
            {
                // add if the processId() is the min rank
                if ( pt_it->processId() > *std::min_element( pt_it->neighborPartitionIds().begin(),pt_it->neighborPartitionIds().end() ) )
                    continue;
            }
        }
        // otherwise we put all the nodes in each file (to have a complete nodeset)

        // warning add 1 to the id in order to be sure that all gmsh id >0
        out << pt_it->id()+1
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
    size_type number_elements= std::distance( eltOnProccess.template get<1>(), eltOnProccess.template get<2>() );

    auto local_numberElements = number_markedfaces+number_elements;
    auto global_numberElements=local_numberElements;

    mpi::all_reduce( this->worldComm(),
                     local_numberElements,
                     global_numberElements,
                     std::plus<size_type>() );

    std::vector<size_type> all_localnumberElements;

    mpi::all_gather( this->worldComm(),
                     local_numberElements,
                     all_localnumberElements );

    size_type indexEltStart = 0;

    for ( int i=0; i<this->worldComm().rank(); ++i )
        indexEltStart+=all_localnumberElements[i];

    return boost::make_tuple( global_numberElements,indexEltStart );
}

template<typename MeshType, int N>
void
ExporterGmsh<MeshType,N>::gmshSaveElementsStart( std::ostream& out, size_type nGlobElt ) const
{
    out << "$Elements\n";

    // write the count of elements
    out << nGlobElt << "\n";
}

template<typename MeshType, int N>
void
ExporterGmsh<MeshType,N>::gmshSaveElementsEnd( std::ostream& out ) const
{
    out << "$EndElements\n";
}

template<typename MeshType, int N>
void
ExporterGmsh<MeshType,N>::gmshSaveElements( std::ostream& out, mesh_ptrtype mesh, size_type indexEltStart ) const
{
    auto eltOnProccess = elements( mesh );
    auto elt_it = eltOnProccess.template get<1>();
    auto elt_en = eltOnProccess.template get<2>();
    std::map<int,int> pids;
    int pid = 0;
    std::for_each( elt_it, elt_en, [&pid,&pids]( typename MeshType::element_type const& e ){ pids[e.id()]=pid++; });

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

    // count the number of faces
    int number_markedfaces = std::distance(allmarkedfaces.template get<1>(), allmarkedfaces.template get<2>());

    for ( ; face_it != face_end; ++face_it )
    {
        // elm-number elm-type number-of-tags < tag > ... node-number-list
        /*
        if(boption(_name="exporter.gmsh.merge") == true)
        {
        */
            out<< elem_number++ <<" ";
            /*
        }
        else
        {
            out<< indexEltStart + face_it->id()+1 <<" ";
        }
        */
        out << ordering_face.type();
        // number-of-tags < tag >

        if ( FEELPP_GMSH_FORMAT_VERSION==std::string( "2.1" ) )
        {
            // out<<" 2 " << face_it->marker().value() << " " << face_it->marker2().value();
            out<<" 3 " << face_it->marker().value() << " " << face_it->marker2().value() << " " << face_it->processId()+1;
        }

        else if ( FEELPP_GMSH_FORMAT_VERSION==std::string( "2.2" ) )
        {
            uint16_type nbTag = 3;
            if ( boption(_name="partition.linear" ) )
            {
                nbTag += 1;
            }
            else
                nbTag += face_it->numberOfPartitions();
            out << " " << nbTag
                << " " << face_it->marker().value()
                << " " << face_it->marker2().value();


            if ( boption(_name="partition.linear" ) )
            {
                out << " " <<  1
                    << " " << pids[face_it->element0().id()]+1;
            }
            else
            {
                out << " " << face_it->numberOfPartitions()
                    << " " << face_it->processId()+1;
                for ( size_type i=0 ; i<face_it->numberOfNeighborPartitions(); ++i )
                    out << " " << -( face_it->neighborPartitionIds()[i]+1 );
            }
        }

        // node-number-list
        for ( uint16_type p=0; p<face_type::numPoints; ++p )
        {
            //DCHECK( mapPointsId.find(face_it->point( ordering_face.fromGmshId( p ) ).id()) != mapPointsId.end() ) << "invalid point id\n";
            //out << " " << mapPointsId.find(face_it->point( ordering_face.fromGmshId( p ) ).id())->second;
            out << " " << face_it->point( ordering_face.fromGmshId( p ) ).id()+1;
        }

        out<<"\n";
    } // faces



    for ( ; elt_it != elt_en; ++elt_it, ++pid )
    {
        /*
        if(boption(_name="exporter.gmsh.merge") == true)
        {
        */
            out << elem_number++ <<" ";
            /*
        }
        else
        {
            out << indexEltStart + number_markedfaces + elt_it->id()+1 <<" ";
        }
    */
        out << ordering.type();

        if ( FEELPP_GMSH_FORMAT_VERSION==std::string( "2.1" ) )
        {
            //out<<" 2 " << elt_it->marker().value() << " " << elt_it->marker2().value();
            out<<" 3 " << elt_it->marker().value() << " " << elt_it->marker2().value() << " " << elt_it->processId()+1;
        }

        else if ( FEELPP_GMSH_FORMAT_VERSION== std::string( "2.2" ) )
        {
            std::vector<int> f;
            uint16_type nbTag = 3;

            if ( boption(_name="partition.linear" ) )
            {
                for ( size_type i=0 ; i< elt_it->nNeighbors(); ++i )
                    if ( elt_it->neighbor(i).first != invalid_size_type_value )
                        f.push_back( elt_it->neighbor(i).first );

                nbTag+=f.size()+1;
            }
            else
                nbTag+=elt_it->numberOfPartitions();
            out << " " << nbTag
                << " " << elt_it->marker().value()
                << " " << elt_it->marker2().value();

            if ( boption(_name="partition.linear" ) )
            {

                out << " " << f.size()+1 << " " << pids[elt_it->id()]+1;
                for( auto i : f )
                    out << " " << -( pids[i]+1 );
            }
            else
            {
                out << " " << elt_it->numberOfPartitions()
                    << " " << elt_it->processId()+1;

                for ( size_type i=0 ; i<elt_it->numberOfNeighborPartitions(); ++i )
                    out << " " << -( elt_it->neighborPartitionIds()[i]+1 );
            }
        }

        for ( uint16_type p=0; p<element_type::numPoints; ++p )
        {
            //DCHECK ( mapPointsId.find(elt_it->point( ordering.fromGmshId( p ) ).id()) != mapPointsId.end() ) << "invalid point id\n";
            //std::cout << "index " << p << " -> " << ordering.fromGmshId(p) << " -> " << elt_it->point( ordering.fromGmshId(p) ).id()+1 << " : " << elt_it->point( ordering.fromGmshId(p) ).node() << "\n";
            //out << " " << mapPointsId.find(elt_it->point( ordering.fromGmshId( p ) ).id())->second;
            out << " " << elt_it->point( ordering.fromGmshId( p ) ).id()+1;
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
ExporterGmsh<MeshType,N>::computeMinMax(step_ptrtype __step, std::map<std::string, std::vector<double> > & minMaxValues) const
{
    typedef typename mesh_type::element_const_iterator element_mesh_const_iterator;

    typedef typename step_type::nodal_scalar_type nodal_scalar_type;
    typedef typename step_type::nodal_scalar_const_iterator nodal_scalar_const_iterator;

    typedef typename step_type::nodal_vector_type nodal_vectorial_type;
    typedef typename step_type::nodal_vector_const_iterator nodal_vectorial_const_iterator;

    typedef typename step_type::element_scalar_type element_scalar_type;

    mesh_ptrtype mesh = __step->mesh();

    uint16_type nLocalDof;
    size_type globaldof;
    typedef typename MeshType::element_type element_type;
    GmshOrdering<element_type> ordering;

    /* relinearize elemnent ids in the case when we want to output all the data in one file */
    auto allmarkedfaces = boundaryfaces( mesh );
    auto face_it = allmarkedfaces.template get<1>();
    auto face_en = allmarkedfaces.template get<2>();
    int number_markedfaces= std::distance( face_it, face_en );

    auto elts = elements( mesh );
    auto elt_it = elts.template get<1>();
    auto elt_en = elts.template get<2>();

    nodal_scalar_const_iterator __varScal = __step->beginNodalScalar();
    nodal_scalar_const_iterator __varScal_end = __step->endNodalScalar();

    for ( ; __varScal!=__varScal_end ; ++__varScal )
    {
        nodal_scalar_type const& __u = __varScal->second;

        // record min-max value for function
        // check if record already exists
        if(minMaxValues.empty() || minMaxValues.find(__varScal->first) == minMaxValues.end())
        {
            minMaxValues[__varScal->first].push_back(__u.min());
            minMaxValues[__varScal->first].push_back(__u.max());
        }
        else
        {
            if(minMaxValues[__varScal->first][0] > __u.min())
            {
                minMaxValues[__varScal->first][0] = __u.min();
            }
            if(minMaxValues[__varScal->first][1] < __u.max())
            {
                minMaxValues[__varScal->first][1] = __u.max();
            }
        }
    }

    nodal_vectorial_const_iterator __varVec = __step->beginNodalVector();
    nodal_vectorial_const_iterator __varVec_end = __step->endNodalVector();

    for ( ; __varVec!=__varVec_end ; ++__varVec )
    {
        nodal_vectorial_type const& __uVec = __varVec->second;

        uint16_type nComponents = __uVec.nComponents;

        element_mesh_const_iterator elt_it;
        element_mesh_const_iterator elt_en;
        boost::tie( boost::tuples::ignore, elt_it, elt_en ) = elements( mesh );
        
        // record min-max value for function
        // check if record already exists
        if(minMaxValues.empty() || minMaxValues.find(__varVec->first) == minMaxValues.end())
        {
            nLocalDof = nodal_vectorial_type::functionspace_type::basis_type::nLocalDof;

            for(int i = 0; i < nLocalDof * 3; i++)
            {
                minMaxValues[__varVec->first].push_back(0.0);
                minMaxValues[__varVec->first].push_back(0.0);
            }
        }
        
        /* need to update min/max for vectorial data */
        /*
        for (; elt_it!=elt_en ; ++elt_it )
        {
            nLocalDof = nodal_vectorial_type::functionspace_type::basis_type::nLocalDof;

            for ( uint16_type l = 0; l < nLocalDof; ++l )
            {
                uint16_type gmsh_l = ordering.fromGmshId( l );

                for ( uint16_type c = 0; c < 3; ++c )
                {
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
        */
    }

    auto __ElmScal = __step->beginElementScalar();
    auto __ElmScal_end = __step->endElementScalar();

    for ( ; __ElmScal!=__ElmScal_end ; ++__ElmScal )
    {
        element_scalar_type const& __u = __ElmScal->second;

        element_mesh_const_iterator elt_it = mesh->beginElement();
        element_mesh_const_iterator elt_en = mesh->endElement();

        if ( !__u.areGlobalValuesUpdated() )
            __u.updateGlobalValues();
        
        // record min-max value for function
        // check if record already exists
        if(minMaxValues.empty() || minMaxValues.find(__ElmScal->first) == minMaxValues.end())
        {
            minMaxValues[__ElmScal->first].push_back(0.0);
            minMaxValues[__ElmScal->first].push_back(0.0);
        }
        
#if 0
        for ( ; elt_it!=elt_en ; ++elt_it )
        {
            globaldof = boost::get<0>( __u.functionSpace()->dof()->localToGlobal( elt_it->id(), 0, 0 ) ); //l,c

            // either use the relinearized version for one file dataset or classic for one file per process
            /*
            if(boption(_name="exporter.gmsh.merge") == true)
            {
            */
                out << elt_pids[elt_it->id()] << " " << /*__u( globaldof)*/__u.container()( globaldof ) << "\n";
            //}
            //else
            //{
                //out << indexEltStart + number_markedfaces+elt_it->id()+1 << " " << [>__u( globaldof)<]__u.container()( globaldof ) << "\n";
            //}
        }
#endif
    }

}


template<typename MeshType, int N>
void
ExporterGmsh<MeshType,N>::gmshSaveElementNodeData( std::ostream& out,
        step_ptrtype __step, size_type indexEltStart) const
{
    typedef typename mesh_type::element_const_iterator element_mesh_const_iterator;

    typedef typename step_type::nodal_scalar_type nodal_scalar_type;
    typedef typename step_type::nodal_scalar_const_iterator nodal_scalar_const_iterator;

    typedef typename step_type::nodal_vector_type nodal_vectorial_type;
    typedef typename step_type::nodal_vector_const_iterator nodal_vectorial_const_iterator;

    typedef typename step_type::element_scalar_type element_scalar_type;

    mesh_ptrtype mesh = __step->mesh();

    uint16_type nLocalDof;
    size_type globaldof;
    typedef typename MeshType::element_type element_type;
    GmshOrdering<element_type> ordering;

    /* relinearize elemnent ids in the case when we want to output all the data in one file */
    auto allmarkedfaces = boundaryfaces( mesh );
    auto face_it = allmarkedfaces.template get<1>();
    auto face_en = allmarkedfaces.template get<2>();
    int number_markedfaces= std::distance( face_it, face_en );
    std::map<int,int> face_pids;
    int pid = indexEltStart + 1;
    //if(boption(_name="exporter.gmsh.merge") == true)
    //{
        std::for_each( face_it, face_en, [&pid, &face_pids]( typename MeshType::face_type const& e ){ face_pids[e.id()]=pid++; });
    //}

    auto elts = elements( mesh );
    auto elt_it = elts.template get<1>();
    auto elt_en = elts.template get<2>();
    std::map<int,int> elt_pids;
    //if(boption(_name="exporter.gmsh.merge") == true)
    //{
        std::for_each( elt_it, elt_en, [&pid,&elt_pids]( typename MeshType::element_type const& e ){ elt_pids[e.id()]=pid++; });
    //}

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
        out << "\"" << this->worldComm().size() << "_" << this->worldComm().rank() << "-" << __varScal->first <<"\"\n";//a scalar node\n";
        out << "1\n";//number of real tag
        out << __step->time() << "\n";//"0.0\n";//the time value (0.0)
        out << "3\n";//number of integer tags:
        out << __step->index()-1 << "\n";//"0\n";//the time step (0; time steps always start at 0)
        out << "1\n";//n-component (1 is scalar) field
        //out << mesh->numElements() << "\n";//number associated nodal values

        element_mesh_const_iterator elt_it;
        element_mesh_const_iterator elt_en;
        boost::tie( boost::tuples::ignore, elt_it, elt_en ) = elements( mesh );

        out << std::distance(elt_it, elt_en) << "\n";

        if ( !__u.areGlobalValuesUpdated() )
            __u.updateGlobalValues();

        for ( ; elt_it!=elt_en ; ++elt_it )
        {
            // either use the relinearized version for one file dataset or classic for one file per process
            /*
            if(boption(_name="exporter.gmsh.merge") == true)
            {
            */
                out << elt_pids[elt_it->id()];
                /*
            }
            else
            {
                out << indexEltStart + number_markedfaces + elt_it->id()+1;
            }
            */
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
#if 0
                if ( ublas::norm_2( boost::get<0>( __u.functionSpace()->dof()->dofPoint( globaldof ) )-elt_it->point( ordering.fromGmshId( l ) ).node() ) > 1e-10 )
                {
                    std::cout << "------------------------------------------------------------\n";
                    std::cout << "invalid dof/mesh points\n";
                    std::cout << "dof global id:" << globaldof << " | local id:" << gmsh_l << "\n";
                    std::cout << "point global id:" <<  elt_it->point( ordering.fromGmshId( l ) ).id() << " | local id:" << gmsh_l << "\n";
                    std::cout << "node dof:  " << boost::get<0>( __u.functionSpace()->dof()->dofPoint( globaldof ) ) << "\n";
                    std::cout << "node element:  " << elt_it->point( ordering.fromGmshId( l ) ).node() << "\n";
                }
#endif
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
        out << "\"" << this->worldComm().size() << "_" << this->worldComm().rank() << "-" << __varVec->first <<"\"\n";//"a vectorial field\n";
        out << "1\n";//number of real tag
        out << __step->time() << "\n";//"0.0\n";//the time value (0.0)
        out << "3\n";//number of integer tags:
        out << __step->index() << "\n";//"0\n";//the time step (0; time steps always start at 0)
        out << "3\n";//n-component (3 is vectorial) field
        //out << mesh->numElements() << "\n";//number associated nodal values

        element_mesh_const_iterator elt_it;
        element_mesh_const_iterator elt_en;
        boost::tie( boost::tuples::ignore, elt_it, elt_en ) = elements( mesh );
        
        out << std::distance(elt_it, elt_en) << "\n";//number associated nodal values

        for (; elt_it!=elt_en ; ++elt_it )
        {
            // either use the relinearized version for one file dataset or classic for one file per process
            /*
            if(boption(_name="exporter.gmsh.merge") == true)
            {
            */
                out << elt_pids[elt_it->id()];
                /*
            }
            else
            {
                out << indexEltStart + number_markedfaces + elt_it->id()+1;
            }
            */
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
        out << "\"" << this->worldComm().size() << "_" << this->worldComm().rank() << "-" << __ElmScal->first <<"\"\n";//a scalar node\n";
        out << "1\n";//number of real tag
        out << __step->time() << "\n";//"0.0\n";//the time value (0.0)
        out << "3\n";//number of integer tags:
        out << __step->index()-1 << "\n";//"0\n";//the time step (0; time steps always start at 0)
        out << "1\n";//n-component (1 is scalar) field
        //out << mesh->numElements() << "\n";//number associated nodal values

        element_mesh_const_iterator elt_it;
        element_mesh_const_iterator elt_en;
        boost::tie( boost::tuples::ignore, elt_it, elt_en ) = elements( mesh );

        out << std::distance(elt_it, elt_en) << "\n";

        if ( !__u.areGlobalValuesUpdated() )
            __u.updateGlobalValues();

        for ( ; elt_it!=elt_en ; ++elt_it )
        {
            globaldof = boost::get<0>( __u.functionSpace()->dof()->localToGlobal( elt_it->id(), 0, 0 ) ); //l,c

            // either use the relinearized version for one file dataset or classic for one file per process
            /*
            if(boption(_name="exporter.gmsh.merge") == true)
            {
            */
                out << elt_pids[elt_it->id()] << " " << /*__u( globaldof)*/__u.container()( globaldof ) << "\n";
            //}
            //else
            //{
                //out << indexEltStart + number_markedfaces+elt_it->id()+1 << " " << [>__u( globaldof)<]__u.container()( globaldof ) << "\n";
            //}
        }

        out << "$EndElementData\n";
    }

}



template<typename MeshType, int N>
template<typename ConvexType>
void
ExporterGmsh<MeshType,N>::gmshSaveOneElementAsMesh( std::string const& filename,
                                                    typename mesh_type::element_type::super const& elt,
                                                    PointSet<ConvexType,typename MeshType::value_type> const& ptset ) const
{
    std::ofstream out( filename.c_str(), std::ios::out );
    gmshSaveFormat( out );

    const uint32_type nPointInPtSet = ptset.nPoints();

    out << "$Nodes\n";
    out << ptset.nPoints()+elt.nPoints() << "\n";//number points

    for ( uint32_type i = 0; i < nPointInPtSet; ++i )
    {
        auto const thepoint = ptset.point(i);
        out << i+1
            << " "  << std::setw( 20 ) << std::setprecision( 16 ) << thepoint(0);
        if ( thepoint.size() >= 2 )
            out << " "  << std::setw( 20 ) << std::setprecision( 16 ) << thepoint(1);
        else
            out << " 0";
        if ( thepoint.size() >= 3 )
            out << " "  << std::setw( 20 ) << std::setprecision( 16 ) << thepoint(2);
        else
            out << " 0";
        out << "\n";
    }


    for ( int i = 0; i < elt.nPoints(); ++i )
    {
        out << nPointInPtSet+i+1
            << " "  << std::setw( 20 ) << std::setprecision( 16 ) << elt.point(i).node()[0];
        if ( mesh_type::nRealDim >= 2 )
            out << " "  << std::setw( 20 ) << std::setprecision( 16 ) << elt.point(i).node()[1];
        else
            out << " 0";
        if ( mesh_type::nRealDim >= 3 )
            out << " "  << std::setw( 20 ) << std::setprecision( 16 ) << elt.point(i).node()[2];
        else
            out << " 0";
        out << "\n";
    }

    out << "$EndNodes\n";

    out << "$Elements\n"
        << nPointInPtSet+1 << "\n";// number element
    for ( uint32_type i = 0; i < nPointInPtSet; ++i )
    {
        out << i+1 << " " << GMSH_ENTITY::GMSH_POINT
            << " 2 0 " << i << " " // add elementary tag as feel id
            << i+1 << "\n";
    }

    typedef typename MeshType::element_type element_type;
    GmshOrdering<element_type> ordering;
    out << nPointInPtSet+1 << " "; // id element
    out << ordering.type();

    if ( FEELPP_GMSH_FORMAT_VERSION==std::string( "2.1" ) )
    {
        out<<" 3 " << elt.marker().value() << " " << elt.marker2().value() << " " << elt.processId()+1;
    }
    else if ( FEELPP_GMSH_FORMAT_VERSION== std::string( "2.2" ) )
    {
        uint16_type nbTag = 3 + elt.numberOfPartitions();
        out << " " << nbTag
            << " " << elt.marker().value()
            << " " << elt.marker2().value()
            << " " << elt.numberOfPartitions()
            << " " << elt.processId()+1;

        for ( size_type i=0 ; i<elt.numberOfNeighborPartitions(); ++i )
            out << " " << -( elt.neighborPartitionIds()[i]+1 );
    }

    for ( uint16_type p=0; p<element_type::numPoints; ++p )
    {
        out << " " << nPointInPtSet +ordering.fromGmshId( p ) + 1;
    }

    out<<"\n";
    out << "$EndElements\n";
    out.close();

}

template<typename MeshType, int N>
template<typename ConvexRefType, typename ConvexPtSetType>
void
ExporterGmsh<MeshType,N>::gmshSaveOneElementAsMesh( std::string const& filename,
                                                    Reference<ConvexRefType,ConvexRefType::nDim,ConvexRefType::nOrder,ConvexRefType::nRealDim >  const& elt,
                                                    PointSet<ConvexPtSetType,typename MeshType::value_type> const& ptset ) const
{
    std::ofstream out( filename.c_str(), std::ios::out );
    gmshSaveFormat( out );

    const uint32_type nPointInPtSet = ptset.nPoints();

    out << "$Nodes\n";
    out << ptset.nPoints()+elt.nPoints() << "\n";//number points

    for ( uint32_type i = 0; i < nPointInPtSet; ++i )
    {
        auto const thepoint = ptset.point(i);
        out << i+1
            << " "  << std::setw( 20 ) << std::setprecision( 16 ) << thepoint(0);
        if ( thepoint.size() >= 2 )
            out << " "  << std::setw( 20 ) << std::setprecision( 16 ) << thepoint(1);
        else
            out << " 0";
        if ( thepoint.size() >= 3 )
            out << " "  << std::setw( 20 ) << std::setprecision( 16 ) << thepoint(2);
        else
            out << " 0";
        out << "\n";
    }

    for ( int i = 0; i < elt.nPoints(); ++i )
    {
        out << nPointInPtSet+i+1
            << " "  << std::setw( 20 ) << std::setprecision( 16 ) << elt.point(i)(0);
        if ( mesh_type::nRealDim >= 2 )
            out << " "  << std::setw( 20 ) << std::setprecision( 16 ) << elt.point(i)(1);
        else
            out << " 0";
        if ( mesh_type::nRealDim >= 3 )
            out << " "  << std::setw( 20 ) << std::setprecision( 16 ) << elt.point(i)(2);
        else
            out << " 0";
        out << "\n";
    }
    out << "$EndNodes\n";


    out << "$Elements\n"
        << nPointInPtSet+1 <<"\n";// number element

    for ( uint32_type i = 0; i < nPointInPtSet; ++i )
    {
        out << i+1 << " " << GMSH_ENTITY::GMSH_POINT
            << " 2 0 " << i << " " // add elementary tag as feel id
            << i+1 << "\n";
    }

    typedef typename MeshType::element_type element_type;
    GmshOrdering<element_type> ordering;
    out << nPointInPtSet+1 << " "; // id element
    out << ordering.type();
    out << " 2 0 0 ";
    for ( uint16_type p=0; p<element_type::numPoints; ++p )
    {
        out << " " << nPointInPtSet + ordering.fromGmshId( p ) + 1;
    }
    out<<"\n";
    out << "$EndElements\n";

    out.close();
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

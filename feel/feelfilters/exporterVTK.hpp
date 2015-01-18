/* -* -mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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
   \file ExporterVTK.hpp
   \author Alexandre Ancel <alexandre.ancel@cemosis.fr>
   \date 2014-11-13
 */
#ifndef __ExporterVTK_H
#define __ExporterVTK_H 1

#if defined(FEELPP_HAS_VTK)

#include <iostream>
#include <fstream>

#include <feel/feelcore/debug.hpp>
#include <feel/feelfilters/exporter.hpp>

#include <cstring>
#include <libxml/parser.h>
#include <libxml/tree.h>

#include <vtkSmartPointer.h>
#include <vtkCellType.h>
#include <vtkIntArray.h>
#include <vtkStringArray.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>

#include <vtkInformation.h>
#include <vtkVertex.h>
#include <vtkLine.h>
#include <vtkPolyLine.h>
#include <vtkQuadraticTriangle.h>
#include <vtkQuad.h>
#include <vtkQuadraticQuad.h>
#include <vtkTetra.h>
#include <vtkQuadraticTetra.h>
#include <vtkHexahedron.h>
#include <vtkQuadraticHexahedron.h>
#include <vtkTriangle.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkXMLMultiBlockDataWriter.h>

/* Only use MPI when we have vtk 5.8+ */
/* features initializing MPI using an external context a missing in 5.8- */
/* but lets aim for the latest major version 6 to reduce the complexity */
#if VTK_MAJOR_VERSION >= 6 && defined(VTK_HAS_PARALLEL)
#include <vtkMPI.h>
#include <vtkMPIController.h>
#include <vtkMPICommunicator.h>
#include <vtkXMLPMultiBlockDataWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>

#if defined(FEELPP_VTK_INSITU_ENABLED)
#include <vtkCPProcessor.h>
#include <vtkCPPipeline.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#endif

#endif


namespace Feel
{
/**
 * \class ExporterVTK
 * \brief Export to VTK format
 *
 * \ingroup Exporter
 * @author Alexandre Ancel
 */
template<typename MeshType, int N>
class ExporterVTK
    :
public Exporter<MeshType,N>
{
public:

    /** @name Typedefs
     */
    //@{

    typedef MeshType mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef Exporter<MeshType,N> super;
    typedef typename mesh_type::value_type value_type;
    typedef typename super::timeset_type timeset_type;
    typedef typename super::timeset_ptrtype timeset_ptrtype;
    typedef typename super::timeset_const_iterator timeset_const_iterator;

    typedef typename timeset_type::step_type step_type;
    typedef typename timeset_type::step_ptrtype step_ptrtype;
    typedef typename timeset_type::step_const_iterator step_const_iterator;

    /* Use the vtkUnstructuredGrid type to store data */
    /* as it handles all the element type needed */
    typedef vtkUnstructuredGrid vtkout_type;
    typedef vtkUnstructuredGridWriter vtkoutwriter_type;

    /* Compute face type from mesh parameters */
    typedef typename 
    /* face type */
    /* if (Mdim == 1) */
    mpl::if_<mpl::equal_to<mpl::int_<MeshType::nDim>, mpl::int_<1> >,
        mpl::identity<vtkVertex>,
        /* if (Mdim == 2) */
        typename mpl::if_<mpl::equal_to<mpl::int_<MeshType::nDim>, mpl::int_<2> >,
            /* if(MShape == SHAPE_TRIANGLE) */    
            mpl::identity<vtkLine>,
            /* if (Mdim == 3) */
            typename mpl::if_<mpl::equal_to<mpl::int_<MeshType::nDim>, mpl::int_<3> >,
                /* if(MShape == SHAPE_TETRA) */    
                typename mpl::if_<mpl::equal_to<mpl::int_<MeshType::Shape>, mpl::size_t<SHAPE_TETRA> >,
                    mpl::identity<vtkTriangle>,
                    mpl::identity<vtkQuad>
                >::type,
                /* We should normally not reach this case */
                /* anyway we set a default vtkTriangle for face type */
                mpl::identity<vtkTriangle>
            >::type
        >::type
    >::type::type vtkface_type;

    /* Compute element type from the parameters */
    typedef typename 
    /* if (Mdim == 1) */
    mpl::if_<mpl::equal_to<mpl::int_<MeshType::nDim>, mpl::int_<1> >,
        mpl::identity<vtkLine>,
        /* if (Mdim == 2) */
        typename mpl::if_<mpl::equal_to<mpl::int_<MeshType::nDim>, mpl::int_<2> >,
            /* if(MShape == SHAPE_TRIANGLE) */    
            typename mpl::if_<mpl::equal_to<mpl::int_<MeshType::Shape>, mpl::size_t<SHAPE_TRIANGLE> >,
                mpl::identity<vtkTriangle>,
                mpl::identity<vtkQuad>
            >::type,
            /* if (Mdim == 3) */
            typename mpl::if_<mpl::equal_to<mpl::int_<MeshType::nDim>, mpl::int_<3> >,
                /* if(MShape == SHAPE_TETRA) */    
                typename mpl::if_<mpl::equal_to<mpl::int_<MeshType::Shape>, mpl::size_t<SHAPE_TETRA> >,
                    mpl::identity<vtkTetra>,
                    mpl::identity<vtkHexahedron>
                >::type,
                /* We should normally not reach this case */
                /* anyway we set a default vtkTetra for face type */
                mpl::identity<vtkTetra>
            >::type
        >::type
    >::type::type vtkelement_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    ExporterVTK( WorldComm const& worldComm = Environment::worldComm() );
    ExporterVTK( std::string const& __p = "default", int freq = 1, WorldComm const& worldComm = Environment::worldComm() );
    ExporterVTK( po::variables_map const& vm, std::string const& exp_prefix = "", WorldComm const& worldComm = Environment::worldComm() ) FEELPP_DEPRECATED;
    ExporterVTK( std::string const& exp_prefix, WorldComm const& worldComm = Environment::worldComm() );

    ExporterVTK( ExporterVTK const & __ex );

    ~ExporterVTK();

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

    void init();

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
    void saveMesh( mesh_ptrtype mesh, vtkSmartPointer<vtkout_type> out ) const;
    template<typename Iterator>
    void saveNodeData( typename timeset_type::step_ptrtype step, Iterator __var, Iterator en, vtkSmartPointer<vtkout_type> out ) const;
    template<typename Iterator>
    void saveElementData( typename timeset_type::step_ptrtype step, Iterator __var, Iterator en, vtkSmartPointer<vtkout_type> out ) const;

    /**
     * As we process the timesteps one by one, we need a way to record each new timestep.
     * To do so, we use a pvd file (Paraview format) that allows use to specify new timesteps
     * using xml syntax.
     */
    int writeTimePVD(std::string xmlFilename, double timestep, std::string dataFilename, int partNo = 0) const;


    /**
     * Build a multi block structure based on the data gathered 
     * on the different processes.
     */
    vtkSmartPointer<vtkMultiBlockDataSet>
        buildMultiBlockDataSet( double time, vtkSmartPointer<vtkout_type> out ) const;

    /**
     * Actual write of the dataset into a file 
     */
    void write( int stepIndex, std::string filename, vtkSmartPointer<vtkMultiBlockDataSet> out) const;

    //@}

private:

    mutable VTKCellType M_face_type;
    mutable VTKCellType M_element_type;

    /* class members for in-situ visualization */
#if VTK_MAJOR_VERSION >= 6 && defined(VTK_HAS_PARALLEL)
    mutable MPI_Comm lComm;
    mutable vtkMPICommunicatorOpaqueComm * opaqueComm;
#if defined(FEELPP_VTK_INSITU_ENABLED)
    mutable vtkSmartPointer<vtkCPProcessor> inSituProcessor;
#endif
#endif

};

} // Feel

//#if !defined( FEELPP_INSTANTIATION_MODE )
# include <feel/feelfilters/exporterVTK_impl.hpp>
//#endif // FEELPP_INSTANTIATION_MODE

#endif // defined(FEELPP_HAS_VTK)

#endif /* __ExporterVTK_H */

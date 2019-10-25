/* -* -mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
#ifndef FEELPP_FILTERS_EXPORTERVTK_HPP
#define FEELPP_FILTERS_EXPORTERVTK_HPP 1

#if defined(FEELPP_HAS_VTK)

#include <iostream>
#include <fstream>

#include <cstring>
#ifdef FEELPP_HAS_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>
#endif

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-W#warnings"
#pragma clang diagnostic ignored "-Winconsistent-missing-override"

#endif
#if defined(__GNUC__) && !(defined(__clang__))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#endif

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

#include <feel/feelfilters/vtkBaseInsituPipeline.hpp>
#endif // FEELPP_VTK_INSITU_ENABLED

#endif // VTK_MAJOR_VERSION >= 6 && defined(VTK_HAS_PARALLEL)

#if defined(__GNUC__) && !(defined(__clang__))
#pragma GCC diagnostic pop
#endif
#if defined(__clang__)
#pragma clang diagnostic pop
#endif


#include <feel/feelfilters/detail/meshcontiguousnumberingmapping.hpp>

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
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    typedef Exporter<MeshType,N> super;
    typedef typename mesh_type::value_type value_type;
    typedef typename super::timeset_type timeset_type;
    typedef typename super::timeset_ptrtype timeset_ptrtype;
    typedef typename super::timeset_const_iterator timeset_const_iterator;

    typedef typename timeset_type::step_type step_type;
    typedef typename timeset_type::step_ptrtype step_ptrtype;
    typedef typename timeset_type::step_const_iterator step_const_iterator;
protected :
    using steps_write_on_disk_type = typename super::steps_write_on_disk_type;
    using mesh_contiguous_numbering_mapping_type = Feel::detail::MeshContiguousNumberingMapping<mesh_type,float>;
    using mesh_contiguous_numbering_mapping_ptrtype = std::shared_ptr<mesh_contiguous_numbering_mapping_type>;
public :

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

    explicit ExporterVTK( worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    ExporterVTK( std::string const& __p = "default", int freq = 1, worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    ExporterVTK( po::variables_map const& vm, std::string const& exp_prefix = "", worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() ) FEELPP_DEPRECATED;
    ExporterVTK( std::string const& exp_prefix, worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );

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

    //@}

    /** @name  Methods
     */
    //@{

    /**
       save the timeset
     */
    void save( steps_write_on_disk_type const& stepsToWriteOnDisk ) const override;

    /**
     * export mesh
     */
    void visit( mesh_type* mesh ) override;

private :
    /**
     * save the \p mesh to the file \p filename
     */
    void saveMesh( timeset_ptrtype __ts, mesh_ptrtype mesh, std::map<std::string,vtkSmartPointer<vtkout_type>> & outs ) const;
    void saveFields( timeset_ptrtype __ts, typename timeset_type::step_ptrtype step, std::map<std::string,vtkSmartPointer<vtkout_type>> & outs ) const;
    template<bool IsNodal,typename Iterator>
    void saveFields( typename timeset_type::step_ptrtype step, mesh_contiguous_numbering_mapping_type const& mp, int part, Iterator __var, Iterator en, vtkSmartPointer<vtkout_type> out ) const;

#ifdef FEELPP_HAS_LIBXML2
    /**
     * As we process the timesteps one by one, we need a way to record each new timestep.
     * To do so, we use a pvd file (Paraview format) that allows use to specify new timesteps
     * using xml syntax.
     */
    int writeTimePVD(std::string xmlFilename, double timestep, std::string dataFilename, int partNo = 0) const;
#endif


    /**
     * Build a multi block structure based on the data gathered
     * on the different processes.
     */
    vtkSmartPointer<vtkMultiBlockDataSet>
    buildMultiBlockDataSet( double time, std::map<std::string,vtkSmartPointer<vtkout_type>> const& outs ) const;

    /**
     * Actual write of the dataset into a file
     */
    void write( int stepIndex, std::string filename, vtkSmartPointer<vtkMultiBlockDataSet> out) const;


    void saveData( vtkSmartPointer<vtkMultiBlockDataSet> mbds, int stepIndex, double time ) const;

    void updateInSituProcessor( vtkSmartPointer<vtkMultiBlockDataSet> mbds, int stepIndex, double time ) const;
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
    mutable std::map<std::string, mesh_contiguous_numbering_mapping_ptrtype > M_cache_mp;
};

} // Feel

//#if !defined( FEELPP_INSTANTIATION_MODE )
# include <feel/feelfilters/exportervtk_impl.hpp>
//#endif // FEELPP_INSTANTIATION_MODE

#endif // defined(FEELPP_HAS_VTK)

#endif /* __ExporterVTK_H */

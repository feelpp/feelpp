/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2014-06-04

  Copyright (C) 2014 Université Joseph Fourier (Grenoble I)

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
   \file thermodyn.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2014-06-04
 */

#include <boost/preprocessor/cat.hpp>

#undef GUARD_FOR_THERMODYNAMICS
#define GUARD_FOR_THERMODYNAMICSBASE 1
#include "thermodynconfig.h"
#undef THERMODYNAMICSBASE_CLASS_NAME
#define THERMODYNAMICSBASE_CLASS_NAME BOOST_PP_CAT(ThermoDynamicsBase,THERMODYNAMICSBASE_NAMECLASS_SPEC)

#if defined( INCLUDE_THERMODYNAMICSBASE_HPP )

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/exporter.hpp>
//#include <feel/feelvf/vf.hpp>
#include <feel/feelts/bdf.hpp>

#include <feel/feelmodels2/feelmodelscore/applibasenumericalsimulationtransitory.hpp>
#include <feel/feelmodels2/feelmodelscore/markermanagement.hpp>
#include <feel/feelmodels2/feelmodelsalg/modelalgebraic.hpp>



namespace Feel
{

namespace FeelModels
{

    class THERMODYNAMICSBASE_CLASS_NAME : public AppliBaseNumericalSimulationTransitory,
                                          public MarkerManagementDirichletBC,
                                          public MarkerManagementNeumannBC

    {
    public:
        typedef AppliBaseNumericalSimulationTransitory super_type;
        typedef THERMODYNAMICSBASE_CLASS_NAME self_type;
        typedef boost::shared_ptr<self_type> self_ptrtype;
        //___________________________________________________________________________________//
        // mesh
        static const uint16_type nDim = THERMODYNAMICS_DIM;
        static const uint16_type nOrderGeo = THERMODYNAMICS_ORDERGEO;
        typedef Simplex<nDim,nOrderGeo,nDim> convex_type;
        typedef Mesh<convex_type> mesh_type;
        typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
        // basis
        static const uint16_type nOrderPoly = THERMODYNAMICS_ORDERPOLY;
        typedef Lagrange<nOrderPoly, Scalar,Continuous,PointSetFekete> basis_temperature_type;
        typedef Lagrange<nOrderPoly, Vectorial,Continuous,PointSetFekete> basis_velocityconvection_type;
        // function space temperature
        typedef FunctionSpace<mesh_type, bases<basis_temperature_type> > space_temperature_type;
        typedef boost::shared_ptr<space_temperature_type> space_temperature_ptrtype;
        typedef typename space_temperature_type::element_type element_temperature_type;
        typedef boost::shared_ptr<element_temperature_type> element_temperature_ptrtype;
        // function space velocity convection
        typedef FunctionSpace<mesh_type, bases<basis_velocityconvection_type> > space_velocityconvection_type;
        typedef boost::shared_ptr<space_velocityconvection_type> space_velocityconvection_ptrtype;
        typedef typename space_velocityconvection_type::element_type element_velocityconvection_type;
        typedef boost::shared_ptr<element_velocityconvection_type> element_velocityconvection_ptrtype;
        // time scheme
        typedef Bdf<space_temperature_type>  bdf_temperature_type;
        typedef boost::shared_ptr<bdf_temperature_type> bdf_temperature_ptrtype;
        // exporter
        typedef Exporter<mesh_type,nOrderGeo> export_type;
        typedef boost::shared_ptr<export_type> export_ptrtype;

        // algebraic solver
        typedef MethodsNum methodsnum_type;
        typedef boost::shared_ptr< methodsnum_type > methodsnum_ptrtype;


        THERMODYNAMICSBASE_CLASS_NAME( bool __isStationary,
                                       std::string __prefix,
                                       WorldComm const& __worldComm,
                                       bool __buildMesh,
                                       std::string __subPrefix,
                                       std::string __appliShortRepository );

        std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"ThermoDynamicsMesh.path"); }
        //___________________________________________________________________________________//
        // mesh, space, element temperature
        mesh_ptrtype const& mesh() const { return M_mesh; }
        double meshSize() const { return M_meshSize; }
        space_temperature_ptrtype const& spaceTemperature() const { return M_Xh; }
        element_temperature_ptrtype const& fieldTemperature() const { return M_fieldTemperature; }
        element_velocityconvection_ptrtype const& fieldVelocityConvection() const { return M_fieldVelocityConvection; }
        bool fieldVelocityConvectionIsUsed() const { return M_fieldVelocityConvectionIsUsed; }
        bool fieldVelocityConvectionIsIncompressible() const { return M_fieldVelocityConvectionIsIncompressible; }
        void setFieldVelocityConvectionIsUsed(bool b) { M_fieldVelocityConvectionIsUsed=b; }
        bool fieldVelocityConvectionIsOperational() const { return (M_fieldVelocityConvection.use_count() > 0); }
        bool fieldVelocityConvectionIsUsedAndOperational() const { return this->fieldVelocityConvectionIsUsed() && this->fieldVelocityConvectionIsOperational(); }
        void setFieldVelocityConvectionIsIncompressible(bool b) { M_fieldVelocityConvectionIsIncompressible=b; }
        //___________________________________________________________________________________//
        // physical parameters
        double thermalConductivity() const { return M_thermalConductivity; }
        double rho() const { return M_rho; } // density
        double heatCapacity() const { return M_heatCapacity; }
        void setThermalConductivity( double d ) { M_thermalConductivity=d; }
        void setRho( double d ) { M_rho=d; }
        void setHeatCapacity( double d ) { M_heatCapacity=d; }
        //___________________________________________________________________________________//
        // algebraic data and solver
        backend_ptrtype const& backend() const { return  M_backend; }
        BlocksBaseVector<double> const& blockVectorSolution() const { return M_blockVectorSolution; }
        BlocksBaseVector<double> & blockVectorSolution() { return M_blockVectorSolution; }
        size_type nLocalDof() const;
        methodsnum_ptrtype const& methodNum() const { return M_methodNum; }
        methodsnum_ptrtype & methodNum() { return M_methodNum; }
        //___________________________________________________________________________________//
        // exporter
        void exportResults() { this->exportResults( this->currentTime() ); }
        void exportResults( double time );
        //___________________________________________________________________________________//
        // time step scheme
        bdf_temperature_ptrtype const& timeStepBdfTemperature() const { return M_bdfTemperature; }
        boost::shared_ptr<TSBase> timeStepBase() { return this->timeStepBdfTemperature(); }
        boost::shared_ptr<TSBase> timeStepBase() const { return this->timeStepBdfTemperature(); }
        void updateBdf();
        void updateTimeStep() { this->updateBdf(); }
        //___________________________________________________________________________________//

        boost::shared_ptr<std::ostringstream> getInfo() const;

        virtual void loadConfigBCFile() = 0;
        virtual void loadConfigMeshFile(std::string const& geofilename) = 0;

        void loadParameterFromOptionsVm();
        void createMesh();
        void createFunctionSpaces();
        void createTimeDiscretisation();
        void createExporters();
        BlocksBaseGraphCSR buildBlockMatrixGraph() const;
        int nBlockMatrixGraph() const { return 1; }
        void init( bool buildMethodNum, methodsnum_type::appli_ptrtype const& app );
        void updateForUseFunctionSpacesVelocityConvection();
        void restartExporters();

        void build();
        void loadMesh( mesh_ptrtype mesh );

        //___________________________________________________________________________________//
        //___________________________________________________________________________________//
        // apply assembly and solver
        void solve();

        virtual void updateLinearPDE( const vector_ptrtype& X, sparse_matrix_ptrtype& A, vector_ptrtype& F, bool buildCstPart,
                                      sparse_matrix_ptrtype& A_extended, bool _BuildExtendedPart,
                                      bool _doClose, bool _doBCStrongDirichlet ) const;
        virtual void updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const = 0;
        virtual void updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const=0;
        virtual void updateSourceTermLinearPDE(vector_ptrtype& F, bool buildCstPart) const =0;

        //___________________________________________________________________________________//
        //___________________________________________________________________________________//
        // update field from expr
        template < typename ExprT >
        void updateFieldVelocityConvection( vf::Expr<ExprT> const& expr )
        {
            if ( M_fieldVelocityConvection )
                this->updateForUseFunctionSpacesVelocityConvection();
            M_fieldVelocityConvection->on(_range=elements(this->mesh()), _expr=expr );
        }


        //private :
    protected :
        double M_meshSize;
        mesh_ptrtype M_mesh;

        space_temperature_ptrtype M_Xh;
        element_temperature_ptrtype M_fieldTemperature;
        bool M_fieldVelocityConvectionIsUsed, M_fieldVelocityConvectionIsIncompressible;
        space_velocityconvection_ptrtype M_XhVelocityConvection;
        element_velocityconvection_ptrtype M_fieldVelocityConvection; // only define with convection effect
        bdf_temperature_ptrtype M_bdfTemperature;

        // physical parameter
        double M_thermalConductivity; // [ W/(m*K) ]
        double M_rho; // density [ kg/(m^3) ]
        double M_heatCapacity; // [ J/(kg*K) ]

        // algebraic data/tools
        backend_ptrtype M_backend;
        methodsnum_ptrtype M_methodNum;
        BlocksBaseVector<double> M_blockVectorSolution;

        export_ptrtype M_exporter;
        bool M_doExportAll, M_doExportVelocityConvection;
    };

} // namespace FeelModels
} // namespace Feel

#endif /* __THERMODYNAMICSBASE_H */

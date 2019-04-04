/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2016-12-12

  Copyright (C) 2016 Feel++ Consortium

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
   \file thermoelectric.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2016-12-12
 */

#ifndef FEELPP_TOOLBOXES_ELECTRIC_HPP
#define FEELPP_TOOLBOXES_ELECTRIC_HPP 1

#include <feel/feelmodels/heat/heat.hpp>
#include <feel/feelmodels/electric/electricpropertiesdescription.hpp>


namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisPotentialType>
class Electric : public ModelNumerical,
                       public MarkerManagementDirichletBC,
                       public MarkerManagementNeumannBC,
                       public MarkerManagementRobinBC,
                       public std::enable_shared_from_this< Electric<ConvexType,BasisPotentialType> >
{

public:
    typedef ModelNumerical super_type;
    typedef Electric<ConvexType,BasisPotentialType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    // mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    typedef mesh_type mesh_electric_type;

    // function space electric-potential
    typedef BasisPotentialType basis_electricpotential_type;
    static const uint16_type nOrderPolyElectricPotential = basis_electricpotential_type::nOrder;
    typedef FunctionSpace<mesh_type, bases<basis_electricpotential_type> > space_electricpotential_type;
    typedef std::shared_ptr<space_electricpotential_type> space_electricpotential_ptrtype;
    typedef typename space_electricpotential_type::element_type element_electricpotential_type;
    typedef std::shared_ptr<element_electricpotential_type> element_electricpotential_ptrtype;
    typedef typename space_electricpotential_type::element_external_storage_type element_electricpotential_external_storage_type;
    // function space electric-field
    typedef Lagrange<nOrderPolyElectricPotential, Vectorial,Discontinuous/*Continuous*/,PointSetFekete> basis_electricfield_type;
    typedef FunctionSpace<mesh_electric_type, bases<basis_electricfield_type> > space_electricfield_type;
    typedef std::shared_ptr<space_electricfield_type> space_electricfield_ptrtype;
    typedef typename space_electricfield_type::element_type element_electricfield_type;
    typedef std::shared_ptr<element_electricfield_type> element_electricfield_ptrtype;

    // mechanical properties desc
    typedef bases<Lagrange<0, Scalar,Discontinuous> > basis_scalar_P0_type;
    typedef FunctionSpace<mesh_type, basis_scalar_P0_type> space_scalar_P0_type;
    typedef std::shared_ptr<space_scalar_P0_type> space_scalar_P0_ptrtype;
    typedef ElectricPropertiesDescription<space_scalar_P0_type> electricproperties_type;
    typedef std::shared_ptr<electricproperties_type> electricproperties_ptrtype;

    // exporter
    typedef Exporter<mesh_type,nOrderGeo> export_type;
    typedef std::shared_ptr<export_type> export_ptrtype;

    // algebraic solver
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef std::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;

    // context for evaluation
    typedef typename space_electricpotential_type::Context context_electricpotential_type;
    typedef std::shared_ptr<context_electricpotential_type> context_electricpotential_ptrtype;


    //___________________________________________________________________________________//
    // constructor
    Electric( std::string const& prefix,
              std::string const& keyword = "electric",
              worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
              std::string const& subPrefix = "",
              ModelBaseRepository const& modelRep = ModelBaseRepository() );
    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"ElectricMesh.path"); }
    std::shared_ptr<std::ostringstream> getInfo() const override;
    void updateInformationObject( pt::ptree & p ) override;

private :
    void loadParameterFromOptionsVm();
    void initMesh();
    void initBoundaryConditions();
    void initPostProcess();
public :
    void setMesh(mesh_ptrtype const& mesh) { M_mesh = mesh; }
    // update for use
    void init( bool buildModelAlgebraicFactory = true );
    BlocksBaseGraphCSR buildBlockMatrixGraph() const override;
    int nBlockMatrixGraph() const;
    void initAlgebraicFactory();

    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );
    void exportFields( double time );
    std::set<std::string> postProcessFieldExported( std::set<std::string> const& ifields, std::string const& prefix = "" ) const;
    bool updateExportedFields( export_ptrtype exporter, std::set<std::string> const& fields, double time );
    void exportMeasures( double time );
    //void setDoExportResults( bool b ) { if (M_exporter) M_exporter->setDoExport( b ); }
    bool hasPostProcessFieldExported( std::string const& key ) const { return M_postProcessFieldExported.find( key ) != M_postProcessFieldExported.end(); }

    void updateParameterValues();
    constexpr auto symbolsExpr() const { return Feel::vf::symbolsExpr( symbolExpr("electric_P",idv(this->fieldElectricPotential()) ) ); }
    //___________________________________________________________________________________//

    mesh_ptrtype const& mesh() const { return M_mesh; }
    elements_reference_wrapper_t<mesh_type> const& rangeMeshElements() const { return M_rangeMeshElements; }

    space_electricpotential_ptrtype const& spaceElectricPotential() const { return M_XhElectricPotential; }
    element_electricpotential_ptrtype const& fieldElectricPotentialPtr() const { return M_fieldElectricPotential; }
    element_electricpotential_type const& fieldElectricPotential() const { return *M_fieldElectricPotential; }

    space_electricfield_ptrtype const& spaceElectricField() const { return M_XhElectricField; }
    element_electricfield_ptrtype const& fieldElectricFieldPtr() const { return M_fieldElectricField; }
    element_electricfield_type const& fieldElectricField() const { return *M_fieldElectricField; }
    element_electricfield_ptrtype const& fieldCurrentDensityPtr() const { return M_fieldCurrentDensity; }
    element_electricfield_type const& fieldCurrentDensity() const { return *M_fieldCurrentDensity; }

    electricproperties_ptrtype const& electricProperties() const { return M_electricProperties; }

    backend_ptrtype const& backend() const { return M_backend; }
    BlocksBaseVector<double> const& blockVectorSolution() const { return M_blockVectorSolution; }
    BlocksBaseVector<double> & blockVectorSolution() { return M_blockVectorSolution; }
    model_algebraic_factory_ptrtype const& algebraicFactory() const { return M_algebraicFactory; }
    model_algebraic_factory_ptrtype & algebraicFactory() { return M_algebraicFactory; }

    //___________________________________________________________________________________//
    // apply assembly and solver
    void solve();

    void updateLinearPDE( DataUpdateLinear & data ) const override;
    void updateLinearPDEDofElimination( DataUpdateLinear & data ) const override;

    void updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const override;
    void updateJacobian( DataUpdateJacobian & data ) const override;
    void updateJacobianDofElimination( DataUpdateJacobian & data ) const override;
    void updateResidual( DataUpdateResidual & data ) const override;
    void updateResidualDofElimination( DataUpdateResidual & data ) const override;


    //___________________________________________________________________________________//
    void updateElectricField();
    void updateCurrentDensity();

    template<typename ExprT>
    void updateCurrentDensity( Expr<ExprT> const& expr, elements_reference_wrapper_t<mesh_type> range )
        {
            M_fieldCurrentDensity->on(_range=range, _expr=expr );
        }
private :
    void updateLinearPDEWeakBC( sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart ) const;
    void updateJacobianWeakBC( element_electricpotential_external_storage_type const& v, sparse_matrix_ptrtype& J, bool buildCstPart ) const;
    void updateResidualWeakBC( element_electricpotential_external_storage_type const& v, vector_ptrtype& R, bool buildCstPart ) const;

private :
    bool M_hasBuildFromMesh, M_isUpdatedForUse;

    mesh_ptrtype M_mesh;
    elements_reference_wrapper_t<mesh_type> M_rangeMeshElements;

    space_electricpotential_ptrtype M_XhElectricPotential;
    element_electricpotential_ptrtype M_fieldElectricPotential;
    space_electricfield_ptrtype M_XhElectricField;
    element_electricfield_ptrtype M_fieldElectricField;
    element_electricfield_ptrtype M_fieldCurrentDensity;

    // physical parameter
    electricproperties_ptrtype M_electricProperties;
    // boundary conditions
    map_scalar_field<2> M_bcDirichlet;
    map_scalar_field<2> M_bcNeumann;
    map_scalar_fields<2> M_bcRobin;
    map_scalar_field<2> M_volumicForcesProperties;

    // algebraic data/tools
    backend_ptrtype M_backend;
    model_algebraic_factory_ptrtype M_algebraicFactory;
    BlocksBaseVector<double> M_blockVectorSolution;
    std::map<std::string,std::set<size_type> > M_dofsWithValueImposed;

    // post-process
    export_ptrtype M_exporter;
    std::set<std::string> M_postProcessFieldExported;
    std::set<std::string> M_postProcessUserFieldExported;


};

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_TOOLBOXES_ELECTRIC_HPP

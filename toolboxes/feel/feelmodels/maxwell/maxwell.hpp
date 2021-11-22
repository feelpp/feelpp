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
   \file maxwell.hpp
   \author Romain Hild <romain.hild@unistra.fr>
   \date 2018-05-03
 */

#ifndef FEELPP_TOOLBOXES_MAXWELL_HPP
#define FEELPP_TOOLBOXES_MAXWELL_HPP 1

#include <boost/mpl/equal.hpp>

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelcore/markermanagement.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feelpoly/raviartthomas.hpp>

#include <feel/feelmodels/maxwell/maxwellpropertiesdescription.hpp>

namespace Feel
{
namespace FeelModels
{

enum class MaxwellPostProcessFieldExported
{
    MagneticPotential = 0, MagneticField, Pid
};

template< typename ConvexType>//, typename BasisPotentialType>
class Maxwell : public ModelNumerical,
                public MarkerManagementDirichletBC,
                public MarkerManagementNeumannBC,
                public MarkerManagementRobinBC,
                public std::enable_shared_from_this< Maxwell<ConvexType>>//,BasisPotentialType> >
{

public:
    typedef ModelNumerical super_type;
    typedef Maxwell<ConvexType>/*,BasisPotentialType>*/ self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    // mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    // function space magnetic-potential
    // typedef BasisPotentialType basis_magneticpotential_type;
    using basis_magneticpotential_type = Nedelec<0, NedelecKind::NED1>;
    static const uint16_type nOrderPolyMagneticPotential = basis_magneticpotential_type::nOrder;
    typedef FunctionSpace<mesh_type, bases<basis_magneticpotential_type> > space_magneticpotential_type;
    typedef std::shared_ptr<space_magneticpotential_type> space_magneticpotential_ptrtype;
    typedef typename space_magneticpotential_type::element_type element_magneticpotential_type;
    typedef std::shared_ptr<element_magneticpotential_type> element_magneticpotential_ptrtype;
    typedef typename space_magneticpotential_type::element_external_storage_type element_magneticpotential_external_storage_type;

    // function space magnetic-field
    // typedef Lagrange<nOrderPolyMagneticPotential-1, Vectorial,Discontinuous/*Continuous*/,PointSetFekete> basis_magneticfield_type;
// #if FEELPP_DIM==3
    using basis_magneticfield_type = typename mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<3> >,
                                                       RaviartThomas<0>,
                                                       Lagrange<0, Scalar, Discontinuous> >::type;
    // using basis_magneticfield_type = RaviartThomas<0>;
// #else
//     using basis_magneticfield_type = Lagrange<0, Scalar, Discontinuous>;
// #endif
    typedef FunctionSpace<mesh_type, bases<basis_magneticfield_type> > space_magneticfield_type;
    typedef std::shared_ptr<space_magneticfield_type> space_magneticfield_ptrtype;
    typedef typename space_magneticfield_type::element_type element_magneticfield_type;
    typedef std::shared_ptr<element_magneticfield_type> element_magneticfield_ptrtype;

    // mechanical properties desc
    typedef bases<Lagrange<0, Scalar,Discontinuous> > basis_scalar_P0_type;
    typedef FunctionSpace<mesh_type, basis_scalar_P0_type> space_scalar_P0_type;
    typedef std::shared_ptr<space_scalar_P0_type> space_scalar_P0_ptrtype;
    typedef MaxwellPropertiesDescription<space_scalar_P0_type> maxwellproperties_type;
    typedef std::shared_ptr<maxwellproperties_type> maxwellproperties_ptrtype;

    // exporter
    typedef Exporter<mesh_type,nOrderGeo> export_type;
    typedef std::shared_ptr<export_type> export_ptrtype;

    // algebraic solver
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef std::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;

    // context for evaluation
    typedef typename space_magneticpotential_type::Context context_magneticpotential_type;
    typedef std::shared_ptr<context_magneticpotential_type> context_magneticpotential_ptrtype;

    using map_dirichlet_field = typename mpl::if_< mpl::equal_to<mpl::int_<nDim>, mpl::int_<3> >,
                                                   map_vector_field<nDim>,
                                                   map_scalar_field<2> >::type;


    //___________________________________________________________________________________//
    // constructor
    Maxwell( std::string const& prefix,
             bool buildMesh = true,
             worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
             std::string const& subPrefix = "",
             ModelBaseRepository const& modelRep = ModelBaseRepository() );

    std::shared_ptr<std::ostringstream> getInfo() const;
private :
    void loadParameterFromOptionsVm();
    void initMesh();
    void initBoundaryConditions();
    template<typename convex = convex_type>
    void initDirichlet(std::enable_if_t<convex::nDim==2>* = nullptr) { this->M_bcDirichlet = this->modelProperties().boundaryConditions().template getScalarFields<nDim>( "magnetic-potential", "Dirichlet" ); }
    template<typename convex = convex_type>
    void initDirichlet(std::enable_if_t<convex::nDim==3>* = nullptr) { this->M_bcDirichlet = this->modelProperties().boundaryConditions().template getVectorFields<nDim>( "magnetic-potential", "Dirichlet" ); }
    void initPostProcess();
public :
    // update for use
    void init( bool buildModelAlgebraicFactory = true );
    BlocksBaseGraphCSR buildBlockMatrixGraph() const;
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

    //___________________________________________________________________________________//

    mesh_ptrtype mesh() const { return super_type::super_model_meshes_type::mesh<mesh_type>( this->keyword() ); }
    void setMesh( mesh_ptrtype const& mesh ) { super_type::super_model_meshes_type::setMesh( this->keyword(), mesh ); }
    elements_reference_wrapper_t<mesh_type> const& rangeMeshElements() const { return M_rangeMeshElements; }

    space_magneticpotential_ptrtype const& spaceMagneticPotential() const { return M_XhMagneticPotential; }
    element_magneticpotential_ptrtype const& fieldMagneticPotentialPtr() const { return M_fieldMagneticPotential; }
    element_magneticpotential_type const& fieldMagneticPotential() const { return *M_fieldMagneticPotential; }

    space_magneticfield_ptrtype const& spaceMagneticField() const { return M_XhMagneticField; }
    element_magneticfield_ptrtype const& fieldMagneticFieldPtr() const { return M_fieldMagneticField; }
    element_magneticfield_type const& fieldMagneticField() const { return *M_fieldMagneticField; }

    maxwellproperties_ptrtype const& maxwellProperties() const { return M_maxwellProperties; }

    //___________________________________________________________________________________//
    // apply assembly and solver
    void solve();

    void updateLinearPDE( DataUpdateLinear & data ) const;
    void updateLinearPDEWeakBC( sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart ) const;
    void updateLinearPDEStrongDirichletBC( sparse_matrix_ptrtype& A, vector_ptrtype& F ) const;

    //___________________________________________________________________________________//
    void updateMagneticField();

    //___________________________________________________________________________________//
    // template<typename T, typename convex = convex_type>
    // auto vcurl(T f, std::enable_if_t<convex::nDim==3>* = nullptr) const -> decltype(curl(f)) { return curl(f); }
    // template<typename T, typename convex = convex_type>
    // auto vcurlt(T f, std::enable_if_t<convex::nDim==3>* = nullptr) const -> decltype(curlt(f)) { return curlt(f); }
    // template<typename T, typename convex = convex_type>
    // auto vcurlv(T f, std::enable_if_t<convex::nDim==3>* = nullptr) const -> decltype(curlv(f)) { return curlv(f); }
    // template<typename T, typename convex = convex_type>
    // auto vcurl(T f, std::enable_if_t<convex::nDim==2>* = nullptr) const -> decltype(curlx(f)) { return curlx(f); }
    // template<typename T, typename convex = convex_type>
    // auto vcurlt(T f, std::enable_if_t<convex::nDim==2>* = nullptr) const -> decltype(curlxt(f)) { return curlxt(f); }
    // template<typename T, typename convex = convex_type>
    // auto vcurlv(T f, std::enable_if_t<convex::nDim==2>* = nullptr) const -> decltype(curlxv(f)) { return curlxv(f); }

private :
    //bool M_hasBuildFromMesh, M_isUpdatedForUse;

    elements_reference_wrapper_t<mesh_type> M_rangeMeshElements;

    space_magneticpotential_ptrtype M_XhMagneticPotential;
    element_magneticpotential_ptrtype M_fieldMagneticPotential;
    space_magneticfield_ptrtype M_XhMagneticField;
    element_magneticfield_ptrtype M_fieldMagneticField;
    // physical parameter
    maxwellproperties_ptrtype M_maxwellProperties;
    // boundary conditions
    map_dirichlet_field M_bcDirichlet;
    map_scalar_field<2> M_bcNeumann;
    map_scalar_fields<2> M_bcRobin;
    map_vector_field<nDim> M_volumicForcesProperties;

    // regularization
    double M_epsilon;

    std::map<std::string,std::set<size_type> > M_dofsWithValueImposed;
    // start dof index fields in matrix (temperature,maxwell-potential,...)
    std::map<std::string,size_type> M_startBlockIndexFieldsInMatrix;

    // post-process
    export_ptrtype M_exporter;
    std::set<std::string> M_postProcessFieldExported;
    std::set<std::string> M_postProcessUserFieldExported;

}; // class Maxwell

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_TOOLBOXES_MAXWELL_HPP

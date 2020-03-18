/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

#ifndef FEELPP_TOOLBOXES_COEFFICIENTFORMPDE_HPP
#define FEELPP_TOOLBOXES_COEFFICIENTFORMPDE_HPP 1


#include <feel/feelmodels/coefficientformpdes/coefficientformpdebase.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/exporter.hpp>
//#include <feel/feelvf/vf.hpp>
#include <feel/feelts/bdf.hpp>


#include <feel/feelmodels/modelcore/stabilizationglsparameterbase.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisUnknownType>
class CoefficientFormPDE : public CoefficientFormPDEBase<ConvexType>,
                           public std::enable_shared_from_this< CoefficientFormPDE<ConvexType,BasisUnknownType> >
{
public:
    typedef CoefficientFormPDEBase<ConvexType> super_type;
    using size_type = typename super_type::size_type;
    typedef CoefficientFormPDE<ConvexType,BasisUnknownType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;
    //___________________________________________________________________________________//
    // mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    // basis
    static const uint16_type nOrderUnknown = BasisUnknownType::nOrder;
    typedef BasisUnknownType basis_unknown_type;
    // function space unknown
    typedef FunctionSpace<mesh_type, bases<basis_unknown_type> > space_unknown_type;
    typedef std::shared_ptr<space_unknown_type> space_unknown_ptrtype;
    typedef typename space_unknown_type::element_type element_unknown_type;
    typedef std::shared_ptr<element_unknown_type> element_unknown_ptrtype;
    typedef typename space_unknown_type::element_external_storage_type element_unknown_external_storage_type;
    // materials properties
    typedef MaterialsProperties<mesh_type> materialsproperties_type;
    typedef std::shared_ptr<materialsproperties_type> materialsproperties_ptrtype;
    // time scheme
    typedef Bdf<space_unknown_type> bdf_unknown_type;
    typedef std::shared_ptr<bdf_unknown_type> bdf_unknown_ptrtype;


    // algebraic solver
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef std::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;

    CoefficientFormPDE( std::string const& prefix,
                        std::string const& keyword = "pde",
                        worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                        std::string const& subPrefix  = "",
                        ModelBaseRepository const& modelRep = ModelBaseRepository() );

    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"CoefficientFormPDEMesh.path"); }

    //___________________________________________________________________________________//
    // mesh, space, element unknown
    space_unknown_ptrtype const& spaceUnknown() const { return M_Xh; }
    element_unknown_ptrtype const& fieldUnknownPtr() const { return M_fieldUnknown; }
    element_unknown_type const& fieldUnknown() const { return *M_fieldUnknown; }
    //___________________________________________________________________________________//
#if 0
    // boundary condition + body forces
    map_scalar_field<2> const& bcDirichlet() const { return M_bcDirichlet; }
    map_scalar_field<2> const& bcNeumann() const { return M_bcNeumann; }
    map_scalar_fields<2> const& bcRobin() const { return M_bcRobin; }
    map_scalar_field<2> const& bodyForces() const { return M_volumicForcesProperties; }
#endif

    //___________________________________________________________________________________//
    // time step scheme
    std::string const& timeStepping() const { return M_timeStepping; }
    bdf_unknown_ptrtype const& timeStepBdfUnknown() const { return M_bdfUnknown; }
    std::shared_ptr<TSBase> timeStepBase() { return this->timeStepBdfUnknown(); }
    std::shared_ptr<TSBase> timeStepBase() const { return this->timeStepBdfUnknown(); }
    void startTimeStep();
    void updateTimeStep();
    //___________________________________________________________________________________//

    std::shared_ptr<std::ostringstream> getInfo() const override;
    void updateInformationObject( pt::ptree & p ) override;


    void init( bool buildModelAlgebraicFactory=true );
    void initAlgebraicFactory();

    //___________________________________________________________________________________//
    // algebraic data and solver
    backend_ptrtype const& backend() const { return  M_backend; }
    BlocksBaseVector<double> const& blockVectorSolution() const { return M_blockVectorSolution; }
    BlocksBaseVector<double> & blockVectorSolution() { return M_blockVectorSolution; }
    size_type nLocalDof() const;
    model_algebraic_factory_ptrtype const& algebraicFactory() const { return M_algebraicFactory; }
    model_algebraic_factory_ptrtype & algebraicFactory() { return M_algebraicFactory; }

    BlocksBaseGraphCSR buildBlockMatrixGraph() const override;
    int nBlockMatrixGraph() const { return 1; }

    void updateParameterValues();
    void setParameterValues( std::map<std::string,double> const& paramValues );

    //___________________________________________________________________________________//
    // apply assembly and solver
    //___________________________________________________________________________________//

    void solve();
private :
    void initFunctionSpaces();
    void initBoundaryConditions();
    void initTimeStep();
    void initInitialConditions();
    void initPostProcess() override;

private :


    space_unknown_ptrtype M_Xh;
    element_unknown_ptrtype M_fieldUnknown;

    // time discretisation
    std::string M_timeStepping;
    bdf_unknown_ptrtype M_bdfUnknown;
    //double M_timeStepThetaValue;
    //vector_ptrtype M_timeStepThetaSchemePreviousContrib;

    // boundary conditions
    map_scalar_field<2> M_bcDirichlet;
    map_scalar_field<2> M_bcNeumann;
    map_scalar_fields<2> M_bcRobin;
    map_scalar_field<2> M_volumicForcesProperties;
    MarkerManagementDirichletBC M_bcDirichletMarkerManagement;
    MarkerManagementNeumannBC M_bcNeumannMarkerManagement;
    MarkerManagementRobinBC M_bcRobinMarkerManagement;

    // algebraic data/tools
    backend_ptrtype M_backend;
    model_algebraic_factory_ptrtype M_algebraicFactory;
    BlocksBaseVector<double> M_blockVectorSolution;
};

} // namespace Feel
} // namespace FeelModels

#endif

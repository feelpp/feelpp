#ifndef _ADVECTIONBASE_HPP
#define _ADVECTIONBASE_HPP 1

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelcore/markermanagement.hpp>
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmodels/modelalg/modelalgebraicfactory.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

/*
VD, TODO :
the stabilization method and the stabilization coefficient might be external options. 
This brings the need of specific options for the advection class. 
In the level set framework for exemple, it is used twice : for the level set and the hamilton jacobi reinit. 
Thus, being able to tune externally the stab method or the stab coeff could be great */

namespace Feel {
namespace FeelModels {

enum class AdvectionStabMethod { NONE=0, GALS, CIP, SUPG, SGS };

template< typename ConvexType, typename BasisAdvectionType >
class AdvectionBase : 
    public ModelNumerical,
    public MarkerManagementDirichletBC,
    public MarkerManagementNeumannBC

{
public :
    typedef ModelNumerical super_type;

    typedef AdvectionBase< ConvexType, BasisAdvectionType > self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    //--------------------------------------------------------------------//
    // Mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    static const uint16_type nRealDim = convex_type::nRealDim;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //--------------------------------------------------------------------//
    // Space advection
    typedef BasisAdvectionType basis_advection_type;
    static const uint16_type nOrder = basis_advection_type::nOrder;
    
    typedef FunctionSpace< mesh_type, bases<basis_advection_type> > space_advection_type;
    typedef boost::shared_ptr<space_advection_type> space_advection_ptrtype;
    
    typedef typename space_advection_type::element_type element_advection_type;
    typedef boost::shared_ptr<element_advection_type> element_advection_ptrtype;

    typedef typename space_advection_type::value_type value_type;
    typedef typename space_advection_type::periodicity_type periodicity_advection_type;

    //--------------------------------------------------------------------//
    // Space advection velocity
    //typedef typename basis_advection_type::vectorial_basis_type basis_vectorial_type;
    //typedef Lagrange<nOrder, Vectorial, basis_advection_continuity_type, basis_advection_pointset_type> basis_vectorial_type;
    typedef Lagrange<nOrder, Vectorial, Continuous, PointSetFekete> basis_vectorial_type;

    typedef bases<basis_vectorial_type> basis_advection_velocity_type;
    typedef FunctionSpace< 
        mesh_type, 
        basis_advection_velocity_type, 
        value_type, 
        periodicity_advection_type > space_advection_velocity_type;
    typedef boost::shared_ptr<space_advection_velocity_type> space_advection_velocity_ptrtype;
    typedef typename space_advection_velocity_type::element_type element_advection_velocity_type;
    typedef boost::shared_ptr<element_advection_velocity_type> element_advection_velocity_ptrtype;

    //--------------------------------------------------------------------//
    // Backend
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    //--------------------------------------------------------------------//
    // Algebraic tools
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef boost::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;
    typedef typename model_algebraic_factory_type::graph_type graph_type;
    typedef typename model_algebraic_factory_type::graph_ptrtype graph_ptrtype;
    typedef typename model_algebraic_factory_type::indexsplit_type indexsplit_type;
    typedef typename model_algebraic_factory_type::indexsplit_ptrtype indexsplit_ptrtype;

    //--------------------------------------------------------------------//
    // Time
    typedef Bdf<space_advection_type> bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;

    //--------------------------------------------------------------------//
    // Exporter
    typedef Exporter<mesh_type, nOrderGeo> exporter_type;
    typedef boost::shared_ptr<exporter_type> exporter_ptrtype;

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//

    //--------------------------------------------------------------------//
    // Constructor
    AdvectionBase( 
            std::string const& prefix,
            WorldComm const& _worldComm = Environment::worldComm(),
            std::string const& subPrefix = "",
            std::string const& rootRepository = ModelBase::rootRepositoryByDefault() );

    AdvectionBase( self_type const& A ) = default;

    //--------------------------------------------------------------------//
    // Initialization
    void build();
    void init( bool buildModelAlgebraicFactory, model_algebraic_factory_type::appli_ptrtype const& app );

    void loadParametersFromOptionsVm();
    void createMesh();
    void createFunctionSpaces();
    void createAlgebraicData();
    void createTimeDiscretization();
    void createExporters();
    
    //--------------------------------------------------------------------//
    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"AdvectionMesh.path"); }

    mesh_ptrtype const& mesh() const { return M_mesh; }

    space_advection_ptrtype const& functionSpace() const { return M_Xh; }

    element_advection_ptrtype & fieldSolutionPtr() { return M_fieldSolution; }
    element_advection_ptrtype const& fieldSolutionPtr() const { return M_fieldSolution; }
    element_advection_type & fieldSolution() { return *M_fieldSolution; }
    element_advection_type const& fieldSolution() const { return *M_fieldSolution; }

    element_advection_velocity_type const& fieldAdvectionVelocity() const { return *M_fieldAdvectionVelocity; }

    //--------------------------------------------------------------------//
    // Algebraic data
    backend_ptrtype backend() { return M_backend; }
    backend_ptrtype const& backend() const { return  M_backend; }
    size_type matrixPattern() const;
    graph_ptrtype buildMatrixGraph() const;
    //indexsplit_ptrtype buildIndexSplit() const;
    model_algebraic_factory_ptrtype algebraicFactory() { return M_algebraicFactory; }
    model_algebraic_factory_ptrtype const& algebraicFactory() const { return M_algebraicFactory; }
    size_type nLocalDof() const;

    //--------------------------------------------------------------------//
    // Time scheme
    bdf_ptrtype timeStepBDF() { return M_bdf; }
    bdf_ptrtype const& timeStepBDF() const { return M_bdf; }
    boost::shared_ptr<TSBase> timeStepBase() { return this->timeStepBDF(); }
    boost::shared_ptr<TSBase> timeStepBase() const { return this->timeStepBDF(); }
    void updateTimeStepBDF();
    void initTimeStep();
    void updateTimeStep() { this->updateTimeStepBDF(); }

    //--------------------------------------------------------------------//
    // Stabilization
    AdvectionStabMethod stabilizationMethod() const { return M_stabMethod; }
    double stabilizationCoefficient() const { return M_stabCoeff; }
    void setStabilizationCoefficient(double coeff) { M_stabCoeff = coeff; }

    //--------------------------------------------------------------------//
    // Algebraic model updates
    // Linear PDE
    void updateLinearPDE( DataUpdateLinear & data ) const;
    void updateLinearPDEStabilization( DataUpdateLinear & data ) const;
    virtual void updateSourceTermLinearPDE(vector_ptrtype& F, bool buildCstPart) const =0;
    virtual void updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const = 0;
    virtual void updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const=0;
    
    void updateBCNeumannLinearPDE( vector_ptrtype& F ) const;
    
    //--------------------------------------------------------------------//
    // Advection velocity update
    template<typename ExprT>
    void updateAdvectionVelocity(vf::Expr<ExprT> const& expr);

    //--------------------------------------------------------------------//
    // Solve
    virtual void solve();

    // Export results
    void exportMeasures( double time );
    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );

protected:

    bool M_isUpdatedForUse;
    //--------------------------------------------------------------------//
    // Mesh
    mesh_ptrtype M_mesh;
    // Advection space
    space_advection_ptrtype M_Xh;
    // Time discretization
    bdf_ptrtype M_bdf;
    //--------------------------------------------------------------------//
    // Algebraic data
    backend_ptrtype M_backend;
    model_algebraic_factory_ptrtype M_algebraicFactory;
    BlocksBaseVector<double> M_blockVectorSolution;
    //--------------------------------------------------------------------//
    // Advection velocity
    space_advection_velocity_ptrtype M_XhAdvectionVelocity;
    element_advection_velocity_ptrtype M_fieldAdvectionVelocity;
    boost::optional<vector_field_expression<nDim,1,2> > M_exprAdvectionVelocity;
    //--------------------------------------------------------------------//
    // Solution
    element_advection_ptrtype M_fieldSolution;
    //--------------------------------------------------------------------//
    // Export
    exporter_ptrtype M_exporter;
    bool M_doExportAll, M_doExportAdvectionVelocity;
    //--------------------------------------------------------------------//
    // Stabilization
    static const std::map<std::string, AdvectionStabMethod> AdvectionStabMethodIdMap;
    AdvectionStabMethod M_stabMethod;
    double M_stabCoeff;
    double M_gamma1;
    
};//AdvectionBase

} // namespace FeelModels
} // namespace Feel


#endif

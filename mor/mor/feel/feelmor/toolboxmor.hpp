#ifndef FEELPP_TOOLBOX_MOR_HPP
#define FEELPP_TOOLBOX_MOR_HPP

#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelmor/modelcrbbase.hpp>
#include <feel/feelmor/crbmodelproperties.hpp>

namespace Feel
{

namespace po = boost::program_options;

FEELPP_EXPORT po::options_description
makeToolboxMorOptions();
FEELPP_EXPORT AboutData
makeToolboxMorAbout( std::string const& str = "opusheat-tb" );

/**
 * @brief Model for ToolboxMor
 * 
 * Should give the assembly process for the rhs and lhs, both in offline and online phase
 * 
 * @tparam MeshType The type of the mesh used by ToolboxMor
 */
template<typename MeshType>
class DeimMorModelBase
{
public:
    using mesh_type = MeshType;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;
    using parameter_type = typename ParameterSpace<>::element_type;
    using backend_type = Backend<double>;
    using vector_ptrtype = typename backend_type::vector_ptrtype;
    using sparse_matrix_ptrtype = typename backend_type::sparse_matrix_ptrtype;
    using deim_function_type = std::function<vector_ptrtype(parameter_type const&)>;
    using mdeim_function_type = std::function<sparse_matrix_ptrtype(parameter_type const&)>;

    /**
     * @brief Should return the function assembling the rhs
     * 
     * @return deim_function_type The function assembling the rhs
     */
    virtual deim_function_type deimFunction() = 0;
    /**
     * @brief Should return the function assembling the lhs
     * 
     * @return mdeim_function_type The function assembling the lhs
     */
    virtual mdeim_function_type mdeimFunction() = 0;
    /**
     * @brief Should return the function assembling the rhs for the online phase
     * 
     * @param mesh the mesh used for the online phase
     * @return deim_function_type the function assembling the rhs
     */
    virtual deim_function_type deimOnlineFunction(mesh_ptrtype const& mesh) = 0;
    /**
     * @brief Should return the function assembling the lhs for the online phase
     * 
     * @param mesh the mesh used for the online phase
     * @return mdeim_function_type the function assembling the lhs
     */
    virtual mdeim_function_type mdeimOnlineFunction(mesh_ptrtype const& mesh) = 0;
};

/**
 * @brief A model for ToolboxMor when the assembly is done by a Toolbox
 * 
 * @tparam ToolboxType The type of the toolbox
 */
template<typename ToolboxType>
class DeimMorModelToolbox : public DeimMorModelBase<typename ToolboxType::mesh_type>
{
    using toolbox_type = ToolboxType;
    using toolbox_ptrtype = std::shared_ptr<toolbox_type>;
    using mesh_type = typename toolbox_type::mesh_type;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;
    using self_type = DeimMorModelToolbox<toolbox_type>;
    using self_ptrtype = std::shared_ptr<self_type>;
    using super_type = DeimMorModelBase<mesh_type>;
    using parameter_type = typename super_type::parameter_type;
    using vector_ptrtype = typename super_type::vector_ptrtype;
    using sparse_matrix_ptrtype = typename super_type::sparse_matrix_ptrtype;
    using deim_function_type = typename super_type::deim_function_type;
    using mdeim_function_type = typename super_type::mdeim_function_type;

    using toolbox_init_function_type = std::function<toolbox_ptrtype(mesh_ptrtype)>;

public:
    /**
    * @brief Construct a new Deim Mor Model Toolbox object for the online phase
    * 
    * @param prefix the prefix used for the toolbox
    */
    DeimMorModelToolbox(std::string const& prefix) :
        M_prefix(prefix)
        {}
    /**
     * @brief Construct a new Deim Mor Model Toolbox object
     * 
     * @param tb The toolbox used for the offline assembly
     */
    DeimMorModelToolbox(toolbox_ptrtype tb) :
        M_prefix(tb->prefix()),
        M_tb(tb)
    {
        M_rhs = M_tb->algebraicFactory()->rhs()->clone();
        M_mat = M_tb->algebraicFactory()->matrix();
    }
    static self_ptrtype New(std::string const& prefix) {return std::make_shared<self_type>(prefix);}
    static self_ptrtype New(toolbox_ptrtype tb) {return std::make_shared<self_type>(tb);}

    deim_function_type deimFunction() override;
    mdeim_function_type mdeimFunction() override;
    deim_function_type deimOnlineFunction(mesh_ptrtype const& mesh) override;
    mdeim_function_type mdeimOnlineFunction(mesh_ptrtype const& mesh) override;

    void setToolboxInitFunction( toolbox_init_function_type fun ) { M_toolboxInitFunction = fun; }

private:
    std::string M_prefix;
    toolbox_ptrtype M_tb;
    vector_ptrtype M_rhs;
    sparse_matrix_ptrtype M_mat;
    toolbox_ptrtype M_tbDeim;
    vector_ptrtype M_rhsDeim;
    sparse_matrix_ptrtype M_matDeim;
    toolbox_ptrtype M_tbMdeim;
    vector_ptrtype M_rhsMdeim;
    sparse_matrix_ptrtype M_matMdeim;

    toolbox_init_function_type M_toolboxInitFunction;
};

/**
 * @brief A class handling the decomposition of a parametrized problem
 * 
 * Given a function to assemble the rhs and one to assemble the lhs, the class uses DEIM to get an affine decomposition
 * 
 * @tparam SpaceType The function space of the solution
 * @tparam Options The options of the ModelCrbBase class
 */
template<typename SpaceType, int Options = 0>
class FEELPP_EXPORT ToolboxMor : public ModelCrbBase< ParameterSpace<>, SpaceType, Options>
{
    using self_type = ToolboxMor<SpaceType, Options>;
    using super_type = ModelCrbBase< ParameterSpace<>, SpaceType, Options>;

  public:

    typedef typename super_type::beta_vector_light_type beta_vector_light_type;
    typedef typename super_type::affine_decomposition_light_type affine_decomposition_light_type;

    typedef typename super_type::element_type element_type;
    typedef typename super_type::element_ptrtype element_ptrtype;
    typedef typename super_type::space_type space_type;
    using mesh_type = typename super_type::mesh_type;
    using mesh_ptrtype = typename super_type::mesh_ptrtype;
    static constexpr uint16_type nDim = space_type::nDim;
    static constexpr bool is_scalar = space_type::is_scalar;

    using super_type::computeBetaQm;
    using parameterspace_type = typename super_type::parameterspace_type;
    using parameter_type = typename super_type::parameter_type;
    using vectorN_type = typename super_type::vectorN_type;
    using beta_vector_type = typename super_type::beta_vector_type;
    using affine_decomposition_type = typename super_type::affine_decomposition_type;
    using sampling_type = typename ParameterSpace<>::sampling_type;
    using sampling_ptrtype = std::shared_ptr<sampling_type>;

    using vector_ptrtype = typename super_type::vector_ptrtype;
    using sparse_matrix_ptrtype = typename super_type::sparse_matrix_ptrtype;

    using deim_type = DEIM<self_type>;
    using deim_ptrtype = std::shared_ptr<deim_type>;
    using mdeim_type = MDEIM<self_type>;
    using mdeim_ptrtype = std::shared_ptr<mdeim_type>;
    using deim_function_type = std::function<vector_ptrtype(parameter_type const&)>;
    using mdeim_function_type = std::function<sparse_matrix_ptrtype(parameter_type const&)>;

    /**
     * @brief Construct a new Toolbox Mor object
     * 
     * @param name Name of the database
     * @param prefix Prefix for options
     */
    explicit ToolboxMor(std::string const& name = "ToolboxMor", std::string const& prefix = "");

    /**
     * @brief Set the Assemble DEIM function
     * 
     * @param fct Function to assemble the rhs during the offline phase
     */
    void setAssembleDEIM(deim_function_type const& fct ) { M_assembleForDEIM = fct; }
    /**
     * @brief Set the Assemble MDEIM function
     * 
     * @param fct Function to assemble the lhs during the offline phase
     */
    void setAssembleMDEIM(mdeim_function_type const& fct ) { M_assembleForMDEIM = fct; }
    /**
     * @brief Returns the reduced mesh used by DEIM
     * 
     * @return mesh_ptrtype The reduced mesh used by DEIM
     */
    mesh_ptrtype getDEIMReducedMesh() { return M_deim->onlineModel()->functionSpace()->mesh(); }
    /**
     * @brief Returns the reduced mesh used by MDEIM
     * 
     * @return mesh_ptrtype The reduced mesh used by MDEIM
     */
    mesh_ptrtype getMDEIMReducedMesh() { return M_mdeim->onlineModel()->functionSpace()->mesh(); }
    /**
     * @brief Set the Online Assemble DEIM function
     * 
     * @param fct Function to assemble the rhs during the online phase
     */
    void setOnlineAssembleDEIM(deim_function_type const& fct ) { M_deim->onlineModel()->setAssembleDEIM(fct); }
    /**
     * @brief Set the Online Assemble MDEIM function
     * 
     * @param fct Function to assemble the lhs during the online phase
     */
    void setOnlineAssembleMDEIM(mdeim_function_type const& fct ) { M_mdeim->onlineModel()->setAssembleMDEIM(fct); }
    /**
     * @brief Initialized the assembly process using a DeimMorModel
     * 
     * @param model The DeimMorModel
     */
    void initToolbox(std::shared_ptr<DeimMorModelBase<mesh_type>> model );
    /**
     * @brief Initialized the online assembly process using a DeimMorModel
     * 
     * @param model The DeimMorModel
     */
    void initOnlineToolbox(std::shared_ptr<DeimMorModelBase<mesh_type>> model);

    /**
     * @brief Initialized the model for the offline phase
     * 
     * Construct the parameter space from a ModelParameters and run the DEIM and MDEIM algorithms
     */
    void initModel() override;
    /**
     * @brief Initialized the model for the online phase
     * 
     * @param ptree json section to load the database
     * @param dbDir database directory
     */
    void setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir ) override;
    /**
     * @brief Initialize the deim online model
     * 
     */
    void initOnlineModel( std::shared_ptr<super_type> const& model ) override;
    /**
     * @brief Need to be called after the initModel to assemble the data
     * 
     */
    void postInitModel();

    sparse_matrix_ptrtype assembleForMDEIM( parameter_type const& mu, int const& tag ) override;
    vector_ptrtype assembleForDEIM( parameter_type const& mu, int const& tag ) override;
    vector_ptrtype assembleOutputMean( parameter_type const& mu, CRBModelOutput& output);
    vector_ptrtype assembleOutputIntegrate( parameter_type const& mu, CRBModelOutput& output);
    vector_ptrtype assembleOutputSensor( parameter_type const& mu, CRBModelOutput& output);
    vector_ptrtype assembleOutputPoint( parameter_type const& mu, CRBModelOutput& output);

    typename super_type::betaqm_type
    computeBetaQm( parameter_type const& mu , double time , bool only_terms_time_dependent=false ) override;
    typename super_type::betaqm_type
    computeBetaQm( parameter_type const& mu ) override;

    double output( int output_index, parameter_type const& mu, element_type &T, bool need_to_solve=false ) override;

    // element_type solve(parameter_type const& mu) override;


    /**
     * @brief Returns the ModelProperties used by ToolboxMor
     * 
     * @return std::shared_ptr<ModelProperties> const& The ModelProperties used by ToolboxMor
     */
    std::shared_ptr<CRBModelProperties> const& modelProperties() const { return M_modelProperties; }
    std::vector<std::string> const& outputDeimName() const { return M_outputDeimName; }

  private :
    void assembleData();

    void updateBetaQ_impl( parameter_type const& mu , double time , bool only_terms_time_dependent );

    template<typename ModelType>
    typename super_type::betaqm_type
    computeBetaQ_impl( parameter_type const& mu , double time , bool only_terms_time_dependent=false,
                       typename std::enable_if< !ModelType::is_time_dependent >::type* = nullptr )
        {
            this->updateBetaQ_impl( mu,time,only_terms_time_dependent);
            return boost::make_tuple( this->M_betaAqm, this->M_betaFqm );
        }

    template<typename ModelType>
    typename super_type::betaqm_type
    computeBetaQ_impl( parameter_type const& mu , double time , bool only_terms_time_dependent=false,
                       typename std::enable_if< ModelType::is_time_dependent >::type* = nullptr )
        {
            this->updateBetaQ_impl( mu,time,only_terms_time_dependent);
            return boost::make_tuple( this->M_betaMqm, this->M_betaAqm, this->M_betaFqm );
        }


    template<typename ModelType>
    typename super_type::betaqm_type
    computeBetaQ_impl( parameter_type const& mu, typename std::enable_if< !ModelType::is_time_dependent >::type* = nullptr )
        {
            this->updateBetaQ_impl( mu,0,false);
            return boost::make_tuple( this->M_betaAqm, this->M_betaFqm );
        }

    template<typename ModelType>
    typename super_type::betaqm_type
    computeBetaQ_impl( parameter_type const& mu, typename std::enable_if< ModelType::is_time_dependent >::type* = nullptr )
        {
            this->updateBetaQ_impl( mu,0,false);
            return boost::make_tuple( this->M_betaMqm, this->M_betaAqm, this->M_betaFqm );
        }

    int Qa();
    int Nl();
    int Ql( int l );
    int mQA( int q );
    int mLQF( int l, int q );
    void resizeQm( bool resizeMat = true );

private:

    std::shared_ptr<CRBModelProperties> M_modelProperties;

    int M_trainsetDeimSize;
    int M_trainsetMdeimSize;

    deim_ptrtype M_deim;
    mdeim_ptrtype M_mdeim;

    deim_function_type M_assembleForDEIM;
    mdeim_function_type M_assembleForMDEIM;

    std::vector<std::string> M_outputName;
    std::map<std::string, int> M_outputDeim;
    std::vector<std::string> M_outputDeimName;

}; // class ToolboxMor

} // namespace Feel

#include <feel/feelmor/toolboxmor_impl.hpp>

#endif

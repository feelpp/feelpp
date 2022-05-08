#ifndef FEELPP_CRBMODELPROPERTIES_HPP
#define FEELPP_CRBMODELPROPERTIES_HPP 1

#include <feel/feelcore/json.hpp>
#include <feel/feelcore/commobject.hpp>
#include <feel/feelmor/crbmodelparameters.hpp>
#include <feel/feelmor/crbmodeloutputs.hpp>

namespace Feel {

/**
 * @brief Class to store properties of a CRB model
 * 
 */
class FEELPP_EXPORT CRBModelProperties : public CommObject
{
public:
    using super = CommObject;
    /**
     * @brief Construct a new CRBModelProperties object
     * 
     * @param directoryLibExpr directory to store the expression library
     * @param world world communicator
     * @param prefix prefix for the option
     * @param vm options set to use
     */
    CRBModelProperties( std::string const& directoryLibExpr = "",
                        worldcomm_ptr_t const& world = Environment::worldCommPtr(),
                        std::string const& prefix="",
                        po::variables_map const& vm = Environment::vm() );

    /**
     * @brief Construct a new CRBModelProperties object
     * 
     * @tparam T a nl::json object
     * @param jarg the json representing the properties
     * @param directoryLibExpr 
     * @param world world communicator
     * @param prefix prefix for the option
     * @param vm options set to use
     */
    template <typename T,std::enable_if_t< std::is_same_v<std::decay_t<T>,nl::json>, bool> = true >
    CRBModelProperties( T && jarg,
                        std::string const& directoryLibExpr = "",
                        worldcomm_ptr_t const& world = Environment::worldCommPtr(),
                        std::string const& prefix="",
                        po::variables_map const& vm = Environment::vm() )
        :
        CRBModelProperties( directoryLibExpr,world,prefix,vm )
    {
        M_jsonData = std::forward<T>( jarg );
        this->setupImpl();
    }

    /**
     * @brief Destroy the CRBModelProperties object
     * 
     */
    virtual ~CRBModelProperties() {}

    //! setup from a collection json filename
    void setup( std::vector<std::string> const& filename );
    //! setup from json filename
    void setup( std::string const& filename ) { this->setup( std::vector<std::string>{ filename } ); }
    //! setup from json object
    template <typename T,std::enable_if_t< std::is_same_v<std::decay_t<T>,nl::json>, bool> = true >
    void setup( T && jarg )
    {
        M_jsonData = std::forward<T>( jarg );
        this->setupImpl();
    }
    //! setup from filename option
    void setupFromFilenameOption( po::variables_map const& vm = Environment::vm() );

    nl::json const& jsonData() const { return M_jsonData; }

    std::string const& name() const {  return M_name; }
    void setName( std::string const& t) { M_name = t; }
    std::string const& shortName() const {  return M_shortname; }
    void setShortName( std::string const& t) { M_shortname = t; }

    std::string const& description() const {  return M_description; }
    void setDescription( std::string const& t) { M_description = t; }

    CRBModelParameters & parameters()  {  return M_params; }
    CRBModelParameters const& parameters() const {  return M_params; }

    CRBModelOutputs & outputs() { return M_outputs; }
    CRBModelOutputs const& outputs() const { return M_outputs; }

  private :
    static nl::json read_json( std::string const& filename, worldcomm_ptr_t const& world );

    void setupImpl();

private:
    std::string M_prefix, M_directoryLibExpr;
    nl::json M_jsonData;
    nl::json M_json_merge_patch;
    nl::json::array_t M_json_patch;

    std::string M_name, M_shortname, M_description;
    CRBModelParameters M_params;
    CRBModelOutputs M_outputs;
};

} // namespace Feel

#endif

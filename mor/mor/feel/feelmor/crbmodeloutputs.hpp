#ifndef FEELPP_CRBCRBMODELOUTPUTS_HPP
#define FEELPP_CRBCRBMODELOUTPUTS_HPP 1

#include <feel/feelcore/commobject.hpp>
#include <feel/feelmodels/modelexpression.hpp>
#include <feel/feelmodels/modelmarkers.hpp>

namespace Feel
{

/**
 * @brief A class to store an output of a CRB model
 * 
 */
class FEELPP_EXPORT CRBModelOutput : public CommObject
{
public:
    using super = CommObject;
    /**
     * @brief Construct a new CRBModelOutput object
     * 
     * @param worldComm world communicator
     */
    CRBModelOutput( worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    CRBModelOutput( CRBModelOutput const& ) = default;
    CRBModelOutput( CRBModelOutput&& ) = default;
    CRBModelOutput& operator=( CRBModelOutput const& ) = default;
    CRBModelOutput& operator=( CRBModelOutput && ) = default;

    /**
     * @brief Construct a new CRBModelOutput object
     * 
     * @param name name of the output
     * @param jarg json section of the output
     * @param worldComm world communicator
     * @param directoryLibExpr directory to store the expression library
     */
    CRBModelOutput( std::string name, nl::json const& jarg,
                    worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                    std::string const& directoryLibExpr = "" );

    std::string name() const { return M_name; }
    std::string type() const { return M_type; }
    std::set<std::string> markers() const { return M_markers; }
    std::vector<double> const& coord() const { return M_coord; }
    double radius() const { return M_radius; }
    int dim() const { return M_dim; }
    std::string getString( std::string const& key ) const;

    template <int M=1,int N=1>
    bool hasExpression() const { return M_expr.hasExpr<M,N>(); }

    template <int M=1,int N=1>
        auto const& expression() const
    {
        if ( !this->hasExpression<M,N>() )
            CHECK( false ) << "no expression defined";
        return M_expr.template expr<M,N>();
    }

    bool isEvaluable() const
    {
        if ( this->type() == "expression" || this->type() == "value" )
            return M_expr.isEvaluable();
        else
            return false;
    }

    auto evaluate() const
    {
        return M_expr.evaluate();
    }
    void setParameterValues( std::map<std::string,double> const& mp )
    {
        M_expr.setParameterValues( mp );
    }

  private:
    nl::json M_p;
    std::string M_directoryLibExpr;

    std::string M_name;
    std::string M_type;
    ModelMarkers M_markers;
    std::vector<double> M_coord;
    double M_radius;
    int M_dim;
    ModelExpression M_expr;

};

/**
 * @brief A class to store the outputs of a CRB model
 * 
 */
class CRBModelOutputs: public std::map<std::string, CRBModelOutput>, public CommObject
{
public:
    using super=CommObject;
    /**
     * @brief Construct a new CRBModelOutputs object
     * 
     * @param world world communicator
     */
    CRBModelOutputs( worldcomm_ptr_t const& world = Environment::worldCommPtr() );
    /**
     * @brief Construct a new CRBModelOutputs object
     * 
     * @param jarg json section of the outputs
     * @param worldComm world communicator
     */
    CRBModelOutputs( nl::json const& jarg, worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    CRBModelOutputs( CRBModelOutputs const& ) = default;
    virtual ~CRBModelOutputs();
    /**
     * @brief set the json section for the outputs
     * 
     * @param jarg json section for the outputs
     */
    void setPTree( nl::json const& jarg );
    /**
     * @brief Set the Directory Lib Expr object
     * 
     * @param directoryLibExpr directory where to store expression libraries
     */
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }
    /**
     * @brief load an output from a json file
     * 
     * @param filename path to the json file
     * @param name name of the output
     * @return CRBModelOutput 
     */
    CRBModelOutput loadOutput( std::string const& filename, std::string const& name);
    /**
     * @brief Get an output from a json section
     * 
     * @param jarg the json section
     * @param name name of the output
     * @return CRBModelOutput 
     */
    CRBModelOutput getOutput( nl::json const& jarg, std::string const& name);

    /**
     * @brief returns all output of type type
     * 
     * @param type type of the outputs
     * @return std::map<std::string, CRBModelOutput> 
     */
    std::map<std::string, CRBModelOutput> ofType(std::string const& type);
    std::map<std::string, CRBModelOutput> ofTypes(std::set<std::string> const& types);

    void setParameterValues( std::map<std::string,double> const& mp );

private:
    void setup();
private:
    nl::json M_p;
    std::string M_directoryLibExpr;
};

} // namespace Feel

#endif

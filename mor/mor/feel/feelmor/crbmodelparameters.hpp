#ifndef FEELPP_CRBCRBMODELPARAMETERS_HPP
#define FEELPP_CRBCRBMODELPARAMETERS_HPP 1

#include <feel/feelcore/commobject.hpp>
#include <feel/feelcore/json.hpp>

namespace Feel {

/**
 * @brief A class to store a parameter of a CRB model
 * 
 */
struct FEELPP_EXPORT CRBModelParameter
{
    CRBModelParameter() = default;
    CRBModelParameter( CRBModelParameter const& ) = default;
    CRBModelParameter( CRBModelParameter&& ) = default;
    CRBModelParameter& operator=( CRBModelParameter const& ) = default;
    CRBModelParameter& operator=( CRBModelParameter && ) = default;

    /**
     * @brief Construct a new CRBModelParameter object
     * 
     * @param name name of the parameter
     * @param min minimum value of the parameter
     * @param max maximum value of the parameter
     * @param value value of the parameter
     * @param samplingSize size of the sampling to use for the parameter
     * @param sampling type of sampling to use for the parameter
     * @param desc description of the parameter
     */
    CRBModelParameter( std::string const& name, double min, double max, double value = 0., int samplingSize = 0, std::string const& sampling = "", std::string const& desc = "" )
        :
        M_name( name ),
        M_min( min ),
        M_max( max ),
        M_value( value ),
        M_samplingSize( samplingSize ),
        M_sampling( sampling ),
        M_desc( desc )
        {}

    std::string const& name() const { return M_name; }
    void setName( std::string const& name ) { M_name = name; }
    double min() const { return M_min; }
    void setMin( double v ) { M_min = v; }
    double max() const { return M_max; }
    void setMax( double v ) { M_max = v; }
    double value() const { return M_value; }
    void setValue( double v ) { M_value = v; }
    int samplingSize() const { return M_samplingSize; }
    void setSamplingSize( int s ) { M_samplingSize = s; } 
    std::string const& sampling() const { return M_sampling; }
    void setSampling( std::string const& sampling ) { M_sampling = sampling; }
    std::string const& description() const { return M_desc; }
    void setDescription( std::string const& desc ) { M_desc = desc; }
 
private:
    std::string M_name;
    double M_min, M_max;
    double M_value;
    int M_samplingSize;
    std::string M_sampling;
    std::string M_desc;

};

/**
 * @brief A class to store the parameters of a CRB model
 * 
 */
class CRBModelParameters: public std::map<std::string,CRBModelParameter>, public CommObject
{
public:
    using super=CommObject;
    /**
     * @brief Construct a new CRBModelParameters object
     * 
     * @param world world communicator
     */
    CRBModelParameters( worldcomm_ptr_t const& world = Environment::worldCommPtr() );
    CRBModelParameters( CRBModelParameters const& ) = default;
    virtual ~CRBModelParameters();
    /**
     * @brief set the json section for the parameters
     * 
     * @param jarg json section for the parameters
     */
    void setPTree( nl::json const& jarg );

private:
    void setup();
private:
    nl::json M_p;
};

} // namespace Feel

#endif

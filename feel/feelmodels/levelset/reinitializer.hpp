#ifndef _REINITIALIZER_HPP
#define _REINITIALIZER_HPP 1

namespace Feel
{

enum class ReinitializerType
{
    FM,
    HJ
};

template<typename FunctionSpaceType>
class Reinitializer
{
public:
    //--------------------------------------------------------------------//
    // Typedefs
    typedef Reinitializer<FunctionSpaceType> reinitializer_type;
    typedef boost::shared_ptr<reinitializer_type> reinitializer_ptrtype;

    typedef FunctionSpaceType functionspace_type;
    typedef typename functionspace_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;
    
    //--------------------------------------------------------------------//
    // Constructor/Destructor
    Reinitializer( std::string const& prefix );
    virtual ~Reinitializer() {}

    //static reinitializer_ptrtype build( std::string const& type, std::string const& prefix="" );
    //--------------------------------------------------------------------//
    std::string const& prefix() const { return M_prefix; } 
    ReinitializerType type() const { return M_reinitializerType; }
    //--------------------------------------------------------------------//
    // Run reintialization
    virtual element_type run( element_type const& phi ) =0;

    element_type operator() ( element_type const& phi ) { return this->run(phi); }

protected:
    ReinitializerType M_reinitializerType;

private:
    std::string M_prefix;
};

template<typename FunctionSpaceType>
Reinitializer<FunctionSpaceType>::Reinitializer( std::string const& prefix )
    : M_prefix( prefix )
{}

} // Feel

#endif

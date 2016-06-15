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
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;

    typedef typename functionspace_type::mesh_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef typename functionspace_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

    typedef typename functionspace_type::periodicity_0_type periodicity_type;
    static const bool is_periodic = functionspace_type::is_periodic;
    
    //--------------------------------------------------------------------//
    // Constructor/Destructor
    Reinitializer( 
            functionspace_ptrtype const& space,
            std::string const& prefix ) : 
        M_space(space),
        M_periodicity(boost::fusion::at_c<0>(space->periodicity())),
        M_prefix(prefix)
        {}
    virtual ~Reinitializer() = default;

    //static reinitializer_ptrtype build( std::string const& type, std::string const& prefix="" );
    //--------------------------------------------------------------------//
    //ReinitializerType type() const { return M_reinitializerType; }
    std::string const& prefix() const { return M_prefix; }

    virtual void loadParametersFromOptionsVm() {}
    //--------------------------------------------------------------------//
    functionspace_ptrtype const& functionSpace() const { return M_space; }
    mesh_ptrtype const& mesh() const { return M_space->mesh(); }
    //--------------------------------------------------------------------//
    // Run reinitialization
    virtual element_type run( element_type const& phi ) =0;

    element_type operator() ( element_type const& phi ) { return this->run(phi); }

protected:
    //ReinitializerType M_reinitializerType;

    functionspace_ptrtype M_space;
    periodicity_type M_periodicity;

    std::string M_prefix;
};

} // Feel

#endif

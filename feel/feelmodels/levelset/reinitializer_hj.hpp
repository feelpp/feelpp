#ifndef _REINITIALIZER_HJ_HPP
#define _REINITIALIZER_HJ_HPP 1

namespace Feel
{

template<typename FunctionSpaceType>
class ReinitializerHJ
: public Reinitializer<FunctionSpaceType>
{
public:
    //--------------------------------------------------------------------//
    // Typedefs
    typedef Reinitializer<FunctionSpaceType> super_type;
    typedef ReinitializerHJ<FunctionSpaceType> self_type;
    typedef boost:shared_ptr<self_type> self_ptrtype;

    typedef FunctionSpaceType functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;

    typedef typename functionspace_type::mesh_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef typename functionspace_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

    typedef typename functionspace_type::periodicity_0_type periodicity_type;
    static const bool is_periodic = functionspace_type::is_periodic;

    //--------------------------------------------------------------------//
    // Hamilton-Jacobi advection
    template<
        typename ConvexType,
        typename BasisAdvectionType,
        typename PeriodicityType = NoPeriodicity
           >
    class AdvectionHJ
        : public AdvectionBase<ConvexType, BasisAdvectionType, PeriodicityType>
        , public boost::enable_shared_from_this< AdvectionHJ<ConvexType, BasisAdvectionType, PeriodicityType> >
    {
    public:
        typedef AdvectionBase<ConvexType, BasisAdvectionType, PeriodicityType> super_type;

        typedef AdvectionHJ<ConvexType, BasisAdvectionType, PeriodicityType> self_type;
        typedef boost::shared_ptr<self_type> self_ptrtype;

        // Constructor
        AdvectionHJ(
                std::string const& prefix,
                );
        
        static self_ptrtype New(
                std::string const& prefix
                );

        // Initialization
        void init( bool buildModelAlgebraicFactory = true );

        // BC and source term assembly
        void updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const;
        void updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const;
        void updateSourceTermLinearPDE(element_advection_ptrtype& fieldSource, bool buildCstPart) const;
        bool hasSourceTerm() const;
    };

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    // Constructor
    ReinitializerHJ( functionspace_ptrtype const& space );
    //--------------------------------------------------------------------//
    // Run reinitialization
    element_type run( element_type const& phi );
};

template<typename FunctionSpaceType>
ReinitializerHJ<FunctionSpaceType>::ReinitializerHJ( functionspace_ptrtype const& space )
    : super_type( space )
{}

template<typename FunctionSpaceType>
typename ReinitializerHJ<FunctionSpaceType>::element_type
ReinitializerHJ<FunctionSpaceType>::run( element_type const& phi )
{
}

} // namespace Feel

#endif

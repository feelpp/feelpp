#ifndef _ADVECTION_HPP
#define _ADVECTION_HPP 1

#include <feel/feelmodels/advection/advectionbase.hpp>

namespace Feel {
namespace FeelModels {

template< typename ConvexType, typename BasisAdvectionType >
class Advection
    : public AdvectionBase<ConvexType, BasisAdvectionType>
    , public boost::enable_shared_from_this< Advection<ConvexType, BasisAdvectionType> >
{
public:
    typedef AdvectionBase<ConvexType, BasisAdvectionType> super_type;

    typedef Advection<ConvexType, BasisAdvectionType> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    typedef typename super_type::space_advection_ptrtype space_advection_ptrtype;

    //--------------------------------------------------------------------//
    // Constructor
    Advection( 
            std::string const& prefix,
            WorldComm const& _worldComm = Environment::worldComm(),
            std::string const& subPrefix = "",
            std::string const& rootRepository = ModelBase::rootRepositoryByDefault() );

    static self_ptrtype New( 
            std::string const& prefix,
            WorldComm const& _worldComm = Environment::worldComm(),
            std::string const& subPrefix = "",
            std::string const& rootRepository = ModelBase::rootRepositoryByDefault() );
    //--------------------------------------------------------------------//
    // Initialization
    void init( bool buildModelAlgebraicFactory = true );
    void loadConfigBCFile();
    //--------------------------------------------------------------------//
    // BC and source term assembly
    void updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const;
    void updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const;
    void updateSourceTermLinearPDE(vector_ptrtype& F, bool buildCstPart) const;

protected:
    // Boundary conditions
    map_scalar_field<2> M_bcDirichlet;
    map_scalar_field<2> M_bcNeumann;
    map_scalar_fields<2> M_bcRobin;

};
    

} // namespace FeelModels
} // namespace Feel

#endif

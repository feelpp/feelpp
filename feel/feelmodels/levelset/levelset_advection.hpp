#ifndef _LEVELSETADVECTION_HPP
#define _LEVELSETADVECTION_HPP 1

#include <feel/feelmodels/advection/advectionbase.hpp>

namespace Feel {
namespace FeelModels {

template<typename ConvexType, typename BasisAdvectionType>
class LevelSetAdvection
    : public AdvectionBase<ConvexType, BasisAdvectionType>
    , public boost::enable_shared_from_this< LevelSetAdvection<ConvexType, BasisAdvectionType> >
{
public:
    typedef AdvectionBase<ConvexType, BasisAdvectionType> super_type;

    typedef LevelSetAdvection<ConvexType, BasisAdvectionType> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    typedef typename super_type::space_advection_ptrtype space_advection_ptrtype;

    //--------------------------------------------------------------------//
    // Constructor
    LevelSetAdvection(
            std::string const& prefix,
            WorldComm const& _worldComm = Environment::worldComm(),
            std::string const& subPrefix = "",
            std::string const& rootRepository = ModelBase::rootRepositoryByDefault() );

    //--------------------------------------------------------------------//
    // Initialization
    void init( bool buildModelAlgebraicFactory = true );

    //--------------------------------------------------------------------//
    // BC and source term assembly
    void updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const;
    void updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const;

    bool hasSourceTerm() const { return false; }

};

template<typename ConvexType, typename BasisAdvectionType>
LevelSetAdvection<ConvexType, BasisAdvectionType>::LevelSetAdvection(
        std::string const& prefix,
        WorldComm const& _worldComm = Environment::worldComm(),
        std::string const& subPrefix = "",
        std::string const& rootRepository = ModelBase::rootRepositoryByDefault() )
: super_type( prefix, worldComm, subPrefix, rootRepository )
{
    this->log("LevelSetAdvection", "constructor", "start" );
    
    this->setFilenameSaveInfo( prefixvm(this->prefix(), "LevelSetAdvection.info") );
    //-----------------------------------------------------------------------------//
    // Set model (advection only)
    this->setModelName( "Advection" );

    this->log("LevelSetAdvection", "constructor", "finish");
}

template<typename ConvexType, typename BasisAdvectionType>
void
LevelSetAdvection<ConvexType, BasisAdvectionType>::init( bool buildModelAlgebraicFactory )
{
    super_type::init( buildModelAlgebraicFactory, this->shared_from_this() );
}
    

}
}

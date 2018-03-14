#ifndef _SURFACETENSIONFORCEMODEL_HPP
#define _SURFACETENSIONFORCEMODEL_HPP 1


#include <feel/feelmodels/multifluid/interfaceforcesmodel.hpp>

namespace Feel {
namespace FeelModels {

template<class LevelSetType>
class SurfaceTensionForceModel
: public virtual InterfaceForcesModel<LevelSetType>
{
    typedef SurfaceTensionForceModel<LevelSetType> self_type;
    typedef InterfaceForcesModel<LevelSetType> super_type;
public:
    typedef typename super_type::levelset_type levelset_type;
    typedef typename super_type::levelset_ptrtype levelset_ptrtype;

    typedef typename super_type::mesh_type mesh_type;

    typedef typename levelset_type::space_levelset_vectorial_type space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;

    typedef typename space_type::element_type element_type;
    typedef typename space_type::element_ptrtype element_ptrtype;

    static const uint16_type Dim = levelset_type::nDim;

    //--------------------------------------------------------------------//
    // Construction
    SurfaceTensionForceModel() = default;
    SurfaceTensionForceModel( SurfaceTensionForceModel const& i ) = default;

    void build( std::string const& prefix, levelset_ptrtype const& ls ) override;

    boost::shared_ptr<std::ostringstream> getInfo() const override;

    void loadParametersFromOptionsVm();
    //--------------------------------------------------------------------//
    double surfaceTensionCoeff() const { return M_surfaceTensionCoeff; }
    void setSurfaceTensionCoeff( double k ) { M_surfaceTensionCoeff = k; }

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
private:
    void updateInterfaceForcesImpl( element_ptrtype & F ) const override;

    //--------------------------------------------------------------------//
    double M_surfaceTensionCoeff;

};

template<typename LevelSetType>
void
SurfaceTensionForceModel<LevelSetType>::build( std::string const& prefix, levelset_ptrtype const& ls )
{
    super_type::build( prefix, ls );
    this->loadParametersFromOptionsVm();
}

template<typename LevelSetType>
boost::shared_ptr<std::ostringstream> 
SurfaceTensionForceModel<LevelSetType>::getInfo() const
{
    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << "Surface tension ("
           << "coeff = " << this->M_surfaceTensionCoeff
           << ")";

    return _ostr;
}

template<typename LevelSetType>
void
SurfaceTensionForceModel<LevelSetType>::loadParametersFromOptionsVm()
{
    M_surfaceTensionCoeff = doption( _name="surface-tension.coeff", _prefix=this->prefix() ); 
}

template<typename LevelSetType>
void
SurfaceTensionForceModel<LevelSetType>::updateInterfaceForcesImpl( element_ptrtype & F ) const
{
    auto N = this->levelset()->N();
    auto K = this->levelset()->K();
    F->on( 
            _range=elements(this->levelset()->mesh()),
            _expr= - this->surfaceTensionCoeff()*idv(K)*idv(N)*idv(this->levelset()->D())
         );
}


} // namespace FeelModels
} // namespace Feel

#endif

#ifndef _HELFRICHFORCEMODEL_HPP
#define _HELFRICHFORCEMODEL_HPP 1

namespace Feel {
namespace FeelModels {

template<class LevelSetType>
class HelfrichForceModel
: public InterfaceForcesModel<LevelSetType>
{
    typedef HelfrichForceModel<LevelSetType> self_type;
    typedef InterfaceForcesModel<LevelSetType> super_type;
public:
    typedef typename super_type::levelset_type levelset_type;
    typedef typename super_type::levelset_ptrtype levelset_ptrtype;

    typedef typename levelset_type::space_levelset_vectorial_type space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;

    typedef typename space_type::element_type element_type;
    typedef typename space_type::element_ptrtype element_ptrtype;

    //--------------------------------------------------------------------//
    // Construction
    HelfrichForceModel( std::string const& prefix, levelset_ptrtype const& ls );
    HelfrichForceModel( HelfrichForceModel const& i ) = default;

    void loadParametersFromOptionsVm();
    //--------------------------------------------------------------------//
    double bendingModulus() const { return M_helfrichBendingModulus; }
    void setBendingModulus( double k ) { M_helfrichBendingModulus = k; }

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
private:
    void updateInterfaceForcesImpl( element_ptrtype & F );
    void addHelfrichForce( element_ptrtype & F, int impl );

    //--------------------------------------------------------------------//
    double M_helfrichBendingModulus;
    int M_forceImpl;

};

template<typename LevelSetType>
HelfrichForceModel<LevelSetType>::HelfrichForceModel(
        std::string const& prefix,
        levelset_ptrtype const& ls
        ) :
    super_type( prefix, ls )
{}

template<typename LevelSetType>
void
HelfrichForceModel<LevelSetType>::loadParametersFromOptionsVm()
{
    M_helfrichBendingModulus = doption( _name="helfrich-bending-modulus", _prefix=this->prefix() ); 
    M_forceImpl = doption( _name="helfrich-force-impl", prefix=this->prefix() );
}

template<typename LevelSetType>
void
HelfrichForceModel<LevelSetType>::updateInterfaceForcesImpl( element_ptrtype & F )
{
    this->addHelfrichForce( F, M_forceImpl );
}

template<typename LevelSetType>
void
HelfrichForceModel<LevelSetType>::addHelfrichForce( element_ptrtype & F, int impl )
{
    switch impl:
    case 0:
    {
        auto phi = this->levelset()->phi();
        auto n_l2 = this->levelset()->N();
        auto k_l2 = this->levelset()->K();
        auto t_int = vf::vec(trans(idv(n_l2)) * vf::oneY(),  - trans(idv(n_l2)) * vf::oneX() );
        auto modgradphi = this->levelset()->smoother()->project( sqrt( gradv( phi ) * trans(gradv(phi)) ) );
        auto AA = this->levelset()->projectorL2Vectorial->project( - idv(k_l2) * idv(k_l2) / 2. * trans(idv(n_l2)) );
        auto BB = this->levelset()->smootherVectorial()->project( trans(t_int) * ( ( gradv(modgradphi) * idv(k_l2) ) * t_int) / max(idv(modgradphi), 0.01) );
        auto Fc_global = this->levelset()->smootherVectorial()->project( this->bendingModulus() * ( divv( AA ) + divv( BB ) ) * gradv( phi ) );
        F += vf::project(
                this->levelset()->functionSpaceVectorial(), 
                elements(this->levelset()->mesh()), 
                idv(Fc_global) * idv(this->levelset()->D()) 
                );
    }
    break;
    default:
        CHECK(false) << "Wrong force implementation\n";
}

} // namespace FeelModels
} // namespace Feel

#endif

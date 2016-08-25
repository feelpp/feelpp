#ifndef _HELFRICHFORCEMODEL_HPP
#define _HELFRICHFORCEMODEL_HPP 1

#include <feel/feelmodels/multifluid/interfaceforcesmodel.hpp>

namespace Feel {
namespace FeelModels {

template<class LevelSetType>
class HelfrichForceModel
: public virtual InterfaceForcesModel<LevelSetType>
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

    static const uint16_type Dim = levelset_type::nDim;

    //--------------------------------------------------------------------//
    // Construction
    HelfrichForceModel() = default;
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
void
HelfrichForceModel<LevelSetType>::loadParametersFromOptionsVm()
{
    M_helfrichBendingModulus = doption( _name="helfrich-bending-modulus", _prefix=this->prefix() ); 
    M_forceImpl = ioption( _name="helfrich-force-impl", _prefix=this->prefix() );
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
    switch (impl)
    {
        case 0:
        {
//#if FEELPP_DIM == 2
            //auto phi = this->levelset()->phi();
            //auto n_l2 = this->levelset()->N();
            //auto k_l2 = this->levelset()->K();
            //auto t_int = vf::vec(trans(idv(n_l2)) * vf::oneY(),  - trans(idv(n_l2)) * vf::oneX() );
            //auto modgradphi = this->levelset()->smoother()->project( sqrt( gradv( phi ) * trans(gradv(phi)) ) );
            //auto AA = this->levelset()->projectorL2Vectorial()->project( - idv(k_l2) * idv(k_l2) / 2. * idv(n_l2) );
            //auto BB = this->levelset()->smootherVectorial()->project( trans( trans(t_int) * ( ( gradv(modgradphi) * idv(k_l2) ) * t_int) / max(idv(modgradphi), 0.01) ) );
            //auto Fc_global = this->levelset()->smootherVectorial()->project( this->bendingModulus() * ( divv( AA ) + divv( BB ) ) * trans(gradv( phi )) );
            //*F += vf::project(
                    //_space=this->levelset()->functionSpaceVectorial(), 
                    //_range=elements(this->levelset()->mesh()), 
                    ////idv(Fc_global) * idv(this->levelset()->D()) 
                    //_expr=trans(gradv(phi))
                    //);
//#else
            auto phi = this->levelset()->phi();
            auto N = this->levelset()->N();
            auto K = this->levelset()->K();
            auto Id = vf::Id<Dim, Dim>();
            auto NxN = idv(N)*trans(idv(N));
            auto modGradPhiK = this->levelset()->smoother()->project( 
                    sqrt( gradv(phi) * trans(gradv(phi)) ) * idv(K)
                    );

            auto Fb_par = this->levelset()->projectorL2Vectorial()->project(
                    - idv(K)*idv(K)/2. * idv(N)
                    );
            auto Fb_ortho = this->levelset()->smootherVectorial()->project(
                    (Id-NxN)*trans(gradv(modGradPhiK)) / sqrt( gradv(phi) * trans(gradv(phi)) )
                    );
            auto Fb_global = this->levelset()->smootherVectorial()->project(
                    this->bendingModulus() * (divv(Fb_par) + divv(Fb_ortho)) * idv(N)
                    );

            *F += vf::project(
                    _space=this->levelset()->functionSpaceVectorial(),
                    _range=elements(this->levelset()->mesh()),
                    _expr=idv(Fb_global) * idv(this->levelset()->D())
                    );

//#endif
        }
        break;
        default:
            CHECK(false) << "Wrong force implementation\n";
    }
}

} // namespace FeelModels
} // namespace Feel

#endif

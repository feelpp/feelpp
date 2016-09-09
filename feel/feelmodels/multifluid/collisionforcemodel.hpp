#ifndef _COLLISIONFORCEMODEL_HPP
#define _COLLISIONFORCEMODEL_HPP 1

#include <feel/feelmodels/multifluid/interfaceforcesmodel.hpp>

namespace Feel {
namespace FeelModels {

template<class LevelSetType>
class CollisionForceModel
: public virtual InterfaceForcesModel<LevelSetType>
{
    typedef CollisionForceModel<LevelSetType> self_type;
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
    CollisionForceModel() = default;
    CollisionForceModel( CollisionForceModel const& i ) = default;

    void loadParametersFromOptionsVm();
    //--------------------------------------------------------------------//
    double collisionCoeff() const { return M_collisionCoeff; }
    void setCollisionCoeff( double k ) { M_collisionCoeff = k; }

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
private:
    void updateInterfaceForcesImpl( element_ptrtype & F );
    void addCollisionForce( element_ptrtype & F, int impl );

    //--------------------------------------------------------------------//
    double M_collisionCoeff;
    int M_forceImpl;

#ifdef DEBUG_COLLISIONFORCEMODEL
    typedef boost::shared_ptr<Exporter<mesh_type, 1>> exporter_ptrtype;
    exporter_ptrtype M_exporter;
#endif
};

template<typename LevelSetType>
void
CollisionForceModel<LevelSetType>::loadParametersFromOptionsVm()
{
    M_collisionCoeff = doption( _name="collision-force-coeff", _prefix=this->prefix() ); 
    M_forceImpl = ioption( _name="collision-force-impl", _prefix=this->prefix() );

#ifdef DEBUG_COLLISIONFORCEMODEL
    M_exporter = Feel::exporter(
            _mesh=this->levelset()->mesh(),
            _name="CollisionForce",
            _path=this->levelset()->exporterPath()
            );
#endif
}

template<typename LevelSetType>
void
CollisionForceModel<LevelSetType>::updateInterfaceForcesImpl( element_ptrtype & F )
{
    this->addCollisionForce( F, M_forceImpl );
}

template<typename LevelSetType>
void
CollisionForceModel<LevelSetType>::addCollisionForce( element_ptrtype & F, int impl )
{
    switch (impl)
    {
        case 0:
        {
            auto phi = this->levelset()->phi();
            auto phi2 = this->levelset()->distanceBetweenLabels(); 
            double epsilon = this->levelset()->thicknessInterface();

            auto gradPhi2 = this->levelset()->projectorL2Vectorial()->project(
                    _expr=trans(gradv(phi2))
                    );

            auto Fc = vf::project(
                    _space=this->levelset()->functionSpaceVectorial(),
                    _range=elements(this->levelset()->mesh()),
                    _expr=this->collisionCoeff() * exp(-idv(phi2)/epsilon) * idv(gradPhi2) / idv(phi2) * idv(this->levelset()->D())
                    );

            *F += Fc;

#ifdef DEBUG_COLLISIONFORCEMODEL
            M_exporter->step(this->levelset()->time())->add(
                    "dirac", "dirac", this->levelset()->D() );
            M_exporter->step(this->levelset()->time())->add(
                    "collision-force", "collision-force", Fc );
            M_exporter->save();
#endif
        }
        break;
        default:
            CHECK(false) << "Wrong force implementation\n";
    }
}

} // namespace FeelModels
} // namespace Feel

#endif

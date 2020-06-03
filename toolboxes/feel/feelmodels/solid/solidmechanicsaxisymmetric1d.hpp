

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisDisplacementType>
class SolidMechanicsAxisymmetric1d : public ModelNumerical,
                                     public ModelPhysics<ConvexType::nRealDim>,
                                     public std::enable_shared_from_this< SolidMechanicsAxisymmetric1d<ConvexType,BasisDisplacementType> >
{
public :
    using super_type = ModelNumerical;
    using size_type = typename super_type::size_type;
    typedef SolidMechanicsAxisymmetric1d<ConvexType,BasisDisplacementType> self_type;
    // TODO static assert only 1d
    typedef std::shared_ptr<self_type> self_ptrtype;
    //___________________________________________________________________________________//
    // mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    static const uint16_type nRealDim = convex_type::nRealDim;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    // basis
    using basis_displacement_type = BasisDisplacementType;
    static const uint16_type nOrderDisplacement = basis_displacement_type::nOrder;
    // function space temperature
    typedef FunctionSpace<mesh_type, bases<basis_displacement_type> > space_displacement_type;
    typedef std::shared_ptr<space_displacement_type> space_displacement_ptrtype;
    typedef typename space_displacement_type::element_type element_displacement_type;
    typedef std::shared_ptr<element_displacement_type> element_displacement_ptrtype;


    SolidMechanicsAxisymmetric1d( std::string const& prefix,
                                  std::string const& keyword,// = "cfpde",
                                  worldcomm_ptr_t const& worldComm,// = Environment::worldCommPtr(),
                                  std::string const& subPrefix,//  = "",
                                  ModelBaseRepository const& modelRep = ModelBaseRepository() )
        :
        super_type( prefix, keyword, worldComm, subPrefix, modelRep ),
        ModelPhysics<nRealDim>( "solid" ),
        ModelBase( prefix, keyword, worldComm, subPrefix, modelRep )
        {}
private :
    //void initMesh();
private :

};

} // namespace FeelModels
} // namespace Feel

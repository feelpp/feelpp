/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

#ifndef FEELPP_TOOLBOXES_COEFFICIENTFORMPDES_HPP
#define FEELPP_TOOLBOXES_COEFFICIENTFORMPDES_HPP 1

#include <feel/feelmodels/coefficientformpdes/coefficientformpde.hpp>


namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename... BasisUnknownType>
class CoefficientFormPDEs : public ModelNumerical,
                            public ModelGenericPDEs<ConvexType::nDim>,
                            public std::enable_shared_from_this< CoefficientFormPDEs<ConvexType,BasisUnknownType...> >
{
    using coefficient_form_pde_base_type = CoefficientFormPDEBase<ConvexType>;
public :
    typedef ModelNumerical super_type;
    using size_type = typename super_type::size_type;

    using self_type = CoefficientFormPDEs<ConvexType,BasisUnknownType...>;

    using mesh_type = typename coefficient_form_pde_base_type::mesh_type;
    using mesh_ptrtype = typename coefficient_form_pde_base_type::mesh_ptrtype;

    // materials properties
    typedef MaterialsProperties<mesh_type> materialsproperties_type;
    typedef std::shared_ptr<materialsproperties_type> materialsproperties_ptrtype;

    // algebraic solver
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef std::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;


    static constexpr auto tuple_type_unknown_basis = hana::to_tuple(hana::tuple_t<BasisUnknownType...>);
private :

    struct traits
    {
        template <typename TheType>
        using remove_hana_type_t = typename std::decay_t<TheType>::type;

        template <typename TheBasisType>
        using coefficient_form_pde_t = CoefficientFormPDE<ConvexType,remove_hana_type_t<TheBasisType>>;

        template <typename... TheType>
        static constexpr auto
        variant_from_tuple( hana::tuple<TheType...> const& t )
            {
                return std::variant<TheType...>{};
            }
    };
public :

    using variant_unknown_basis_type = std::decay_t<decltype(traits::variant_from_tuple(tuple_type_unknown_basis))>;

    CoefficientFormPDEs( std::string const& prefix,
                         std::string const& keyword = "pdes",
                         worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                         std::string const& subPrefix  = "",
                         ModelBaseRepository const& modelRep = ModelBaseRepository() );

    void init( bool buildModelAlgebraicFactory=true );
    void initAlgebraicFactory();

    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"mesh.path"); }

    //___________________________________________________________________________________//
    // algebraic data and solver
    backend_ptrtype const& backend() const { return  M_backend; }
    BlocksBaseVector<double> const& blockVectorSolution() const { return M_blockVectorSolution; }
    BlocksBaseVector<double> & blockVectorSolution() { return M_blockVectorSolution; }
    size_type nLocalDof() const;
    model_algebraic_factory_ptrtype const& algebraicFactory() const { return M_algebraicFactory; }
    model_algebraic_factory_ptrtype & algebraicFactory() { return M_algebraicFactory; }

    BlocksBaseGraphCSR buildBlockMatrixGraph() const override;
    //int nBlockMatrixGraph() const { return 1; }

    static
    std::string const& unknowBasisTag( variant_unknown_basis_type const& vb );
private :
    void initMesh();
    void initMaterialProperties();

private :

    static const std::vector<std::string> S_unknownBasisTags;

    mesh_ptrtype M_mesh;

    // physical parameters
    materialsproperties_ptrtype M_materialsProperties;

    // algebraic data/tools
    backend_ptrtype M_backend;
    model_algebraic_factory_ptrtype M_algebraicFactory;
    BlocksBaseVector<double> M_blockVectorSolution;

    std::vector<std::shared_ptr<coefficient_form_pde_base_type>> M_coefficientFormPDEs;

};

} // namespace Feel
} // namespace FeelModels

#endif


#ifndef FEELPP_TOOLBOXES_MESH_METRICMESHADAPTATION_H
#define FEELPP_TOOLBOXES_MESH_METRICMESHADAPTATION_H 1


namespace Feel
{
namespace FeelModels
{

template<typename ConvexType>
class MetricMeshAdaptation
{
public :
    typedef ConvexType convex_type;
    typedef Mesh< convex_type > mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    typedef bases<Lagrange<1, Vectorial> > basis_vectorial_type;
    typedef FunctionSpace< mesh_type, basis_vectorial_type > space_vectorial_type;
    typedef std::shared_ptr<space_vectorial_type> space_vectorial_ptrtype;
    typedef typename space_vectorial_type::element_type element_vectorial_type;
    typedef  std::shared_ptr<element_vectorial_type> element_vectorial_ptrtype;

    typedef typename space_vectorial_type::component_functionspace_type space_scalar_type;
    typedef std::shared_ptr<space_scalar_type> space_scalar_ptrtype;
    typedef typename space_scalar_type::element_type element_scalar_type;
    typedef std::shared_ptr<element_scalar_type> element_scalar_ptrtype;

    MetricMeshAdaptation( space_vectorial_ptrtype space, std::string const& prefix );

    void init();
    void update( Expr<GinacExVF<2>> const& e );
    void update( element_scalar_type const& u );

    template < typename ExprT >
    void updateScalarMetric(  Expr<ExprT> const& e )
        {
            M_scalarMetric->on(_range=elements(M_scalarMetric->mesh()),_expr=e);
        }

    mesh_ptrtype mesh() const { return M_vectorialSpace->mesh(); }
    space_vectorial_ptrtype vectorialSpace() const { return M_vectorialSpace; }
    space_scalar_ptrtype scalarSpace() const { return M_vectorialSpace->compSpace(); }
    element_scalar_type const& scalarMetric() const { return *M_scalarMetric; }

    template < typename ExprT >
    element_scalar_ptrtype feProjection( Expr<ExprT> const& e, std::string const& type = "nodal" ) const
        {
            CHECK( type == "nodal" ) << "only nodal projection";
            auto u = this->scalarSpace()->elementPtr();
            u->on(_range=elements(this->scalarSpace()->mesh()),_expr=e);
            return u;
        }
private :

    space_vectorial_ptrtype M_vectorialSpace;
    element_scalar_ptrtype M_scalarMetric;

};


} // namespace FeelModels
} // namespace Feel

#endif

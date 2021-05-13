#include <feel/feelmodels/hdg/mixedpoisson.cpp>

namespace Feel {
namespace FeelModels {

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateHDG & data ) const
{
    condensed_matrix_ptr_t<value_type>& A = data.matrix();
    condensed_vector_ptr_t<value_type>& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string scst = buildCstPart ? "(build cst part)" : "build non cst part)";
    this->log("MixedPoisson", "updateLinearPDE", "start"+scst);

    auto ps = this->spaceProduct();
    auto bbf = blockform2( ps, A );
    auto u = this->fieldFlux();
    auto p = this->fieldPotential();
    auto phat = this->fieldTrace();

    auto tau_constant = cst(M_tauCst);
    auto sc_param = M_useSC ? 0.5 : 1.0;

    // -(p,div(v))_Omega
    bbf( 0_c, 1_c ) += integrate(_range=elements(support(M_Wh)),_expr=-(idt(p)*div(u)));

    // <phat,v.n>_Gamma\Gamma_I
    bbf( 0_c, 2_c ) += integrate(_range=internalfaces(support(M_Wh)),
                                 _expr=( idt(phat)*(leftface(normal(u))+rightface(normal(u)))) );
    bbf( 0_c, 2_c ) += integrate(_range=M_gammaMinusIntegral,
                                 _expr=idt(phat)*normal(u));

    // (div(j),q)_Omega
    bbf( 1_c, 0_c ) += integrate(_range=elements(support(M_Wh)), _expr=id(p)*divt(u));


    // <tau p, w>_Gamma
    bbf( 1_c, 1_c ) += integrate(_range=internalfaces(support(M_Wh)),
                                 _expr=tau_constant *
                                 ( leftfacet( idt(p))*leftface(id(p)) +
                                   rightfacet( idt(p))*rightface(id(p) )));
    bbf( 1_c, 1_c ) += integrate(_range=boundaryfaces(support(M_Wh)),
                                 _expr=tau_constant * id(p)*idt(p));


    // <-tau phat, w>_Gamma\Gamma_I
    bbf( 1_c, 2_c ) += integrate(_range=internalfaces(support(M_Wh)),
                                 _expr=-tau_constant * idt(phat) *
                                 ( leftface( id(p) )+
                                   rightface( id(p) )));
    bbf( 1_c, 2_c ) += integrate(_range=M_gammaMinusIntegral,
                                 _expr=-tau_constant * idt(phat) * id(p) );


    // <j.n,mu>_Omega/Gamma
    bbf( 2_c, 0_c ) += integrate(_range=internalfaces(support(M_Wh)),
                                 _expr=( id(phat)*(leftfacet(normalt(u))+
                                                   rightfacet(normalt(u))) ) );

    // <tau p, mu>_Omega/Gamma
    bbf( 2_c, 1_c ) += integrate(_range=internalfaces(support(M_Wh)),
                                 _expr=tau_constant * id(phat) * ( leftfacet( idt(p) )+
                                                                   rightfacet( idt(p) )));

    // <-tau phat, mu>_Omega/Gamma
    bbf( 2_c, 2_c ) += integrate(_range=internalfaces(support(M_Wh)),
                                 _expr=-sc_param*tau_constant * idt(phat) * id(phat) );

    this->updateLinearPDE( data, this->modelContext() );
}

}
}

#include <feel/feelmodels/hdg/stokes.cpp>

namespace Feel {
namespace FeelModels {

STOKES_CLASS_TEMPLATE_DECLARATIONS
void
STOKES_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateHDG & data ) const
{
    condensed_matrix_ptr_t<value_type>& A = data.matrix();
    condensed_vector_ptr_t<value_type>& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string scst = buildCstPart ? "(build cst part)" : "build non cst part)";
    this->log("Stokes", "updateLinearPDE", "start"+scst);

    auto ps = this->spaceProduct();
    auto bbf = blockform2( ps, A );
    auto t = this->fieldTensor();
    auto u = this->fieldFlux();
    auto p = this->fieldPotential();
    auto uhat = this->fieldTrace();

    auto tau_constant = cst(M_tauCst);
    auto sc_param = M_useSC ? 0.5 : 1.0;

    // (delta,gamma)_Omega
    bbf( 0_c, 0_c ) +=  integrate(_range=elements(support(M_Wh)),
                                  _expr=inner(idt(t),id(t)) );
    // (u,div(delta))_Omega
    bbf( 0_c, 1_c ) += integrate(_range=elements(support(M_Wh)),
                               _expr=(trans(idt(u))*div(t)) );

    // <uhat,sn>_Gamma\Gamma_I
    bbf( 0_c, 3_c) += integrate(_range=internalfaces(support(M_Wh)),
                                _expr=-(trans(idt(uhat))*leftface((id(t)*N()))+
                                        trans(idt(uhat))*rightface((id(t)*N())) ));
    bbf( 0_c, 3_c) += integrate(_range=boundaryfaces(support(M_Wh)),
                                _expr=-trans(idt(uhat))*(id(t)*N()) );


    // (mu v, tn)_Omega
    bbf( 1_c, 0_c) += integrate(_range=elements(support(M_Wh)),
                                _expr=2*mu*inner(idt(t),grad(u)));
    bbf( 1_c, 0_c) += integrate(_range=internalfaces(support(M_Wh)),
                                _expr=-2*mu*(trans(leftface(id(u)))*leftfacet((idt(t)*N()))+
                                             trans(rightface(id(u)))*rightfacet((idt(t)*N())) ));
    bbf( 1_c, 0_c) += integrate(_range=boundaryfaces(support(M_Wh)),
                                _expr=-2*mu*trans(id(u))*(idt(t)*N()) ); 


    // <mu tau u, v>_Gamma
    bff( 1_c, 1_c) += integrate(_range=boundaryfaces(support(M_Wh)),
                              _expr=mu*tau_constant*trans(id(v))*idt(u) );
    bff( 1_c, 1_c) += integrate(_range=internalfaces(support(M_Wh)),
                              _expr=mu*tau_constant*(trans(leftface(id(v)))*leftfacet(idt(u))+
                                                     trans(rightface(id(v)))*rightfacet(idt(u))) );

#if 0
    if( !this->isStationary() )
    {
        // (1/delta_t p, w)_Omega  [only if it is not stationary]
        auto coeff = this->timeStepBdfPotential()->polyDerivCoefficient(0);
        bbf( 1_c, 1_c) += integrate(_range=elements(support(M_Wh)),
                                    _expr=coeff*inner(idt(p), id(p)) );
    }
#endif
    // <-p, v>_Omega
    bbf( 1_c, 2_c) += integrate(_range=elements(support(M_Wh)),
                                _expr=(-1)*idt(p)*div(v) );
    bbf( 1_c, 2_c) += integrate(_range=internalfaces(support(M_Wh)),
                                _expr=inner(jumpt(idt(p)),id(v)) );
    bbf( 1_c, 2_c) += integrate(_range=boundaryfaces(support(M_Wh)),
                                _expr=normal(v)*idt(p) );

    // <-mu tau uhat, v>_Omega/Gamma
    bbf( 1_c, 3_c) += integrate(_range=internalfaces(support(M_Wh)),
                                _expr=-mu*tau_constant*(trans(leftface(id(v)))*idt(uhat)+
                                                        trans(rightface(id(v)))*idt(uhat)) );
    bbf( 1_c, 3_c) += integrate(_range=boundaryfaces(support(M_Wh)),
                                _expr=-mu*tau_constant*(trans(id(v))*idt(uhat)) );

    // <-u, grad(q)>_Omega
    bbf( 2_c, 1_c) += integrate(_range=elements(support(M_Wh)),
                                _expr=(-1)*grad(q)*idt(u) );
    
    // <uhat, qn>_Gamma
    bbf( 2_c, 3_c) += integrate(_range=internalfaces(support(M_Wh)),
                                _expr=inner(jump(id(q)),idt(uhat)) );
    bbf( 2_c, 3_c) += integrate(_range=boundaryfaces(support(M_Wh)),
                                _expr=id(q)*(trans(idt(uhat))*N()) );
#if 0
    if ( bc_only_dirichlet )
    {
        a( 2_c, 4_c) += integrate(_range=elements(mesh),
                                  _expr=idt(qm)*id(q) );       
        a( 4_c, 2_c) += integrate(_range=elements(mesh),
                                  _expr=id(qm)*idt(q) ); 
    }
    else
    {
        a( 4_c, 4_c) += integrate(_range=elements(mesh),
                              _expr=id(qm)*idt(qm) );
    }
#endif

    // (-2 mu tn, m)_Gamma
    bbf( 3_c, 0_c) += integrate(_range=internalfaces(support(M_Wh)),
                                _expr=-2*mu*(trans(id(m))*leftfacet(idt(delta)*N())+
                                           trans(id(m))*rightfacet(idt(delta)*N())) );
    bbf( 3_c, 0_c ) += integrate(_range = markedfaces( support(M_Wh), "Neumann" ),
                                 _expr = -2*mu * trans(id(m))*(idt(delta)*N() ));
 
    // (mu tau u, m)_Gamma
    bbf( 3_c, 1_c) += integrate(_range=internalfaces(support(M_Wh)),
                                _expr=mu*tau_constant*(trans(leftfacet(idt(u)))*id(m)+
                                                   trans(rightfacet(idt(u)))*id(m)) );
    bbf( 3_c, 1_c ) += integrate( _range = markedfaces( support(M_Wh), "Neumann" ),
                                  _expr = mu * tau_constant * ( trans( idt( u ) ) * id( m ) ) );

    // (pn, m)_Gamma
    bbf( 3_c, 2_c) += integrate(_range=internalfaces(support(M_Wh)),
                                _expr=inner(id(m),jumpt(idt(p))) ); //normalt(m)*(rightface(id(p))-leftface(id(p))) );
    bbf( 3_c, 2_c) += integrate(_range=markedfaces(support(M_Wh),"Neumann"),
                                _expr=idt(p)*trans(id(m))*N() );

    // (-mu tau uhat, m)_Gamma
    bbf( 3_c, 3_c) += integrate(_range=internalfaces(support(M_Wh)),                  
                                _expr=-sc_param*mu*tau_constant*trans(idt(uhat))*id(m) );
    bbf( 3_c, 3_c) += integrate(_range=markedfaces(support(M_Wh),"Dirichlet"),                  
                                _expr=trans(idt(uhat))*id(m) );
    bbf( 3_c, 3_c ) += integrate(_range = markedfaces( support(M_Wh), "Neumann" ),
                                 _expr = -mu*tau_constant*trans( idt( uhat ) ) * id( m ) );

    this->updateLinearPDE( data, this->modelContext() );
}

STOKES_CLASS_TEMPLATE_DECLARATIONS
void
STOKES_CLASS_TEMPLATE_TYPE::updatePostPDE( DataUpdateHDG & data ) const
{
    this->updatePostPDE( data, this->modelContext() );
}

}
}

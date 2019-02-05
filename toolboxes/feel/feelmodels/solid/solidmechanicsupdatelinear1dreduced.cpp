#include <feel/feelmodels/solid/solidmechanics.hpp>


namespace Feel
{
namespace FeelModels
{

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateLinearGeneralisedStringGeneralisedAlpha( DataUpdateLinear & data ) const
{
#if (SOLIDMECHANICS_DIM==2)
    this->log( "SolidMechanics","updateLinearGeneralisedStringGeneralisedAlpha","start");

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool _buildCstPart = data.buildCstPart();
    bool _doBCStrongDirichlet = data.doBCStrongDirichlet();


    bool BuildNonCstPart = !_buildCstPart;
    bool BuildCstPart = _buildCstPart;

    //auto mesh = M_mesh_1dReduced;
    auto M_rangeMeshElements1dReduced = elements(M_mesh_1dReduced);
    auto Xh1 = M_Xh_1dReduced;
    auto u = Xh1->element(), v = Xh1->element();

    double epp=this->thickness1dReduced();//0.1;
    double k=2.5;
    double G=1e5;
    double gammav=0.01;

    double E=this->mechanicalProperties()->cstYoungModulus();//M_youngmodulus;
    double mu=this->mechanicalProperties()->cstCoeffPoisson();//M_coeffpoisson;
    double R0=this->radius1dReduced();//0.5;
    double rho=this->mechanicalProperties()->cstRho();
    bool robinOutletCoef = 1./math::sqrt(k*G/rho);
    //---------------------------------------------------------------------------------------//

    //double rho_s=0.8;
    //double alpha_f=M_genAlpha_alpha_f;//1./(1.+rho_s);
    //double alpha_m=M_genAlpha_alpha_m;//(2.-rho_s)/(1.+rho_s);
    //double alpha_vel = M_genAlpha_alpha_f;//0.5*(3-rho_s)/(1+rho_s);

    //double gamma=0.5+alpha_m-alpha_f;
    //double beta=0.25*(1+alpha_m-alpha_f)*(1+alpha_m-alpha_f);

    auto deltaT =  M_newmark_displ_1dReduced->timeStep();
    auto const& buzz1 = M_newmark_displ_1dReduced->previousUnknown();

    //---------------------------------------------------------------------------------------//

    auto rowStartInMatrix = this->rowStartInMatrix();
    auto colStartInMatrix = this->colStartInMatrix();
    auto rowStartInVector = this->rowStartInVector();
    auto bilinearForm1dreduced = form2( _test=Xh1,_trial=Xh1,_matrix=A,
                                        _rowstart=rowStartInMatrix,
                                        _colstart=colStartInMatrix );
    auto linearForm1dreduced = form1( _test=Xh1, _vector=F,
                                      _rowstart=rowStartInVector );

    //---------------------------------------------------------------------------------------//
    // acceleration term
    if (!this->isStationary())
    {
        if (BuildCstPart)
        {
            bilinearForm1dreduced +=
                integrate( _range=M_rangeMeshElements1dReduced,
                           _expr=this->timeStepNewmark1dReduced()->polySecondDerivCoefficient()*rho*epp*idt(u)*id(v) );
        }

        if (BuildNonCstPart)
        {
            linearForm1dreduced +=
                integrate( _range=M_rangeMeshElements1dReduced,
                           _expr=rho*epp*idv(M_newmark_displ_1dReduced->polyDeriv() )*id(v) );
        }
    }
    //---------------------------------------------------------------------------------------//
    // reaction term
    if (BuildCstPart)
    {
        bilinearForm1dreduced +=
            integrate( _range=M_rangeMeshElements1dReduced,
                       _expr=(E*epp/((1-mu*mu)*R0*R0))*idt(u)*id(v) );
    }
    //---------------------------------------------------------------------------------------//
    // diffusion term
    if (BuildCstPart)
    {
        bilinearForm1dreduced +=
            integrate( _range=M_rangeMeshElements1dReduced,
                       _expr=k*G*epp*dxt(u)*dx(v) );
    }

    //---------------------------------------------------------------------------------------//
    // viscoealstic term
    if (BuildCstPart)
    {
        bilinearForm1dreduced +=
            integrate( _range=M_rangeMeshElements1dReduced,
                       _expr=this->timeStepNewmark1dReduced()->polyFirstDerivCoefficient()*gammav*dxt(u)*dx(v) );
    }
    if (BuildNonCstPart)
    {
        linearForm1dreduced +=
            integrate( _range=M_rangeMeshElements1dReduced,
                       _expr=gammav*dxv(this->timeStepNewmark1dReduced()->polyFirstDeriv())*dx(v) );
    }
    //---------------------------------------------------------------------------------------//
    // source term
#if 0
    if (BuildNonCstPart)
    {
        linearForm1dreduced +=
            integrate( _range=M_rangeMeshElements1dReduced,
                       _expr=idv(*M_stress_1dReduced)*id(v) );
    }
#endif
    //---------------------------------------------------------------------------------------//
    // dirichlet bc
    if (this->hasMarkerDirichletBCelimination() && BuildNonCstPart && _doBCStrongDirichlet)
    {
        this->updateBCDirichletStrongLinearPDE( A,F );
    }

    //---------------------------------------------------------------------------------------//

    this->log( "SolidMechanics","updateLinearGeneralisedStringGeneralisedAlpha","finish");
#endif
}

//-------------------------------------------------------------------------------------//



} // FeelModels

} // Feel


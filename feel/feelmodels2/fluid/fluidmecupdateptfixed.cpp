
#include "fluidmec.hpp"

#include <fsi/fsialg/functionSup.cpp>

#define STAB_FIC 0


namespace Feel
{
namespace FSI
{

    void
    FLUIDMECHANICS_CLASS_NAME::updatePtFixe( const vector_ptrtype& Xold, sparse_matrix_ptrtype& A , vector_ptrtype& F,
                                             bool BuildCstPart, bool _doClose, bool _doBCStrongDirichlet) const
    {
#if defined(FSI_FLUID_BUILD_PTFIXED_CODE)
        using namespace Feel::vf;

        if (this->verbose()) std::cout << "[FluidMechanics] updatePtFixe start\n";

        auto mesh = this->mesh();
        auto Xh = this->functionSpace();

        element_fluid_type U( Xh, "u" );
        element_fluid_type V( Xh, "v" );
        element_fluid_0_type u = U_fluid->element<0>();
        element_fluid_0_type v = V.element<0>();
        element_fluid_1_type p = U_fluid->element<1>();
        element_fluid_1_type q = V.element<1>();
#if STAB_FIC
        element_fluid_2_type c = U.element<2>();
        element_fluid_3_type pp = U.element<3>();
#endif

        element_fluid_type Uold( Xh, "uold" );Uold = *Xold;
        element_fluid_0_type uold = Uold.element<0>();
        element_fluid_1_type pold = Uold.element<1>();

        //--------------------------------------------------------------------------------------------------//

        //Deformations tensor (trial)
        auto deft = sym(gradt(u));

        //Identity Matrix
#if (FLUIDMECHANICS_DIM==2)
        auto Id = mat<2,2>(cst(1.),cst(0.),
                               cst(0.),cst(1.) );
#else
        auto Id = mat<3,3>(cst(1.),cst(0.),cst(0.),
                           cst(0.),cst(1.),cst(0.),
                           cst(0.),cst(0.),cst(1.) );
#endif

        // Strain tensor (trial)
        auto Sigmat = -idt(p)*Id + 2*idv(*M_P0Mu)*deft;
        //volume force
#if defined(FLUIDMECHANICS_VOLUME_FORCE)
        auto f = FLUIDMECHANICS_VOLUME_FORCE(this->shared_from_this()) ;
#endif
        //boundaries conditions
        auto bcDef = FLUIDMECHANICS_BC(this->shared_from_this());

        //--------------------------------------------------------------------------------------------------//

        form2( Xh,Xh, A, _init=true );

        //convective term
        if (M_isMoveDomain)
            form2( Xh,Xh, A )  +=
                integrate (_range=elements(mesh),
                           _expr=idv(*M_P0Rho)*trans(gradt(u)*val(idv(uold) - idv(this->meshVelocity())))*id(v) );
        else
            form2( Xh,Xh, A )  +=
                integrate (_range=elements(mesh),
                           _expr=idv(*M_P0Rho)*trans(gradt(u)*idv(uold))*id(v));

        // sigma : grad(v) on Omega
        form2( Xh, Xh, A ) +=
            integrate( _range=elements(mesh),
                       _expr=trace(Sigmat*trans(grad(v))) );

        // incompressibility term
        form2( Xh, Xh, A ) +=
            integrate( _range=elements(mesh),
                       _expr= divt(u)*id(q) );

        //add this term if we use lagrange multiplier
        this->addLagrangeMult(mpl::bool_<useLagMult>(),A);

        /*
        //terme en plus pour div=0 (voir these goncalo)
        form2( Xh, Xh, A ) +=
        integrate( elements( mesh ),// _Q<3*uOrder-1+3*(GeoOrder-1)>(),
        divv(uold)*trans(idt(u))*id(v) );
        */

        // volume force
#if defined(FLUIDMECHANICS_VOLUME_FORCE)
        form1( Xh, F, _init=true ) =
            integrate( _range=elements(mesh),
                       _expr=trans(f)*id(v) );
#endif
        if (M_haveSourceAdded)
            form1( Xh, F ) +=
                integrate( _range=elements(mesh),
                           _expr=trans(idv(*M_SourceAdded))*id(v) );

        //neumann condition
        ForEachBC( bcDef,cl::neumann_scal,
                   form1( Xh, F ) +=
                   integrate( _range=markedfaces(mesh,PhysicalName),
                              _expr= trans(Expression*N())*id(v) )  );

        // weak formulation of the boundaries conditions
        if (M_weakCL)
            {
                ForEachBC( bcDef, cl::dirichlet_vec,
                           form2( Xh, Xh, A ) +=
                           integrate( _range=markedfaces(mesh,PhysicalName),
                                      _expr= -trans(Sigmat*N())*id(v)
                                      + penalbc*trans(idt(u))*id(v)/hFace() ) );
                ForEachBC( bcDef, cl::dirichlet_vec,
                           form1( Xh, F ) +=
                           integrate( _range=markedfaces(mesh,PhysicalName),
                                      _expr= penalbc*trans(Expression)*id(v)/hFace() ) );

                ForEachBC( bcDef, cl::paroi_mobile,
                           form2( Xh, Xh, A ) +=
                           integrate( _range=markedfaces(mesh,PhysicalName),
                                      _expr= -trans(Sigmat*N())*id(v)
                                      + penalbc*trans(idt(u))*id(v)/hFace() ) );

                ForEachBC( bcDef, cl::paroi_mobile,
                           form1( Xh, F ) +=
                           integrate( _range=markedfaces(mesh,PhysicalName),
                                      _expr= penalbc*trans(idv(this->meshVelocity2()))*id(v)/hFace() ) );
            }

        //--------------------------------------------------------------------------------------------------//

        auto P = Id-N()*trans(N());
        double gammaN = this->application()->vm()[prefixvm(this->prefix(),"bc-slip-gammaN")].as<double>();
        double gammaTau =  this->application()->vm()[prefixvm(this->prefix(),"bc-slip-gammaTau")].as<double>();
        auto Beta = M_bdf_fluid->poly();
        //auto beta = Beta.element<0>();
        auto beta = vf::project( Beta.element<0>().functionSpace(), boundaryfaces(Beta.element<0>().mesh()), idv(*M_P0Rho)*idv(Beta.element<0>()) );
        auto Cn = gammaN*max(abs(trans(idv(beta))*N()),idv(M_P0Mu)/vf::h());
        auto Ctau = gammaTau*idv(M_P0Mu)/vf::h() + max( -trans(idv(beta))*N(),cst(0.) );

        ForEachBC( bcDef, cl::slip,
                   form2( Xh, Xh, A ) +=
                   integrate( _range=markedfaces(mesh,PhysicalName),
                              _expr=val(Cn)*(trans(idt(u))*N())*(trans(id(v))*N())+
                              val(Ctau)*trans(idt(u))*id(v)
                              //+ trans(idt(p)*Id*N())*id(v)
                              //- trans(id(v))*N()* trans(2*idv(*M_P0Mu)*deft*N())*N()
                       ) );

        //transients terms
        if (!this->isStationary())
            {
                form2( Xh, Xh, A ) +=
                    integrate( _range=elements(mesh),
                               _expr= idv(*M_P0Rho)*trans(idt(u))*id(v)*M_bdf_fluid->polyDerivCoefficient(0) );

                auto Buzz = M_bdf_fluid->polyDeriv();
                auto buzz = Buzz.element<0>();

                form1( Xh, F ) +=
                    integrate( _range=elements(mesh),
                               _expr= idv(*M_P0Rho)*trans(idv(buzz))*id(v) );
            }

        //--------------------------------------------------------------------------------------------------//

        A->close();
        F->close();

        // strong formulation of the boundaries conditions
        if (!M_weakCL)
            {
                ForEachBC( bcDef, cl::paroi_mobile,
                           form2( Xh, Xh, A ) +=
                           on( _range=markedfaces(mesh, PhysicalName),
                               _element=u,
                               _rhs=F,
                               _expr=idv(this->meshVelocity2()) ) );
                ForEachBC( bcDef, cl::dirichlet_vec,
                           form2( _test=Xh, _trial=Xh, _matrix=A ) +=
                           on( _range=markedfaces(mesh, PhysicalName),
                               _element=u,
                               _rhs=F,
                               _expr=Expression ) );
            }

        //--------------------------------------------------------------------------------------------------//

        if (this->verbose()) std::cout << "[FluidMechanics] updatePtFixe finish\n";
#endif
    } // updatePtFixe

} // end namespace FSI
} // end namespace Feel

//template class Feel::FSI::FluidMechanics<FLUIDMECHANICS_DIM,FLUIDMECHANICS_ORDER_VELOCITY,Feel::Simplex>;


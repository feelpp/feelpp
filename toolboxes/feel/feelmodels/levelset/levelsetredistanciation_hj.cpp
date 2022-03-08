//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file levelsetredistanciation_hj.cpp
//! @author Thibaut Metivet <thibaut.metivet@inria.fr>
//! @date 15 Jan 2020
//! @copyright 2020 Feel++ Consortium
//!

#include <feel/feelmodels/levelset/levelsetredistanciation_hj.hpp>

namespace Feel
{

#define LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_DECLARATIONS \
   template<typename FunctionSpaceType>
#define LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_TYPE \
    LevelSetRedistanciationHJ<FunctionSpaceType>

LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_DECLARATIONS
LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_TYPE::LevelSetRedistanciationHJ( 
        functionspace_ptrtype const& space,
        std::string const& prefix )
: 
    super_type( space, prefix ),
    M_advectionHJ( new advectionhj_type( prefix, prefix, space->worldCommPtr() ) ),
    M_nGlobalIter(1)
{
    this->loadParametersFromOptionsVm();

    //M_advectionHJ = advectionhj_type::New( prefix );
    M_advectionHJ->setModelName( "Advection" );
    M_advectionHJ->setTimeInitial(0.);
    M_advectionHJ->setTimeFinal( M_timeStep * M_maxIterations );
    M_advectionHJ->setTimeStep( M_timeStep );
    M_advectionHJ->setFunctionSpace( space );
    M_advectionHJ->init();

    M_functionSpaceP1Vec = functionspace_P1v_type::New( space->mesh() );
    M_projectorL2Vec = projector(M_functionSpaceP1Vec, M_functionSpaceP1Vec, 
            backend(_name=prefixvm("levelset", "projector-l2-vectorial")) );
    
    //if( this->M_useVolumeConstraint )
        M_functionSpaceP0 = functionspace_P0_type::New( space->mesh() );
}

LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_TYPE::loadParametersFromOptionsVm()
{
    M_tolerance = doption( _name="tol", _prefix=this->prefix() );
    M_timeStep = doption( _name="time-step", _prefix=this->prefix() );
    M_maxIterations = ioption( _name="max-iter", _prefix=this->prefix() );
    M_thicknessHeaviside = doption( _name="thickness-heaviside", _prefix=this->prefix() );
    M_useVolumeConstraint = boption( _name="keep-volume", _prefix=this->prefix() );
}

LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_TYPE::element_type
LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_TYPE::run( element_type const& phi ) const
{
    auto mesh = M_advectionHJ->mesh();
    auto space = M_advectionHJ->functionSpace();
    // Set initial value
    M_advectionHJ->fieldSolution() = phi;
    //auto gradPhi = M_projectorL2Vec->project( _expr=trans(gradv(phi)) );
    ////auto modGradPhi0 = vf::project(
            ////_space=this->functionSpaceP0(),
            ////_range=elements(mesh),
            ////_expr=sqrt( gradv(phi) * trans(gradv(phi)) )
            ////);
    //auto modGradPhi = vf::project(
            //_space=this->functionSpace(),
            //_range=elements(mesh),
            //_expr=sqrt( trans(idv(gradPhi)) * idv(gradPhi) )
            //);
    //M_advectionHJ->fieldSolutionPtr()->on(
            //_range=elements(mesh),
            //_expr=idv(phi)/idv(modGradPhi)
            //);

    //M_advectionHJ->timeStepBase()->setTimeFinal( M_advectionHJ->timeStep() * M_maxIterations );
    M_advectionHJ->initTimeStep();
    // Init convergence test utilities
    double err_dist;
    double err_disto = integrate(
        _range=elements(this->mesh()),
        _expr=vf::abs( sqrt(gradv(phi)*trans(gradv(phi))) - 1.0 )
            ).evaluate()(0,0);

    double rateChangePhiL2;
    double rateChangePhiL2o = 1.;
    double relativeRateChangePhiL2;

    //for(int i = 0; i < this->maxIterations(); ++i)
    for( int i = 0; !M_advectionHJ->timeStepBase()->isFinished(); M_advectionHJ->updateTimeStep(), ++i )
    {
        LOG(INFO) << "iter redist : " << i << std::endl;

        auto const& phi_redist = M_advectionHJ->fieldSolutionPtr();
        auto const& phi_redisto = M_advectionHJ->timeStepBDF()->unknowns()[1];

        auto gradPhiReinit = M_projectorL2Vec->project( _expr=trans(gradv(phi_redist)) );
#if 1
        //auto phi_sign = idv(phi_redist) / vf::max( vf::sqrt(gradv(phi_redist)*trans(gradv(phi_redist))), 0.92 );
        auto phi_sign = idv(phi_redist) / vf::sqrt( trans(idv(gradPhiReinit))*idv(gradPhiReinit) );
        auto H_redist = vf::project(
                _space=space,
                _range=elements(mesh),
                _expr=vf::abs(
                    ( phi_sign < -M_thicknessHeaviside )*vf::constant(0.0)
                    +
                    ( phi_sign >= -M_thicknessHeaviside && phi_sign <= M_thicknessHeaviside )*
                    0.5*(1 + phi_sign/M_thicknessHeaviside + 1/M_PI*vf::sin( M_PI*phi_sign/M_thicknessHeaviside ) )
                    +
                    ( phi_sign > M_thicknessHeaviside )*vf::constant(1.0) )
                );
        auto Sign = 2*( idv(H_redist)-0.5 );
#else
        auto gradPhi = M_projectorL2Vec->project( _expr=trans(gradv(phi)) );
        //auto psi = idv(phi) / vf::sqrt(gradv(phi)*trans(gradv(phi)));
        auto psi = idv(phi) / vf::sqrt( trans(idv(gradPhi))*idv(gradPhi) );
        auto H_expr = vf::chi( psi<-M_thicknessHeaviside )*vf::constant(0.0)
            +
            vf::chi( psi>=-M_thicknessHeaviside )*vf::chi( psi<=M_thicknessHeaviside )*
            0.5*(1 + psi/M_thicknessHeaviside + 1/M_PI*vf::sin( M_PI*psi/M_thicknessHeaviside ) )
            +
            vf::chi(psi>M_thicknessHeaviside)*vf::constant(1.0);
        auto Sign = 2. * ( H_expr - 0.5 );
        //double h = mesh->hAverage();
        //auto Sign = psi / vf::sqrt( psi*psi + h*h );
#endif

        // Update advection source
        M_advectionHJ->updateSourceAdded( Sign );
        // Update advection velocity
        //auto beta = Sign * trans(gradv(phi_redist)) / vf::max( vf::sqrt(gradv(phi_redist)*trans(gradv(phi_redist))), 0.92 );
        if( this->M_useVolumeConstraint )
        {
            auto it_elt = mesh->beginOrderedElement();
            auto en_elt = mesh->endOrderedElement();

            const rank_type pid = mesh->worldCommElements().localRank();
            const int ndofv = functionspace_type::fe_type::nDof;

            double thickness = M_thicknessHeaviside;
            typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_type elements_reference_wrapper_type;
            typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype elements_reference_wrapper_ptrtype;
            elements_reference_wrapper_ptrtype interfaceElts( new elements_reference_wrapper_type );

            for (; it_elt!=en_elt; it_elt++)
            {
                auto const& elt = boost::unwrap_ref( *it_elt );
                if ( elt.processId() != pid )
                    continue;
                bool mark_elt = true;
                for (int j=0; j<ndofv; j++)
                {
                    if ( std::abs( phi_redist->localToGlobal(elt.id(), j, 0) ) > thickness )
                    {
                        mark_elt = false;
                        break; //don't need to do the others dof
                    }
                }
                if( mark_elt )
                    interfaceElts->push_back( boost::cref(elt) );
            }

            elements_reference_wrapper_t<mesh_type> interfaceElements = 
                boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                        interfaceElts->begin(),
                        interfaceElts->end(),
                        interfaceElts
                        );

            //auto modGradPhiReinit = vf::sqrt(gradv(phi_redist)*trans(gradv(phi_redist)));
            auto modGradPhiReinit = vf::sqrt( trans(idv(gradPhiReinit))*idv(gradPhiReinit) );
            auto Delta = Feel::FeelModels::levelsetDelta( _element=*phi_redist, _thickness=M_thicknessHeaviside );
            auto spaceP0 = this->functionSpaceP0();
            auto LambdaNum = integrate(
                    _range=elements(mesh),
                    _expr=Delta * Sign * ( 1. - modGradPhiReinit )
                    ).broken( spaceP0 );
            auto LambdaDen = integrate(
                    _range=elements(mesh),
                    _expr=Delta * Delta * modGradPhiReinit
                    ).broken( spaceP0 );
            auto Lambda = this->functionSpace()->element();
            Lambda.on( _range=interfaceElements, _expr=idv(LambdaNum) / idv(LambdaDen) );
            //auto Lambda = idv(LambdaNum) / idv(LambdaDen);
            auto beta = ( Sign + idv(Lambda)*Delta ) * idv(gradPhiReinit) / modGradPhiReinit;
            M_advectionHJ->updateAdvectionVelocity( beta );
        }
        else
        {
            //auto beta = Sign * trans(gradv(phi_redist)) / vf::sqrt(gradv(phi_redist)*trans(gradv(phi_redist)));
            auto beta = Sign * idv(gradPhiReinit) / vf::sqrt( trans(idv(gradPhiReinit))*idv(gradPhiReinit) );
            M_advectionHJ->updateAdvectionVelocity( beta );
        }

        // Solve
        M_advectionHJ->solve();
        M_advectionHJ->exportResults(M_nGlobalIter);

        // Test convergence
        err_dist = integrate(
                _range=elements(this->mesh()),
                _expr=vf::abs( sqrt(gradv(phi_redist)*trans(gradv(phi_redist))) - 1.0 )
                ).evaluate()(0,0);
        double delta_err = std::abs(err_disto - err_dist) / err_disto;
        LOG(INFO)<<"delta err = "<<delta_err<<"\n";

        rateChangePhiL2 = std::sqrt( integrate(
                                         _range=elements(mesh),
                                         _expr=(idv(phi_redist) - idv(phi_redisto))*(idv(phi_redist) - idv(phi_redisto))
                ).evaluate()(0,0) );
        relativeRateChangePhiL2 = std::abs( rateChangePhiL2 - rateChangePhiL2o) / rateChangePhiL2o;

        ++M_nGlobalIter;

        if ( (i!=0) && (relativeRateChangePhiL2 <= M_tolerance ) )
        {
            LOG(INFO) << "stop Hamilton Jacobi iterations because tolerance threshold has been reached\n";
            LOG(INFO) << "relative rate of change = " << relativeRateChangePhiL2 << std::endl;
            break;
        }
        
        // Update "old" variables
        err_disto = err_dist;
        rateChangePhiL2o = rateChangePhiL2;
    }

    //Feel::cout << "redist done in " << __iter - start_iter << " iter\n";
    Feel::cout << "final relative rate of change = " << relativeRateChangePhiL2 << std::endl;

    return M_advectionHJ->fieldSolution();
}

//--------------------------------------------------------------------//
//--------------------------------------------------------------------//
//--------------------------------------------------------------------//
// Hamilton-Jacobi advection
#define ADVECTIONHJ_CLASS_TEMPLATE_DECLARATIONS \
    template<typename SpaceType>
#define ADVECTIONHJ_CLASS_TEMPLATE_TYPE \
    AdvectionHJ<SpaceType>

LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_DECLARATIONS
ADVECTIONHJ_CLASS_TEMPLATE_DECLARATIONS
LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_TYPE::ADVECTIONHJ_CLASS_TEMPLATE_TYPE::AdvectionHJ(
        std::string const& prefix )
    : super_type( prefix )
{
    this->setModelName( "Advection" );
}

LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_DECLARATIONS
ADVECTIONHJ_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_TYPE::template ADVECTIONHJ_CLASS_TEMPLATE_TYPE::self_ptrtype
LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_TYPE::ADVECTIONHJ_CLASS_TEMPLATE_TYPE::New(
        std::string const& prefix )
{
    return std::make_shared<self_type>( prefix );
}

LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_DECLARATIONS
ADVECTIONHJ_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_TYPE::ADVECTIONHJ_CLASS_TEMPLATE_TYPE::init(
        functionspace_ptrtype const& space,
        bool buildModelAlgebraicFactory )
{
    this->setFunctionSpace( space );
    //super_type::init( buildModelAlgebraicFactory, this->shared_from_this() );
    super_type::init( buildModelAlgebraicFactory );
}

LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_DECLARATIONS
ADVECTIONHJ_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_TYPE::ADVECTIONHJ_CLASS_TEMPLATE_TYPE::updateWeakBCLinearPDE(
        sparse_matrix_ptrtype& A, 
        vector_ptrtype& F,
        bool buildCstPart) const
{
}

LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_DECLARATIONS
ADVECTIONHJ_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_TYPE::ADVECTIONHJ_CLASS_TEMPLATE_TYPE::updateBCStrongDirichletLinearPDE(
        sparse_matrix_ptrtype& A, 
        vector_ptrtype& F) const
{}

#undef ADVECTIONHJ_CLASS_TEMPLATE_DECLARATIONS
#undef ADVECTIONHJ_CLASS_TEMPLATE_TYPE

#undef LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_DECLARATIONS
#undef LEVELSETREDISTANCIATIONHJ_CLASS_TEMPLATE_TYPE

} // namespace Feel

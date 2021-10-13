/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Thibaut Metivet <thibaut.metivet@univ-grenoble-alpes.fr>
 Date: 2016-05-20

 Copyright (C) 2016 Université Joseph Fourier (Grenoble I)

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 3.0 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
/**
 \file reinitializer_hj.hpp
 \author Thibaut Metivet <thibaut.metivet@univ-grenoble-alpes.fr>
 \date 2016-05-20
 */
#ifndef _REINITIALIZER_HJ_HPP
#define _REINITIALIZER_HJ_HPP 1

#include <feel/feelmodels/levelset/reinitializer.hpp>
#include <feel/feelmodels/advection/advection.hpp>
#include <feel/feelmodels/levelset/levelsetdeltaexpr.hpp>

namespace Feel
{

template<typename FunctionSpaceType>
class ReinitializerHJ
: public Reinitializer<FunctionSpaceType>
{
public:
    //--------------------------------------------------------------------//
    // Typedefs
    typedef Reinitializer<FunctionSpaceType> super_type;
    typedef ReinitializerHJ<FunctionSpaceType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    typedef FunctionSpaceType functionspace_type;
    typedef std::shared_ptr<functionspace_type> functionspace_ptrtype;

    typedef typename functionspace_type::mesh_type mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    typedef typename functionspace_type::element_type element_type;
    typedef std::shared_ptr<element_type> element_ptrtype;

    typedef typename functionspace_type::periodicity_0_type periodicity_type;
    static const bool is_periodic = functionspace_type::is_periodic;

    typedef FunctionSpace<mesh_type, bases<Lagrange<1,Vectorial,Continuous>>> functionspace_P1v_type;
    typedef std::shared_ptr<functionspace_P1v_type> functionspace_P1v_ptrtype;

    typedef Projector<functionspace_P1v_type, functionspace_P1v_type> projectorL2_vectorial_type;
    typedef std::shared_ptr<projectorL2_vectorial_type> projectorL2_vectorial_ptrtype;

    typedef FunctionSpace<mesh_type, bases<Lagrange<0,Scalar,Discontinuous> > > functionspace_P0_type;
    typedef std::shared_ptr<functionspace_P0_type> functionspace_P0_ptrtype;

    //--------------------------------------------------------------------//
    // Hamilton-Jacobi advection
    template<typename SpaceType>
    class AdvectionHJ
        : public Feel::FeelModels::AdvDiffReac<SpaceType/*,TODO*/>
        , public std::enable_shared_from_this< AdvectionHJ<SpaceType> >
    {
    public:
        typedef Feel::FeelModels::AdvDiffReac<SpaceType> super_type;

        typedef AdvectionHJ<SpaceType> self_type;
        typedef std::shared_ptr<self_type> self_ptrtype;

        typedef SpaceType functionspace_type;
        typedef std::shared_ptr<functionspace_type> functionspace_ptrtype;

        // Constructor
        AdvectionHJ(
                std::string const& prefix = ""
                );
        
        static self_ptrtype New(
                std::string const& prefix = ""
                );

        // Initialization
        void init( functionspace_ptrtype const& space, bool buildModelAlgebraicFactory = true );

        // BC and source term assembly
        void updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const;
        void updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const;
        bool hasSourceTerm() const { return false; }
    };

    //typedef AdvectionHJ<FunctionSpaceType> advectionhj_type;
    typedef Feel::FeelModels::AdvDiffReac<FunctionSpaceType> advectionhj_type;
    typedef std::shared_ptr<advectionhj_type> advectionhj_ptrtype;

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    // Constructor
    ReinitializerHJ( 
            functionspace_ptrtype const& space,
            std::string const& prefix = "" );
    //--------------------------------------------------------------------//
    // Parameters
    void loadParametersFromOptionsVm();

    double tolerance() const { return M_tolerance; }
    void setTolerance( double tol ) { M_tolerance = tol; }
    
    double timeStep() const { return M_advectionHJ->timeStep(); }
    void setTimeStep( double dt ) { M_advectionHJ->setTimeStep( dt ); }

    int maxIterations() const { return M_maxIterations; }
    void setMaxIterations( int max ) { M_maxIterations = max; }

    double thicknessHeaviside() const { return M_thicknessHeaviside; }
    void setThicknessHeaviside( double eps ) { M_thicknessHeaviside = eps; }

    functionspace_P0_ptrtype functionSpaceP0() const { return M_functionSpaceP0; }
    //--------------------------------------------------------------------//
    // Run reinitialization
    element_type run( element_type const& phi );

private:
    advectionhj_ptrtype M_advectionHJ;

    double M_tolerance;
    double M_timeStep;
    int M_maxIterations;
    double M_thicknessHeaviside;

    int M_nGlobalIter;

    functionspace_P1v_ptrtype M_functionSpaceP1Vec;
    projectorL2_vectorial_ptrtype M_projectorL2Vec;

    bool M_useVolumeConstraint;
    functionspace_P0_ptrtype M_functionSpaceP0;
};

#define REINITIALIZERHJ_CLASS_TEMPLATE_DECLARATIONS \
   template<typename FunctionSpaceType>
#define REINITIALIZERHJ_CLASS_TEMPLATE_TYPE \
    ReinitializerHJ<FunctionSpaceType>

REINITIALIZERHJ_CLASS_TEMPLATE_DECLARATIONS
REINITIALIZERHJ_CLASS_TEMPLATE_TYPE::ReinitializerHJ( 
        functionspace_ptrtype const& space,
        std::string const& prefix )
: 
    super_type( space, prefix ),
    M_advectionHJ( new advectionhj_type( prefix, space->worldCommPtr() ) ),
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

REINITIALIZERHJ_CLASS_TEMPLATE_DECLARATIONS
void
REINITIALIZERHJ_CLASS_TEMPLATE_TYPE::loadParametersFromOptionsVm()
{
    M_tolerance = doption( _name="tol", _prefix=this->prefix() );
    M_timeStep = doption( _name="time-step", _prefix=this->prefix() );
    M_maxIterations = ioption( _name="max-iter", _prefix=this->prefix() );
    M_thicknessHeaviside = doption( _name="thickness-heaviside", _prefix=this->prefix() );
    M_useVolumeConstraint = boption( _name="keep-volume", _prefix=this->prefix() );
}

REINITIALIZERHJ_CLASS_TEMPLATE_DECLARATIONS
typename REINITIALIZERHJ_CLASS_TEMPLATE_TYPE::element_type
REINITIALIZERHJ_CLASS_TEMPLATE_TYPE::run( element_type const& phi )
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
            elements(this->mesh()),
            vf::abs( sqrt(gradv(phi)*trans(gradv(phi))) - 1.0 )
            ).evaluate()(0,0);

    double rateChangePhiL2;
    double rateChangePhiL2o = 1.;
    double relativeRateChangePhiL2;

    //for(int i = 0; i < this->maxIterations(); ++i)
    for( int i = 0; !M_advectionHJ->timeStepBase()->isFinished(); M_advectionHJ->updateTimeStep(), ++i )
    {
        LOG(INFO) << "iter reinit : " << i << std::endl;

        auto const& phi_reinit = M_advectionHJ->fieldSolutionPtr();
        auto const& phi_reinito = M_advectionHJ->timeStepBDF()->unknowns()[1];

        auto gradPhiReinit = M_projectorL2Vec->project( _expr=trans(gradv(phi_reinit)) );
#if 1
        //auto phi_sign = idv(phi_reinit) / vf::max( vf::sqrt(gradv(phi_reinit)*trans(gradv(phi_reinit))), 0.92 );
        auto phi_sign = idv(phi_reinit) / vf::sqrt( trans(idv(gradPhiReinit))*idv(gradPhiReinit) );
        auto H_reinit = vf::project(
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
        auto Sign = 2*( idv(H_reinit)-0.5 );
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
        //auto beta = Sign * trans(gradv(phi_reinit)) / vf::max( vf::sqrt(gradv(phi_reinit)*trans(gradv(phi_reinit))), 0.92 );
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
                    if ( std::abs( phi_reinit->localToGlobal(elt.id(), j, 0) ) > thickness )
                    {
                        mark_elt = false;
                        break; //don't need to do the others dof
                    }
                }
                if( mark_elt )
                    interfaceElts->push_back( boost::cref(elt) );
            }

            elements_reference_wrapper_t<mesh_type> interfaceElements = boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                    interfaceElts->begin(),
                    interfaceElts->end(),
                    interfaceElts
                    );

            //auto modGradPhiReinit = vf::sqrt(gradv(phi_reinit)*trans(gradv(phi_reinit)));
            auto modGradPhiReinit = vf::sqrt( trans(idv(gradPhiReinit))*idv(gradPhiReinit) );
            auto Delta = Feel::FeelModels::levelsetDelta( *phi_reinit, M_thicknessHeaviside );
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
            //auto beta = Sign * trans(gradv(phi_reinit)) / vf::sqrt(gradv(phi_reinit)*trans(gradv(phi_reinit)));
            auto beta = Sign * idv(gradPhiReinit) / vf::sqrt( trans(idv(gradPhiReinit))*idv(gradPhiReinit) );
            M_advectionHJ->updateAdvectionVelocity( beta );
        }

        // Solve
        M_advectionHJ->solve();
        M_advectionHJ->exportResults(M_nGlobalIter);

        // Test convergence
        err_dist = integrate(
                elements(this->mesh()),
                vf::abs( sqrt(gradv(phi_reinit)*trans(gradv(phi_reinit))) - 1.0 )
                ).evaluate()(0,0);
        double delta_err = std::abs(err_disto - err_dist) / err_disto;
        LOG(INFO)<<"delta err = "<<delta_err<<"\n";

        rateChangePhiL2 = std::sqrt( integrate(
                elements(mesh),
                (idv(phi_reinit) - idv(phi_reinito))*(idv(phi_reinit) - idv(phi_reinito))
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

    //Feel::cout << "reinit done in " << __iter - start_iter << " iter\n";
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

REINITIALIZERHJ_CLASS_TEMPLATE_DECLARATIONS
ADVECTIONHJ_CLASS_TEMPLATE_DECLARATIONS
REINITIALIZERHJ_CLASS_TEMPLATE_TYPE::ADVECTIONHJ_CLASS_TEMPLATE_TYPE::AdvectionHJ(
        std::string const& prefix )
    : super_type( prefix )
{
    this->setModelName( "Advection" );
}

REINITIALIZERHJ_CLASS_TEMPLATE_DECLARATIONS
ADVECTIONHJ_CLASS_TEMPLATE_DECLARATIONS
typename REINITIALIZERHJ_CLASS_TEMPLATE_TYPE::template ADVECTIONHJ_CLASS_TEMPLATE_TYPE::self_ptrtype
REINITIALIZERHJ_CLASS_TEMPLATE_TYPE::ADVECTIONHJ_CLASS_TEMPLATE_TYPE::New(
        std::string const& prefix )
{
    return std::make_shared<self_type>( prefix );
}

REINITIALIZERHJ_CLASS_TEMPLATE_DECLARATIONS
ADVECTIONHJ_CLASS_TEMPLATE_DECLARATIONS
void
REINITIALIZERHJ_CLASS_TEMPLATE_TYPE::ADVECTIONHJ_CLASS_TEMPLATE_TYPE::init(
        functionspace_ptrtype const& space,
        bool buildModelAlgebraicFactory )
{
    this->setFunctionSpace( space );
    //super_type::init( buildModelAlgebraicFactory, this->shared_from_this() );
    super_type::init( buildModelAlgebraicFactory );
}

REINITIALIZERHJ_CLASS_TEMPLATE_DECLARATIONS
ADVECTIONHJ_CLASS_TEMPLATE_DECLARATIONS
void
REINITIALIZERHJ_CLASS_TEMPLATE_TYPE::ADVECTIONHJ_CLASS_TEMPLATE_TYPE::updateWeakBCLinearPDE(
        sparse_matrix_ptrtype& A, 
        vector_ptrtype& F,
        bool buildCstPart) const
{
}

REINITIALIZERHJ_CLASS_TEMPLATE_DECLARATIONS
ADVECTIONHJ_CLASS_TEMPLATE_DECLARATIONS
void
REINITIALIZERHJ_CLASS_TEMPLATE_TYPE::ADVECTIONHJ_CLASS_TEMPLATE_TYPE::updateBCStrongDirichletLinearPDE(
        sparse_matrix_ptrtype& A, 
        vector_ptrtype& F) const
{}

#undef ADVECTIONHJ_CLASS_TEMPLATE_DECLARATIONS
#undef ADVECTIONHJ_CLASS_TEMPLATE_TYPE

#undef REINITIALIZERHJ_CLASS_TEMPLATE_DECLARATIONS
#undef REINITIALIZERHJ_CLASS_TEMPLATE_TYPE

} // namespace Feel

#endif

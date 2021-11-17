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
//! @file
//! @author Lorenzo Sala <sala@unistra.fr>
//! @date 12 Feb 2019
//! @copyright 2019 Feel++ Consortium
//!


#ifndef FEELPP_POROELASTIC_HPP
#define FEELPP_POROELASTIC_HPP

#include <feel/feelmodels/hdg/mixedpoisson.hpp>
#include <feel/feelmodels/hdg/mixedelasticity.hpp>

namespace Feel
{

namespace FeelModels
{

inline
po::options_description
makeMixedPoissonElasticityOptions( std::string prefix = "hdg.poroelasticity" )
{
    po::options_description mpOptions( "Mixed Poisson Elasticity HDG options");
    mpOptions.add ( mixedpoisson_options("hdg.poisson") );
    mpOptions.add ( makeMixedElasticityOptions("hdg.elasticity") );

    return mpOptions;
}

inline po::options_description
makeMixedPoissonElasticityLibOptions( std::string prefix = "hdg.poroelasticity" )
{
    po::options_description mpLibOptions( "Mixed Poisson Elasticity HDG Lib options");
    // if ( !prefix.empty() )
    //     mpLibOptions.add( backend_options( prefix ) );
    return mpLibOptions;
}


template <int Dim, int Order,int G_Order=1, int E_Order=4>
class MixedPoissonElasticity
{

public:
    typedef double value_type;
    using convex_type = Simplex<Dim, G_Order>;
    using mp_type = MixedPoisson<convex_type, Order>;
    typedef MixedElasticity<Dim,Order,G_Order,E_Order> me_type;
    typedef std::shared_ptr<mp_type> mp_ptrtype;
    typedef std::shared_ptr<me_type> me_ptrtype;
    typedef MixedPoissonElasticity<Dim,Order,G_Order,E_Order> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;
    typedef typename mp_type::mesh_type mesh_type;
    typedef typename mp_type::mesh_ptrtype mesh_ptrtype;

    using sparse_matrix_type = backend_type::sparse_matrix_type;
    using sparse_matrix_ptrtype = backend_type::sparse_matrix_ptrtype;
    using vector_type = backend_type::vector_type;
    using vector_ptrtype = backend_type::vector_ptrtype;

    using Vh_tEL =  Pdhms_type<mesh_type,Order>;
    using Wh_tEL =  Pdhv_type<mesh_type,Order>;
    using Vh_tPOI = Pdhv_type<mesh_type,Order>;
    using Wh_tPOI = Pdh_type<mesh_type,Order>;


    using op_interp_ptrtypeEL = std::shared_ptr<OperatorInterpolation<Wh_tEL, Pdhv_type<mesh_type,Order>>>;
    using opv_interp_ptrtypeEL = std::shared_ptr<OperatorInterpolation<Vh_tEL, Pdhms_type<mesh_type,Order>>>;

    using op_interp_ptrtypePOI = std::shared_ptr<OperatorInterpolation<Wh_tPOI, Pdh_type<mesh_type,Order>>>;
    using opv_interp_ptrtypePOI = std::shared_ptr<OperatorInterpolation<Vh_tPOI, Pdhv_type<mesh_type,Order>>>;


    //! Model properties type
    using model_prop_type = ModelProperties;
    using model_prop_ptrtype = std::shared_ptr<model_prop_type>;

    typedef Exporter<mesh_type,1> exporter_type;
    typedef std::shared_ptr <exporter_type> exporter_ptrtype;



private:

    mesh_ptrtype M_mesh;

    mp_ptrtype M_PoissonModel;
    me_ptrtype M_ElasticityModel;

    exporter_ptrtype M_exporter;

public:
    MixedPoissonElasticity(mesh_ptrtype _meshPoisson, mesh_ptrtype _meshElasticity, mesh_ptrtype meshCommon, mesh_ptrtype meshVisu = nullptr)
        {
            M_mesh = meshCommon;


            M_PoissonModel = mp_type::New(_prefix="hdg.poisson");
            M_PoissonModel->setMesh( _meshPoisson );
            M_PoissonModel->init();
            M_PoissonModel->algebraicFactory()->addFunctionLinearAssembly(std::bind( &self_type::poissonAssembly, std::ref(*this), std::placeholders::_1 ));

            M_ElasticityModel = me_type::New("hdg.elasticity");
            M_ElasticityModel->init( _meshElasticity , meshVisu );


            this->initExporter( meshVisu );
        }

    mesh_ptrtype mesh() const { return M_mesh; }

    void poissonAssembly( ModelAlgebraic::DataUpdateLinear & data );
    // void assembleF_Poisson();               // this is the assembleF of MixedPoisson

    void assembleF_Elasticity();            // this is the assembleF of MixedElasticity

    // void solvePoisson() { this->assembleF_Poisson(); M_PoissonModel->solve(); }
    // void solveElasticity() { this->assembleF_Elasticity(); M_ElasticityModel->solve(); }

    void run( op_interp_ptrtypeEL Idh_el = nullptr, opv_interp_ptrtypeEL Idhv_el = nullptr,
              op_interp_ptrtypePOI Idh_poi = nullptr, opv_interp_ptrtypePOI Idhv_poi = nullptr );

    // Exporter methods
    exporter_ptrtype GetExporter() { return M_exporter; }
    void initExporter( mesh_ptrtype meshVisu = nullptr );
    void exportResults ( double Time, mesh_ptrtype mesh = nullptr,  op_interp_ptrtypeEL Idh_el = nullptr, opv_interp_ptrtypeEL Idhv_el = nullptr,
                         op_interp_ptrtypePOI Idh_poi = nullptr, opv_interp_ptrtypePOI Idhv_poi = nullptr  ) ;

    void exportResults( mesh_ptrtype mesh = nullptr, op_interp_ptrtypeEL Idh_el = nullptr, opv_interp_ptrtypeEL Idhv_el = nullptr,
                        op_interp_ptrtypePOI Idh_poi = nullptr, opv_interp_ptrtypePOI Idhv_poi = nullptr  )
        {
            this->exportResults (M_PoissonModel->currentTime(), mesh , Idh_el, Idhv_el, Idh_poi, Idhv_poi);
            M_exporter -> save();
        }

}; // end class declaration

template <int Dim, int Order,int G_Order, int E_Order>
void
MixedPoissonElasticity<Dim,Order,G_Order,E_Order>::poissonAssembly( ModelAlgebraic::DataUpdateLinear & data )
{
    // build each time, not just for the constant part
    bool buildCstPart = data.buildCstPart();
    if( buildCstPart )
        return;
    // retrieve matrix and vector already assemble
    auto F = std::dynamic_pointer_cast<condensed_vector_t<typename mp_type::value_type>>(data.rhs());
    auto ps = M_PoissonModel->spaceProduct();
    auto blf = blockform1( ps, F);
    auto w = M_PoissonModel->spacePotential()->element();
    auto dt = M_PoissonModel->timeStepBdfPotential()->timeStep();
    auto disp = M_ElasticityModel-> potentialField();
    auto dispm = M_ElasticityModel->timeStepNM()->previousUnknown();

    blf(1_c) += integrate( _range=elements(M_mesh),
                           _expr=-(div(disp)-div(dispm))*id(w)/dt );
}

// template <int Dim, int Order,int G_Order, int E_Order>
// void
// MixedPoissonElasticity<Dim,Order,G_Order,E_Order>::assembleF_Poisson()
// {
//     auto ps = M_PoissonModel->getPS();
//     auto F = M_PoissonModel->getF();

//     auto blf = blockform1 (*ps, F);
//     auto w = M_PoissonModel->potentialSpace()->element();
//     auto dt = M_PoissonModel->timeStepBDF()->timeStep();

//     // - <d/dt div(u),w>
//     blf(1_c) += integrate( _range=elements( M_mesh ),
//                            _expr= -( div(M_ElasticityModel-> potentialField())-div(M_ElasticityModel->timeStepNM()->previousUnknown()) ) * id(w) / dt );

// }

template <int Dim, int Order,int G_Order, int E_Order>
void
MixedPoissonElasticity<Dim,Order,G_Order,E_Order>::assembleF_Elasticity()
{

    auto ps = product( M_ElasticityModel->fluxSpace(), M_ElasticityModel->potentialSpace(), M_ElasticityModel->traceSpace() );

    auto blf = blockform1 ( ps, M_ElasticityModel->getF() );
    auto v = M_ElasticityModel->fluxSpace()->element( "v" );
    auto m = M_ElasticityModel->traceSpace()->element( "m" );

    auto pressure = M_PoissonModel->fieldPotential();

    // - < pI , v>
    blf( 0_c ) += integrate( _range=elements( M_mesh ),
                             _expr= - inner( idv(pressure)*eye<Dim,Dim>(),  id(v)) );

    // Adding extra-term for special Neumann
    // - <pI n, m>
    std::string marker = "";

    for( auto const& pairMat : M_ElasticityModel->modelProperties().materials() )
    {
        auto material = pairMat.second;
        marker = material.getString("special_neumann");
    }

    if (!marker.empty())
        blf( 2_c ) += integrate( _range=markedfaces(M_mesh,marker),
                                 _expr= - inner( idv(pressure)*eye<Dim,Dim>() * N(), id(m)) );

}

template <int Dim, int Order,int G_Order, int E_Order>
void
MixedPoissonElasticity<Dim,Order,G_Order,E_Order>::run( op_interp_ptrtypeEL Idh_el, opv_interp_ptrtypeEL Idhv_el,
                                                        op_interp_ptrtypePOI Idh_poi, opv_interp_ptrtypePOI Idhv_poi  )
{
    if (M_PoissonModel->isStationary() || M_ElasticityModel->isStationary() )
    {
        Feel::cout << std::endl << "ERROR: this model has to be unsteady." << std::endl << std::endl;
        return;
    }


    auto t_init = M_PoissonModel->timeStepBdfPotential()->timeInitial();
    auto dt = M_PoissonModel->timeStepBdfPotential()->timeStep();
    auto t_fin = M_PoissonModel->timeStepBdfPotential()->timeFinal();

    // M_PoissonModel->assembleCstPart();

    for (; !M_PoissonModel->timeStepBase()->isFinished() && !M_ElasticityModel->timeStepBase()->isFinished() ; M_PoissonModel->updateTimeStep() )
    {
        Feel::cout << "===============================================" << std::endl;
        Feel::cout << "time simulation: " << M_PoissonModel->time() << "s \n";
        Feel::cout << "===============================================" << std::endl;

        // Elasticity problem
        M_ElasticityModel->assembleCst();
        M_ElasticityModel->assembleNonCst();
        this->assembleF_Elasticity();
        M_ElasticityModel->solve();
        M_ElasticityModel->exportResults( M_ElasticityModel->mesh(), Idh_el, Idhv_el );

        // Poisson problem
        // M_PoissonModel->assembleNonCstPart();
        // M_PoissonModel->assembleAll();
        // this->assembleF_Poisson();
        M_PoissonModel->solve();
        M_PoissonModel->exportResults();

        // Exporter
        // this->exportResults( mesh, Idh_el, Idhv_el, Idh_poi, Idhv_poi );

        // update
        M_ElasticityModel->updateTimeStep();
    }

}

// EXPORTER
template <int Dim, int Order,int G_Order, int E_Order>
void
MixedPoissonElasticity<Dim,Order,G_Order,E_Order>::initExporter( mesh_ptrtype meshVisu )
{
    std::string geoExportType="static"; //change_coords_only, change, static
    M_exporter = exporter ( _mesh=meshVisu?meshVisu:this->mesh(),
                            _name="Export",
                            _geo=geoExportType,
                            _path=M_PoissonModel->exporterPath() );

}

template <int Dim, int Order,int G_Order, int E_Order>
void
MixedPoissonElasticity<Dim,Order,G_Order,E_Order>::exportResults ( double time, mesh_ptrtype mesh,
                                                                   op_interp_ptrtypeEL Idh_el, opv_interp_ptrtypeEL Idhv_el,
                                                                   op_interp_ptrtypePOI Idh_poi, opv_interp_ptrtypePOI Idhv_poi )
{
    Feel::cout << "[MixedPoissonElasticity]:[exportResults]:[start]" << std::endl;


    if ( M_exporter->exporterGeometry() != EXPORTER_GEOMETRY_STATIC && mesh  )
    {
        LOG(INFO) << "exporting on visualisation mesh at time " << time;
        M_exporter->step( time )->setMesh( mesh );
    }

    // Export computed solutions for elasticity
    auto prefix = M_ElasticityModel->prefix();
    auto postProcess = M_ElasticityModel->modelProperties().postProcess();
    auto BC_Elasticity = M_ElasticityModel->modelProperties().boundaryConditions();
    auto itField = postProcess.find( "Fields");
    if ( itField != postProcess.end() )
    {
        for ( auto const& field : (*itField).second )
        {
            if ( field == "strain" )
            {
                auto M_up = M_ElasticityModel->fluxField();
                LOG(INFO) << "exporting strain at time " << time;
                M_exporter->step(time)->add(prefixvm(prefix, "strain"), Idhv_el?(*Idhv_el)(M_up):M_up );
            }
            else if ( field == "displacement" )
            {
                auto M_pp = M_ElasticityModel->potentialField();
                LOG(INFO) << "exporting displacement at time " << time;
                M_exporter->step(time)->add(prefixvm(prefix, "displacement"),Idh_el?(*Idh_el)(M_pp): M_pp ) ;

                auto itField = BC_Elasticity.find("ExactSolution");
                if ( itField != BC_Elasticity.end() )
                {
                    auto mapField = (*itField).second;
                    auto itType = mapField.find( "u_exact" );
                    if (itType != mapField.end() )
                    {
                        for (auto const& exAtMarker : (*itType).second )
                        {
                            if (exAtMarker.isExpression() )
                            {
                                auto u_exact = expr<Dim,1> (exAtMarker.expression());
                                if ( !M_ElasticityModel->isStationary() )
                                    u_exact.setParameterValues( { {"t", time } } );

                                auto export_uEX = project( _space=M_ElasticityModel->potentialSpace(), _range=elements( M_ElasticityModel->mesh() ), _expr=u_exact);
                                M_exporter->step(time)->add(prefixvm(prefix, "u_exact"), Idh_el?(*Idh_el)( export_uEX): export_uEX );

                                auto l2err_u = normL2( _range=elements(M_ElasticityModel->mesh()), _expr=u_exact - idv(M_ElasticityModel->potentialField()) );
                                auto l2norm_uex = normL2( _range=elements(M_ElasticityModel->mesh()), _expr=u_exact );
                                if (l2norm_uex < 1)
                                    l2norm_uex = 1.0;
                                Feel::cout << "----- Computed Errors -----" << std::endl;
                                Feel::cout << "||u-u_ex||_L2=\t" << l2err_u/l2norm_uex << std::endl;
                                // Export the errors
                                M_exporter -> step( time )->add(prefixvm(prefix, "u_error_L2"), l2err_u/l2norm_uex );
                            }
                        }
                    }
                }
            }
        }
    }

    // Export computed results for poisson
    prefix = M_PoissonModel->prefix();
    auto postProcess2 = M_PoissonModel->modelProperties().postProcess();
    auto BC_Poisson = M_PoissonModel->modelProperties().boundaryConditions();
    auto itField2 = postProcess2.find( "Fields");
    if ( itField2 != postProcess2.end() )
    {
        for ( auto const& field : (*itField2).second )
        {
            if ( field == "flux" )
            {
                auto M_up = M_PoissonModel->fieldFlux();
                LOG(INFO) << "exporting flux at time " << time;
                M_exporter->step( time )->add(prefixvm(prefix, "flux"), Idhv_poi?(*Idhv_poi)( M_up):M_up );
            }
            else if ( field == "potential" )
            {
                auto M_pp = M_PoissonModel->fieldPotential();
                LOG(INFO) << "exporting potential at time " << time;
                M_exporter->step( time )->add(prefixvm(prefix, "pressure"), Idh_poi?(*Idh_poi)(M_pp):M_pp );

                auto itField = BC_Poisson.find("Exact solution");
                if ( itField != BC_Poisson.end() )
                {
                    auto mapField = (*itField).second;
                    auto itType = mapField.find( "p_exact" );
                    if (itType != mapField.end() )
                    {
                        for (auto const& exAtMarker : (*itType).second )
                        {
                            if (exAtMarker.isExpression() )
                            {
                                auto p_exact = expr(exAtMarker.expression() );
                                if ( !M_PoissonModel->isStationary() )
                                    p_exact.setParameterValues( { {"t", time } } );
                                M_exporter->step( time )->add(prefixvm(prefix, "p_exact"),
                                                              project( _space=M_PoissonModel->spacePotential(),
                                                                       _range=elements(M_PoissonModel->mesh()),
                                                                       _expr=p_exact) );

                                auto l2err_p = normL2( _range=elements(M_PoissonModel->mesh()), _expr=p_exact - idv(M_PoissonModel->fieldPotential()) );
                                auto l2norm_pex = normL2( _range=elements(M_PoissonModel->mesh()), _expr=p_exact );
                                if (l2norm_pex < 1)
                                    l2norm_pex = 1.0;
                                Feel::cout << "||p-p_ex||_L2=\t" << l2err_p/l2norm_pex << std::endl;
                                Feel::cout << "---------------------------" << std::endl;
                                // Export the errors
                                M_exporter -> step( time )->add(prefixvm(prefix, "p_error_L2"),
                                                                l2err_p/l2norm_pex );
                            }
                        }
                    }
                }
            }
        }
    }

    Feel::cout << "[MixedPoissonElasticity]:[exportResults]:[finish]" << std::endl;
}




} // end namespace FeelModels

} // end namespace Feel

#endif


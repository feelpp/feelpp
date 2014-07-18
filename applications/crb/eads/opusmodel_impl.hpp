/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-12-10

  Copyright (C) 2008-2010 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file opusmodel.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-12-10
 */
#if !defined(OPUSMODEL_IMPL_HPP_)
#define OPUSMODEL_IMPL_HPP_ 1

#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>

#include <feel/feel.hpp>
#include <feel/feelvf/print.hpp>

#include <opusmodel.hpp>
#include <opusmodelthermal.hpp>
#include <opusmodelfluidpoiseuille.hpp>
#include <opusmodelfluidoseen.hpp>

namespace Feel
{
/**
 * \addtogroup Models
 * \\@{
 */
template<int OrderU, int OrderP, int OrderT>
OpusModel<OrderU,OrderP,OrderT>::OpusModel( OpusModel const& om )
    :
    super( om ),
    M_force_rebuild( false ),
    M_mesh( om.M_mesh ),
    M_thermal( om.M_thermal )
{}

template<int OrderU, int OrderP, int OrderT>
OpusModel<OrderU,OrderP,OrderT>::OpusModel( po::variables_map const& vm )
    :
    super( vm ),
    mu_file( vm["mufile"].as<std::string>() ),
    M_dt( vm["bdf.time-step"].as<double>() ),
    M_meshsize( vm["hsize"].as<double>() ),
    M_force_rebuild( false ),
    M_mesh( new mesh_type ),
    M_mesh_air( new mesh_type ),
    M_mesh_line( new mesh12_type ),
    M_mesh_cross_section_2( new mesh12_type ),
    M_thermal (),
    M_exporter(),
    M_exporter_fluid()
{
    //this->init();
}

template<int OrderU, int OrderP, int OrderT>
OpusModel<OrderU,OrderP,OrderT>::OpusModel()
    :
    super(),
    M_dt( 0.1 ),
    M_meshsize( 0.2 ),
    M_force_rebuild( false ),
    M_mesh( new mesh_type ),
    M_mesh_air( new mesh_type ),
    M_mesh_line( new mesh12_type ),
    M_mesh_cross_section_2( new mesh12_type ),
    M_thermal (),
    M_exporter(),
    M_exporter_fluid()
{
    //this->init();
}
template<int OrderU, int OrderP, int OrderT>
void
OpusModel<OrderU,OrderP,OrderT>::init()
{
    M_exporter = boost::shared_ptr<export_type> ( Exporter<mesh_type,1>::New( "ensight", "opus" ) );
    M_exporter->setPrefix( "opus" );
    M_exporter_fluid = boost::shared_ptr<export_type> ( Exporter<mesh_type,1>::New( "ensight", "fluid" ) );
    M_exporter_fluid->setPrefix( "fluid" );

    M_mesh = createGMSHMesh( _mesh=new mesh_type,
                             _desc =  this->data()->createMesh( M_meshsize ),
                             _update = MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                             _force_rebuild = M_force_rebuild );
    LOG(INFO) << "Imported mesh thermal\n";
    M_mesh_air = createGMSHMesh( _mesh=new mesh_type,
                                 _desc =  this->data()->createMeshAir( M_meshsize ),
                                 _update = MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                 _force_rebuild = M_force_rebuild );
    LOG(INFO) << "Imported mesh air\n";
    M_mesh_line = createGMSHMesh( _mesh=new mesh12_type,
                                  _desc =  this->data()->createMeshLine( 1 ),
                                  _update = MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                  _force_rebuild = M_force_rebuild );
    LOG(INFO) << "Imported mesh line\n";
    M_mesh_cross_section_2 = createGMSHMesh( _mesh=new mesh12_type,
                             _desc =  this->data()->createMeshCrossSection2( 0.2 ),
                             _update = MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                             _force_rebuild = M_force_rebuild );
    LOG(INFO) << "[init] Imported mesh cross section 2\n";

    M_P1h = p1_functionspace_type::New( M_mesh_line );
    LOG(INFO) << "[init] P1 mesh\n";
    M_P0h = p0_space_type::New( M_mesh );
    LOG(INFO) << "[init] P0 mesh\n";
    typedef typename node<double>::type node_type;

    node_type period( 2 );
    period[0]=this->data()->component( "PCB" ).e()+this->data()->component( "AIR" ).e();
    period[1]=0;
    LOG(INFO) << "[init] period=" << period[0] << "," << period[1] << "\n";
    //M_Th = temp_functionspace_type::New( _mesh=M_mesh,
    //                                     _periodicity=Periodic<1,2,value_type>( period ) );
    M_Th = temp_functionspace_type::New( _mesh=M_mesh,
                                         _periodicity=periodicity( Periodic<>( 1,2, period ) ) );
    LOG(INFO) << "[init] M_Th init done\n";
    M_grad_Th = grad_temp_functionspace_type::New( _mesh=M_mesh );
    LOG(INFO) << "[init] M_grad_Th init done\n";
    M_Xh = fluid_functionspace_type::New( M_mesh );
    LOG(INFO) << "[init] M_Xh init done\n";
    //M_Xh = fluid_functionspace_type::New( M_mesh_air );
    //M_Xh = oseen_functionspace_type::New( M_mesh->extract(  ) );

    LOG(INFO) << "Generated function space\n";
    LOG(INFO) << " o        number of elements :  " << M_mesh->numElements() << "\n";
    LOG(INFO) << " o          number of points :  " << M_mesh->numPoints() << "\n";
    LOG(INFO) << " o number of local dof in Th :  " << M_Th->nLocalDof() << "\n";
    LOG(INFO) << " o       number of dof in Th :  " << M_Th->nDof() << "\n";
    LOG(INFO) << " o       number of dof in Th :  " << M_Th->dof()->nDof() << "\n";



    M_thermal = thermal_operator_ptrtype( new thermal_operator_type( this->vm(), M_Th ) );
    LOG(INFO) << "Generated thermal operator\n";
    M_fluid = fluid_operator_ptrtype( new fluid_operator_type( this->vm(), M_Xh ) );

    LOG(INFO) << "Generated fluid operator\n";

    using namespace vf;

    LOG(INFO) << "[OpusModel::OpusModel] start bdf\n";
    M_temp_bdf = bdf( _space=M_Th, _vm=this->vm(), _name="temperature" );
    M_temp_bdf->print();
    LOG(INFO) << "[OpusModel::OpusModel] temp bdf done\n";
    M_fluid_bdf = bdf( _space=M_Xh, _vm=this->vm(), _name="fluid" );
    M_fluid_bdf->print();
    LOG(INFO) << "[OpusModel::OpusModel] bdf stops\n";
}

template<int OrderU, int OrderP, int OrderT>
OpusModel<OrderU,OrderP,OrderT>::~OpusModel()
{}

template<int OrderU, int OrderP, int OrderT>
void
OpusModel<OrderU,OrderP,OrderT>::run ( const double * X, unsigned long N,
                                       double * Y, unsigned long P )
{
    LOG(INFO) << "[OpusModel::run] input/output relationship\n";

    for ( int i = 0; i < N; ++i )
        LOG(INFO) << "[OpusModel::run] X[" << i << "]=" << X[i] << "\n";

    this->data()->component( "IC1" ).setK( X[0] );
    this->data()->component( "IC2" ).setK( X[0] );
    this->data()->component( "AIR" ).setFlowRate( X[1] );
    this->data()->component( "IC1" ).setQ( X[2] );
    this->data()->component( "IC2" ).setQ( X[2] );

    LOG(INFO) << "[OpusModel::run] parameters set\n";

    // check if the mesh size or e_a have been changed since last run, if yes
    // then the geometry and mesh need to be rebuilt
    if ( ( math::abs( M_meshsize - X[5] ) > 1e-5 ) ||
            ( math::abs( this->data()->component( "AIR" ).e() - X[4] ) > 1e-10 ) )
        M_force_rebuild = true;

    else
        M_force_rebuild = false;

    M_force_rebuild = true;
    this->data()->component( "AIR" ).setE( X[4] );
    M_meshsize = X[5];

    this->data()->print();

    LOG(INFO) << "[OpusModel::run] parameters print\n";

    LOG(INFO) << "[OpusModel::run] start init\n";
    this->init();

    LOG(INFO) << "[OpusModel::run] init doned\n";
    M_thermal->setData( this->data() );
    M_fluid->setData( this->data() );
    M_thermal->setThermalConductance( X[3] );
    M_fluid->setFluidFlowRate( X[1] );

    LOG(INFO) << "[OpusModel::run] parameters set\n";
    this->data()->print();
    LOG(INFO) << "[OpusModel::run] run\n";
    this->run();
    Y[0]=s1;
    Y[1]=s2;
    LOG(INFO) << "[OpusModel::run] run done, set outputs\n";

    for ( int i = 0; i < P; ++i )
        LOG(INFO) << "[OpusModel::run] Y[" << i << "]=" << Y[i] << "\n";

}

template<int OrderU, int OrderP, int OrderT>
void
OpusModel<OrderU,OrderP,OrderT>::run()
{
    LOG(INFO) << "[OpusModel::run] starts\n";
    using namespace vf;
    domains = p0_element_ptrtype( new p0_element_type( M_P0h, "domains" ) );
    *domains = vf::project( M_P0h, elements( M_P0h->mesh() ),
                            chi( emarker() == M_Th->mesh()->markerName( "PCB" ) )*M_Th->mesh()->markerName( "PCB" )+
                            chi( emarker() == M_Th->mesh()->markerName( "AIR123" ) )*M_Th->mesh()->markerName( "AIR" )+
                            chi( emarker() == M_Th->mesh()->markerName( "AIR4" ) )*M_Th->mesh()->markerName( "AIR" )+
                            chi( emarker() == M_Th->mesh()->markerName( "IC1" ) )*M_Th->mesh()->markerName( "IC1" ) +
                            chi( emarker() == M_Th->mesh()->markerName( "IC2" ) )*M_Th->mesh()->markerName( "IC2" ) );

    k = p0_element_ptrtype( new p0_element_type( M_P0h, "k" ) );
    *k = vf::project( M_P0h, elements( M_P0h->mesh() ),
                      chi( emarker() == M_Th->mesh()->markerName( "PCB" ) )*this->data()->component( "PCB" ).k()+
                      chi( emarker() == M_Th->mesh()->markerName( "AIR123" ) )*this->data()->component( "AIR" ).k()+
                      chi( emarker() == M_Th->mesh()->markerName( "AIR4" ) )*this->data()->component( "AIR" ).k()+
                      chi( emarker() == M_Th->mesh()->markerName( "IC1" ) )*this->data()->component( "IC1" ).k()+
                      chi( emarker() == M_Th->mesh()->markerName( "IC2" ) )*this->data()->component( "IC2" ).k() );
    rhoC = p0_element_ptrtype( new p0_element_type( M_P0h, "rhoC" ) );
    *rhoC = vf::project( M_P0h, elements( M_P0h->mesh() ),
                         chi( emarker() == M_Th->mesh()->markerName( "PCB" ) )*this->data()->component( "PCB" ).rhoC()+
                         chi( emarker() == M_Th->mesh()->markerName( "AIR123" ) )*this->data()->component( "AIR" ).rhoC()+
                         chi( emarker() == M_Th->mesh()->markerName( "AIR4" ) )*this->data()->component( "AIR" ).rhoC()+
                         chi( emarker() == M_Th->mesh()->markerName( "IC1" ) )*this->data()->component( "IC1" ).rhoC() +
                         chi( emarker() == M_Th->mesh()->markerName( "IC2" ) )*this->data()->component( "IC2" ).rhoC() );

    Q = p0_element_ptrtype( new p0_element_type( M_P0h, "Q" ) );
    *Q = vf::project( M_P0h, elements( M_P0h->mesh() ),
                      chi( emarker() == M_Th->mesh()->markerName( "IC1" ) )*this->data()->component( "IC1" ).Q() +
                      chi( emarker() == M_Th->mesh()->markerName( "IC2" ) )*this->data()->component( "IC2" ).Q() );
    LOG(INFO) << "[OpusModel::OpusModel] P0 functions allocated\n";



    std::ofstream os6( "Outputs.dat" );
    os6.precision( 10 );
    os6.width( 15 );
    os6.setf( std::ios::right );

    double surf1 = integrate( markedelements( M_mesh,M_mesh->markerName( "IC2" ) ),constant( 1.0 ),_Q<0>() ).evaluate()( 0,0 );
    double len2 = ( integrate( markedfaces( M_mesh,M_mesh->markerName( "Gamma_3_AIR3" ) ),constant( 1.0 ),_Q<0>() ).evaluate()( 0,0 ) +
                    integrate( markedfaces( M_mesh,M_mesh->markerName( "Gamma_3_AIR4" ) ),constant( 1.0 ),_Q<0>() ).evaluate()( 0,0 ) );
    std::ofstream outputs( "outputs.dat" );
    std::ostringstream os ;
    LOG(INFO) << "output file set\n";
    temp_element_type T( M_Th, "temperature" );

    temp_element_type Temperature( M_Th, "temperature" );

    boost::timer ti;

    T = vf::project( M_Th, elements( M_Th->mesh() ), constant( M_thermal->T0() ) );
    fluid_element_type U( M_Xh, "fluid" );
    fluid_element_0_type u = U.template element<0>();
    fluid_element_1_type p = U.template element<1>();

    double flow_rate = this->vm()["fluid.flow-rate"].template as<double>();
    double e_AIR = this->data()->component( "AIR" ).e();
    double e_PCB = this->data()->component( "PCB" ).e();
    double e_IC = this->data()->component( "IC1" ).e();
    //double L_IC = this->data()->component("IC1").h();

    LOG(INFO) << "[opusmodel] flow_rate = " << flow_rate << "\n";
    LOG(INFO) << "[opusmodel] e_AIR = " << e_AIR << "\n";
    LOG(INFO) << "[opusmodel] e_PCB = " << e_PCB << "\n";
    LOG(INFO) << "[opusmodel] e_IC = " << e_IC << "\n";
    double time = M_temp_bdf->timeInitial();
    auto chi_AIR = chi( Px() >= e_PCB+e_IC );
    auto ft = constant( 1.0-( !this->data()->isSteady() )*math::exp( -time/3.0 ) );
    auto vy = ft*constant( 3. )/( 2.*( e_AIR-e_IC ) )*flow_rate*( 1.-vf::pow( ( Px()-( ( e_AIR+e_IC )/2+e_PCB ) )/( ( e_AIR-e_IC )/2 ),2 ) );

    u = vf::project( _space=M_Xh->template functionSpace<0>(), _expr=vec( constant( 0. ),vy ), _range=markedelements( M_Xh->mesh(), "AIR4" ) );
    p = vf::project( _space=M_Xh->template functionSpace<1>(), _expr=constant( 0. ) , _range=markedelements( M_Xh->mesh(), "AIR4" ));


    LOG(INFO) << "fluid and temperature fields set\n";

    if ( !this->data()->isSteady() )
    {
        if ( M_temp_bdf->timeInitial() > 0.0 )
        {
            T = M_temp_bdf->unknown( 0 );
        }

        else
        {
            M_temp_bdf->initialize( T );
        }

        if ( M_fluid_bdf->timeInitial() > 0.0 )
        {
            U = M_fluid_bdf->unknown( 0 );
        }

        else
        {
            M_fluid_bdf->initialize( U );
        }
    }

    else
    {
        M_temp_bdf->setSteady();
        M_fluid_bdf->setSteady();
    }

    LOG(INFO) << "[initialization] done in " << ti.elapsed() << "\n";

    for ( M_temp_bdf->start(), M_fluid_bdf->start() ;
            ( M_temp_bdf->isFinished() == false ) && ( M_fluid_bdf->isFinished() == false );
            M_temp_bdf->next(), M_fluid_bdf->next() )
    {
        LOG(INFO) << "============================================================\n";
        LOG(INFO) << "time(T): " << M_temp_bdf->time() << "s, iteration: " << M_temp_bdf->iteration() << " order:"  << M_temp_bdf->timeOrder() << "\n";
        LOG(INFO) << "time(U): " << M_fluid_bdf->time() << "s, iteration: " << M_fluid_bdf->iteration() << " order:"  << M_fluid_bdf->timeOrder() << "\n";
        std::cout << "============================================================\n";
        std::cout << "time(T): " << M_temp_bdf->time() << "s, iteration: " << M_temp_bdf->iteration() << " order:"  << M_temp_bdf->timeOrder() << "\n";
        std::cout << "time(U): " << M_fluid_bdf->time() << "s, iteration: " << M_fluid_bdf->iteration() << " order:"  << M_fluid_bdf->timeOrder() << "\n";

        // Fluide
        ti.restart();
        M_fluid->update( M_fluid_bdf->time() );
        LOG(INFO) << "[fluid] update done in " << ti.elapsed() << "\n";
        ti.restart();
        M_fluid->solve( U );
        //M_exporter->step(time)->setMesh( U.functionSpace()->mesh() );
        M_exporter_fluid->step( M_fluid_bdf->time() )->setMesh( U.functionSpace()->mesh() );
        M_exporter_fluid->step( M_fluid_bdf->time() )->add( "Velocity",  U.template element<0>() );
        M_exporter_fluid->step( M_fluid_bdf->time() )->add( "Pressure",  U.template element<1>() );
        M_exporter_fluid->save();

        LOG(INFO) << "[fluid] solve done in " << ti.elapsed() << "\n";

        // Thermal
        ti.restart();

        double thetime = M_temp_bdf->time();

        M_thermal->update( thetime,
                           ( !this->data()->isSteady() )*idv( *rhoC )*M_temp_bdf->polyDerivCoefficient( 0 ), // mass
                           ( idv( *k ) ), // diff
                           ( idv( *rhoC )*idv( U.template element<0>() ) ), // conv
                           //( print(idv(*Q)*(1.0-vf::exp(-cst_ref(thetime))),"Q=") + // source term
                           ( print( idv( *Q ),"Q=" ) + // source term
                             ( !this->data()->isSteady() )*print( idv( *rhoC )*print( idv( M_temp_bdf->polyDeriv() ),"Tn=" ),"rhoc*Tn" ) ) // bdf contrib
                         );

        LOG(INFO) << "[thermal] update done in " << ti.elapsed() << "\n";
        ti.restart();
        M_thermal->solve( T );

        LOG(INFO) << "[thermal] solve done in " << ti.elapsed() << "\n";

        // export results
        ti.restart();
        this->exportResults( thetime, T, U );

        LOG(INFO) << "[export] export done in " << ti.elapsed() << "\n";

        ti.restart();
        //double surf1 = integrate( markedelements(M_mesh,M_mesh->markerName( "IC2" )), _Q<0>(),constant(1.0)).evaluate()(0,0);
        double surf1_exact = this->data()->component( "IC2" ).e()*this->data()->component( "IC2" ).h();

        if ( math::abs( surf1_exact - surf1 ) > 1e-10 )
        {
            LOG(INFO) << "[s1] Invalid IC surface computation\n";
            LOG(INFO) << "[s1] surface = " << surf1 << "\n";
            LOG(INFO) << "[s1] surface(exact) = " << surf1_exact << "\n";
        }

        s1 = integrate( markedelements( M_mesh,M_mesh->markerName( "IC2" ) ), idv( T ),_Q<OrderT>() ).evaluate()( 0,0 )/surf1;
        //double len2 = integrate( markedfaces(M_mesh,M_mesh->markerName( "Gamma_3_AIR4" )), _Q<0>(),constant(1.0)).evaluate()(0,0);
        double len2_exact = this->data()->component( "AIR" ).e();

        if ( math::abs( len2_exact - len2 ) > 1e-10 )
        {
            LOG(INFO) << "[s2] Invalid Gamma_3_AIR4 length computation\n";
            LOG(INFO) << "[s2] length = " << len2 << "\n";
            LOG(INFO) << "[s2] length(exact) = " << len2_exact << "\n";
        }

        s2 = ( integrate( markedfaces( M_mesh,M_mesh->markerName( "Gamma_3_AIR3" ) ), idv( T ),_Q<OrderT>() ).evaluate()( 0,0 )+
               integrate( markedfaces( M_mesh,M_mesh->markerName( "Gamma_3_AIR4" ) ), idv( T ),_Q<OrderT>() ).evaluate()( 0,0 ) )/len2;
        outputs.precision( 10 );
        outputs.setf( std::ios_base::scientific, std::ios_base::floatfield );
        outputs.width( 15 );
        outputs << s1 << " " << s2 << std::endl;
        LOG(INFO) << "[s1,s2] postprocess done in " << ti.elapsed() << "\n";
        LOG(INFO) << "s1=" << s1 << "\n";
        LOG(INFO) << "s2=" << s2 << "\n";
        std::cout << "s1=" << s1 << " s2=" << s2 << "\n";
        os6 << M_temp_bdf->time() << " "  << s1 << "  "  << s2 << std::endl;

        AUTO( N_IC_PCB,vec( constant( -1. ),constant( 0. ) ) );
        double meas_PCB = integrate( markedfaces( M_mesh, M_mesh->markerName( "Gamma_IC1_PCB" ) ),constant( 1.0 ),_Q<0>() ).evaluate()( 0,0 ) ;
        double mean_jump_1 = integrate( markedfaces( M_mesh, M_mesh->markerName( "Gamma_IC1_PCB" ) ),
                                        trans( jumpv( idv( T ) ) )*N_IC_PCB,_Q<OrderT>() ).evaluate()( 0,0 )/meas_PCB;
        double mean_jump_2 = integrate( markedfaces( M_mesh, M_mesh->markerName( "Gamma_IC2_PCB" ) ),
                                        trans( jumpv( idv( T ) ) )*N_IC_PCB,_Q<OrderT>() ).evaluate()( 0,0 )/meas_PCB;

        LOG(INFO) <<  "meas(Gamma_IC1_PCB) = " << meas_PCB << "\n";
        LOG(INFO) <<  "mean([[T]],IC1) = " << mean_jump_1 << "\n";
        LOG(INFO) <<  "mean([[T]],IC2) = " << mean_jump_2 << "\n";

        LOG(INFO) << "flux between IC1 and PCB = " << integrate( markedfaces( M_mesh, M_mesh->markerName( "Gamma_IC1_PCB" ) ),
                jumpv( idv( *k )*gradv( T ) ),_Q<OrderT-1>() ).evaluate()( 0,0 ) << "\n";
        LOG(INFO) << "flux between IC2 and PCB = " << integrate( markedfaces( M_mesh, M_mesh->markerName( "Gamma_IC2_PCB" ) ),
                jumpv( idv( *k )*gradv( T ) ),_Q<OrderT-1>() ).evaluate()( 0,0 ) << "\n";
        LOG(INFO) << "[k]_{IC2 and PCB} = " << integrate( markedfaces( M_mesh, M_mesh->markerName( "Gamma_IC2_PCB" ) ),
                jumpv( idv( *k ) )*N_IC_PCB,_Q<0>() ).evaluate()( 0,0 ) << "\n";


        ti.restart();
        M_temp_bdf->shiftRight( T );
        M_fluid_bdf->shiftRight( U );
        LOG(INFO) << "[bdf] shifRight done in " << ti.elapsed() << "\n";
        LOG(INFO) << "time spent in iteration = " << M_temp_bdf->realTimePerIteration() << "s\n";
    }

}


template<int OrderU, int OrderP, int OrderT>
void
OpusModel<OrderU,OrderP,OrderT>::exportResults( double time, temp_element_type& T, fluid_element_type& U, bool force_export  )
{

    std::ostringstream osstr ;

    int j = time;
    osstr<<j;

    if ( force_export || this->data()->doExport() )
    {
        //M_exporter->step(time)->setMesh( U.functionSpace()->mesh() );
        M_exporter->step( time )->setMesh( T.functionSpace()->mesh() );
        M_exporter->step( time )->add( "Domains", *domains );
        M_exporter->step( time )->add( "k", *k );
        M_exporter->step( time )->add( "rhoC", *rhoC );
        M_exporter->step( time )->add( "Q", *Q );
        M_exporter->step( time )->add( "Temperature", T );
        M_exporter->step( time )->add( "TVelocity",  U.template element<0>() );
        M_exporter->step( time )->add( "TPressure",  U.template element<1>() );
        M_exporter->save();
    }
}

/** \\@} */
}

#endif

//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
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
//! @author <you>
//! @date 15 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
//!
#define FEELPP_INSTANTIATE_BIOTSAVART_THERMOELECTRIC
#include <feel/feelcrb/crbplugin.hpp>


#include "biotsavart.hpp"
#include <thermoelectric-linear.hpp>

namespace Feel {

template<typename te_rb_model_type>
BiotSavartCRB<te_rb_model_type>::BiotSavartCRB()
    : super_type("biotsavart"),
      M_mode(( CRBModelMode )ioption(_name="biotsavart.run-mode" ) )
{
    std::string mshFileName = Environment::expand(soption("gmsh.filename"));
    bool repart = boption("biotsavart.repart");

    auto mesh = loadMesh( _mesh=new mesh_type,
                          _filename=mshFileName,
                          _rebuild_partitions=repart,
                          _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES);

    std::vector<std::string> conductor = vsoption("biotsavart.conductor");
    std::vector<std::string> mgn = vsoption("biotsavart.mgn");
    if( conductor.empty() || mgn.empty() )
    {
        LOG(FATAL) << "Marker for conductor or mgn empty! "
                   << "Use mesh.conductor and mesh.mgn to pass the markers";
    }
    std::list<std::string> conductorList(conductor.begin(), conductor.end());
    M_meshCond = createSubmesh(mesh, markedelements(mesh, conductorList));
    std::list<std::string> mgnList(mgn.begin(), mgn.end());
    M_meshMgn = createSubmesh(mesh, markedelements(mesh, mgnList));

    tic();
    M_teCrbModel = std::make_shared<te_rb_model_type>(M_meshCond);
    M_crbModel = std::make_shared<crb_model_type>(M_teCrbModel, crb::stage::offline);
    M_crb = crb_type::New("biotsavart_crb", M_crbModel);
    toc("constructor + eim");

    M_Xh = M_teCrbModel->functionSpace();
    M_XhMgn = vec_space_type::New( M_meshMgn);

    this->setupCommunicatorsBS();
}

template<typename te_rb_model_type>
void BiotSavartCRB<te_rb_model_type>::initModel()
{}

template<typename te_rb_model_type>
void BiotSavartCRB<te_rb_model_type>::runBS()
{
    M_mu = this->paramFromOption();
    if( boption("biotsavart.compute-offline") )
    {
        this->offline();
        if( boption("biotsavart.compute-online") )
            this->online(M_mu);
    }

    if( boption("biotsavart.compute-fe") )
        this->computeFE(M_mu);

    this->exportResults();
}

template<typename te_rb_model_type>
void BiotSavartCRB<te_rb_model_type>::offline()
{
    tic();

    tic();
    M_crb->offline();
    toc("CRB offline");

    this->setDimensions();

    this->loadIntegrals();

    toc("offline");
}

template<typename te_rb_model_type>
void BiotSavartCRB<te_rb_model_type>::setupCommunicatorsBS()
{
    tic();
    M_dofMgn.clear();
    M_commsC1M.clear();
    // Identify repartition of dofs (Omega_C, Omega_M)
    std::vector<int> isIn( 2, 0 );
    if( M_Xh->nLocalDof()  > 0 ) // current proc has dofs in Omega_C
        isIn[0] = 1;
    if( M_XhMgn->nLocalDof() > 0 ) // current proc has dofs in Omega_M
        isIn[1] = 1;

    std::vector<int> isInGlob( 2*Environment::worldComm().size() );
    mpi::all_gather( Environment::worldComm(), isIn.data(), 2, isInGlob );

    std::vector<int> pC,pM,pC1M;
    for(int i=0; i<Environment::worldComm().size(); i++)
    {
        if( isInGlob[2*i] )
            pC.push_back( i );
        if( isInGlob[2*i + 1] )
            pM.push_back( i ); //Check if pM has at least one element for each procs
    }

    // commC = communicator with all procs which in pC
    auto groupC = Environment::worldComm().group().include(pC.begin(), pC.end());
    mpi::communicator commC(Environment::worldComm(), groupC);
    // commM = communicator with all procs which in pM
    auto groupM = Environment::worldComm().group().include(pM.begin(), pM.end());
    mpi::communicator commM(Environment::worldComm(), groupM);

    dof_points_type dofM; //dofs

    for( int i=0; i<pM.size(); i++)
    {
        pC1M.clear();
        pC1M.resize( pC.size() );
        pC1M = pC;
        auto find = std::find( pC1M.begin(), pC1M.end(), pM[i] );
        if( find != pC1M.end() )
            pC1M.erase(find);  //No double definition
        pC1M.insert(pC1M.begin(), pM[i] );

        auto groupC1M = Environment::worldComm().group().include(pC1M.begin(), pC1M.end());
        mpi::communicator commC1M(Environment::worldComm(), groupC1M);
        M_commsC1M.push_back( commC1M );
        M_dofMgn.insert(std::make_pair(i,dofM));
    }

    // dofM = map<int, dof_points_type>
    int dofSize;
    for(int i=0; i<M_commsC1M.size(); i++)
    {
        std::cout << "Proc " << Environment::worldComm().globalRank()
                  << " makes " << i << "th group loop" << std::endl;
        if( M_commsC1M[i] ) // Current communicator
        {
            if(M_commsC1M[i].rank() == 0 ) // Proc of rank 0 is the Omega_M part
            {
                for ( size_type dof_id = 0; dof_id < M_XhMgn->nLocalDofWithGhost() ; ++dof_id )
                {
                    auto dofpoint = M_XhMgn->dof()->dofPoint(dof_id);
                    M_dofMgn.at(i).push_back( dofpoint );
                }
                dofSize = M_dofMgn.at(i).size();

                std::cout << "Proc : " << Environment::worldComm().globalRank()
                          << " has " << dofSize << " dofs in Omega_M" << std::endl;
            }

            mpi::broadcast( M_commsC1M[i], dofSize, 0);
            M_dofMgn.at(i).resize( dofSize );
            mpi::broadcast( M_commsC1M[i], M_dofMgn.at(i).data(), dofSize, 0);
        }
    }
    toc("setup comms");
}

template<typename te_rb_model_type>
void BiotSavartCRB<te_rb_model_type>::setDimensions()
{
    int mMax = M_teCrbModel->mMaxSigma(); // eim_sigma->mMax()
    int nMax = M_crb->dimension();
    int M = ioption("biotsavart.eim-dimension");
    int N = ioption("biotsavart.crb-dimension");
    M_N = N < nMax ? N : nMax;
    M_M = M < mMax ? M : mMax;
    Feel::cout << "Computing Biot&Savart using "
               << M_N << " reduced basis and "
               << M_M << " EIM basis." << std::endl;
}

template<typename te_rb_model_type>
void BiotSavartCRB<te_rb_model_type>::loadIntegrals()
{
    if ( boption("biotsavart.rebuild-database") )
    {
        this->computeIntegrals();
        return;
    }

    fs::path dbDir(M_crb->dbDirectory());
    fs::path bs(soption("biotsavart.path-to-database"));
    fs::path bsDbDir = dbDir / bs;
    fs::path filename(boost::str( boost::format("BS_integrals_MND_p%1%.db")
                                  % Environment::rank() ));
    fs::path db = bsDbDir / filename;
    if ( !fs::exists( db ) )
    {
        this->computeIntegrals();
        return;
    }

    fs::ifstream ifs( db );
    if ( ifs )
    {
        int M, N;
        tic();
        boost::archive::binary_iarchive ia( ifs );
        ia >> M;
        ia >> N;
        ia >> M_intMND;

        if ( M > M_M )
            M_intMND.resize(M_M);

        for ( int m = 0; m < M_intMND.size(); ++m )
        {
            if ( M_intMND[m].size() > 0 )
            {
                if ( N > M_N )
                    M_intMND[m].conservativeResize(Eigen::NoChange, M_N);
            }
        }

        if ( M < M_M || N < M_N )
            this->computeIntegrals(M, N);

        toc("load integrals");
        return;
    }
    else
    {
        this->computeIntegrals();
        return;
    }
}

template<typename te_rb_model_type>
void BiotSavartCRB<te_rb_model_type>::saveIntegrals()
{
    tic();
    fs::path dbDir(M_crb->dbDirectory());
    fs::path bs(soption("biotsavart.path-to-database"));
    fs::path bsDbDir = dbDir / bs;
    fs::path filename(boost::str( boost::format("BS_integrals_MND_p%1%.db")
                                  % Environment::rank() ));
    fs::path db = bsDbDir / filename;
    if ( Environment::isMasterRank() )
    {
        if ( !fs::exists( bsDbDir ) )
            fs::create_directories( bsDbDir );
    }
    Environment::worldComm().barrier();

    fs::ofstream ofs( db );
    if ( ofs )
    {
        boost::archive::binary_oarchive oa( ofs );
        oa << M_M;
        oa << M_N;
        oa << M_intMND;
    }
    toc("save integrals");
}

template<typename te_rb_model_type>
void BiotSavartCRB<te_rb_model_type>::computeIntegrals(int M, int N)
{
    tic();
    auto coeff = 1/(4*M_PI);
    auto mu0 = 4*M_PI*1e-7; //SI unit : H.m-1 = m.kg.s-2.A-2

    Feel::cout << "Computing " << (M_M-M)*(M_N-N) << " integrals for "
               << M_XhMgn->nDof() << " dofs distributed on " << M_dofMgn.size()
               << " procs on a domain with:"<< std::endl;
    std::cout << "Proc[" << Environment::rank() << "] " << M_Xh->nLocalDof()
              << " dofs" << std::endl;

    std::vector<double> intLocD, intSumD;
    if ( M == 0 )
        M_intMND = std::vector<Eigen::MatrixXd>(M_M);
    else
        M_intMND.resize(M_M);

    for( int i = 0; i < M_dofMgn.size(); ++i)
    {
        int dofSize = M_dofMgn.at(i).size();
        intLocD.resize(dofSize);
        intSumD.resize(dofSize);

        for( int m = 0; m < M_M; ++m)
        {
            if( M_commsC1M[i].rank() == 0 )
            {
                if ( N == 0 )
                    M_intMND[m] = Eigen::MatrixXd::Zero(dofSize, M_N);
                else
                    M_intMND[m].conservativeResize(Eigen::NoChange, M_N);
            }

            // retrieve qm
            auto qm = M_teCrbModel->eimSigmaQ(m);

            for( int n = 0; n < M_N; ++n)
            {
                tic();
                // if already computed, skip
                if ( m < M && n < N )
                    continue;

                // retrieve xi_n
                element_type xiVT_n = M_crbModel->rBFunctionSpace()->primalBasisElement(n);
                auto xi_n = xiVT_n.template element<0>();

                // we are in the communicator and we have some dof in Omega_cond
                if( M_commsC1M[i] && M_Xh->nLocalDof() > 0 )
                {
                    std::vector<Eigen::Matrix<double,3,1>> coords( dofSize/3 );
                    // fill vector of coordinates
                    for(int d = 0; d < dofSize; d+= 3)
                    {
                        auto dof_coord = M_dofMgn.at(i)[d].template get<0>();
                        Eigen::Matrix<double,3,1> coord;
                        coord << dof_coord[0],dof_coord[1],dof_coord[2];
                        coords[d/3] = coord;
                    }

                    auto dist = inner( _e1v-P(),_e1v-P(),
                                       mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>());
                    auto mgnField
                        = integrate(_range=elements( M_Xh->mesh() ),
                                    _expr=-mu0*coeff*idv(qm)*cross(trans(gradv(xi_n)),
                                                                   _e1v-P() )/(dist*dist*dist),
                                    _quad=_Q<1>() ).template evaluate(coords);
                    for( int d = 0; d < dofSize; ++d)
                    {
                        auto dofComp = M_dofMgn.at(i)[d].template get<2>();
                        intLocD[d] = mgnField[d/3](dofComp,0);
                    }
                }
                if( M_commsC1M[i] )
                    boost::mpi::reduce(M_commsC1M[i], intLocD.data(),
                                       dofSize, intSumD.data(), std::plus<double>(), 0);
                if( M_commsC1M[i].rank() == 0 )
                    M_intMND[m].col(n) = Eigen::Map<Eigen::VectorXd>(intSumD.data(), dofSize);
                double time = toc("Done computing integral", false );
                Feel::cout << "[Done computing integralN" << n << "M" << m << "] Time : "
                           << time << "s" << std::endl;
            }
        }
    }
    toc("compute integrals");
    this->saveIntegrals();
}

template<typename te_rb_model_type>
typename BiotSavartCRB<te_rb_model_type>::parameter_type
BiotSavartCRB<te_rb_model_type>::paramFromOption()
{
    std::vector<double> muOpt = vdoption("biotsavart.param");
    auto paramSpace = M_crbModel->parameterSpace();
    int paramSize = paramSpace->dimension();
    if ( muOpt.size() != paramSize )
    {
        Feel::cout << "ERROR: The number of parameter from 'biotsavart.param' is different than from the parameter space: "
                   << paramSize << ". Stopping here!" << std::endl;
        boost::mpi::environment::abort(1);
    }
    auto muMin = paramSpace->min();
    auto muMax = paramSpace->max();
    for ( int i = 0; i < paramSize; ++i)
    {
        if( muOpt[i] < muMin(i) || muOpt[i] > muMax(i) )
        {
            Feel::cout << "WARNING: The " << i << "-th parameter is not in the right range! "
                       << "Using min: " << muMin(i) << " instead." << std::endl;
            muOpt[i] = muMin(i);
        }
    }
    auto mu = paramSpace->element();
    for ( int i = 0; i < paramSize; ++i)
        mu(i) = muOpt[i];
    return mu;
}

template<typename te_rb_model_type>
void BiotSavartCRB<te_rb_model_type>::online( parameter_type & mu )
{
    tic();

    this->computeUn(mu, M_N);
    this->computeB(M_uN);

    toc("online");
}

template<typename te_rb_model_type>
void BiotSavartCRB<te_rb_model_type>::computeUn( parameter_type & mu, int N )
{
    tic();
    double online_tol = doption(_name="crb.online-tolerance");
    vectorN_type time_crb;
    auto o = M_crb->run( mu, time_crb, online_tol, N);
    toc("crb run");

    auto solutions = o.template get<2>();
    M_uN = solutions.template get<0>()[0];
    M_VT = M_crb->expansion( M_uN, M_uN.size() );
    M_betaMu = M_teCrbModel->eimSigmaBeta(mu);
}

template<typename te_rb_model_type>
void BiotSavartCRB<te_rb_model_type>::computeB( vectorN_type & uN )
{
    tic();
    Eigen::VectorXd results;

    M_B = M_XhMgn->element();

    for( int i = 0; i < M_dofMgn.size(); ++i)
    {
        int dofSize = M_dofMgn.at(i).size();
        results = Eigen::VectorXd::Zero(dofSize);
        if( M_commsC1M[i].rank() == 0 )
        {
            for( int m = 0; m < M_M; ++m)
            {
                auto betaM = M_betaMu(m);
                results += betaM*(M_intMND[m]*uN);
            }
            for( int d = 0; d < dofSize; ++d)
                M_B.set(M_dofMgn.at(i)[d].template get<1>(), results(d));
        }
    }
    M_B.close();
    toc("compute B");
}

template<typename te_rb_model_type>
void BiotSavartCRB<te_rb_model_type>::computeFE( parameter_type & mu )
{
    tic();

    tic();
    auto Vh = current_space_type::New( M_meshCond );
    M_j = Vh->element();
    M_VTFe = M_teCrbModel->solve(mu);
    M_teCrbModel->template computeTruthCurrentDensity( M_j, mu, M_VTFe);
    toc("compute j FE");

    tic();
    M_BFe = M_XhMgn->element();
    auto coeff = 1/(4*M_PI);
    auto mu0 = 4*M_PI*1e-7; //SI unit : H.m-1 = m.kg.s-2.A-2

    Feel::cout << "Computing integrals for "
               << M_XhMgn->nDof() << " dofs distributed on " << M_dofMgn.size()
               << " procs on a domain with:"<< std::endl;
    std::cout << "Proc[" << Environment::rank() << "] " << M_Xh->nLocalDof()
              << " dofs" << std::endl;

    std::vector<double> intLocD, intSumD;
    for(int i = 0; i < M_dofMgn.size(); ++i)
    {
        if( M_commsC1M[i] )
        {
            int dofSize = M_dofMgn.at(i).size();
            intLocD.resize( dofSize );
            intSumD.resize( dofSize );

            tic();
            if( M_Xh->nLocalDof() > 0 )
            {
                std::vector<Eigen::Matrix<double,3,1>> coords( dofSize/3 );
                for( int d = 0; d < dofSize; d+= 3 )
                {
                    auto dof_coord = M_dofMgn.at(i)[d].template get<0>();
                    Eigen::Matrix<double,3,1> coord;
                    coord << dof_coord[0],dof_coord[1],dof_coord[2];
                    coords[d/3] = coord;
                }

                auto dist = inner( _e1v-P(),_e1v-P(),
                                   mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
                auto mgnField = integrate(_range = elements( M_Xh->mesh() ),
                                          _expr=mu0*coeff*cross(idv(M_j), _e1v-P())/(dist*dist*dist),
                                          _quad=_Q<1>()
                    ).template evaluate(coords);

                for( int d = 0; d < dofSize; d++)
                {
                    auto dofComp = M_dofMgn.at(i)[d].template get<2>();
                    intLocD[d] = mgnField[d/3](dofComp,0);
                }
            }
            toc("integral_computation");

            tic();
            mpi::reduce(M_commsC1M[i], intLocD.data(), dofSize, intSumD.data(), std::plus<double>(), 0);
            toc("reduce_contributions_of_integrals");

            if( M_commsC1M[i].rank() == 0 )
            {
                for(int d=0; d<dofSize; d++)
                    M_BFe.set(M_dofMgn.at(i)[d].template get<1>(), intSumD[d]);
            }
        }
        Environment::worldComm().barrier();
    }

    M_BFe.close();
    toc("compute BiotSavart FE");

    toc("compute FE");
}

template<typename te_rb_model_type>
void BiotSavartCRB<te_rb_model_type>::exportResults()
{
    exporter_ptrtype eC( exporter_type::New( "conductor") );
    eC->setMesh(M_Xh->mesh());
    exporter_ptrtype eM( exporter_type::New( "mgn") );
    eM->setMesh(M_XhMgn->mesh());

    if(boption("biotsavart.compute-offline") && boption("biotsavart.compute-online") )
    {
        auto WN = M_crb->wn();
        element_type VT = M_crb->expansion( M_uN, M_uN.size() );
        auto V = VT.template element<0>();
        auto T = VT.template element<1>();

        auto Jh = current_space_type::New( M_meshCond );
        auto j = Jh->element();
        for( int n = 0; n < M_uN.size(); ++n )
        {
            element_type xiVT_n = M_crbModel->rBFunctionSpace()->primalBasisElement(n);
            auto xi_n = xiVT_n.template element<0>();
            for( int m = 0; m < M_betaMu.size(); ++m )
            {
                auto qm = M_teCrbModel->eimSigmaQ(m);
                j += vf::project( _space=Jh, _range=elements(M_Xh->mesh()), _expr=M_betaMu(m)*M_uN(n)*idv(qm)*trans(gradv(xi_n)) );
            }
        }

        eC->add( "V", V);
        eC->add( "T", T);
        eC->add( "j", j);
        eM->add( "B", M_B);
    }
    if( boption("biotsavart.compute-fe") )
    {
        eC->add( "j_fe", M_j);
        eM->add( "B_fe", M_BFe );
    }

    eC->save();
    eM->save();
}

template<typename te_rb_model_type>
typename BiotSavartCRB<te_rb_model_type>::value_type
BiotSavartCRB<te_rb_model_type>::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve )
{
    
}

using biotsavartcrbthermoelectric = BiotSavartCRB<ThermoElectric>;

template class FEELPP_EXPORT BiotSavartCRB<ThermoElectric>;

FEELPP_CRB_PLUGIN( biotsavartcrbthermoelectric, biotsavartcrbthermoelectric )
//FEELPP_CRB_PLUGIN( BiotSavartCRB<NLThermoelectric>, biotsavartcrbnlthermoelectric )

}

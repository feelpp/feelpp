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
#define FEELPP_INSTANTIATE_BIOTSAVART_ALPHA_ELECTRIC
#include <feel/feelcrb/crbplugin.hpp>


#include "biotsavartalphaelectric.hpp"
#include <alphaelectric.hpp>

namespace Feel {

template<typename te_rb_model_type>
BiotSavartAlphaElectroCRB<te_rb_model_type>::BiotSavartAlphaElectroCRB(crb::stage stage)
    : super_type("biotsavart")
{
    // if( stage == crb::stage::online )
    // {
    //     M_teCrbModel = boost::make_shared<te_rb_model_type>();
    //     M_crbModel = boost::make_shared<crb_model_type>(M_teCrbModel, stage);
    //     M_crb = boost::make_shared<crb_type>("biotsavartalphaelectro_crb", M_crbModel);
    //     M_crb->loadDBLast(crb::last::modified, crb::load::all);
    //     M_XhCond = M_teCrbModel->functionSpace();
    // }
}

template<typename te_rb_model_type>
void BiotSavartAlphaElectroCRB<te_rb_model_type>::initModel()
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
    M_teCrbModel = boost::make_shared<te_rb_model_type>(M_meshCond);
    M_crbModel = boost::make_shared<crb_model_type>(M_teCrbModel, crb::stage::offline);
    M_crb = boost::make_shared<crb_type>("biotsavartalphaelectro_crb", M_crbModel);
    toc("constructor + eim");
    tic();
    M_crb->offline();
    toc("rb construction");
    Feel::cout << tc::green << "Construction of DEIM, MDEIM and CRB finished!!"
               << tc::reset <<std::endl;

    auto Pset = this->parameterSpace()->sampling();
    int N = ioption("biotsavart.trainset-deim-size");
    std::string supersampling = (boost::format("DmuDeimBS-N%1%") % N ).str();
    std::ifstream file( supersampling);
    bool all_proc_same_sampling=true;
    if( ! file )
    {
        Pset->randomize( N , all_proc_same_sampling , supersampling );
        Pset->writeOnFile( supersampling );
    }
    else
    {
        Pset->clear();
        Pset->readFromFile(supersampling);
    }

    M_XhCond = M_teCrbModel->functionSpace();
    this->setFunctionSpaces( space_type::New( M_meshMgn) );
    Feel::cout << "Conductor nDof = " << M_XhCond->nDof()
               << " Box nDof = " << this->Xh->nDof() << std::endl;
    this->setupCommunicatorsBS();

    M_deim = deim( _model=boost::dynamic_pointer_cast<self_type>(this->shared_from_this()),
                   _sampling=Pset, _prefix="bs");
    this->M_deim->run();
    Feel::cout << tc::green << "Construction of BiotSavart DEIM finished!!"
               << tc::reset << std::endl;
    auto online_model = M_deim->M_online_model;
    online_model->setupCRB(M_crb);
}

template<typename te_rb_model_type>
void BiotSavartAlphaElectroCRB<te_rb_model_type>::setupCRB( crb_ptrtype crb )
{
    M_crb = crb;
    M_crbModel = crb->model();
    M_teCrbModel = M_crbModel->model();
    M_XhCond = M_teCrbModel->functionSpace();
    Feel::cout << "Online model : Conductor  nDof = " << M_XhCond->nDof()
               << " Box nDof = " << this->Xh->nDof() << std::endl;
    this->setupCommunicatorsBS();
}

template<typename te_rb_model_type>
void BiotSavartAlphaElectroCRB<te_rb_model_type>::runBS()
{
    M_mu = this->paramFromOption();
    if( boption("biotsavart.compute-offline") )
    {
        if( boption("biotsavart.compute-online") )
            this->online(M_mu);
    }

    if( boption("biotsavart.compute-fe") )
        this->computeFE(M_mu);

    if( boption("biotsavart.compute-offline")
        && boption("biotsavart.compute-online")
        && boption("biotsavart.compute-fe") )
        this->computeErrors();

    this->exportResults(M_mu);
}

template<typename te_rb_model_type>
typename BiotSavartAlphaElectroCRB<te_rb_model_type>::vector_ptrtype
BiotSavartAlphaElectroCRB<te_rb_model_type>::assembleForDEIM( parameter_type const& mu )
{
    tic();
    double online_tol = doption(_name="crb.online-tolerance");
    vectorN_type time_crb;
    auto o = M_crb->run( mu, time_crb, online_tol);

    auto solutions = o.template get<2>();
    M_uN = solutions.template get<0>()[0];
    this->expand();

    vector_ptrtype Bvec = backend()->newVector( this->Xh );
    auto B = this->Xh->element( Bvec );

    auto alphaExpr = expr(M_teCrbModel->alpha(mu));
    auto alphaPrimeExpr = expr(M_teCrbModel->alphaPrime(mu));

    auto psi = vec( cos(alphaExpr)*Px() + sin(alphaExpr)*Py(),
                    -sin(alphaExpr)*Px() + cos(alphaExpr)*Py(),
                    Pz() );
    auto Jinv = mat<3,3>( cos(alphaExpr), -sin(alphaExpr), -alphaPrimeExpr*Py(),
                          sin(alphaExpr), cos(alphaExpr), alphaPrimeExpr*Px(),
                          cst(0.), cst(0.), cst(1.) );

    auto coeff = 1/(4*M_PI);
    auto mu0 = 4*M_PI*1e-4; //SI unit : H.m-1 = m.kg.s-2.A-2
    auto sigma = M_teCrbModel->sigma();

    std::vector<double> intLocD, intSumD;
    for(int i = 0; i < M_dofMgn.size(); ++i)
    {
        if( M_commsC1M[i] )
        {
            int dofSize = M_dofMgn.at(i).size();
            intLocD.resize( dofSize );
            intSumD.resize( dofSize );

            if( M_XhCond->nLocalDof() > 0 )
            {
                std::vector<Eigen::Matrix<double,3,1>> coords( dofSize/3 );
                for( int d = 0; d < dofSize; d+= 3 )
                {
                    auto dofCoord = M_dofMgn.at(i)[d].template get<0>();
                    Eigen::Matrix<double,3,1> coord;
                    coord << dofCoord[0], dofCoord[1], dofCoord[2];
                    coords[d/3] = coord;
                }

                auto dist = inner( _e1v-psi, _e1v-psi,
                                   mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
                auto mgnField = integrate(_range=elements( M_XhCond->mesh() ),
                                          _expr=-mu0*coeff*sigma*cross(trans(gradv(M_V)*Jinv),
                                                                       _e1v-psi)/(dist*dist*dist),
                                          _quad=_Q<1>()
                                          ).template evaluate(coords);

                for( int d = 0; d < dofSize; ++d )
                {
                    auto dofComp = M_dofMgn.at(i)[d].template get<2>();
                    intLocD[d] = mgnField[d/3](dofComp,0);
                }
            }

            mpi::reduce(M_commsC1M[i], intLocD.data(), dofSize, intSumD.data(), std::plus<double>(), 0);

            if( M_commsC1M[i].rank() == 0 )
            {
                for(int d=0; d<dofSize; d++)
                    B.set(M_dofMgn.at(i)[d].template get<1>(), intSumD[d]);
            }
        }
        Environment::worldComm().barrier();
    }

    B.close();

    toc("assembleBiotSavart");

    return Bvec;
}

template<typename te_rb_model_type>
void BiotSavartAlphaElectroCRB<te_rb_model_type>::setupCommunicatorsBS()
{
    tic();
    M_dofMgn.clear();
    M_commsC1M.clear();
    // Identify repartition of dofs (Omega_C, Omega_M)
    std::vector<int> isIn( 2, 0 );
    if( M_XhCond->nLocalDof()  > 0 ) // current proc has dofs in Omega_C
        isIn[0] = 1;
    if( this->Xh->nLocalDof() > 0 ) // current proc has dofs in Omega_M
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
        // std::cout << "Proc " << Environment::worldComm().globalRank()
        //           << " makes " << i << "th group loop" << std::endl;
        if( M_commsC1M[i] ) // Current communicator
        {
            if(M_commsC1M[i].rank() == 0 ) // Proc of rank 0 is the Omega_M part
            {
                for ( size_type dof_id = 0; dof_id < this->Xh->nLocalDofWithGhost() ; ++dof_id )
                {
                    auto dofpoint = this->Xh->dof()->dofPoint(dof_id);
                    M_dofMgn.at(i).push_back( dofpoint );
                }
                dofSize = M_dofMgn.at(i).size();

                // std::cout << "Proc : " << Environment::worldComm().globalRank()
                //           << " has " << dofSize << " dofs in Omega_M" << std::endl;
            }

            mpi::broadcast( M_commsC1M[i], dofSize, 0);
            M_dofMgn.at(i).resize( dofSize );
            mpi::broadcast( M_commsC1M[i], M_dofMgn.at(i).data(), dofSize, 0);
        }
    }
    toc("setup comms");
}

template<typename te_rb_model_type>
typename BiotSavartAlphaElectroCRB<te_rb_model_type>::parameter_type
BiotSavartAlphaElectroCRB<te_rb_model_type>::paramFromOption()
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
typename BiotSavartAlphaElectroCRB<te_rb_model_type>::parameter_type
BiotSavartAlphaElectroCRB<te_rb_model_type>::param0()
{
    std::vector<double> x(this->nbParameters(), 0.);
    auto mu = this->newParameter();
    mu = Eigen::Map<Eigen::VectorXd const>( x.data(), mu.size());
    return mu;
}

template<typename te_rb_model_type>
void BiotSavartAlphaElectroCRB<te_rb_model_type>::online( parameter_type & mu, int M )
{
    tic();

    tic();
    double online_tol = doption(_name="crb.online-tolerance");
    vectorN_type time_crb;
    auto o = M_crb->run( mu, time_crb, online_tol, M_N);
    toc("crb run");

    auto solutions = o.template get<2>();
    M_uN = solutions.template get<0>()[0];
    this->expand();

    auto beta = M_deim->beta(mu, M);
    auto q = M_deim->q();
    vector_ptrtype Bvec = backend()->newVector( this->Xh );
    Bvec->zero();
    for( int m = 0; m < beta.size(); ++m )
        Bvec->add( beta(m), q[m]);

    Bvec->close();
    M_B = this->Xh->element( Bvec );

    toc("online");
}

template<typename te_rb_model_type>
void BiotSavartAlphaElectroCRB<te_rb_model_type>::expand(int N )
{
    M_V = M_crb->expansion( M_uN, N);
}

template<typename te_rb_model_type>
void BiotSavartAlphaElectroCRB<te_rb_model_type>::computeFE( parameter_type & mu )
{
    tic();

    tic();
    M_VFe = M_teCrbModel->solve(mu);
    toc("compute j FE");

    auto alphaExpr = expr(M_teCrbModel->alpha(mu));
    auto alphaPrimeExpr = expr(M_teCrbModel->alphaPrime(mu));

    auto psi = vec( cos(alphaExpr)*Px() + sin(alphaExpr)*Py(),
                    -sin(alphaExpr)*Px() + cos(alphaExpr)*Py(),
                    Pz() );
    auto Jinv = mat<3,3>( cos(alphaExpr), -sin(alphaExpr), -alphaPrimeExpr*Py(),
                          sin(alphaExpr), cos(alphaExpr), alphaPrimeExpr*Px(),
                          cst(0.), cst(0.), cst(1.) );
    tic();
    M_BFe = this->Xh->element();
    auto coeff = 1/(4*M_PI);
    auto mu0 = 4*M_PI*1e-4; //SI unit : H.m-1 = m.kg.s-2.A-2
    auto sigma = M_teCrbModel->sigma();

    Feel::cout << "Computing integrals for "
               << this->Xh->nDof() << " dofs distributed on " << M_dofMgn.size()
               << " procs on a domain with:"<< std::endl;
    std::cout << "Proc[" << Environment::rank() << "] " << M_XhCond->nLocalDof()
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
            if( M_XhCond->nLocalDof() > 0 )
            {
                std::vector<Eigen::Matrix<double,3,1>> coords( dofSize/3 );
                for( int d = 0; d < dofSize; d+= 3 )
                {
                    auto dofCoord = M_dofMgn.at(i)[d].template get<0>();
                    Eigen::Matrix<double,3,1> coord;
                    coord << dofCoord[0], dofCoord[1], dofCoord[2];
                    coords[d/3] = coord;
                }

                auto dist = inner( _e1v-psi,_e1v-psi,
                                   mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
                auto mgnField = integrate(_range = elements( M_XhCond->mesh() ),
                                          _expr=-mu0*coeff*sigma*cross(trans(gradv(M_VFe)*Jinv),
                                                                       _e1v-psi)/(dist*dist*dist),
                                          _quad=_Q<1>()
                                          ).template evaluate(coords);

                for( int d = 0; d < dofSize; ++d )
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
typename BiotSavartAlphaElectroCRB<te_rb_model_type>::cond_element_type
BiotSavartAlphaElectroCRB<te_rb_model_type>::alpha( parameter_type const& mu )
{
    auto alpha = expr(M_teCrbModel->alpha(mu));
    auto a = M_XhCond->element(alpha);
    return a;
}

template<typename te_rb_model_type>
std::vector<double> BiotSavartAlphaElectroCRB<te_rb_model_type>::computeErrors()
{
    std::vector<double> err(4,0);
    double normVFe = normL2( elements(M_meshCond), idv(M_VFe) );
    double normBFe = normL2( elements(M_meshMgn), idv(M_BFe) );
    err[0] = normL2( elements(M_meshCond), idv(M_V)-idv(M_VFe) );
    err[1] = err[0]/normVFe;
    err[2] = normL2( elements(M_meshMgn), idv(M_B)-idv(M_BFe) );
    err[3] = err[2]/normBFe;
    Feel::cout << "ErrV: " << err[0] << " relative: " << err[1] << std::endl
               << "ErrB: " << err[2] << " relative: " << err[3] << std::endl;
    return err;
}

template<typename te_rb_model_type>
void BiotSavartAlphaElectroCRB<te_rb_model_type>::exportResults(parameter_type const& mu)
{
    exporter_ptrtype eC( exporter_type::New( "conductor") );
    eC->setMesh(M_XhCond->mesh());
    exporter_ptrtype eM( exporter_type::New( "mgn") );
    eM->setMesh(this->Xh->mesh());

    auto a = this->alpha(mu);
    eC->add("alpha", a);

    if( boption("biotsavart.compute-offline") && boption("biotsavart.compute-online") )
    {
        // auto WN = M_crb->wn();
        // cond_element_type V = M_crb->expansion( M_uN, M_uN.size() );

        // auto Jh = current_space_type::New( M_meshCond );
        // auto j = Jh->element();
        // for( int n = 0; n < M_uN.size(); ++n )
        // {
        //     element_type xiVT_n = M_crbModel->rBFunctionSpace()->primalBasisElement(n);
        //     auto xi_n = xiVT_n.template element<0>();
        //     for( int m = 0; m < M_betaMu.size(); ++m )
        //     {
        //         auto qm = M_teCrbModel->eimSigmaQ(m);
        //         j += vf::project( _space=Jh, _range=elements(M_XhCond->mesh()), _expr=M_betaMu(m)*M_uN(n)*idv(qm)*trans(gradv(xi_n)) );
        //     }
        // }

        eC->add( "V", M_V);
        // eC->add( "j", j);
        eM->add( "B", M_B);
    }
    if( boption("biotsavart.compute-fe") )
    {
        eC->add( "V_fe", M_VFe);
        // eC->add( "j_fe", M_j);
        eM->add( "B_fe", M_BFe );
    }

    eC->save();
    eM->save();
}

template<typename te_rb_model_type>
double BiotSavartAlphaElectroCRB<te_rb_model_type>::homogeneity( element_type& B)
{
    auto m = minmax( _range=elements(M_meshMgn), _pset=_Q<2>(),
                     _expr=abs(idv(B)(2,0)) );
    return m.max()/m.min() - 1;
}

template<typename te_rb_model_type>
typename BiotSavartAlphaElectroCRB<te_rb_model_type>::value_type
BiotSavartAlphaElectroCRB<te_rb_model_type>::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve )
{
    return 0;
}

using alphabiotsavartcrbelectric = BiotSavartAlphaElectroCRB<AlphaElectric>;

template class FEELPP_EXPORT BiotSavartAlphaElectroCRB<AlphaElectric>;

FEELPP_CRB_PLUGIN( alphabiotsavartcrbelectric, "alphabiotsavartcrbelectric")

}

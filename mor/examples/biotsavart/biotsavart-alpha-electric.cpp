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
#include <feel/feelmor/crbplugin.hpp>


#include "biotsavart-alpha-electric.hpp"
#include <electric-alpha.hpp>

namespace Feel {

template<typename te_rb_model_type>
BiotSavartAlphaElectricCRB<te_rb_model_type>::BiotSavartAlphaElectricCRB(crb::stage stage)
    : super_type("biotsavart"),
      M_propertyPath(Environment::expand(soption("biotsavart.filename"))),
      M_repart(boption("biotsavart.repart")),
      M_computeFe(boption("biotsavart.compute-fe")),
      M_computeOffline(boption("biotsavart.compute-offline")),
      M_computeOnline(boption("biotsavart.compute-online")),
      M_rebuildDb(boption("biotsavart.rebuild-database")),
      M_pathToDb(soption("biotsavart.path-to-database")),
      M_trainsetDeimSize(ioption("biotsavart.trainset-deim-size"))
{
    // if( stage == crb::stage::online )
    // {
    //     M_teCrbModel = std::make_shared<te_rb_model_type>();
    //     M_crbModel = std::make_shared<crb_model_type>(M_teCrbModel, stage);
    //     M_crb = std::make_shared<crb_type>("biotsavartalphaelectro_crb", M_crbModel);
    //     M_crb->loadDBLast(crb::last::modified, crb::load::all);
    //     M_XhCond = M_teCrbModel->functionSpace();
    // }
}

template<typename te_rb_model_type>
void BiotSavartAlphaElectricCRB<te_rb_model_type>::initModel()
{
    M_modelProps = std::make_shared<prop_type>(M_propertyPath);
    M_materials = M_modelProps->materials().materialWithPhysic("magnetic");

    std::string mshFileName = Environment::expand(soption("gmsh.filename"));

    auto M_mesh = loadMesh( _mesh=new mesh_type,
                            _filename=mshFileName,
                            _rebuild_partitions=M_repart,
                            _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES);

    tic();
    M_teCrbModel = std::make_shared<te_rb_model_type>(M_mesh);
    M_crbModel = std::make_shared<crb_model_type>(M_teCrbModel, crb::stage::offline);
    M_crb = crb_type::New("alphaelectric", M_crbModel, crb::stage::offline);
    toc("constructor + eim");

    tic();
    M_crb->offline();
    toc("rb construction");
    Feel::cout << tc::green << "Construction of DEIM, MDEIM and CRB finished!!"
               << tc::reset <<std::endl;

    auto Pset = this->parameterSpace()->sampling();
    int N = M_trainsetDeimSize;
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
    std::vector<std::string> range;
    for( auto const& mp : M_materials )
        range.push_back(mp.first);
    auto domain = markedelements(M_mesh, range);
    this->setFunctionSpaces(space_type::New( _mesh=M_mesh, _range=domain ) );
    Feel::cout << "Conductor nDof = " << M_XhCond->nDof()
               << " Box nDof = " << this->Xh->nDof() << std::endl;
    this->setupCommunicatorsBS();

    M_deim = Feel::deim( _model=std::dynamic_pointer_cast<self_type>(this->shared_from_this()),
                         _sampling=Pset, _prefix="bs");
    M_deim->run();
    Feel::cout << tc::green << "Construction of BiotSavart DEIM finished!!"
               << tc::reset << std::endl;
    auto online_model = M_deim->onlineModel();
    online_model->setupCRB(M_crb);
}

template<typename te_rb_model_type>
void BiotSavartAlphaElectricCRB<te_rb_model_type>::setupCRB( crb_ptrtype crb )
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
void BiotSavartAlphaElectricCRB<te_rb_model_type>::runBS()
{
    M_mu = this->paramFromProperties();
    if( M_computeOffline )
    {
        if( M_computeOnline )
            this->online(M_mu);
    }

    if( M_computeFe )
        this->computeFE(M_mu);

    if( M_computeOffline && M_computeOnline && M_computeFe )
        this->computeErrors();

    this->exportResults(M_mu);
}

template<typename te_rb_model_type>
typename BiotSavartAlphaElectricCRB<te_rb_model_type>::vector_ptrtype
BiotSavartAlphaElectricCRB<te_rb_model_type>::assembleForDEIM( parameter_type const& mu, int const& tag )
{
    tic();
    double online_tol = doption(_name="crb.online-tolerance");
    vectorN_type time_crb;
    auto o = M_crb->run( mu, time_crb, online_tol);

    auto solutions = o.template get<2>();
    M_uN = solutions.template get<0>()[0];
    this->expand();

    auto mesh = M_XhCond->mesh();
    vector_ptrtype Bvec = backend()->newVector( this->Xh );
    auto B = this->Xh->element( Bvec );

    auto coeff = 1/(4*M_PI);
    auto mu0 = 4*M_PI*1e-4; //SI unit : H.m-1 = m.kg.s-2.A-2

    std::vector<double> intLocD, intSumD;
    for(int i = 0; i < M_dofMgn.size(); ++i)
    {
        if( M_commsC1M[i] )
        {
            int dofSize = M_dofMgn.at(i).size();
            intLocD.resize( dofSize, 0 );
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

                std::vector<std::vector<Eigen::Matrix<double,3,1> > > mgnFields(M_teCrbModel->materials().size());
                int j = 0;
                for( auto const& material : M_teCrbModel->materialsWithGeo() )
                {
                    auto sigma = material.second.getDouble("sigma");

                    auto alphaExpr = expr(M_teCrbModel->alpha(mu, material.second));
                    auto alphaPrimeExpr = expr(M_teCrbModel->alphaPrime(mu, material.second));

                    auto Jinv = mat<3,3>( cos(alphaExpr), -sin(alphaExpr), -alphaPrimeExpr*Py(),
                                          sin(alphaExpr), cos(alphaExpr), alphaPrimeExpr*Px(),
                                          cst(0.), cst(0.), cst(1.) );

                    auto psi = vec( cos(alphaExpr)*Px() + sin(alphaExpr)*Py(),
                                    -sin(alphaExpr)*Px() + cos(alphaExpr)*Py(),
                                    Pz() );
                    auto dist = inner( _e1v-psi, _e1v-psi,
                                       mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );

                    mgnFields[j++] = integrate(_range=markedelements( mesh, material.first ),
                                               _expr=-mu0*coeff*sigma*cross(trans(gradv(M_V)*Jinv),
                                                                            _e1v-psi)/(dist*dist*dist),
                                               _quad=_Q<1>()
                                               ).template evaluate(coords);
                }
                for( auto const& material : M_teCrbModel->materialsWithoutGeo() )
                {
                    auto sigma = material.second.getDouble("sigma");
                    auto dist = inner( _e1v-P(), _e1v-P(),
                                       mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
                    mgnFields[j++] = integrate(_range=markedelements( mesh, material.first ),
                                               _expr=-mu0*coeff*sigma*cross(trans(gradv(M_V)),
                                                                            _e1v-P())/(dist*dist*dist),
                                               _quad=_Q<1>()
                                               ).template evaluate(coords);
                }

                for( int d = 0; d < dofSize; ++d )
                {
                    auto dofComp = M_dofMgn.at(i)[d].template get<2>();
                    int j = 0;
                    for( auto const& mat : M_teCrbModel->materials() )
                        intLocD[d] += mgnFields[j][d/3](dofComp,0);
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
void BiotSavartAlphaElectricCRB<te_rb_model_type>::setupCommunicatorsBS()
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
typename BiotSavartAlphaElectricCRB<te_rb_model_type>::parameter_type
BiotSavartAlphaElectricCRB<te_rb_model_type>::paramFromOption()
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
typename BiotSavartAlphaElectricCRB<te_rb_model_type>::parameter_type
BiotSavartAlphaElectricCRB<te_rb_model_type>::param0()
{
    std::vector<double> x(this->nbParameters(), 0.);
    auto mu = this->newParameter();
    mu = Eigen::Map<Eigen::VectorXd const>( x.data(), mu.size());
    return mu;
}

template<typename te_rb_model_type>
typename BiotSavartAlphaElectricCRB<te_rb_model_type>::parameter_type
BiotSavartAlphaElectricCRB<te_rb_model_type>::paramFromProperties() const
{
    auto mu = this->newParameter();
    int i = 0;
    auto parameters = M_modelProps->parameters();
    for( auto const& parameterPair : parameters )
    {
        if( parameterPair.second.hasMinMax() )
        {
            mu(i++) = parameterPair.second.value();
        }
    }
    return mu;
}

template<typename te_rb_model_type>
void BiotSavartAlphaElectricCRB<te_rb_model_type>::online( parameter_type & mu, int M )
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
void BiotSavartAlphaElectricCRB<te_rb_model_type>::expand(int N )
{
    M_V = M_crb->expansion( M_uN, N);
}

template<typename te_rb_model_type>
void BiotSavartAlphaElectricCRB<te_rb_model_type>::computeFE( parameter_type & mu )
{
    tic();

    tic();
    M_VFe = M_teCrbModel->solve(mu);
    toc("compute j FE");

    tic();
    auto mesh = Xh->mesh();
    M_BFe = this->Xh->element();
    auto coeff = 1/(4*M_PI);
    auto mu0 = 4*M_PI*1e-4; //SI unit : H.m-1 = m.kg.s-2.A-2

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
            intLocD.resize( dofSize, 0 );
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

                std::vector<std::vector<Eigen::Matrix<double,3,1> > > mgnFields(M_teCrbModel->materials().size());
                int j = 0;
                for( auto const& material : M_teCrbModel->materialsWithGeo() )
                {
                    auto sigma = material.second.getDouble("sigma");

                    auto alphaExpr = expr(M_teCrbModel->alpha(mu, material.second));
                    auto alphaPrimeExpr = expr(M_teCrbModel->alphaPrime(mu, material.second));

                    auto Jinv = mat<3,3>( cos(alphaExpr), -sin(alphaExpr), -alphaPrimeExpr*Py(),
                                          sin(alphaExpr), cos(alphaExpr), alphaPrimeExpr*Px(),
                                          cst(0.), cst(0.), cst(1.) );

                    auto psi = vec( cos(alphaExpr)*Px() + sin(alphaExpr)*Py(),
                                    -sin(alphaExpr)*Px() + cos(alphaExpr)*Py(),
                                    Pz() );
                    auto dist = inner( _e1v-psi,_e1v-psi,
                                       mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
                    mgnFields[j++] = integrate(_range=markedelements( mesh, material.first ),
                                               _expr=-mu0*coeff*sigma*cross(trans(gradv(M_V)*Jinv),
                                                                            _e1v-psi)/(dist*dist*dist),
                                               _quad=_Q<1>()
                                               ).template evaluate(coords);
                }
                for( auto const& material : M_teCrbModel->materialsWithoutGeo() )
                {
                    auto sigma = material.second.getDouble("sigma");
                    auto dist = inner( _e1v-P(), _e1v-P(),
                                       mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
                    mgnFields[j++] = integrate(_range=markedelements( mesh, material.first ),
                                               _expr=-mu0*coeff*sigma*cross(trans(gradv(M_V)),
                                                                            _e1v-P())/(dist*dist*dist),
                                               _quad=_Q<1>()
                                               ).template evaluate(coords);
                }

                for( int d = 0; d < dofSize; ++d )
                {
                    auto dofComp = M_dofMgn.at(i)[d].template get<2>();
                    int j = 0;
                    for( auto const& mat : M_teCrbModel->materials() )
                        intLocD[d] += mgnFields[j][d/3](dofComp,0);
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
typename BiotSavartAlphaElectricCRB<te_rb_model_type>::cond_element_type
BiotSavartAlphaElectricCRB<te_rb_model_type>::alpha( parameter_type const& mu )
{
    auto a = M_XhCond->element();
    auto mesh = M_XhCond->mesh();
    for( auto const& material : M_teCrbModel->materials() )
        a += project( _space=M_XhCond, _range=markedelements(mesh, material.first),
                      _expr=expr(M_teCrbModel->alpha(mu, material.second)) );
    return a;
}

template<typename te_rb_model_type>
std::vector<double> BiotSavartAlphaElectricCRB<te_rb_model_type>::computeErrors()
{
    std::vector<double> err(4,0);
    auto rangeV = M_VFe.functionSpace()->dof()->meshSupport()->rangeElements();
    auto rangeB = M_BFe.functionSpace()->dof()->meshSupport()->rangeElements();
    double normVFe = normL2( rangeV, idv(M_VFe) );
    double normBFe = normL2( rangeB, idv(M_BFe) );
    err[0] = normL2( rangeV, idv(M_V)-idv(M_VFe) );
    err[1] = err[0]/normVFe;
    err[2] = normL2( rangeB, idv(M_B)-idv(M_BFe) );
    err[3] = err[2]/normBFe;
    Feel::cout << "ErrV: " << err[0] << " relative: " << err[1] << std::endl
               << "ErrB: " << err[2] << " relative: " << err[3] << std::endl;
    return err;
}

template<typename te_rb_model_type>
void BiotSavartAlphaElectricCRB<te_rb_model_type>::exportResults(parameter_type const& mu)
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
double BiotSavartAlphaElectricCRB<te_rb_model_type>::homogeneity( element_type& B)
{
    auto range = B.functionSpace()->dof()->meshSupport()->rangeElements();
    auto m = minmax( _range=range, _pset=_Q<2>(),
                     _expr=abs(idv(B)(2,0)) );
    return m.max()/m.min() - 1;
}

template<typename te_rb_model_type>
typename BiotSavartAlphaElectricCRB<te_rb_model_type>::value_type
BiotSavartAlphaElectricCRB<te_rb_model_type>::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve )
{
    return 0;
}

using alphabiotsavartcrbelectric = BiotSavartAlphaElectricCRB<AlphaElectric>;

template class FEELPP_EXPORT BiotSavartAlphaElectricCRB<AlphaElectric>;

FEELPP_CRB_PLUGIN( alphabiotsavartcrbelectric, alphabiotsavartcrbelectric)

}

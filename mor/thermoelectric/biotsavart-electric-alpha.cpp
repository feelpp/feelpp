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
#include <feel/feelmesh/concatenate.hpp>

#include "biotsavart-electric-alpha.hpp"
#include "electric-alpha.hpp"

namespace Feel {

template<typename te_rb_model_type>
po::options_description
BiotSavartAlphaElectricCRB<te_rb_model_type>::makeOptions( std::string const& prefix )
{
    po::options_description opt("BiotSavart options");
    opt.add_options()
        ( "biotsavart.filename", po::value<std::string>()->default_value("biotsavart.json"),
          "json file")
        ( "biotsavart.trainset-deim-size", po::value<int>()->default_value(10),
          "size of the trainset for DEIM of BiotSavart" )
        ( "biotsavart.db.base", po::value<std::string>()->default_value("alphabiotsavart"), "basename for crb db")
        ( "biotsavart.use-rb-in-deim", po::value<bool>()->default_value(true), "" )
        ( "biotsavart.verbose", po::value<int>()->default_value(0), "" )
        ( "biotsavart.use-eq", po::value<bool>()->default_value(false), "use empirical quadrature" )
        ;
    opt.add(deimOptions("bs"));
    opt.add(crbOptions(prefix));
    opt.add(crbSEROptions(prefix));
    opt.add(eq_options());
    opt.add(AlphaElectric::makeOptions("electric"));
    return opt;
}

template<typename te_rb_model_type>
typename BiotSavartAlphaElectricCRB<te_rb_model_type>::self_ptrtype
BiotSavartAlphaElectricCRB<te_rb_model_type>::New(crb::stage stage, std::string const& prefix)
{
    return std::make_shared<self_type>(stage, prefix);
}

template<typename te_rb_model_type>
BiotSavartAlphaElectricCRB<te_rb_model_type>::BiotSavartAlphaElectricCRB(std::string const& prefix)
    : BiotSavartAlphaElectricCRB<te_rb_model_type>(crb::stage::online, prefix)
{}

template<typename te_rb_model_type>
BiotSavartAlphaElectricCRB<te_rb_model_type>::BiotSavartAlphaElectricCRB(crb::stage stage, std::string const& prefix)
    : super_type(soption("biotsavart.db.base"), Environment::worldCommPtr(), prefix),
      M_propertyPath(Environment::expand(soption("biotsavart.filename"))),
      M_trainsetDeimSize(ioption("biotsavart.trainset-deim-size")),
      M_dbBasename(soption("biotsavart.db.base")),
      M_useRbInDeim(boption("biotsavart.use-rb-in-deim")),
      M_verbose(ioption("biotsavart.verbose")),
      M_useEQ(boption("biotsavart.use-eq"))
{
}

template<typename te_rb_model_type>
void BiotSavartAlphaElectricCRB<te_rb_model_type>::initModel()
{
    M_modelProps = std::make_shared<prop_type>(M_propertyPath);
    M_modelProps->enableBoundaryConditions2();
    M_materials = M_modelProps->materials().materialWithPhysic("magnetic");

    M_mesh = loadMesh( _mesh=new mesh_type,
                       _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES);

    tic();
    M_teCrbModel = std::make_shared<te_rb_model_type>(M_mesh, "electric");
    M_crbModel = std::make_shared<crb_model_type>(M_teCrbModel, crb::stage::offline);
    M_crb = crb_type::New(M_teCrbModel->dbBasename(), M_crbModel, crb::stage::offline);
    toc("constructor + eim", M_verbose > 0);

    tic();
    M_crb->offline();
    M_N = M_crb->dimension();
    toc("rb construction", M_verbose > 0);
    Feel::cout << tc::green << "Construction of DEIM, MDEIM and CRB finished!!"
               << tc::reset <<std::endl;

    auto Pset = this->parameterSpace()->sampling();
    int N = M_trainsetDeimSize;
    std::string supersampling = (boost::format("DmuDeimBS-P%1%-N%2%") % this->parameterSpace()->dimension() % N ).str();
    std::ifstream file( supersampling);
    if( ! file )
    {
        bool all_proc_same_sampling=true;
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

    auto d = Feel::deim( _model=std::dynamic_pointer_cast<self_type>(this->shared_from_this()),
                         _sampling=Pset, _prefix="bs");
    this->addDeim(d);
    this->deim()->run();
    Feel::cout << tc::green << "Construction of BiotSavart DEIM finished!!"
               << tc::reset << std::endl;

    auto onlineModel = d->onlineModel();
    onlineModel->setupCRB(M_crb);
    onlineModel->setIndices(d->indexR());
    if( M_useEQ )
    {
        onlineModel->setEmpiricalQuadrature();
        onlineModel->offlineEq();
    }
}

template<typename te_rb_model_type>
void BiotSavartAlphaElectricCRB<te_rb_model_type>::setupCRB( crb_ptrtype crb )
{
    M_crb = crb;
    M_crbModel = crb->model();
    M_teCrbModel = M_crbModel->model();
    M_XhCond = M_teCrbModel->functionSpace();
    M_N = M_crb->dimension();
    Feel::cout << "Online model : Conductor  nDof = " << M_XhCond->nDof()
               << " Box nDof = " << this->Xh->nDof() << std::endl;
    this->setupCommunicatorsBS();
}

template<typename te_rb_model_type>
typename BiotSavartAlphaElectricCRB<te_rb_model_type>::vector_ptrtype
BiotSavartAlphaElectricCRB<te_rb_model_type>::assembleForDEIM( parameter_type const& mu, int const& tag )
{
    tic();
    if( this->isOnlineModel() && this->M_useEQ )
        return this->assembleWithEQ(mu);

    if( this->isOnlineModel() || M_useRbInDeim )
    {
        int timeSteps = 1;
        std::vector<vectorN_type> uNs(timeSteps, vectorN_type(M_N));
        std::vector<vectorN_type> uNolds(timeSteps, vectorN_type(M_N));
        std::vector<double> outputs(timeSteps, 0);
        M_crb->fixedPointPrimal(M_N, mu, uNs, uNolds, outputs);
        M_uN = uNs[0];
        this->expandV();
    }
    else
    {
        M_V = M_teCrbModel->solve(mu);
    }

    auto mesh = M_XhCond->mesh();
    vector_ptrtype Bvec = backend()->newVector( this->Xh );
    auto B = this->Xh->element( Bvec );

    auto coeff = 1/(4*M_PI);
    auto mu0 = 4*M_PI*1e-7; //SI unit : H.m-1 = m.kg.s-2.A-2

    std::vector<double> intLocD, intSumD;
    for(int i = 0; i < M_dofMgn.size(); ++i)
    {
        if( M_commsC1M[i] )
        {
            int dofSize;
            int coordSize;
            if( !this->isOnlineModel() )
            {
                dofSize = M_dofMgn.at(i).size();
                coordSize = dofSize/3;
            }
            else
            {
                dofSize = M_indexR.size();
                coordSize = dofSize;
            }

            intLocD.resize( dofSize, 0 );
            intSumD.resize( dofSize );

            if( M_XhCond->nLocalDof() > 0 )
            {
                std::vector<Eigen::Matrix<double,3,1>> coords( coordSize );
                for( int d = 0; d < coordSize; ++d )
                {
                    node_type dofCoord;
                    if( !this->isOnlineModel() )
                        dofCoord = M_dofMgn.at(i)[d*3].template get<0>();
                    else
                        dofCoord = this->Xh->dof()->dofPoint(M_indexR[d]).template get<0>();
                    Eigen::Matrix<double,3,1> coord;
                    coord << dofCoord[0], dofCoord[1], dofCoord[2];
                    coords[d] = coord;
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
                    uint16_type dofComp;
                    int j = 0;
                    if( !this->isOnlineModel() )
                    {
                        dofComp = M_dofMgn.at(i)[d].template get<2>();
                        for( auto const& mat : M_teCrbModel->materials() )
                            intLocD[d] += mgnFields[j++][d/3](dofComp,0);
                    }
                    else
                    {
                        dofComp = this->Xh->dof()->dofPoint(M_indexR[d]).template get<2>();
                        for( auto const& mat : M_teCrbModel->materials() )
                            intLocD[d] += mgnFields[j++][d](dofComp,0);
                    }
                }
            }

            mpi::reduce(M_commsC1M[i], intLocD.data(), dofSize, intSumD.data(), std::plus<double>(), 0);

            if( M_commsC1M[i].rank() == 0 )
            {
                for(int d=0; d<dofSize; d++)
                {
                    if( !this->isOnlineModel() )
                        B.set(M_dofMgn.at(i)[d].template get<1>(), intSumD[d]);
                    else
                        B.set(this->Xh->dof()->dofPoint(M_indexR[d]).template get<1>(), intSumD[d]);
                }
            }
        }
        Environment::worldComm().barrier();
    }

    B.close();

    toc("assembleBiotSavart", M_verbose > 0);

    return Bvec;
}

template<typename te_rb_model_type>
typename BiotSavartAlphaElectricCRB<te_rb_model_type>::vector_ptrtype
BiotSavartAlphaElectricCRB<te_rb_model_type>::assembleWithEQ( parameter_type const& mu )
{
    tic();
    auto mesh = M_XhCond->mesh();
    vector_ptrtype Bvec = backend()->newVector( this->Xh );
    auto B = this->Xh->element( Bvec );

    std::vector<double> intLocD, intSumD;
    for(int i = 0; i < M_dofMgn.size(); ++i)
    {
        if( M_commsC1M[i] )
        {
            int dofSize = M_indexR.size();

            intLocD.resize( dofSize, 0 );
            intSumD.resize( dofSize, 0 );

            if( M_XhCond->nLocalDof() > 0 )
            {
                int i = 0;
                if( !M_teCrbModel->materialsWithoutGeo().empty() )
                {
                    auto eq = M_eqs[i++];
                    for( int m = 0; m < eq->nbExpressions(); ++m )
                        intLocD[m] += eq->evaluate(mu, m);
                }
                for( auto const& material : M_teCrbModel->materialsWithGeo() )
                {
                    auto eq = M_eqs[i++];
                    for( int m = 0; m < eq->nbExpressions(); ++m )
                        intLocD[m] += eq->evaluate(mu, m);
                }
            }

            mpi::reduce(M_commsC1M[i], intLocD.data(), dofSize, intSumD.data(), std::plus<double>(), 0);
            if( M_commsC1M[i].rank() == 0 )
            {
                for(int d=0; d<dofSize; d++)
                    B.set(this->Xh->dof()->dofPoint(M_indexR[d]).template get<1>(), intSumD[d]);
            }
        }
        Environment::worldComm().barrier();
    }


    toc("assembleWithEQ", M_verbose > 0);

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
            }

            mpi::broadcast( M_commsC1M[i], dofSize, 0);
            M_dofMgn.at(i).resize( dofSize );
            mpi::broadcast( M_commsC1M[i], M_dofMgn.at(i).data(), dofSize, 0);
        }
    }
    toc("setup comms", M_verbose > 0);
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
void BiotSavartAlphaElectricCRB<te_rb_model_type>::online( parameter_type const& mu, int M )
{
    tic();

    tic();
    int timeSteps = 1;
    std::vector<vectorN_type> uNs(timeSteps, vectorN_type(M_N)), uNolds(timeSteps, vectorN_type(M_N));
    std::vector<double> outputs(timeSteps, 0);
    M_crb->fixedPointPrimal(M_N, mu, uNs, uNolds, outputs);
    M_uN = uNs[0];
    this->expandV();
    toc("crb run", M_verbose > 1);

    auto beta = this->deim()->beta(mu, M);
    // Feel::cout << "coefficient for m=" << M << ":\n";
    // for(int m = 0; m < beta.size(); ++m)
    //     Feel::cout << "beta[" << m << "]=" << beta(m) << std::endl;
    auto q = this->deim()->q();
    vector_ptrtype Bvec = backend()->newVector( this->Xh );
    Bvec->zero();
    for( int m = 0; m < beta.size(); ++m )
        Bvec->add( beta(m), q[m]);

    Bvec->close();
    M_B = this->Xh->element( Bvec );

    toc("online", M_verbose > 0);
}

template<typename te_rb_model_type>
void BiotSavartAlphaElectricCRB<te_rb_model_type>::expandV(int N )
{
    M_V = M_crb->expansion( M_uN, N);
}

template<typename te_rb_model_type>
void BiotSavartAlphaElectricCRB<te_rb_model_type>::computeVFE( parameter_type const& mu )
{
    tic();
    M_VFe = M_teCrbModel->solve(mu);
    toc("compute j FE", M_verbose > 1);
}

template<typename te_rb_model_type>
void BiotSavartAlphaElectricCRB<te_rb_model_type>::computeVRB( parameter_type const& mu, int N )
{
    tic();
    if( N < 1 || N > M_N )
        N = M_N;
    int timeSteps = 1;
    std::vector<vectorN_type> uNs(timeSteps, vectorN_type(N)), uNolds(timeSteps, vectorN_type(N));
    std::vector<double> outputs(timeSteps, 0);
    M_crb->fixedPointPrimal(N, mu, uNs, uNolds, outputs);
    M_uN = uNs[0];
    toc("crb run", M_verbose > 1);
}

template<typename te_rb_model_type>
void BiotSavartAlphaElectricCRB<te_rb_model_type>::computeFE( parameter_type const& mu )
{
    tic();
    tic();
    computeVFE(mu);
    toc("compute V FE", M_verbose > 1);
    tic();
    M_BFe = computeB(mu, M_VFe);
    toc("compute B FE", M_verbose > 1);
    toc("compute FE", M_verbose > 0);
}

template<typename te_rb_model_type>
void BiotSavartAlphaElectricCRB<te_rb_model_type>::computeRB( parameter_type const& mu, int N )
{
    if( N < 1 || N > M_N )
        N = M_N;
    tic();
    computeVRB(mu, N);
    expandV(N);
    tic();
    M_B = computeB(mu, M_V);
    toc("compute B RB", M_verbose > 1);
    toc("compute RB", M_verbose > 0);
}

template<typename te_rb_model_type>
typename BiotSavartAlphaElectricCRB<te_rb_model_type>::element_type
BiotSavartAlphaElectricCRB<te_rb_model_type>::computeB( parameter_type const& mu, cond_element_type const& V )
{
    tic();
    auto mesh = M_XhCond->mesh();
    auto B = this->Xh->element();
    auto coeff = 1/(4*M_PI);
    auto mu0 = 4*M_PI*1e-7; //SI unit : H.m-1 = m.kg.s-2.A-2

    if( M_verbose > 1 )
        Feel::cout << "Computing integrals for "
                   << this->Xh->nDof() << " dofs distributed on " << M_dofMgn.size()
                   << " procs on a domain with:"<< std::endl
                   << "Proc[" << Environment::rank() << "] " << M_XhCond->nLocalDof()
                   << " dofs" << std::endl;

    std::vector<double> intLocD, intSumD;
    for(int i = 0; i < M_dofMgn.size(); ++i)
    {
        if( M_commsC1M[i] )
        {
            int dofSize;
            int coordSize;
            if( !this->isOnlineModel() )
            {
                dofSize = M_dofMgn.at(i).size();
                coordSize = dofSize/3;
            }
            else
            {
                dofSize = M_indexR.size();
                coordSize = dofSize;
            }
            intLocD.resize( dofSize, 0 );
            intSumD.resize( dofSize );

            tic();
            if( M_XhCond->nLocalDof() > 0 )
            {
                std::vector<Eigen::Matrix<double,3,1>> coords( coordSize );
                for( int d = 0; d < coordSize; ++d )
                {
                    node_type dofCoord;
                    if( !this->isOnlineModel() )
                        dofCoord = M_dofMgn.at(i)[d*3].template get<0>();
                    else
                        dofCoord = this->Xh->dof()->dofPoint(M_indexR[d]).template get<0>();
                    Eigen::Matrix<double,3,1> coord;
                    coord << dofCoord[0], dofCoord[1], dofCoord[2];
                    coords[d] = coord;
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
                                               _expr=-mu0*coeff*sigma*cross(trans(gradv(V)*Jinv),
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
                                               _expr=-mu0*coeff*sigma*cross(trans(gradv(V)),
                                                                            _e1v-P())/(dist*dist*dist),
                                               _quad=_Q<1>()
                                               ).template evaluate(coords);
                }

                for( int d = 0; d < dofSize; ++d )
                {
                    uint16_type dofComp;
                    int j = 0;
                    if( !this->isOnlineModel() )
                    {
                        dofComp = M_dofMgn.at(i)[d].template get<2>();
                        for( auto const& mat : M_teCrbModel->materials() )
                            intLocD[d] += mgnFields[j++][d/3](dofComp,0);
                    }
                    else
                    {
                        dofComp = this->Xh->dof()->dofPoint(M_indexR[d]).template get<2>();
                        for( auto const& mat : M_teCrbModel->materials() )
                            intLocD[d] += mgnFields[j++][d](dofComp,0);
                    }
                }
            }
            toc("integral_computation", M_verbose > 2);

            tic();
            mpi::reduce(M_commsC1M[i], intLocD.data(), dofSize, intSumD.data(), std::plus<double>(), 0);
            toc("reduce_contributions_of_integrals", M_verbose > 2);

            if( M_commsC1M[i].rank() == 0 )
            {
                for(int d=0; d<dofSize; d++)
                {
                    if( !this->isOnlineModel() )
                        B.set(M_dofMgn.at(i)[d].template get<1>(), intSumD[d]);
                    else
                        B.set(this->Xh->dof()->dofPoint(M_indexR[d]).template get<1>(), intSumD[d]);
                }
            }
        }
        Environment::worldComm().barrier();
    }

    B.close();
    toc("compute BiotSavart", M_verbose > 1);
    return B;
}

template<typename te_rb_model_type>
typename BiotSavartAlphaElectricCRB<te_rb_model_type>::cond_element_type
BiotSavartAlphaElectricCRB<te_rb_model_type>::alpha( parameter_type const& mu )
{
    auto a = M_XhCond->element();
    auto mesh = M_XhCond->mesh();
    for( auto const& material : M_teCrbModel->materialsWithGeo() )
        a += project( _space=M_XhCond, _range=markedelements(mesh, material.first),
                      _expr=expr(M_teCrbModel->alpha(mu, material.second)) );
    return a;
}

template<typename te_rb_model_type>
double BiotSavartAlphaElectricCRB<te_rb_model_type>::homogeneity( element_type const& B)
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

template<typename te_rb_model_type>
void
BiotSavartAlphaElectricCRB<te_rb_model_type>::setIndices( std::vector<int> const& index )
{
    M_indexR = index;
}

template<typename te_rb_model_type>
void
BiotSavartAlphaElectricCRB<te_rb_model_type>::setEmpiricalQuadrature()
{
    auto mesh = this->M_XhCond->mesh();
    M_V = this->M_XhCond->element();
    auto mu = this->newParameter();
    auto Dmu = this->parameterSpace();

    auto coeff = 1/(4*M_PI);
    auto mu0 = 4*M_PI*1e-7; //SI unit : H.m-1 = m.kg.s-2.A-2

    auto matWithoutGeo = M_teCrbModel->materialsWithoutGeo();
    if( !matWithoutGeo.empty() )
    {
        auto it = matWithoutGeo.begin();
        auto r = markedelements(mesh, it->first);
        for( ++it; it != matWithoutGeo.end(); ++it )
            r = concatenate(r, markedelements(mesh, it->first) );

        M_eqs.push_back(std::make_shared<eq_type>(r, Dmu));

        auto sigma = matWithoutGeo.begin()->second.getDouble("sigma");
        for(int m = 0; m < M_indexR.size(); ++m )
        {
            auto dofCoord = this->Xh->dof()->dofPoint(M_indexR[m]).template get<0>();
            auto dofComp = this->Xh->dof()->dofPoint(M_indexR[m]).template get<2>();
            auto coord = vec(cst(dofCoord[0]),cst(dofCoord[1]),cst(dofCoord[2]));
            auto dist = inner( coord-P(), coord-P(), mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
            auto ex = -mu0*coeff*sigma*cross(trans(gradv(M_V)), coord-P())/(dist*dist*dist);
            using fn_t = std::function<bool(parameter_type const&, cond_element_type&)>;
            fn_t lambda = [this](parameter_type const& mu, cond_element_type& v) mutable -> bool {
                              Feel::cout << "lambda" << std::endl;
                              this->computeVRB(mu);
                              this->expandV();
                              Feel::cout << "lambda end" << std::endl;
                              return true;
                          };
            M_eqs.back()->addExpression(ex,mu,M_V,lambda, dofComp);
        }
        Feel::cout << "add EQ without geo" << std::endl;
    }

    auto matWithGeo = M_teCrbModel->materialsWithGeo();
    for( auto const& matp : matWithGeo )
    {
        auto marker = matp.first;
        auto material = matp.second;
        auto r = markedelements(mesh, marker);

        M_eqs.push_back(std::make_shared<eq_type>(r, Dmu));

        auto sigma = material.getDouble("sigma");
        auto paramNames = material.getVecString("params");
        std::vector<Expr<Feel::vf::Cst<boost::reference_wrapper<double> > >> paramRefs;
        for(int i = 0; i < paramNames.size(); ++i )
            paramRefs.push_back(cst_ref(mu.parameterNamed(paramNames[i])));

        auto alphaStr = M_teCrbModel->alphaRef(material);
        auto alphaExpr = expr(alphaStr, paramNames, paramRefs);
        auto alphaPrimeStr = M_teCrbModel->alphaPrimeRef(material);
        auto alphaPrimeExpr = expr(alphaPrimeStr, paramNames, paramRefs);

        auto Jinv = mat<3,3>( cos(alphaExpr), -sin(alphaExpr), -alphaPrimeExpr*Py(),
                              sin(alphaExpr), cos(alphaExpr), alphaPrimeExpr*Px(),
                              cst(0.), cst(0.), cst(1.) );

        auto psi = vec( cos(alphaExpr)*Px() + sin(alphaExpr)*Py(),
                        -sin(alphaExpr)*Px() + cos(alphaExpr)*Py(),
                        Pz() );

        for(int m = 0; m < M_indexR.size(); ++m )
        {
            auto dofCoord = this->Xh->dof()->dofPoint(M_indexR[m]).template get<0>();
            auto dofComp = this->Xh->dof()->dofPoint(M_indexR[m]).template get<2>();
            auto coord = vec(cst(dofCoord[0]),cst(dofCoord[1]),cst(dofCoord[2]));
            auto dist = inner( coord-psi,coord-psi,
                               mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
            auto ex = -mu0*coeff*sigma*cross(trans(gradv(M_V)*Jinv), coord-psi)/(dist*dist*dist);

            using fn_t = std::function<bool(parameter_type const&, cond_element_type&)>;
            fn_t lambda = [this](parameter_type const& mu, cond_element_type& v) mutable -> bool {
                              this->computeVRB(mu);
                              this->expandV();
                              v = this->potential();
                              return true;
                          };
            M_eqs.back()->addExpression(ex,mu,M_V,lambda,dofComp);
        }
        Feel::cout << "add EQ with geo for material " << marker << std::endl;
    }
}

template<typename te_rb_model_type>
void
BiotSavartAlphaElectricCRB<te_rb_model_type>::offlineEq()
{
    int err = 0;
    for( auto const& eq : M_eqs )
        err = eq->offline();
    if( err )
        Feel::cout << "error in offline EQ" << std::endl;
}

using alphabiotsavartcrbelectric = BiotSavartAlphaElectricCRB<AlphaElectric>;

template class FEELPP_EXPORT BiotSavartAlphaElectricCRB<AlphaElectric>;

FEELPP_CRB_PLUGIN( alphabiotsavartcrbelectric, alphabiotsavartcrbelectric)

}

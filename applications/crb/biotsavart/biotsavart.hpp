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

#ifndef FEELPP_BIOTSAVART_HPP
#define FEELPP_BIOTSAVART_HPP

#include <feel/feel.hpp>
#include <feel/feelcrb/crb.hpp>
#include <feel/feelcrb/eim.hpp>
#include <feel/feelcrb/crbmodel.hpp>
#include <feel/feelcrb/modelcrbbase.hpp>
#include <feel/feelfilters/exporter.hpp>

namespace Feel
{

#if 0
/*
 * base class for thermoelectric model for biotsavart
 * need also:
 * - a constructor taking a mesh
 * - a function makeThermoElectricCRBOptions() returning the options
 */
class ThermoElectricBase
{
public:
    virtual int mMaxSigma() = 0;
    virtual auto eimSigmaQ(int m) = 0;
    virtual Eigen::VectorXd eimSigmaBeta( ParameterSpace<5> mu ) = 0;
    virtual template<typename space_type> computeTruthCurrentDensity( ParameterSpace<5> mu ) = 0;
};
#endif

FEELPP_EXPORT po::options_description biotsavartOptions()
{
    po::options_description opt("BiotSavart options");
    opt.add_options()
        ( "biotsavart.conductor", po::value<std::vector<std::string> >()->multitoken(),
          "marker for the conductor" )
        ( "biotsavart.mgn", po::value<std::vector<std::string> >()->multitoken(),
          "marker for the magnetic part" )
        ( "biotsavart.repart", po::value<bool>()->default_value( false ),
          "repartition the mesh" )
        ( "biotsavart.run-mode", po::value<int>()->default_value( 2 ),
          "execution mode: pfem=0, scm=1, crb=2, scm_online=3, crb_online=4" )
        ( "biotsavart.param", po::value<std::vector<double> >()->multitoken(),
          "parameter to evaluate" )
        ( "biotsavart.compute-fe", po::value<bool>()->default_value(false),
          "compute the finite element version of Biot Savart" )
        ( "biotsavart.rebuild-database", po::value<bool>()->default_value(false),
          "rebuild the integrals or not" )
        ( "biotsavart.path-to-database", po::value<std::string>()->default_value("BiotSavart"),
          "path to the database" )
        ( "biotsavart.crb-dimension", po::value<int>()->default_value(-1),
          "number of reduced basis to use" )
        ( "biotsavart.eim-dimension", po::value<int>()->default_value(-1),
          "number of eim basis to use" )
        ;
    return opt;
}

template<typename te_rb_model_type>
FEELPP_EXPORT class BiotSavartCRB
    : public ModelCrbBase<typename te_rb_model_type::parameter_space_type,
                          typename te_rb_model_type::function_space_type,
                          NonLinear, // BiotSavart ??
                          typename te_rb_model_type::eim_definition_type >
{
    using super_type = ModelCrbBase<typename te_rb_model_type::parameter_space_type,
                                    typename te_rb_model_type::function_space_type,
                                    NonLinear,
                                    typename te_rb_model_type::eim_definition_type >;

    using te_rb_model_ptrtype = boost::shared_ptr<te_rb_model_type>;
    using crb_model_type = CRBModel<te_rb_model_type>;
    using crb_model_ptrtype = boost::shared_ptr<crb_model_type>;
    using crb_type = CRB<crb_model_type>;
    using crb_ptrtype = boost::shared_ptr<crb_type>;

    using value_type = typename te_rb_model_type::value_type;
    using param_space_type = typename te_rb_model_type::parameterspace_type;
    using parameter_type = typename param_space_type::element_type;
    using vectorN_type = Eigen::VectorXd;
    using eigen_vector_type = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using beta_vector_type = typename te_rb_model_type::beta_vector_type;

    using mesh_type = typename te_rb_model_type::mesh_type;
    using mesh_ptrtype = boost::shared_ptr<mesh_type>;
    using exporter_type = Exporter<mesh_type>;
    using exporter_ptrtype = boost::shared_ptr<exporter_type>;

    using space_type = typename te_rb_model_type::space_type;
    using space_ptrtype = boost::shared_ptr<space_type>;
    using element_type = typename space_type::element_type;
    using element_ptrtype = boost::shared_ptr<element_type>;
    using V_space_type = typename space_type::template sub_functionspace<0>::type;
    using V_space_ptrtype = boost::shared_ptr<V_space_type>;
    using V_element_type = typename V_space_type::element_type;
    using V_element_ptrtype = boost::shared_ptr<V_element_type>;
    using T_space_type = typename space_type::template sub_functionspace<1>::type;
    using T_space_ptrtype = boost::shared_ptr<T_space_type>;
    using T_element_type = typename T_space_type::element_type;
    using T_element_ptrtype = boost::shared_ptr<T_element_type>;

    using vec_fct_type = Lagrange<1, Vectorial>;
    using vec_basis_type = bases<vec_fct_type>;
    using vec_space_type = FunctionSpace<mesh_type, vec_basis_type>;
    using vec_space_ptrtype = boost::shared_ptr<vec_space_type>;
    using vec_element_type = typename vec_space_type::element_type;
    using vec_element_ptrtype = boost::shared_ptr<vec_element_type>;

    using disc_vec_fct_type = Lagrange<0, Vectorial, Discontinuous>;
    using disc_vec_basis_type = bases<disc_vec_fct_type>;
    using disc_vec_space_type = FunctionSpace<mesh_type, disc_vec_basis_type>;
    using disc_vec_space_ptrtype = boost::shared_ptr<disc_vec_space_type>;
    using disc_vec_element_type = typename disc_vec_space_type::element_type;
    using disc_vec_element_ptrtype = boost::shared_ptr<disc_vec_element_type>;

    using dof_point_type = boost::tuple<node_type, size_type, uint16_type >;
    using dof_points_type = typename std::vector<dof_point_type>;


public:
    BiotSavartCRB();

    void initModel();
    // setupSpecificityModel
    value_type
    output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);

    parameter_type paramFromOption();
    void runBS();
    void offline();
    void setupCommunicatorsBS();
    void setDimensions();
    void loadIntegrals();
    void saveIntegrals();
    void computeIntegrals(int M = 0, int N = 0);
    void online( parameter_type & mu );
    void computeB( vectorN_type & uN, eigen_vector_type & betaMu );
    void computeFE( parameter_type & mu );
    void exportResults();

    parameter_type newParameter() { return M_teCrbModel->newParameter(); }
    int nbParameters() const { return M_crbModel->parameterSpace()->dimension();  }
    auto parameterSpace() const { return M_crbModel->parameterSpace(); }
    element_type potentialTemperature() const { return M_VT; }
    vec_element_type magneticFlux() const { return M_B; }
    mesh_ptrtype meshCond() const { return M_meshCond; }
    mesh_ptrtype meshMgn() const { return M_meshMgn; }

protected:
    CRBModelMode M_mode;
    te_rb_model_ptrtype M_teCrbModel;
    crb_model_ptrtype M_crbModel;
    crb_ptrtype M_crb;

    std::vector<Eigen::MatrixXd> M_intMND;

    mesh_ptrtype M_meshCond;
    mesh_ptrtype M_meshMgn;
    space_ptrtype M_Xh;
    vec_space_ptrtype M_XhMgn;
    element_type M_VT;
    vec_element_type M_B;
    disc_vec_element_type M_j;
    vec_element_type M_BFe;
    vectorN_type M_uN;
    eigen_vector_type M_betaMu;
    parameter_type M_mu;
    int M_N;
    int M_M;

    std::vector< mpi::communicator > M_commsC1M;
    std::map<int, dof_points_type> M_dofMgn;

}; // class BiotSavartCRB

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
    M_teCrbModel = boost::make_shared<te_rb_model_type>(M_meshCond);
    M_crbModel = boost::make_shared<crb_model_type>(M_teCrbModel);
    M_crb = boost::make_shared<crb_type>("biotsavart_crb", M_crbModel);
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
    this->offline();
    M_mu = this->paramFromOption();
    this->online(M_mu);

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
                // if already computed, skip
                if ( m < M && n < N )
                    continue;

                // retrieve xi_n
                element_type xiVT_n = M_crbModel->rBFunctionSpace()->primalBasisElement(n);
                auto xi_n = xiVT_n.template element<0>();

                // we are in the communicator and we have some dof in Omega_cond
                if( M_commsC1M[i] && M_Xh->nLocalDof() > 0 )
                {
                    std::vector<Eigen::Matrix<double,3,1>> coords( dofSize );
                    // fill vector of coordinates
                    for(int d = 0; d < dofSize; ++d)
                    {
                        auto dof_coord = M_dofMgn.at(i)[d].template get<0>();
                        Eigen::Matrix<double,3,1> coord;
                        coord << dof_coord[0],dof_coord[1],dof_coord[2];
                        coords[d] = coord;
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
                        intLocD[d] = mgnField[d](dofComp,0);
                    }
                }
                if( M_commsC1M[i] )
                    boost::mpi::reduce(M_commsC1M[i], intLocD.data(),
                                       dofSize, intSumD.data(), std::plus<double>(), 0);
                if( M_commsC1M[i].rank() == 0 )
                    M_intMND[m].col(n) = Eigen::Map<Eigen::VectorXd>(intSumD.data(), dofSize);
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

    tic();
    double online_tol = doption(_name="crb.online-tolerance");
    vectorN_type time_crb;
    auto o = M_crb->run( mu, time_crb, online_tol, M_N);
    toc("crb run");

    auto solutions = o.template get<2>();
    M_uN = solutions.template get<0>()[0];
    auto WN = M_crb->wn();
    M_VT = M_crb->expansion( M_uN, M_uN.size(), WN );
    M_betaMu = M_teCrbModel->eimSigmaBeta(mu);

    this->computeB(M_uN, M_betaMu);

    toc("online");
}

template<typename te_rb_model_type>
void BiotSavartCRB<te_rb_model_type>::computeB( vectorN_type & uN, eigen_vector_type & betaMu )
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
                auto betaM = betaMu(m);
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
    M_j = M_teCrbModel->template computeTruthCurrentDensity<disc_vec_space_type>(mu);
    toc("compute j FE");

    tic();
    M_BFe = M_XhMgn->element();
    auto coeff = 1/(4*M_PI);
    auto mu0 = 4*M_PI*1e-7; //SI unit : H.m-1 = m.kg.s-2.A-2

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
                std::vector<Eigen::Matrix<double,3,1>> coords( dofSize );
                for(int d=0; d<dofSize; d++)
                {
                    auto dof_coord = M_dofMgn.at(i)[d].template get<0>();
                    Eigen::Matrix<double,3,1> coord;
                    coord << dof_coord[0],dof_coord[1],dof_coord[2];
                    coords[d] = coord;
                }

                auto dist = inner( _e1v-P(),_e1v-P(),
                                   mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
                auto mgnField = integrate(_range = elements( M_Xh->mesh() ),
                                          _expr=mu0*coeff*cross(idv(M_j), _e1v-P())/(dist*dist*dist),
                                          _quad=_Q<1>()
                    ).template evaluate(coords);

                for(int d=0; d<dofSize; d++)
                {
                    auto dofComp = M_dofMgn.at(i)[d].template get<2>();
                    intLocD[d] = mgnField[d](dofComp,0);
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

    auto err = normL2( elements(M_meshMgn), idv(M_B)-idv(M_BFe) );
    Feel::cout << "dimension CRB: " << M_N
               << " dimension EIM: " << M_M
               << " ||B_{rb} - B_{fe}||_L2 = " << err << std::endl;
}

template<typename te_rb_model_type>
void BiotSavartCRB<te_rb_model_type>::exportResults()
{
    auto WN = M_crb->wn();
    element_type VT = M_crb->expansion( M_uN, M_uN.size(), WN );
    auto V = VT.template element<0>();
    auto T = VT.template element<1>();

    auto Jh = disc_vec_space_type::New(M_Xh->mesh() );
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

    exporter_ptrtype eC( exporter_type::New( "conductor") );
    eC->setMesh(M_Xh->mesh());
    eC->add( "V", V);
    eC->add( "T", T);
    eC->add( "j", j);
    if( boption("biotsavart.compute-fe") )
        eC->add( "j_fe", M_j);
    eC->save();
    exporter_ptrtype eM( exporter_type::New( "mgn") );
    eM->setMesh(M_XhMgn->mesh());
    eM->add( "B", M_B);
    if( boption("biotsavart.compute-fe") )
        eM->add( "B_fe", M_BFe );
    eM->save();
}

} // namespace Feel

#endif

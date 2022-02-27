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

#ifndef BIOTSAVART_HPP
#define BIOTSAVART_HPP 1

#include <feel/feelmor/ser.hpp>
#include <feel/feelmor/crbsaddlepoint.hpp>
#include <feel/feelmor/crbmodelsaddlepoint.hpp>
#include <feel/feelmor/reducedbasisspace.hpp>
#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feelmor/crbplugin.hpp>

#include "thermoelectric-nonlinear.hpp"

namespace Feel
{
template <typename ThermoElectricModelT>
class BiotSavart : public CRBModelBase,
                   public boost::enable_shared_from_this<BiotSavart<ThermoElectricModelT> >
{
public:
    using self_type = BiotSavart;
    using self_ptrtype = boost::shared_ptr<self_type>;

    using parameterspace_type = ParameterSpaceX;
    using parameterspace_ptrtype = boost::shared_ptr<parameterspace_type>;
    using parameter_type = parameterspace_type::element_type;

    using thermoelectric_model_type = ThermoElectricModelT;
    using thermoelectric_model_ptrtype = boost::shared_ptr<thermoelectric_model_type>;
    using crbmodel_type = CRBModelSaddlePoint<thermoelectric_model_type>;
    using crbmodel_ptrtype = boost::shared_ptr<crbmodel_type>;
    using crb_type = CRBSaddlePoint<crbmodel_type>;
    using crb_ptrtype = boost::shared_ptr<crb_type>;
    using ser_type = SER<crb_type>;
    using ser_ptrtype = boost::shared_ptr<ser_type>;

    using prop_type = ModelProperties;
    using prop_ptrtype = boost::shared_ptr<prop_type>;
    using mat_type = ModelMaterial;
    using map_mat_type = std::map<std::string, mat_type>;

    /* needed for CRBPlugin */
    using mesh_type = typename thermoelectric_model_type::mesh_type;
    using mesh_ptrtype = boost::shared_ptr<mesh_type>;
    using vec_fct_type = Lagrange<1, Vectorial>;
    using vec_basis_type = bases<vec_fct_type>;
    using functionspace_type = FunctionSpace<mesh_type, vec_basis_type>;
    using space_type = functionspace_type;
    using functionspace_ptrtype = boost::shared_ptr<functionspace_type>;
    using element_type = typename functionspace_type::element_type;
    using element_ptrtype = boost::shared_ptr<element_type>;

    using rbfunctionspace_type = ReducedBasisSpace<self_type>;
    using rbfunctionspace_ptrtype = boost::shared_ptr<rbfunctionspace_type>;
    /* end needed  for CRBPlugin */

    using vt_space_type = typename thermoelectric_model_type::space_type;
    using vt_space_ptrtype = boost::shared_ptr<vt_space_type>;
    using vt_element_type = typename vt_space_type::element_type;
    using vt_element_ptrtype = boost::shared_ptr<vt_element_type>;
    using vt_rbfunctionspace_type = typename thermoelectric_model_type::rbfunctionspace_type;
    using vt_rbfunctionspace_ptrtype = boost::shared_ptr<vt_rbfunctionspace_type>;

    using dof_point_type = boost::tuple<node_type, size_type, uint16_type >;
    using dof_points_type = typename std::vector<dof_point_type>;

private:
    thermoelectric_model_ptrtype M_teModel;
    crbmodel_ptrtype M_crbModel;
    crb_ptrtype M_crb;
    ser_ptrtype M_ser;

    uuids::uuid M_uuid;

    parameterspace_ptrtype M_Dmu;
    int M_N;

    std::string M_propertyPath;
    prop_ptrtype M_modelProps;
    map_mat_type M_elecMaterials;
    map_mat_type M_magnMaterials;

    mesh_ptrtype M_mesh;
    functionspace_ptrtype M_Xh;
    rbfunctionspace_ptrtype M_XN;
    vt_space_ptrtype M_XhCond;

    std::vector< mpi::communicator > M_commsC1M;
    std::map<int, dof_points_type> M_dofMgn;

public:
    static po::options_description makeOptions( std::string const& prefix="biotsavart" );
protected:
    BiotSavart( crb::stage stage = crb::stage::online, std::string const& prefix="biotsavart" );
public:
    static self_ptrtype New( crb::stage stage = crb::stage::online, std::string const& prefix="biotsavart" );
    virtual ~BiotSavart() {}
    void init( std::string const& prefix="biotsavart" );

    mesh_ptrtype const& mesh() const { return M_mesh; }
    functionspace_ptrtype const& functionSpace() const { return M_Xh; } /* needed for CRBPlugin */
    rbfunctionspace_ptrtype const& rBFunctionSpace() const { return M_XN; } /* needed for CRBPlugin */
    typename functionspace_type::mesh_support_vector_type functionspaceMeshSupport( mesh_ptrtype const& mesh ) const; /* needed for CRBPlugin */
    parameterspace_ptrtype parameterSpace() const { return M_Dmu; } /* needed for CRBPlugin */
    parameter_type paramFromProperties() const { return M_teModel->paramFromProperties(); }

    /* needed for CRBModelBase */
    bool isSteady() const override { return true; }
    int numberOfTimeStep() const override { return 1; }
    double timeStep() const override { return 1e30; }
    double timeInitial() const override { return 0.; }
    double timeFinal() const override { return 1e30; }
    /* end needed for CRBModelBase */

    bool useMonolithicRbSpace() const { return true; }
    uuids::uuid uuid() const { return M_teModel->uuid(); }

    int Qb();
    int mMaxB();

    void crbOffline();
    int computeIntegrals( int n, int m );
    void setupCommunicators();
    vectorN_type computeCoeffVT( parameter_type const& mu );
    vectorN_type computeCoeff( parameter_type const& mu );
    vectorN_type computeCoeff( parameter_type const& mu, vectorN_type const& vtN );
    vt_element_type expansionVT( vectorN_type const& vtN, int N );

    int N() const { return M_N; }
    void setN( int N ) { M_N = N; }
};

template <typename ThermoElectricModelT>
po::options_description
BiotSavart<ThermoElectricModelT>::makeOptions( std::string const& prefix )
{
    po::options_description bsOptions("BiotSavart options" );
    bsOptions.add_options()
        ( "biotsavart.rebuild-database", po::value<bool>()->default_value(false), "rebuild the database" )
    //     ( prefixvm( prefix, "filename").c_str(), po::value<std::string>()->default_value("biotsavart.json"), "model file for biotsavart" )
        ;
    po::options_description teOptions = thermoelectric_model_type::makeOptions(prefix);
    bsOptions.add( teOptions );
    return bsOptions;
}

template <typename ThermoElectricModelT>
typename BiotSavart<ThermoElectricModelT>::self_ptrtype
BiotSavart<ThermoElectricModelT>::New( crb::stage stage, std::string const& prefix )
{
    auto bs = boost::shared_ptr<self_type>(new self_type(stage, prefix));
    if( stage == crb::stage::offline )
        bs->init( prefix );
    return bs;
}

template <typename ThermoElectricModelT>
BiotSavart<ThermoElectricModelT>::BiotSavart( crb::stage stage, std::string const& prefix )
    : M_XN( new rbfunctionspace_type(Environment::worldComm()) )
{
    tic();
    M_teModel = boost::make_shared<thermoelectric_model_type>(prefix);
    M_crbModel = boost::make_shared<crbmodel_type>( M_teModel, stage);
    M_crb = crb_type::New( M_teModel->modelName(), M_crbModel, stage);
    M_ser = boost::make_shared<ser_type>( M_crb, M_crbModel );
    M_Dmu = M_crb->Dmu();
    M_N = 0;
    toc("constructor + eim");
}

template <typename ThermoElectricModelT>
void
BiotSavart<ThermoElectricModelT>::init(std::string const& prefix)
{
    M_propertyPath = Environment::expand( soption(prefixvm( prefix,"filename").c_str()) );
    M_modelProps = boost::make_shared<prop_type>(M_propertyPath);
    M_elecMaterials = M_modelProps->materials().materialWithPhysic("electric");
    M_magnMaterials = M_modelProps->materials().materialWithPhysic("magneto");
    M_mesh = M_teModel->mesh();
    M_XhCond = M_teModel->functionSpace();
    std::vector<std::string> magnRange;
    for( auto const& mat : M_magnMaterials )
        magnRange.push_back( mat.first );
    M_Xh = functionspace_type::New( _mesh=M_mesh, _range=markedelements(M_mesh, magnRange) );
    M_XN->setModel( this->shared_from_this() );
    this->setupCommunicators();
}

template <typename ThermoElectricModelT>
typename BiotSavart<ThermoElectricModelT>::functionspace_type::mesh_support_vector_type
BiotSavart<ThermoElectricModelT>::functionspaceMeshSupport( mesh_ptrtype const& mesh ) const
{
    std::vector<std::string> magnRange;
    for( auto const& mat : M_magnMaterials )
        magnRange.push_back( mat.first );
    auto supp = std::make_shared<MeshSupport<mesh_type>>(mesh, markedelements(mesh, magnRange ) );
    return fusion::make_vector(supp);
}

template <typename ThermoElectricModelT>
void
BiotSavart<ThermoElectricModelT>::crbOffline()
{
    M_ser->run();
}

template <typename ThermoElectricModelT>
void
BiotSavart<ThermoElectricModelT>::setupCommunicators()
{
    tic();
    M_dofMgn.clear();
    M_commsC1M.clear();
    // Identify repartition of dofs (Omega_C, Omega_M)
    std::vector<int> isIn( 2, 0 );
    if( M_XhCond->nLocalDof()  > 0 ) // current proc has dofs in Omega_C
        isIn[0] = 1;
    if( M_Xh->nLocalDof() > 0 ) // current proc has dofs in Omega_M
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
        M_dofMgn.insert(std::make_pair(i, dofM));
    }

    // dofM = map<int, dof_points_type>
    int dofSize;
    for(int i=0; i<M_commsC1M.size(); i++)
    {
        Feel::cout << "Proc " << Environment::worldComm().globalRank()
                   << " makes " << i << "th group loop" << std::endl;
        if( M_commsC1M[i] ) // Current communicator
        {
            if(M_commsC1M[i].rank() == 0 ) // Proc of rank 0 is the Omega_M part
            {
                for ( size_type dof_id = 0; dof_id < M_Xh->nLocalDofWithGhost() ; ++dof_id )
                {
                    auto dofpoint = M_Xh->dof()->dofPoint(dof_id);
                    M_dofMgn.at(i).push_back( dofpoint );
                }
                dofSize = M_dofMgn.at(i).size();

                Feel::cout << "Proc : " << Environment::worldComm().globalRank()
                           << " has " << dofSize << " dofs in Omega_M" << std::endl;
            }

            mpi::broadcast( M_commsC1M[i], dofSize, 0);
            M_dofMgn.at(i).resize( dofSize );
            mpi::broadcast( M_commsC1M[i], M_dofMgn.at(i).data(), dofSize, 0);
        }
    }
    toc("setup comms");
}

template <typename ThermoElectricModelT>
int
BiotSavart<ThermoElectricModelT>::Qb()
{
    return M_crb->dimension();
}

template <typename ThermoElectricModelT>
int
BiotSavart<ThermoElectricModelT>::mMaxB()
{
    return std::accumulate( M_elecMaterials.begin(), M_elecMaterials.end(), 0,
                            [&](double a, auto const& matp) {
                                return a+M_teModel->mMaxSigma(matp.first);
                            } );
}

template <typename ThermoElectricModelT>
int
BiotSavart<ThermoElectricModelT>::computeIntegrals( int n, int m)
{
    tic();
    double coeff = /*M_modelProps->unit() == "mm" ? 1e-4 :*/ 1e-7;

    auto matIt = M_elecMaterials.begin();
    while( m >= M_teModel->mMaxSigma( matIt->first ) )
    {
        m -= M_teModel->mMaxSigma( matIt->first );
        matIt++;
    }
    auto mat = matIt->first;
    Feel::cout << "[" << n << "," << m << "] mat: " << mat << std::endl;

    std::vector<double> intLocD, intSumD;

    auto q_m = M_teModel->eimSigmaQ(mat, m);
    auto xi_n = M_crbModel->rBFunctionSpace()->template rbFunctionSpace<0>()->primalBasisElement(n);
    auto bqm = M_Xh->element();

    for( int i = 0; i < M_dofMgn.size(); ++i)
    {
        int dofSize = M_dofMgn.at(i).size();
        intLocD.resize(dofSize);
        intSumD.resize(dofSize);

        // we are in the communicator and we have some dof in Omega_cond
        if( M_commsC1M[i] && M_XhCond->nLocalDof() > 0 )
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
                = integrate(_range=markedelements( M_mesh, mat),
                            _expr=-coeff*idv(q_m)*cross(trans(gradv(xi_n)),
                                                        _e1v-P() )/(dist*dist*dist),
                            _quad=_Q<1>() ).template evaluate(coords);

            // fill vector with vector field components
            for( int d = 0; d < dofSize; ++d)
            {
                auto dofComp = M_dofMgn.at(i)[d].template get<2>();
                intLocD[d] = mgnField[d/3](dofComp,0);
            }
            // sum all contribution
            if( M_commsC1M[i] )
                boost::mpi::reduce(M_commsC1M[i], intLocD.data(),
                                   dofSize, intSumD.data(), std::plus<double>(), 0);

            if( M_commsC1M[i].rank() == 0 )
            {
                for( int d = 0; d < dofSize; ++d)
                    bqm.set(M_dofMgn.at(i)[d].template get<1>(), intSumD[d]);
            }
            bqm.close();
        } // comm
    } // i
    M_XN->addPrimalBasisElement( bqm );
    M_N++;
    return M_N;
}

template <typename ThermoElectricModelT>
vectorN_type
BiotSavart<ThermoElectricModelT>::computeCoeffVT( parameter_type const& mu )
{
    std::vector<vectorN_type> vtN(1), vtNdu(1);
    std::vector<double> outputs(1);
    M_crb->fixedPointPrimal( M_crb->dimension(), mu, vtN, vtNdu, outputs );
    return vtN[0];
}

template <typename ThermoElectricModelT>
vectorN_type
BiotSavart<ThermoElectricModelT>::computeCoeff( parameter_type const& mu )
{
    auto vtN = this->computeCoeffVT( mu );
    return this->computeBetaQm( mu, vtN );
}

template <typename ThermoElectricModelT>
vectorN_type
BiotSavart<ThermoElectricModelT>::computeCoeff( parameter_type const& mu, vectorN_type const& vtN )
{
    vectorN_type beta(M_N);
    int index = 0;
    for( auto const& mat : M_elecMaterials )
    {
        auto betaSigma = M_teModel->eimSigmaBeta(mat.first, mu, vtN);
        for( int m = 0; m < M_teModel->mMaxSigma( mat.first ); ++m)
        {
            for( int n = 0; n < M_crb->dimension(); ++n )
            {
                beta(index++) = betaSigma(m)*vtN(n);
            }
        }
    }
    return beta;
}

template <typename ThermoElectricModelT>
typename BiotSavart<ThermoElectricModelT>::vt_element_type
BiotSavart<ThermoElectricModelT>::expansionVT( vectorN_type const& vtN, int N )
{
    return M_crb->expansion(vtN, N);
}

FEELPP_BIOTSAVART_PLUGIN( ThermoElectric, BiotSavart, biotsavart )

} // namespace
#endif

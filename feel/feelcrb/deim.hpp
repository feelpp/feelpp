/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-05-02

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
   \file deim.hpp
   \author JB Wahl
   \date 2017-01-27
 */
#ifndef _FEELPP_DEIM_HPP
#define _FEELPP_DEIM_HPP 1

#include <limits>
#include <numeric>
#include <string>
#include <iostream>
#include <fstream>

#include <boost/ref.hpp>
#include <boost/next_prior.hpp>
#include <boost/type_traits.hpp>
#include <boost/tuple/tuple.hpp>

#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelalg/topetsc.hpp>
#include <Eigen/Core>
#include <feel/feelcrb/crbdb.hpp>
#include <feel/feelcrb/crb.hpp>
#include <feel/feelcrb/crbmodel.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>
#include <feel/feelfilters/loadmesh.hpp>

namespace Feel
{

/**
 * \brief Base class for DEIM algorithm
 *
 * Considering a parametrized Tensor \f$ T(\mu)\f$, which can be
 * either a vector or a matrix, DEIM builds an approximation of affine
 * decomposition for \f$ T(\mu)\f$ under the form : \f[
 * T(\mu)=\sum_{m=M}\beta^m(\mu)T^m \f] where the tensors \f$ T^m\f$ do not
 * depend of \f$ \mu\f$ anymore. If such an affine decomposition
 * exists for the tensor \f$ T(\mu)\f$, the approximation returned by
 * DEIM will be exact, provided that \f$ M\f$ is large enough.
 */
template <typename ModelType, typename TensorType,
          template <class T> class CRBType, template <class T> class CRBModelType>
class DEIMBase : public CRBDB
{
public :
    typedef  CRBDB super;

    typedef ModelType model_type;
    typedef boost::shared_ptr<model_type> model_ptrtype;
    typedef typename model_type::value_type value_type;

    typedef CRBModelType<model_type> crbmodel_type;
    typedef boost::shared_ptr<crbmodel_type> crbmodel_ptrtype;
    typedef CRBType<crbmodel_type> crb_type;
    typedef boost::shared_ptr<crb_type> crb_ptrtype;


    typedef typename ModelType::parameterspace_type parameterspace_type;
    typedef typename ModelType::parameterspace_ptrtype parameterspace_ptrtype;
    typedef typename ModelType::parameter_type parameter_type;
    typedef typename parameterspace_type::sampling_type sampling_type;
    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;

    typedef TensorType tensor_type;
    typedef boost::shared_ptr<tensor_type> tensor_ptrtype;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    typedef std::function<tensor_ptrtype ( parameter_type const& mu)> assemble_function_type;

    typedef Eigen::MatrixXd matrixN_type;
    typedef Eigen::VectorXd vectorN_type;

    typedef boost::tuple<parameter_type,double> bestfit_type;

    static const bool is_matrix = std::is_same<tensor_type,sparse_matrix_type>::value;
    typedef typename mpl::if_< mpl::bool_<is_matrix>,
                               std::pair<int,int>, int >::type indice_type;
    typedef typename mpl::if_< mpl::bool_<is_matrix>,
                               std::pair<std::set<int>,std::set<int>>, std::set<int> >::type ptsset_type;

    typedef boost::tuple<double,int> vectormax_type;
    typedef boost::tuple<double,std::pair<int,int>> matrixmax_type;

    typedef typename model_type::space_type space_type;
    typedef typename model_type::space_ptrtype space_ptrtype;
    typedef typename model_type::element_type element_type;
    typedef typename model_type::mesh_type mesh_type;
    typedef typename model_type::mesh_ptrtype mesh_ptrtype;

    //! Default Constructor
    DEIMBase()
    {}

    /**
     * \brief Constructor of DEIM class \p Dmu a pointer to the
     * considered parameter space \p sampling the sampling on which
     * the greedy algorithm will iterate. Can be null and so a default
     * sampling will be build, considering the options : (int)
     * deim.default-sampling-size and (string)
     * deim.default-sampling-mode
     */
    DEIMBase( model_ptrtype model, sampling_ptrtype sampling, std::string prefix="",
              WorldComm const& worldComm = Environment::worldComm() ) :
        super ( ( boost::format( "%1%DEIM-%2%" ) %(is_matrix ? "M":"") %prefix  ).str(),
                "deim",
                model->uuid(),
                worldComm ),
        M_model( model ),
        mesh( model->functionSpace()->mesh() ),
        M_parameter_space( model->parameterSpace() ),
        M_trainset( sampling ),
        M_M(0),
        M_tol(1e-8),
        M_prefix( prefix ),
        M_rebuild( boption( prefixvm( M_prefix, "deim.rebuild-db") ) ),
        M_nl_assembly(false)
    {
        using Feel::cout;

        if ( !M_rebuild )
        {
            if ( this->loadDB() )
                cout<<"DEIM : Database loaded with " << M_M << " basis functions\n";
            else
                cout <<"DEIM : No Database loaded : start greedy algorithm from beginning\n";
        }
        else
            cout << "DEIM : option deim.rebuild-db=true : start greedy algorithm from beginning\n";


        if ( !M_trainset )
            M_trainset = M_parameter_space->sampling();
        if ( M_trainset->empty() )
        {
            int sampling_size = ioption( prefixvm( M_prefix, "deim.default-sampling-size" ) );
            std::string file_name = ( boost::format("deim_trainset_%1%") % sampling_size ).str();
            std::string sampling_mode = soption( prefixvm( M_prefix, "deim.default-sampling-mode" ) );
            std::ifstream file ( file_name );
            if( ! file )
            {
                M_trainset->sample( sampling_size, sampling_mode, true, file_name );
            }
            else
            {
                M_trainset->clear();
                M_trainset->readFromFile(file_name);
            }
            cout << "DEIM sampling created with " << sampling_size << " points.\n";
        }

        M_online_model = model_ptrtype( new model_type("deim-online", true) );
    }

    //! Destructor
    virtual ~DEIMBase()
    {}

    /**
     * The assemble function has to be provided by the model. This
     * function will assemble the tensor \f$ T(\mu) \f$.
     * \p mu the considered parameter
     * \return the tensor \f$ T(\mu) \f$ as vector_ptrtype
     */
    virtual tensor_ptrtype assemble( parameter_type const& mu, bool online=false )=0;
    virtual tensor_ptrtype assemble( parameter_type const& mu, element_type const& u, bool online=false )=0;

    /**
     * Run the greedy algotithm and build the affine decomposition.
     * Greedy algorithm stops when the affine decomposition is exact
     * or when the maximum number of terms is the decomposition is
     * reached. This maximum is defined by the option
     * deim.dimension-max
     */
    void run()
    {

        using Feel::cout;
        tic();
        int mMax = ioption(  prefixvm( M_prefix, "deim.dimension-max" ) );
        double error=0;
        auto mu = M_trainset->max().template get<0>();

        if ( M_M==0 )
        {
            cout <<"===========================================\n";
            cout << "DEIM : Start algorithm with mu="<<muString(mu)<<std::endl;

            tic();
            addNewVector(mu);
            toc("Add new vector in DEIM basis");
        }
        if ( M_M<mMax )
        {
            auto best_fit = computeBestFit();
            error = best_fit.template get<1>();
            mu = best_fit.template get<0>();
            cout << "DEIM : Current error : "<< error
                 <<", tolerance : " << M_tol <<std::endl;
            cout <<"===========================================\n";
        }

        while( M_M<mMax && error>M_tol )
        {
            cout << "DEIM : Construction of basis "<<M_M+1<<"/"<<mMax<<", with mu="<<muString(mu)<<std::endl;

            tic();
            addNewVector(mu);
            toc("Add new vector in DEIM basis");

            if ( M_M<mMax )
            {
                auto best_fit = computeBestFit();
                error = best_fit.template get<1>();
                mu = best_fit.template get<0>();
                cout << "DEIM : Current error : "<< error
                           <<", tolerance : " << M_tol <<std::endl;
                cout <<"===========================================\n";
            }

        }

        cout << "DEIM : Stopping greedy algorithm. Number of basis function : "<<M_M<<std::endl;

        toc("DEIM : Total Time");
    }

    //! \return the \f$ \beta^m(\mu)\f$ for a specific parameter \p mu
    vectorN_type beta( parameter_type const& mu )
    {
        return computeCoefficient( mu );
    }

    vectorN_type beta( parameter_type const& mu, element_type const& u )
    {
        return computeCoefficient( mu, u );
    }


    //! \return the tensors \f$ T^m\f$ of the affine decomposition
    std::vector<tensor_ptrtype> q() const
    {
        return M_bases;
    }
    tensor_ptrtype q( int m ) const
    {
        return M_bases[m];
    }


    //! \return the current size of the affine decompostion
    int size() const
    {
        return M_M;
    }

    //! save the database
    void saveDB() override;
    //! load the database
    bool loadDB() override;
    void loadDB( std::string const& filename, crb::load l ) override {}


protected :
    //! add a new Tensor in the base, evaluated for parameter \p mu
    void addNewVector( parameter_type const& mu )
    {
        LOG(INFO) << "DEIM : addNewVector() start with "<<muString(mu);
        M_mus.push_back(mu);
        tensor_ptrtype Phi = residual( mu );

        auto vec_max = vectorMaxAbs( Phi );
        auto i = vec_max.template get<1>();
        double max = vec_max.template get<0>();

        M_M++;

        Phi->scale( 1./max );
        M_bases.push_back( Phi );
        M_index.push_back(i);

        M_B.conservativeResize(M_M,M_M);
        // update last row of M_B
        for ( int j = 0; j<M_M; j++ )
            M_B(M_M-1, j) = evaluate( M_bases[j], M_index[M_M-1]);

        //update last col of M_B
        for ( int i=0; i<M_M-1; i++ )
            M_B(i, M_M-1) = evaluate( M_bases[M_M-1], M_index[i] );

        updateSubMesh();
        //this->saveDB();
        LOG(INFO) << "DEIM : addNewVector() end";
    }

    //! \return the value of the component \p index of the vector \p V
    double evaluate( vector_ptrtype V, int const& index, bool seq=false )
    {
        if ( seq )
            return V->operator()( index );

        double value=0;
        int proc_number = V->map().procOnGlobalCluster(index);

        if ( Environment::worldComm().globalRank()==proc_number )
            value = V->operator()( index - V->map().firstDofGlobalCluster() );

        boost::mpi::broadcast( Environment::worldComm(), value, proc_number );
        return value;
    }

    //! \return the value of the entry \p idx of the matrix \p M
    double evaluate( sparse_matrix_ptrtype M, std::pair<int,int> const& idx, bool seq=false )
    {
        int i = idx.first;
        int j =idx.second;

        if ( seq )
            return M->operator() ( i,j );

        double value=0;
        int proc_number = M->mapRow().procOnGlobalCluster(i);

        if ( !M->closed() )
            M->close();
        if ( Environment::worldComm().globalRank()==proc_number )
            value = M->operator() ( i - M->mapRow().firstDofGlobalCluster(),
                                    j - M->mapCol().firstDofGlobalCluster());

        boost::mpi::broadcast( Environment::worldComm(), value, proc_number );
        return value;
    }

    /**
     * Evaluate the maximum of the components of \p V in absolute
     * value.
     *
     * \return a couple (max,index) where max is the maximum in
     * absolute value and index is the index of the component where
     * the max is reached
     */
    vectormax_type vectorMaxAbs( vector_ptrtype V )
    {
        auto newV = V->clone();
        *newV = *V;
        newV->abs();
        int index = 0;
        double max = newV->maxWithIndex( &index );

        double max_eval=math::abs(evaluate(V,index));
        DCHECK( max_eval==max ) << "evaluation of V(i)=" <<max_eval<<", maximum="<<max <<std::endl;
        return boost::make_tuple( max, index );
    }

    /**
     * Evaluate the maximum of the entry of \p M in absolute
     * value.
     *
     * \return a couple (max,index) where max is the maximum in
     * absolute value and index is the couple of index of the entry
     * where the max is reached
     */
    matrixmax_type vectorMaxAbs( sparse_matrix_ptrtype M )
    {
        DCHECK( M->isInitialized() ) << "MatrixPetsc<> not initialized";
        if ( !M->closed() )
            M->close();

        auto V = backend()->newVector( M->mapRowPtr() );
        PetscInt idx[M->mapRow().nLocalDof()];

        int ierr=MatGetRowMaxAbs( toPETSc(M)->mat(), toPETSc(V)->vec(), idx );
        CHKERRABORT( M->comm(),ierr );

        int i_row = 0;
        int i_col = 0;
        PetscReal val=0;
        VecMax(toPETSc(V)->vec(), &i_row, &val);
        double max = static_cast<Real>( val );

        int proc_number = V->map().procOnGlobalCluster(i_row);
        if ( Environment::worldComm().globalRank()==proc_number )
            i_col = idx[i_row - V->map().firstDofGlobalCluster()];
        boost::mpi::broadcast( Environment::worldComm(), i_col, proc_number );

        std::pair<int,int> index (i_row,i_col);
        double max_eval=math::abs(evaluate(M,index));
        DCHECK( max_eval==max ) << "evaluation of M(i,j)=" <<max_eval<<", maximum="<<max <<std::endl;

        return boost::make_tuple( max, index );
    }

    //! \return the beta coefficient for parameter \p mu
    vectorN_type computeCoefficient( parameter_type const& mu )
    {
        tensor_ptrtype T = assemble( mu, true );
        return computeCoefficient( T );
    }
    vectorN_type computeCoefficient( parameter_type const& mu, element_type const& u )
    {
        tensor_ptrtype T = assemble( mu, u, true );
        return computeCoefficient( T );
    }


    //! Compute the beta coefficients for a assembled tensor \p T
    vectorN_type computeCoefficient( tensor_ptrtype T, bool online=true )
    {
        vectorN_type rhs (M_M);
        vectorN_type coeff (M_M);
        if ( M_M>0 )
        {
            for ( int i=0; i<M_M; i++ )
                rhs(i) = evaluate( T, online ? M_indexR[i]:M_index[i], online );
            coeff = M_B.lu().solve( rhs );
        }
        return coeff;
    }

    /**
     * Choose the best parameter to build the next basis vector.
     *
     * For each parameter mu in sampling, commpute the residual \f$
     * R_M(\mu)=T(\mu)-\hat T_M(\mu) \f$ where \f$ \hat T_M(\mu) \f$
     * is the DEIM approximation of size \f$ M\f$.
     * \return a couple (mu_max, max). mu_max=\f$ argmax_{\mu\in S}R(\mu)\f$
     * and max=\f$ max_{\mu\in S}R(\mu)\f$
     */
    bestfit_type computeBestFit()
    {
        tic();
        LOG(INFO) << "DEIM : computeBestFit() start";
        double max=0;
        auto mu_max = M_trainset->max().template get<0>();
        for ( auto const& mu : *M_trainset )
        {
            tensor_ptrtype T = residual( mu );

            double norm = T->linftyNorm();
            if ( norm>max )
            {
                max = norm;
                mu_max = mu;
            }
        }

        LOG(INFO) << "DEIM : computeBestFit() end";
        toc("DEIM : compute best fit");

        return boost::make_tuple( mu_max, max );
    }

    /**
     * Evavaluate the residual \f$ R_M(\mu)=T(\mu)-\hat T_M(\mu) \f$
     * where \f$ \hat T_M(\mu) \f$ is the DEIM approximation of size
     * \f$ M\f$.
     * \return a shared pointer on the Residual.
     */
    tensor_ptrtype residual( parameter_type const& mu )
    {
        LOG(INFO) << "DEIM : residual() start with "<< muString(mu);
        tensor_ptrtype T;

        T = this->assemble( mu );

        double norm = T->linftyNorm();
        vectorN_type coeff = computeCoefficient(T, false);

        auto newT = copyTensor( T );

        for ( int i=0; i<M_M; i++ )
            add( newT, -coeff(i), M_bases[i] );
        LOG(INFO) << "DEIM : residual() end";
        newT->scale( 1./norm );
        return newT;
    }

    //! evaluate V= V+a*vec
    void add( vector_ptrtype V, double const& a, vector_ptrtype vec )
    {
        V->add( a, vec );
    }

    //! evaluate M= M+a*mat
    void add( sparse_matrix_ptrtype M, double const& a, sparse_matrix_ptrtype mat )
    {
        M->addMatrix( a, mat );
    }

    // \return a shared pointer on a copy of \p V
    vector_ptrtype copyTensor( vector_ptrtype V )
    {
        vector_ptrtype newV = V->clone();
        *newV = *V;
        return newV;
    }

    // \return a shared pointer on a copy of \p M
    sparse_matrix_ptrtype copyTensor( sparse_matrix_ptrtype M )
    {
        sparse_matrix_ptrtype newM = backend()->newMatrix( M->mapColPtr(),
                                                           M->mapRowPtr(),
                                                           M->graph() );
        *newM=*M;
        return newM;
    }

    std::string muString( parameter_type const& mu )
    {
        std::ostringstream mu_str;
        mu_str << "["<<mu[0];
        for ( int i=1; i<mu.size(); i++ )
            mu_str <<","<< mu[i];
        mu_str <<"]";
        return mu_str.str();
    }


    virtual void updateSubMesh()
    {}


    friend class boost::serialization::access;

    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const;

    template<class Archive>
    void load( Archive & ar, const unsigned int version ) ;

    BOOST_SERIALIZATION_SPLIT_MEMBER()

protected :
    model_ptrtype M_model, M_online_model;
    mesh_ptrtype mesh;
    parameterspace_ptrtype M_parameter_space;
    sampling_ptrtype M_trainset;
    int M_M;
    double M_tol;
    matrixN_type M_B;
    std::vector< tensor_ptrtype > M_bases;
    std::vector<indice_type> M_index, M_indexR, M_ldofs;
    std::vector<ptsset_type> M_ptsset;
    std::vector<parameter_type> M_mus;
    std::set<int> M_elts_ids;

    std::string M_prefix;
    bool M_rebuild, M_nl_assembly;
};

template <typename ModelType, typename TensorType,
          template <class T> class CRBType, template <class T> class CRBModelType>
template<class Archive>
void
DEIMBase<ModelType,TensorType,CRBType,CRBModelType>::load( Archive & ar, const unsigned int version )
{
    ar & boost::serialization::base_object<super>( *this );
    ar & BOOST_SERIALIZATION_NVP( M_parameter_space );
    ar & BOOST_SERIALIZATION_NVP( M_M );
    ar & BOOST_SERIALIZATION_NVP( M_B );
    ar & BOOST_SERIALIZATION_NVP( M_index );

    for( int m=0; m<M_M; m++ )
        ar & BOOST_SERIALIZATION_NVP( M_bases[m] );

    DCHECK( M_bases.size()==M_M )<<"Wrong size : M_bases.size()="<< M_bases.size() <<", M_M=" <<M_M<<std::endl;
    for ( int i=0; i<M_M; i++)
        DCHECK( M_bases[i] )<<"Null ptr at i="<<i<<std::endl;
}

template <typename ModelType, typename TensorType,
          template <class T> class CRBType, template <class T> class CRBModelType>
template<class Archive>
void
DEIMBase<ModelType,TensorType,CRBType,CRBModelType>::save( Archive & ar, const unsigned int version ) const
{
#if 0
    ar & boost::serialization::base_object<super>( *this );
    ar & BOOST_SERIALIZATION_NVP( M_parameter_space );
    ar & BOOST_SERIALIZATION_NVP( M_M );
    ar & BOOST_SERIALIZATION_NVP( M_B );
    ar & BOOST_SERIALIZATION_NVP( M_index );
    for( int m=0; m<M_M; m++ )
        ar & BOOST_SERIALIZATION_NVP( M_bases[m] );
#endif
}


template <typename ModelType, typename TensorType,
          template <class T> class CRBType, template <class T> class CRBModelType>
void
DEIMBase<ModelType,TensorType,CRBType,CRBModelType>::saveDB()
{
#if 0
    Feel::cout<< "DEIM : saving DB\n";
    fs::ofstream ofs( this->dbLocalPath() / this->dbFilename() );
    if ( ofs )
    {
        //boost::archive::text_oarchive oa( ofs );
        boost::archive::binary_oarchive oa( ofs );
        // write class instance to archive
        oa << *this;
        // archive and stream closed when destructors are called
    }
#endif
}


template <typename ModelType, typename TensorType,
          template <class T> class CRBType, template <class T> class CRBModelType>
bool
DEIMBase<ModelType,TensorType,CRBType,CRBModelType>::loadDB()
{
    return false;
    if( this->isDBLoaded() )
    {
        return true;
    }


    fs::path db = this->lookForDB();
    if ( db.empty()  )
    {
        return false;
    }

    fs::ifstream ifs( db );
    if ( ifs )
    {
        //boost::archive::text_iarchive ia( ifs );
        boost::archive::binary_iarchive ia( ifs );
        //write class instance to archive
        ia >> *this;
        this->setIsLoaded( true );
        return true;
    }

    return false;
}



template <typename ModelType,
          template <class T> class CRBType=CRB,
          template <class T> class CRBModelType=CRBModel>
class DEIM :
        public DEIMBase<ModelType, typename Backend<typename ModelType::value_type>::vector_type,
                        CRBType, CRBModelType>
{
public :
    typedef DEIMBase<ModelType, typename Backend<typename ModelType::value_type>::vector_type,CRBType,CRBModelType> super_type;
    typedef typename super_type::parameter_type parameter_type;
    typedef typename super_type::parameterspace_ptrtype parameterspace_ptrtype;
    typedef typename super_type::sampling_ptrtype sampling_ptrtype;
    typedef typename super_type::model_ptrtype model_ptrtype;
    typedef typename super_type::tensor_ptrtype vector_ptrtype;
    typedef typename super_type::element_type element_type;
    typedef typename super_type::mesh_type mesh_type;
    typedef typename super_type::space_type space_type;

    DEIM() :
        super_type()
    {}

    DEIM( model_ptrtype model, sampling_ptrtype sampling=nullptr, std::string prefix="" ) :
        super_type( model, sampling, prefix )
    {
        auto mu = this->M_parameter_space->element();
        auto T = this->assemble(mu);
        if (!T)
        {
            auto u = this->M_model->functionSpace()->element();
            auto Tnl = this->assemble(mu,u);
            CHECK( Tnl ) << "You want to use DEIM but you did not implement assmbleForDEIM functions\n";
            this->M_nl_assembly=true;
        }
    }

    ~DEIM()
    {}

    vector_ptrtype assemble( parameter_type const& mu, bool online=false )
    {
        if ( this->M_nl_assembly )
        {
            CHECK(!online) << "Call of online nl assembly with no solution u\n";
            auto u = this->M_model->solve(mu);
            return this->M_model->assembleForDEIM(mu,u);
        }
        if (online)
            return this->M_online_model->assembleForDEIM(mu);
          return this->M_model->assembleForDEIM(mu);
    }
    vector_ptrtype assemble( parameter_type const& mu, element_type const& u, bool online=false )
    {
        CHECK(this->M_nl_assembly) << "You called nl coefficient for DEIM but you implemented assembleForDEIM(mu)\n";
        if ( online )
            return this->M_online_model->assembleForDEIM(mu,u);
        return this->M_model->assembleForDEIM(mu,u);
    }

private :
    void updateSubMesh()
    {
        // Last added index
        int index = this->M_index.back();
        auto Xh = this->M_model->functionSpace();
        auto mesh = Xh->mesh();
        int proc_n = Xh->dof()->procOnGlobalCluster(index);

        // recover the elements which share this index
        std::set<int> pts_ids;
        int ldof=-1;

        if ( Environment::worldComm().globalRank()==proc_n )
        {
            for ( auto const& dof : Xh->dof()->globalDof(index-Xh->dof()->firstDofGlobalCluster()) )
            {
                this->M_elts_ids.insert( dof.second.elementId() );
                if (ldof==-1)
                {
                    ldof = dof.second.localDof();
                    auto elt = mesh->element( dof.second.elementId() );
                    for ( int p=0; p<elt.nPoints(); p++ )
                        pts_ids.insert( elt.point(p).id()+1 );
                }
            }
        }
        else
        {
            auto e = Xh->dof()->searchGlobalProcessDof( index );
            if ( e.template get<0>() )
            {
                auto process_dof = e.template get<1>();
                for ( auto const& dof : Xh->dof()->globalDof( process_dof) )
                    this->M_elts_ids.insert( dof.second.elementId() );
            }
        }

        // store the information in all proc !
        boost::mpi::broadcast( Environment::worldComm(), pts_ids, proc_n );
        boost::mpi::broadcast( Environment::worldComm(), ldof, proc_n );
        this->M_ldofs.push_back( ldof );
        this->M_ptsset.push_back( pts_ids );

        // create new submesh with the new elements and reread it in sequential
        auto submesh = createSubmesh( mesh, idelements(mesh,this->M_elts_ids.begin(), this->M_elts_ids.end()) );
        saveGMSHMesh( _mesh=submesh, _filename="deim-submesh.msh" );
        auto seqmesh = loadMesh( _mesh=new mesh_type,
                                 _filename="deim-submesh.msh", _worldcomm= Environment::worldCommSeq() );
        auto Rh = space_type::New( seqmesh,
                                   _worldscomm=std::vector<WorldComm>(space_type::nSpaces,Environment::worldCommSeq()) );

        // create map between points id and elements id
        std::map<std::set<int>,int> elts_map;
        for ( auto const& eltWrap : elements(seqmesh) )
        {
            auto const& elt = unwrap_ref( eltWrap );
            std::set<int> pts_id;
            for ( int p=0; p<elt.nPoints(); p++ )
                pts_id.insert( elt.point(p).id() );
            elts_map[pts_id] = elt.id();
        }

        // on each proc : store the new indexR corresponding to the reduced space
        this->M_indexR.clear();
        for ( int i=0; i<this->M_ptsset.size(); i++ )
        {
            auto map_it = elts_map.find( this->M_ptsset[i] );
            CHECK( map_it!=elts_map.end() ) <<"DEIM : elt id not found in map, on proc : "
                                            <<Environment::worldComm().globalRank() << std::endl;
            int elt_idR = map_it->second;
            this->M_indexR.push_back( Rh->dof()->localToGlobalId( elt_idR, this->M_ldofs[i] ) );
        }

        this->M_online_model->setFunctionSpaces( Rh );
    }


};


template <typename ModelType,
          template <class T> class CRBType=CRB,
          template <class T> class CRBModelType=CRBModel>
class MDEIM :
        public DEIMBase<ModelType,
                        typename Backend<typename ModelType::value_type>::sparse_matrix_type,
                        CRBType, CRBModelType>
{
public :
    typedef DEIMBase<ModelType, typename Backend<typename ModelType::value_type>::sparse_matrix_type,CRBType,CRBModelType> super_type;

    typedef typename super_type::parameter_type parameter_type;
    typedef typename super_type::parameterspace_ptrtype parameterspace_ptrtype;
    typedef typename super_type::sampling_ptrtype sampling_ptrtype;
    typedef typename super_type::model_ptrtype model_ptrtype;
    typedef typename super_type::tensor_ptrtype vector_ptrtype;
    typedef typename super_type::element_type element_type;
    typedef typename super_type::mesh_type mesh_type;
    typedef typename super_type::space_type space_type;


    MDEIM() :
        super_type()
    {}

    MDEIM( model_ptrtype model, sampling_ptrtype sampling=nullptr, std::string prefix="" ) :
        super_type( model, sampling, prefix )
    {
        auto mu = this->M_parameter_space->element();
        auto T = this->assemble(mu);
        if (!T)
        {
            auto u = this->M_model->functionSpace()->element();
            auto Tnl = this->assemble(mu,u);
            CHECK( Tnl ) << "You want to use MDEIM but you did not implement assmbleForMDEIM functions\n";
            this->M_nl_assembly=true;
        }
    }

    ~MDEIM()
    {}

    vector_ptrtype assemble( parameter_type const& mu, bool online=false )
    {
        if ( this->M_nl_assembly )
        {
            CHECK(!online) << "Call of online nl assembly with no solution u\n";
            auto u = this->M_model->solve(mu);
            return this->M_model->assembleForMDEIM(mu,u);
        }
        if (online)
            return this->M_online_model->assembleForMDEIM(mu);
         return this->M_model->assembleForMDEIM(mu);
    }
    vector_ptrtype assemble( parameter_type const& mu, element_type const& u, bool online=false )
    {
        CHECK(this->M_nl_assembly) << "You called nl coefficient for MDEIM but you implemented assembleForMDEIM(mu)\n";
        if ( online )
            return this->M_online_model->assembleForMDEIM(mu,u);
        return this->M_model->assembleForMDEIM(mu,u);
    }


private :
    void updateSubMesh()
    {
        auto index = this->M_index.back();
        int i1 = index.first;
        int i2 = index.second;
        auto Xh = this->M_model->functionSpace();
        auto mesh = Xh->mesh();
        bool is_same = (i1==i2);

        int proc_n1 = Xh->dof()->procOnGlobalCluster(i1);
        int proc_n2 = Xh->dof()->procOnGlobalCluster(i2);

        std::set<int> pts_ids1, pts_ids2;
        int ldof1=-1;
        int ldof2=-1;


        if ( Environment::worldComm().globalRank()==proc_n1 )
        {
            for ( auto const& dof : Xh->dof()->globalDof(i1-Xh->dof()->firstDofGlobalCluster()) )
            {
                this->M_elts_ids.insert( dof.second.elementId() );
                if ( ldof1==-1 )
                {
                    ldof1 = dof.second.localDof();
                    auto elt = mesh->element(dof.second.elementId() );
                    for ( int p=0; p<elt.nPoints(); p++ )
                        pts_ids1.insert( elt.point(p).id()+1 );
                }
            }
        }
        else
        {
            auto e = Xh->dof()->searchGlobalProcessDof( i1 );
            if ( e.template get<0>() )
            {
                auto process_dof = e.template get<1>();
                for ( auto const& dof : Xh->dof()->globalDof( process_dof) )
                    this->M_elts_ids.insert( dof.second.elementId() );
            }
        }

        boost::mpi::broadcast( Environment::worldComm(), pts_ids1, proc_n1 );
        boost::mpi::broadcast( Environment::worldComm(), ldof1, proc_n1 );

        if ( is_same )
        {
            pts_ids2 = pts_ids1;
            ldof2 = ldof1;
        }
        else
        {
            if ( Environment::worldComm().globalRank()==proc_n2 )
            {
                for ( auto const& dof : Xh->dof()->globalDof(i2-Xh->dof()->firstDofGlobalCluster()) )
                {
                    this->M_elts_ids.insert( dof.second.elementId() );
                    if ( ldof2==-1 )
                    {
                        ldof2 = dof.second.localDof();
                        auto elt = mesh->element(dof.second.elementId() );
                        for ( int p=0; p<elt.nPoints(); p++ )
                            pts_ids2.insert( elt.point(p).id()+1 );
                    }
                }
            }
            else
            {
                auto e = Xh->dof()->searchGlobalProcessDof( i2 );
                if ( e.template get<0>() )
                {
                    auto process_dof = e.template get<1>();
                    for ( auto const& dof : Xh->dof()->globalDof( process_dof) )
                        this->M_elts_ids.insert( dof.second.elementId() );
                }
            }
            boost::mpi::broadcast( Environment::worldComm(), pts_ids2, proc_n2 );
            boost::mpi::broadcast( Environment::worldComm(), ldof2, proc_n2 );
        }

        this->M_ldofs.push_back( std::make_pair(ldof1,ldof2) );
        this->M_ptsset.push_back( std::make_pair(pts_ids1,pts_ids2) );

        // create new submesh with the new elements and reread it in sequential
        auto submesh = createSubmesh( mesh, idelements(mesh,this->M_elts_ids.begin(), this->M_elts_ids.end()) );
        saveGMSHMesh( _mesh=submesh, _filename="mdeim-submesh.msh" );
        auto seqmesh = loadMesh( _mesh=new mesh_type,
                                 _filename="mdeim-submesh.msh", _worldcomm= Environment::worldCommSeq() );
        auto Rh = space_type::New( seqmesh,
                                   _worldscomm=std::vector<WorldComm>(space_type::nSpaces,Environment::worldCommSeq()) );

        // create map between points id and elements id
        std::map<std::set<int>,int> elts_map;
        for ( auto const& eltWrap : elements(seqmesh) )
        {
            auto const& elt = unwrap_ref( eltWrap );
            std::set<int> pts_id;
            for ( int p=0; p<elt.nPoints(); p++ )
                pts_id.insert( elt.point(p).id() );
            elts_map[pts_id] = elt.id();
        }


        // on each proc : store the new indexR corresponding to the reduced space
        this->M_indexR.clear();
        for ( int i=0; i<this->M_ptsset.size(); i++ )
        {
            auto ptsset1 = this->M_ptsset[i].first;
            auto ptsset2 = this->M_ptsset[i].second;

            auto map_it = elts_map.find( ptsset1 );
            CHECK( map_it!=elts_map.end() ) << "MDEIM : elt id not found in map, on proc : "
                                            << Environment::worldComm().globalRank() << std::endl;

            int elt_idR1 = map_it->second;
            int elt_idR2 = elt_idR1;
            if ( ptsset1!=ptsset2 )
            {
                map_it = elts_map.find( ptsset2 );
                CHECK( map_it!=elts_map.end() ) << "MDEIM : elt id not found in map, on proc : "
                                                << Environment::worldComm().globalRank() << std::endl;
                elt_idR2 = map_it->second;
            }

            auto indexR1 = Rh->dof()->localToGlobalId( elt_idR1, this->M_ldofs[i].first );
            auto indexR2 = Rh->dof()->localToGlobalId( elt_idR2, this->M_ldofs[i].second );

            this->M_indexR.push_back( std::make_pair(indexR1, indexR2) );
        }

        this->M_online_model->setFunctionSpaces( Rh );

    }


};


po::options_description deimOptions( std::string const& prefix ="");

} //namespace Feel

#if 0
namespace boost
{
namespace serialization
{
template< typename M,typename T >
struct version< Feel::DEIMBase<M,T> >
{
    typedef mpl::int_<1> type;
    typedef mpl::integral_c_tag tag;
    static const unsigned int value = version::type::value;
};
template<typename M, typename T> const unsigned int version<Feel::DEIMBase<M,T> >::value;
}
}
#endif

#endif

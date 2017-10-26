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

class ModelCrbBaseBase;

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
template <typename ParameterSpaceType, typename SpaceType, typename TensorType>
class DEIMBase : public CRBDB
{
public :
    typedef  CRBDB super;

    typedef typename TensorType::value_type value_type;

    typedef ParameterSpaceType parameterspace_type;
    typedef typename boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef typename parameterspace_type::element_type parameter_type;
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

    typedef SpaceType space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename space_type::mesh_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef std::map<int,tensor_ptrtype> solutionsmap_type;

    typedef CRBBase<space_type,parameterspace_type> crb_type;
    typedef boost::shared_ptr<crb_type> crb_ptrtype;

    //! Default Constructor
    DEIMBase()
    {}

    /**
     * \brief Constructor of DEIM class
     */
    DEIMBase(  space_ptrtype Xh,
               parameterspace_ptrtype Dmu,
               sampling_ptrtype sampling,
               uuids::uuid const& uid,
               std::string const& modelname,
               std::string prefix="",
               WorldComm const& worldComm = Environment::worldComm() ) :
        super ( ( boost::format( "%1%DEIM%2%" ) %(is_matrix ? "M":"") %prefix  ).str(),
                "deim",
                uid,
                worldComm ),
        M_prefix( prefix ),
        M_Dmu( Dmu ),
        M_trainset( sampling ),
        M_M(0),
        M_user_max( ioption(  prefixvm( M_prefix, "deim.dimension-max" ) ) ),
        M_n_rb( 0 ),
        M_tol( doption( prefixvm( M_prefix, "deim.greedy.rtol") ) ),
        M_Atol( doption( prefixvm( M_prefix, "deim.greedy.atol") ) ),
        M_max_value( -1 ),
        M_rebuild( boption( prefixvm( M_prefix, "deim.rebuild-db") ) ),
        M_nl_assembly(false),
        M_store_tensors( false ),
        M_write_nl_solutions( boption( prefixvm( M_prefix, "deim.elements.write") ) ),
        M_write_nl_directory( soption(prefixvm( M_prefix, "deim.elements.directory") ) ),
        M_optimized_online( boption( prefixvm( this->M_prefix, "deim.optimized-online") ) ),
        M_crb_built( false ),
        M_offline_step( true ),
        M_restart( true ),
        M_use_ser( ioption(_name="ser.eim-frequency") || ioption(_name="ser.rb-frequency") ),
        M_ser_use_rb( false ),
        M_online_model_updated( false )
    {
        using Feel::cout;
        this->setDBDirectory( modelname,uid );
        this->addDBSubDirectory( "deim"+M_prefix );

        if ( !M_rebuild )
        {
            if ( this->loadDB() )
                cout<<"DEIM : Database loaded with " << M_M << " basis functions\n";
            else
            {
                cout <<"DEIM : No Database loaded : start greedy algorithm from beginning\n";
                M_rebuild=true;
            }
        }
        else
            cout << "DEIM : option deim.rebuild-db=true : start greedy algorithm from beginning\n";
    }

    //! Destructor
    virtual ~DEIMBase()
    {}

    std::string name( bool lower=false )
    {
        std::string name = ( boost::format( "%1%DEIM%2%" ) %(is_matrix ? "M":"") %M_prefix  ).str();
        if ( lower )
            return algorithm::to_lower_copy(name);
        return name;
    }

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
    void offline() { this->run(); }

    void run()
    {
        using Feel::cout;

        if ( M_write_nl_solutions )
        {
            if ( this->worldComm().isMasterRank() )
            {
                boost::filesystem::path dir( M_write_nl_directory );
                if ( boost::filesystem::exists(dir) && boption( prefixvm( M_prefix, "deim.elements.clean-directory") ) )
                {
                    boost::filesystem::remove_all(dir);
                    boost::filesystem::create_directory(dir);
                }
                else if ( !boost::filesystem::exists(dir) )
                {
                    boost::filesystem::create_directory(dir);
                }
            }
        }
        this->worldComm().barrier();

        if ( !M_rebuild )
            return reAssembleFromDb();

        tic();

        if ( this->worldComm().isMasterRank() )
        {
            if ( !fs::exists( this->dbLocalPath() ) )
                fs::create_directories( this->dbLocalPath() );
        }


        if ( !M_trainset )
            M_trainset = M_Dmu->sampling();
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
            cout << "DEIM sampling created\n";
        }
        cout << "DEIM Offline sampling size = "<< M_trainset->size()<<std::endl;

        int sampling_size = M_trainset->size();
        if ( M_user_max>sampling_size )
        {
            cout << "DEIM : Sampling size (="<< sampling_size
                 << ") smaller than deim.dimension-max (=" << M_user_max
                 << "), dimension max is now " << sampling_size << std::endl;
            M_user_max = sampling_size;
        }


        int mMax = M_user_max;
        double error=0;
        double r_error=0;

        if ( M_restart )
        {
            if ( M_M==0 )
            {
                auto mu = M_trainset->max().template get<0>();
                cout <<"===========================================\n";
                cout << name() +" : Start algorithm with mu="<< mu.toString() <<std::endl;

                tic();
                addNewVector(mu);
                toc("Add new vector in DEIM basis");
            }
        }
        if ( M_use_ser )
        {
            if ( M_restart )
                mMax=1;
            else
            {
                mMax = M_M + M_ser_frequency;
                if ( mMax>M_user_max)
                {
                    cout << this->name() + " : max number of basis reached\n";
                    this->setOfflineStep(false);
                }
            }
        }

        while( M_M<mMax && offlineStep() )
        {
            auto best_fit = computeBestFit();
            error = best_fit.template get<1>();
            auto mu = best_fit.template get<0>();

            if ( M_max_value!=0 )
                r_error = error/M_max_value;

            cout << this->name() + " : Current max error="<<error <<", Atol="<< M_Atol
                 << ", relative max error="<< r_error <<", Rtol="<< M_tol
                 <<", for mu="<< mu.toString() <<std::endl;
            cout <<"===========================================\n";

            if ( error<M_Atol || r_error<M_tol )
            {
                cout << this->name() + " : Tolerance reached !\n";
                this->setOfflineStep(false);
                break;
            }

            cout << this->name() + " : Construction of basis "<<M_M+1<<"/"<<mMax<<", with mu="<<mu.toString()<<std::endl;

            addNewVector(mu);
        }
        M_solutions.clear();

        this->saveDB();

        cout <<"===========================================\n";
        cout << this->name() + " : Stopping greedy algorithm. Number of basis function : "<<M_M<<std::endl;

        toc(this->name() + " : Offline Total Time");
    }

    void reAssembleFromDb()
    {
        tic();
        Feel::cout << name() +" : Start reassambling "<<M_M<<" basis vectors\n";
        M_M=0;
        for ( auto const& mu : M_mus )
        {
            Feel::cout << name()+" : reassemble for mu=" << mu.toString()<<std::endl;
            tensor_ptrtype Phi = residual( mu );

            auto vec_max = vectorMaxAbs( Phi );
            auto i = vec_max.template get<1>();
            double max = vec_max.template get<0>();
            if ( M_max_value==-1 )
                M_max_value=max;
            Phi->scale( 1./max );
            M_bases.push_back( Phi );
            M_M++;
        }
        M_solutions.clear();
        cout <<"===========================================\n";
        cout << this->name() + " : Stopping greedy algorithm. Number of basis function : "<<M_M<<std::endl;
        toc( name() + " : Reassemble "+std::to_string(M_M)+" basis" );
    }

    //! \return the \f$ \beta^m(\mu)\f$ for a specific parameter \p mu
    vectorN_type beta( parameter_type const& mu, int M = -1 )
    {
        return computeCoefficient( mu, M );
    }

    vectorN_type beta( parameter_type const& mu, element_type const& u, int M = -1 )
    {
        // in this case the given element_type u is //
        // so we cannot use the online model to assemble
        // this function should not be called online !
        return computeCoefficient( mu, u, false, M );
    }
    vectorN_type beta( parameter_type const& mu, vectorN_type urb, int M=-1 )
    {
        return computeCoefficient( mu, deimExpansion(urb), true, M );
    }

    //! \Return the tensors \f$ T^m\f$ of the affine decomposition
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

    void setRB( crb_ptrtype rb )
    {
        M_crb = rb;
        M_crb_built=true;
    }

    virtual void updateRb( std::vector<element_type> const& wn )=0;

    void setOfflineStep( bool b ) { M_offline_step = b; }
    bool offlineStep() const { return M_offline_step; }

    void setRestart( bool b ) { M_restart = b; }
    bool restart() const { return M_restart; }

    void setSerFrequency( int freq ) { M_ser_frequency=freq; }
    int serFrequency() const { return M_ser_frequency; }

    void setSerUseRB( bool use_rb ) { M_ser_use_rb=use_rb; }

    //! save the database
    void saveDB() override;
    //! load the database
    bool loadDB() override;
    //!
    //! loadDB from \p filename with load strately \p l
    //!
    void loadDB( std::string const& filename, crb::load l ) override {}


protected :
    //! add a new Tensor in the base, evaluated for parameter \p mu
    void addNewVector( parameter_type const& mu )
    {
        tic();
        LOG(INFO) << this->name() + " : addNewVector() start with "<<mu.toString();
        tensor_ptrtype Phi = residual( mu );

        auto vec_max = vectorMaxAbs( Phi );
        auto i = vec_max.template get<1>();
        double max = vec_max.template get<0>();


        if ( M_max_value==-1 )
            M_max_value=max;
        M_M++;

        Phi->scale( 1./max );
        M_bases.push_back( Phi );
        M_index.push_back(i);
        M_mus.push_back(mu);

        M_B.conservativeResize(M_M,M_M);
        // update last row of M_B
        for ( int j = 0; j<M_M; j++ )
            M_B(M_M-1, j) = evaluate( M_bases[j], M_index[M_M-1]);

        //update last col of M_B
        for ( int i=0; i<M_M-1; i++ )
            M_B(i, M_M-1) = evaluate( M_bases[M_M-1], M_index[i] );

        if ( M_optimized_online )
            updateSubMesh();
        //this->saveDB();
        LOG(INFO) << this->name() + " : addNewVector() end";
        toc( this->name() +" : Add new vector" );
    }

    //! \return the value of the component \p index of the vector \p V
    double evaluate( vector_ptrtype V, int const& index, bool seq=false )
    {
        if ( seq )
            return V->operator()( index );

        V->close();
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

        M->close();
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
    vectorN_type computeCoefficient( parameter_type const& mu, int M = -1 )
    {
        tensor_ptrtype T = assemble( mu, true );
        return computeCoefficient( T, true, M );
    }
    vectorN_type computeCoefficient( parameter_type const& mu, element_type const& u, bool online=true, int M = -1 )
    {
        tensor_ptrtype T = assemble( mu, u, online );
        return computeCoefficient( T, online, M );
    }


    //! Compute the beta coefficients for a assembled tensor \p T
    vectorN_type computeCoefficient( tensor_ptrtype T, bool online=true, int M = -1 )
    {
        bool optimized = online && this->M_optimized_online;
        if( (M < 0) || (M > M_M) )
            M = M_M;
        vectorN_type rhs (M);
        vectorN_type coeff (M);
        if ( M > 0 )
        {
            for ( int i=0; i<M; i++ )
                rhs(i) = evaluate( T, optimized ? M_indexR[i]:M_index[i], optimized );
            coeff = M_B.block(0,0,M,M).fullPivLu().solve( rhs );
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
        LOG(INFO) << this->name() + " : computeBestFit() start";
        double max=0;
        auto mu_max = M_trainset->max().template get<0>();
        for ( auto const& mu : *M_trainset )
        {
            tensor_ptrtype T = residual( mu );
            auto vec_max = vectorMaxAbs(T);
            auto norm = vec_max.template get<0>();

            if ( norm>max )
            {
                max = norm;
                mu_max = mu;
            }
        }
        LOG(INFO) << this->name() + " : computeBestFit() end";
        toc(this->name() + " : compute best fit");

        return boost::make_tuple( mu_max, max );
    }

    /**
     * Evaluate the residual \f$ R_M(\mu)=T(\mu)-\hat T_M(\mu) \f$
     * where \f$ \hat T_M(\mu) \f$ is the DEIM approximation of size
     * \f$ M\f$.
     * \return a shared pointer on the Residual.
     */
    tensor_ptrtype residual( parameter_type const& mu )
    {
        LOG(INFO) << this->name() + " : residual() start with "<< mu.toString();
        tensor_ptrtype T;

        if( !M_store_tensors || M_use_ser || !M_solutions[mu.key()] )
        {
            T = this->assemble( mu );
            T->close();
            if( M_store_tensors && !M_use_ser )
            {
                LOG(INFO)<< this->name() + " : tensor stored in memory for mu="
                         << mu.toString()<<" / "<<mu.key();
                M_solutions[mu.key()] = copyTensor(T);
            }
        }
        else
        {
            LOG(INFO) << this->name() + " : tensor read in memory for mu="
                      << mu.toString()<<" / "<<mu.key();
            T = M_solutions[mu.key()];
        }

#if 0 // produce Deadlock in // with vectors.
        double norm = T->linftyNorm();
#else
        auto vec_max = vectorMaxAbs( T );
        double norm = vec_max.template get<0>();
#endif

        vectorN_type coeff = computeCoefficient(T, false);

        auto newT = copyTensor( T );
        for ( int i=0; i<M_M; i++ )
            add( newT, -coeff(i), M_bases[i] );
        newT->scale( 1./norm );

        LOG(INFO) << this->name() + " : residual() end";
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
        //vector_ptrtype newV = V->clone();
        vector_ptrtype newV = backend()->newVector( V->mapPtr() );
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


    virtual void updateSubMesh()=0;

    virtual element_type deimExpansion( vectorN_type const& urb )=0;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & __ar, const unsigned int __version );

protected :
    std::string M_prefix;
    parameterspace_ptrtype M_Dmu;
    std::vector<parameter_type> M_mus;

    sampling_ptrtype M_trainset;
    int M_M, M_user_max, M_n_rb;
    double M_tol, M_Atol, M_max_value;
    matrixN_type M_B;
    std::vector< tensor_ptrtype > M_bases;

    std::vector<indice_type> M_index, M_indexR, M_ldofs;
    std::vector<ptsset_type> M_ptsset;
    std::set<int> M_elts_ids;

    solutionsmap_type M_solutions;

    bool M_rebuild, M_nl_assembly, M_store_tensors, M_write_nl_solutions;
    std::string M_write_nl_directory;
    bool M_optimized_online; // to be removed when optimized is ok

    crb_ptrtype M_crb;
    bool M_crb_built,M_offline_step, M_restart,M_use_ser,M_ser_use_rb;
    int M_ser_frequency;

    std::vector<element_type> M_rb_basis;
    bool M_online_model_updated;
    space_ptrtype Rh;
};



template <typename ModelType,
          typename TensorType>
class DEIMModel :
        public DEIMBase<typename ModelType::parameterspace_type, typename ModelType::space_type,
                        TensorType>
{
public :
    typedef DEIMBase<typename ModelType::parameterspace_type, typename ModelType::space_type, TensorType> super_type;
    typedef ModelType model_type;
    typedef boost::shared_ptr<model_type> model_ptrtype;
    typedef typename model_type::parameterspace_type parameterspace_type;
    typedef typename super_type::parameterspace_ptrtype parameterspace_ptrtype;
    typedef typename super_type::parameter_type parameter_type;
    typedef typename super_type::sampling_ptrtype sampling_ptrtype;
    typedef typename super_type::tensor_ptrtype tensor_ptrtype;
    typedef typename super_type::element_type element_type;
    typedef typename super_type::mesh_type mesh_type;
    typedef typename super_type::space_type space_type;


    DEIMModel() :
        super_type()
    {}

    DEIMModel( model_ptrtype model, sampling_ptrtype sampling=nullptr, std::string prefix="" ) :
        super_type( model->functionSpace(),
                    model->parameterSpace(),
                    sampling,
                    model->uuid(),
                    model->modelName(),
                    prefix ),
        M_model( model )
    {
        if ( this->M_optimized_online )
        {
            this->M_online_model = model_ptrtype( new model_type() );
            this->M_online_model->setModelOnlineDeim( prefixvm(this->M_prefix,"deim-online") );
        }
        else
            this->M_online_model = M_model;

        if ( Rh )
            M_online_model->setFunctionSpaces( Rh );
    }

    virtual ~DEIMModel()
    {}

    void init()
    {
        if ( this->M_rebuild )
        {
            auto mu = this->M_trainset->max().template get<0>();
            auto T = this->assemble(mu);
            if (!T)
            {
                this->M_nl_assembly=true;
                auto u = M_model->functionSpace()->element();
                auto Tnl = this->assemble(mu,u);
                CHECK( Tnl ) << "You want to use DEIM but you did not implement assmbleForDEIM functions\n";
            }
        }
    }

    tensor_ptrtype assemble( parameter_type const& mu, bool online=false ) override
    {
        if ( this->M_nl_assembly )
        {
            if ( online )
                Feel::cout << this->name() + " : WARNING : Call of online nl assembly with no solution u\n";
            //CHECK(!online) << "Call of online nl assembly with no solution u\n";

            auto u = M_model->functionSpace()->element();
            bool need_solve = true;

            if ( this->M_ser_use_rb && this->M_crb && !online )
            {
                std::vector<vectorN_type> uN, uNdu, uNold, uNduold;
                auto o = this->M_crb->lb( this->M_crb->dimension(), mu, uN, uNdu , uNold, uNduold );
                int size = uN.size();

                u = this->M_crb->expansion( uN[size-1], this->M_crb->dimension(), false );
            }
            else
            {
                if ( this->M_ser_use_rb && !this->M_crb )
                    Feel::cout <<this->name() + " WARNING : Suppose to use crb expansion with no crb class ! u will be computed using model->solve\n";

                if ( M_write_nl_solutions )
                {
                    need_solve = !u.load( _path=M_write_nl_directory,
                                          _suffix=std::to_string(mu.key()), _type="hdf5" );
                    if ( need_solve )
                        LOG(INFO) << this->name() + " : Unable to load nl solution in direcotry "
                                  << M_write_nl_directory << ", for parameter : " << mu.toString()
                                  <<" / " << mu.key()<< ". Solve function will be called.";
                    else
                        LOG(INFO) << this->name() + " : NL solution loaded in direcotry "
                                  << M_write_nl_directory << ", for parameter : " << mu.toString()
                                  <<" / " << mu.key();
                }

                if ( need_solve )
                {
                    LOG(INFO) << this->name() + " : calling solve function for parameter " << mu.toString()
                              <<" / " << mu.key();
                    u = M_model->solve(mu);

                    if ( M_write_nl_solutions )
                    {
                        LOG(INFO) << this->name() + " : Wrting solution on disk in directory "
                                  << M_write_nl_directory << ", for parameter : " << mu.toString()
                                  <<" / " << mu.key();
                        u.save( _path=M_write_nl_directory,
                                _suffix=std::to_string(mu.key()), _type="hdf5" );
                    }
                }

            }

            return modelAssemble(mu,u);
        }
        return modelAssemble(mu,online);
    }

    tensor_ptrtype assemble( parameter_type const& mu, element_type const& u, bool online=false ) override
    {
        CHECK(this->M_nl_assembly) << "You called nl coefficient for DEIM but you implemented assembleForDEIM(mu)\n";
        return modelAssemble(mu,u,online);
    }

    void updateRb( std::vector<element_type> const& wn ) override
    {
        if ( this->M_optimized_online )
        {
            // if the interpolation space has been updated we have to reproject
            // all the rb basis vectors on the new interpolation space
            if ( this->M_online_model_updated )
            {
                for ( int i=0; i<M_rb_basis.size(); i++ )
                {
                    M_rb_basis[i] = Rh->element();
                    M_rb_basis[i].on( elements(Rh->mesh()), idv(wn[i]) );
                }

                this->M_online_model_updated=false;
            }

            // we now project the new basis vector on the interpolation space and we stock them.
            for ( int i=M_rb_basis.size(); i<wn.size(); i++ )
            {
                M_rb_basis.push_back( Rh->element() );
                M_rb_basis[i].on( elements(Rh->mesh()), idv(wn[i]) );
            }
        }
        else // to be removed when optimized is fixed
            for ( int i=M_rb_basis.size(); i<wn.size(); i++ )
                M_rb_basis.push_back( wn[i] );

        this->M_n_rb = M_rb_basis.size();
    }

    void gatherContributions( element_type& u )
    {

    }

    element_type deimExpansion( vectorN_type const& urb ) override
    {
        int N = urb.size();
        FEELPP_ASSERT( N <= M_rb_basis.size() )( N )( M_rb_basis.size() ).error( "invalid expansion size ( N and M_rb_basis ) ");
        return Feel::expansion( M_rb_basis, urb, N );
    }

    model_ptrtype & onlineModel()
    {
        return M_online_model;
    }

protected :
    virtual tensor_ptrtype modelAssemble( parameter_type const& mu, bool online=false )=0;
    virtual tensor_ptrtype modelAssemble( parameter_type const& mu, element_type const& u, bool online=false )=0;

protected :
    model_ptrtype M_model, M_online_model;

    using super_type::M_write_nl_solutions;
    using super_type::M_write_nl_directory;
    using super_type::M_rb_basis;
    using super_type::Rh;
};


template <typename ModelType,
          typename TestSpace=typename ModelType::space_type>
class DEIM :
        public DEIMModel<ModelType,typename Backend<typename ModelType::value_type>::vector_type>
{
public :
    typedef DEIMModel<ModelType,typename Backend<typename ModelType::value_type>::vector_type>  super_type;

    typedef ModelType model_type;
    typedef boost::shared_ptr<model_type> model_ptrtype;
    typedef typename super_type::parameter_type parameter_type;
    typedef typename super_type::parameterspace_ptrtype parameterspace_ptrtype;
    typedef typename super_type::sampling_ptrtype sampling_ptrtype;
    typedef typename super_type::tensor_ptrtype vector_ptrtype;
    typedef typename super_type::element_type element_type;
    typedef typename super_type::mesh_type mesh_type;
    typedef typename super_type::space_type space_type;
    typedef TestSpace testspace_type;
    typedef boost::shared_ptr<testspace_type> testspace_ptrtype;

    DEIM() :
        super_type()
    {}

    DEIM( model_ptrtype model, testspace_ptrtype test_space,
          sampling_ptrtype sampling=nullptr, std::string prefix="" ) :
        super_type( model, sampling, prefix ),
        Teh( test_space )
    {
        this->M_store_tensors = boption( prefixvm( this->M_prefix, "deim.store-vectors") );
        this->init();
    }

    ~DEIM()
    {}

private :
    vector_ptrtype modelAssemble( parameter_type const& mu, bool online=false ) override
    {
        if ( online )
            return this->M_online_model->assembleForDEIM(mu);
       return this->M_model->assembleForDEIM(mu);
    }

    vector_ptrtype modelAssemble( parameter_type const& mu, element_type const& u, bool online=false ) override
    {
        if ( online )
            return this->M_online_model->assembleForDEIMnl(mu,u);
        return this->M_model->assembleForDEIMnl(mu,u);
    }


private :
    virtual void updateSubMesh() override
    {
        // Last added index
        auto index = this->M_index.back();
        //auto Xh = this->M_model->functionSpace();
        auto mesh = Teh->mesh();
        int proc_n = Teh->dof()->procOnGlobalCluster(index);

        // recover the elements which share this index
        std::set<int> pts_ids;
        int ldof=-1;

        if ( Environment::worldComm().globalRank()==proc_n )
        {
            for ( auto const& dof : Teh->dof()->globalDof(index-Teh->dof()->firstDofGlobalCluster()))
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
            auto e = Teh->dof()->searchGlobalProcessDof( index );
            if ( e.template get<0>() )
            {
                auto process_dof = e.template get<1>();
                for ( auto const& dof : Teh->dof()->globalDof( process_dof) )
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
        saveGMSHMesh( _mesh=submesh, _filename=this->name(true)+"-submesh.msh" );
        auto seqmesh = loadMesh( _mesh=new mesh_type,
                                 _filename=this->name(true)+"-submesh.msh",
                                 _worldcomm= Environment::worldCommSeq() );
        Rh = space_type::New( seqmesh,
                              _worldscomm=std::vector<WorldComm>(space_type::nSpaces,Environment::worldCommSeq()) );
        testspace_ptrtype RTeh = testspace_type::New( seqmesh,
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
            CHECK( map_it!=elts_map.end() ) <<this->name() + " : elt id not found in map, on proc : "
                                            <<Environment::worldComm().globalRank() << std::endl;
            int elt_idR = map_it->second;
            this->M_indexR.push_back( RTeh->dof()->localToGlobalId( elt_idR, this->M_ldofs[i] ) );
        }

        this->M_online_model->setFunctionSpaces( Rh );

        this->M_online_model_updated = true;
    }

private :
    testspace_ptrtype Teh;
    using super_type::Rh;
};


    template <typename ModelType>
class MDEIM :
        public DEIMModel<ModelType,
                         typename Backend<typename ModelType::value_type>::sparse_matrix_type>
{
public :
    typedef DEIMModel<ModelType, typename Backend<typename ModelType::value_type>::sparse_matrix_type> super_type;

    typedef ModelType model_type;
    typedef boost::shared_ptr<model_type> model_ptrtype;
    typedef typename super_type::parameter_type parameter_type;
    typedef typename super_type::parameterspace_ptrtype parameterspace_ptrtype;
    typedef typename super_type::sampling_ptrtype sampling_ptrtype;
    typedef typename super_type::tensor_ptrtype sparse_matrix_ptrtype;
    typedef typename super_type::element_type element_type;
    typedef typename super_type::mesh_type mesh_type;
    typedef typename super_type::space_type space_type;

    MDEIM() :
        super_type()
    {}

    MDEIM( model_ptrtype model, sampling_ptrtype sampling=nullptr, std::string prefix="" ) :
        super_type( model, sampling, prefix )
    {
        this->M_store_tensors = boption( prefixvm( this->M_prefix, "deim.store-matrices") );
        this->init();
    }

    ~MDEIM()
    {}

private :
    sparse_matrix_ptrtype modelAssemble( parameter_type const& mu, bool online=false ) override
    {
        if ( online )
            return this->M_online_model->assembleForMDEIM(mu);
        return this->M_model->assembleForMDEIM(mu);
    }

    sparse_matrix_ptrtype modelAssemble( parameter_type const& mu, element_type const& u, bool online=false ) override
    {
        if ( online )
            return this->M_online_model->assembleForMDEIMnl(mu,u);
        return this->M_model->assembleForMDEIMnl(mu,u);
    }

private :
    void updateSubMesh() override
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
        saveGMSHMesh( _mesh=submesh, _filename=this->name(true)+"-submesh.msh" );
        auto seqmesh = loadMesh( _mesh=new mesh_type,
                                 _filename=this->name(true)+"-submesh.msh",
                                 _worldcomm= Environment::worldCommSeq() );
        Rh = space_type::New( seqmesh,
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
            CHECK( map_it!=elts_map.end() ) << this->name()+" : elt id not found in map, on proc : "
                                            << Environment::worldComm().globalRank() << std::endl;

            int elt_idR1 = map_it->second;
            int elt_idR2 = elt_idR1;
            if ( ptsset1!=ptsset2 )
            {
                map_it = elts_map.find( ptsset2 );
                CHECK( map_it!=elts_map.end() ) << this->name()+" : elt id not found in map, on proc : "
                                                << Environment::worldComm().globalRank() << std::endl;
                elt_idR2 = map_it->second;
            }

            auto indexR1 = Rh->dof()->localToGlobalId( elt_idR1, this->M_ldofs[i].first );
            auto indexR2 = Rh->dof()->localToGlobalId( elt_idR2, this->M_ldofs[i].second );

            this->M_indexR.push_back( std::make_pair(indexR1, indexR2) );
        }

        this->M_online_model->setFunctionSpaces( Rh );

        this->M_online_model_updated = true;
    }

private :
    using super_type::Rh;
};



template <typename ParameterSpaceType, typename SpaceType, typename TensorType>
template<class Archive>
void
DEIMBase<ParameterSpaceType,SpaceType,TensorType>::serialize(Archive & __ar, const unsigned int __version )
{
    __ar & BOOST_SERIALIZATION_NVP( M_M );
    DVLOG(2) << "M_M saved/loaded\n";
    __ar & BOOST_SERIALIZATION_NVP( M_B );
    DVLOG(2) << "M_B saved/loaded\n";
    __ar & BOOST_SERIALIZATION_NVP( M_nl_assembly );
    DVLOG(2) << "M_nl_assembly saved/loaded\n";
    __ar & BOOST_SERIALIZATION_NVP( M_indexR );
    DVLOG(2) << "M_indexR saved/loaded\n";
    __ar & BOOST_SERIALIZATION_NVP( M_index );
    DVLOG(2) << "M_index saved/loaded\n";
    __ar & BOOST_SERIALIZATION_NVP( M_n_rb );
    DVLOG(2) << "M_n_rb saved/loaded\n";
    __ar & BOOST_SERIALIZATION_NVP( M_mus );
    DVLOG(2) << "M_mus saved/loaded\n";

    if ( Archive::is_loading::value )
    {
        if ( M_rb_basis.size()==0 && Rh )
        {
            for(int i = 0; i < M_n_rb; i++ )
                M_rb_basis.push_back( Rh->element() );
            for( int i = 0; i < M_n_rb; ++ i )
                __ar & BOOST_SERIALIZATION_NVP( M_rb_basis[i] );
            DVLOG(2) << "M_rb_basis saved/loaded\n";
        }
    }
    else
    {
        for( int i = 0; i < M_n_rb; ++ i )
            __ar & BOOST_SERIALIZATION_NVP( M_rb_basis[i] );
    }

}


template <typename ParameterSpaceType, typename SpaceType, typename TensorType>
void
DEIMBase<ParameterSpaceType,SpaceType,TensorType>::saveDB()
{
    if ( this->worldComm().isMasterRank() )
    {
        if ( Rh )
        {
            fs::path mesh_name( this->name(true)+"-submesh.msh" );
            auto filename = this->dbLocalPath() / mesh_name;
            saveGMSHMesh( _mesh=Rh->mesh(), _filename=filename.string() );
        }

        fs::ofstream ofs( this->dbLocalPath() / this->dbFilename() );
        boost::archive::binary_oarchive oa( ofs );
        oa << *this;
    }
}


template <typename ParameterSpaceType, typename SpaceType, typename TensorType>
bool
DEIMBase<ParameterSpaceType,SpaceType,TensorType>::loadDB()
{
    if ( Rh )
    {
        return true;
    }

    fs::path db = this->lookForDB();
    if ( db.empty() || !fs::exists( db ) )
    {
        return false;
    }

    fs::ifstream ifs( db );
    if ( ifs )
    {
        fs::path mesh_name( this->name(true)+"-submesh.msh" );
        auto filename = this->dbLocalPath() / mesh_name;
        if ( fs::exists( filename ))
        {
            auto seqmesh = loadMesh( _mesh=new mesh_type,
                                     _filename=filename.string(),
                                     _worldcomm= Environment::worldCommSeq() );
            Rh = space_type::New( seqmesh,
                                  _worldscomm=std::vector<WorldComm>(space_type::nSpaces,Environment::worldCommSeq()) );
        }

        boost::archive::binary_iarchive ia( ifs );
        ia >> *this;
        this->setIsLoaded( true );

        return true;
    }

    return false;
}

namespace detail
{
template <typename Args>
struct compute_deim_return
{
    typedef typename boost::remove_reference<typename boost::remove_pointer<typename parameter::binding<Args, tag::model>::type>::type>::type::element_type model1_type;
    typedef typename boost::remove_const<typename boost::remove_pointer<model1_type>::type>::type model_type;
    typedef typename boost::remove_reference<typename parameter::binding<Args, tag::test, typename model_type::space_ptrtype>::type>::type::element_type test_type;

    typedef DEIM<model_type> type;
    typedef boost::shared_ptr<type> ptrtype;
};

template <typename Args>
struct compute_mdeim_return
{
    typedef typename boost::remove_reference<typename boost::remove_pointer<typename parameter::binding<Args, tag::model>::type>::type>::type::element_type model1_type;
    typedef typename boost::remove_const<typename boost::remove_pointer<model1_type>::type>::type model_type;

    typedef MDEIM<model_type> type;
    typedef boost::shared_ptr<type> ptrtype;
};

}

BOOST_PARAMETER_FUNCTION(
                         ( typename Feel::detail::compute_deim_return<Args>::ptrtype ), // 1. return type
                         deim,                        // 2. name of the function template
                         tag,                                        // 3. namespace of tag types
                         ( required
                           ( in_out(model),          * )
                           ) // required
                         ( optional
                           ( sampling, *, nullptr )
                           ( prefix, *( boost::is_convertible<mpl::_,std::string> ), "" )
                           ( test, *, model->functionSpace() )
                           ) // optionnal
                         )
{
    typedef typename Feel::detail::compute_deim_return<Args>::type deim_type;
    return boost::make_shared<deim_type>( model,test, sampling, prefix );
}


BOOST_PARAMETER_FUNCTION(
                         ( typename Feel::detail::compute_mdeim_return<Args>::ptrtype ), // 1. return type
                         mdeim,                        // 2. name of the function template
                         tag,                                        // 3. namespace of tag types
                         ( required
                           ( in_out(model),          * )
                           ) // required
                         ( optional
                           ( sampling, *, nullptr )
                           ( prefix, *( boost::is_convertible<mpl::_,std::string> ), "" )
                           ) // optionnal
                         )
{
    typedef typename Feel::detail::compute_mdeim_return<Args>::type mdeim_type;
    return boost::make_shared<mdeim_type>( model, sampling, prefix );
}

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

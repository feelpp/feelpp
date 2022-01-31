/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): JB WAHL <wahl.jb@gmail.com>
       Date: 2018-03-27

  Copyright (C) 2018 Feel++ Consortium

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
   \file deimbase.hpp
   \author JB Wahl
   \date 2018-03-27
*/

#ifndef _FEELPP_DEIMBASE_HPP
#define _FEELPP_DEIMBASE_HPP 1

#include <limits>
#include <numeric>
#include <string>
#include <iostream>
#include <fstream>

#include <Eigen/Core>

#include <boost/ref.hpp>
#include <boost/next_prior.hpp>
#include <boost/type_traits.hpp>
#include <boost/tuple/tuple.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feelalg/topetsc.hpp>

#include <feel/feelmor/parameterspace.hpp>
#include <feel/feelmor/crbdb.hpp>
#include <feel/feelmor/crb.hpp>
#include <feel/feelmor/crbmodel.hpp>

#include <feel/feelfilters/savegmshmesh.hpp>
#include <feel/feelfilters/loadmesh.hpp>

#include <feel/feeldiscr/reducedbasisspace.hpp>
#include <feel/feeldiscr/twospacesmap.hpp>



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
    typedef CRBDB super;
    typedef DEIMBase<ParameterSpaceType,SpaceType,TensorType> self_type;

public :
    typedef typename TensorType::value_type value_type;

    typedef ParameterSpaceType parameterspace_type;
    typedef typename std::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef typename parameterspace_type::element_type parameter_type;
    typedef typename parameterspace_type::sampling_type sampling_type;
    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;

    typedef TensorType tensor_type;
    typedef std::shared_ptr<tensor_type> tensor_ptrtype;

    typedef Backend<value_type> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;
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
    typedef std::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef std::shared_ptr<element_type> element_ptrtype;
    typedef typename space_type::mesh_type mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    typedef ReducedBasisSpace<space_type> rbspace_type;
    typedef std::shared_ptr<rbspace_type> rbspace_ptrtype;

    typedef std::map<int,tensor_ptrtype> solutionsmap_type;

    typedef CRBBase<space_type,parameterspace_type> crb_type;
    typedef std::shared_ptr<crb_type> crb_ptrtype;

    typedef TwoSpacesMap<space_type> spaces_map_type;
    typedef std::shared_ptr<spaces_map_type> spaces_map_ptrtype;


    static const int n_spaces = space_type::nSpaces;
    template< int T >
    using sub_space = typename space_type::template sub_functionspace_type<T>;
    typedef typename mpl::range_c< int, 0, n_spaces > rangespace_type;
    template < typename  T >
    struct SubElementVec
    {
        typedef std::vector<typename sub_space<T::value>::element_ptrtype> type;
    };
    typedef typename mpl::transform< rangespace_type, SubElementVec<mpl::_1>, mpl::back_inserter<fusion::vector<> > >::type block_rb_basis_type;
    template <int T>
    using sub_rb_type = typename fusion::result_of::at_c<block_rb_basis_type,T>::type;


    //! Default Constructor
    DEIMBase()
        {}

protected :
    //! Constructor from command line options
    DEIMBase(  space_ptrtype Xh, parameterspace_ptrtype Dmu, sampling_ptrtype sampling,
               uuids::uuid const& uid, std::string const& modelname,
               std::string prefix, std::string const& dbfilename, std::string const& dbdirectory,
               worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
               std::string const& prefixModel = "" );

public :
    //! Destructor
    virtual ~DEIMBase()
        {}

    //! the name of this instance
    std::string name( bool lower=false )
        {
            std::string name = ( boost::format( "%1%DEIM%2%%3%" ) %(is_matrix ? "M":"") %M_prefix %tag() ).str();
            if ( lower )
                return algorithm::to_lower_copy(name);
            return name;
        }

    virtual int tag()=0;
    /**
     * Run the greedy algotithm and build the DEIM basis.
     * Greedy algorithm stops when the affine decomposition is exact
     * or when the maximum number of terms is the decomposition is
     * reached. This maximum is defined by the option
     * deim.dimension-max
     */
    void offline() { this->run(); }
    void run();

    //! \return the \f$ \beta^m(\mu)\f$ for a specific parameter \p mu
    vectorN_type beta( parameter_type const& mu, int M = -1 )
        {
            if ( mu!=M_last_mu || ( ((M < 0) || (M > M_M)) && M_last_beta.size() != M_M ) || M_last_beta.size() != M )
            {
                M_last_mu=mu;
                M_last_beta=computeCoefficient( mu, M );
            }
            return M_last_beta;
            //return computeCoefficient( mu, M );
        }

    //! \return the \f$ \beta^m(\mu)\f$ for a specific parameter \p mu and a RB solution \p urb
    vectorN_type beta( parameter_type const& mu, vectorN_type urb, int M=-1 )
        { return computeCoefficient( mu, deimExpansion(urb), true, M ); }

    //! \return the \f$ \beta^m(\mu)\f$ for a specific parameter \p mu and a FE solution \p u
    vectorN_type beta( parameter_type const& mu, element_type const& u, int M = -1 )
        {
            auto ur = Rh->element();
            M_map->project( ur, u );
            return computeCoefficient( mu, ur, true, M );
        }

    //! \Return the basis tensors \f$ T^m\f$ of the affine approximation
    std::vector<tensor_ptrtype> q() const
        { return M_bases; }
    //! \Return the m-th basis tensor \f$ T^m\f$ of the affine approximation
    tensor_ptrtype q( int m ) const
        { return M_bases[m]; }

    //! \return the current size of the affine decompostion
    int size() const { return M_M; }

    /**
     * Update the RB stored in DEIM. This RB is projected on the interpolation mesh
     * \p XN the reduced basis space containing the RB
     */
    virtual void updateRb( rbspace_ptrtype const& XN, std::vector<std::vector<int>> const& subN )=0;

    //! Set the value of M_offline_step. Used with SER
    void setOfflineStep( bool b ) { M_offline_step = b; }
    //! Get the value of M_offline_step. Used with SER
    bool offlineStep() const { return M_offline_step; }

    //! Set the value of M_restart. Used with SER
    void setRestart( bool b ) { M_restart = b; }
    //! Get the value of M_restart. Used with SER
    bool restart() const { return M_restart; }

    //! Set the value of M_ser_frequency. Used with SER
    void setSerFrequency( int freq ) { M_ser_frequency=freq; }
    //! Get the value of M_ser_frequency. Used with SER
    int serFrequency() const { return M_ser_frequency; }

    void setSerUseRB( bool use_rb ) { M_ser_use_rb=use_rb; }
    void setRB( crb_ptrtype rb )
        {
            M_crb = rb;
            M_crb_built=true;
        }



    //! save the database
    void saveDB() override;
    //! load the database
    bool loadDB() override;
    //!
    //! loadDB from \p filename with load strately \p l
    //!
    void loadDB( std::string const& filename, crb::load l ) override {}

    std::string const& prefixModel() const { return M_prefixModel; }
    void setPrefixModel( std::string const& prefix ) { M_prefixModel = prefix; }

    std::vector<indice_type> const& index() const { return M_index; }
    std::vector<indice_type> const& indexR() const { return M_indexR; }
    std::vector<parameter_type> const& mus() const { return M_mus; }

protected :
    /**
     * The assemble function calls the function provided by the model. This
     * function will assemble the tensor \f$ T(\mu) \f$.
     * \p mu the considered parameter
     * \p online =true if the assembly is made on the online model
     * \return the tensor \f$ T(\mu) \f$ as vector_ptrtype
     */
    virtual tensor_ptrtype assemble( parameter_type const& mu, bool online=false, bool force_fem=false )=0;
    virtual tensor_ptrtype assemble( parameter_type const& mu, element_type const& u, bool online=false )=0;

    /**
     * Reassemble the DEIM basis from the database
     * This method is called when the option deim.rebuild-database=false
     */
    void reAssembleFromDb();

    /**
     * \return the T-th projected reduced basis.
     * Used for non-linear problem with block structure
     */
    template<int T>
    sub_rb_type<T>& subRb() { return fusion::at_c<T>( M_block_rb ); }

    //! add a new Tensor in the base, evaluated for parameter \p mu
    void addNewVector( parameter_type const& mu );

    //! \return the value of the component \p index of the vector \p V
    double evaluate( vector_ptrtype V, int const& index, bool seq=false );

    //! \return the value of the entry \p idx of the matrix \p M
    double evaluate( sparse_matrix_ptrtype M, std::pair<int,int> const& idx, bool seq=false );

    /**
     * Evaluate the maximum of the components of \p V in absolute
     * value.
     *
     * \return a couple (max,index) where max is the maximum in
     * absolute value and index is the index of the component where
     * the max is reached
     */
    vectormax_type vectorMaxAbs( vector_ptrtype V );

    /**
     * Evaluate the maximum of the entry of \p M in absolute
     * value.
     *
     * \return a couple (max,index) where max is the maximum in
     * absolute value and index is the couple of index of the entry
     * where the max is reached
     */
    matrixmax_type vectorMaxAbs( sparse_matrix_ptrtype M );

    /**
     * Choose the best parameter to build the next basis vector.
     *
     * For each parameter mu in sampling, commpute the residual \f$
     * R_M(\mu)=T(\mu)-\hat T_M(\mu) \f$ where \f$ \hat T_M(\mu) \f$
     * is the DEIM approximation of size \f$ M\f$.
     * \return a couple (mu_max, max). mu_max=\f$ argmax_{\mu\in S}R(\mu)\f$
     * and max=\f$ max_{\mu\in S}R(\mu)\f$
     */
    bestfit_type computeBestFit();


    /**
     * Evaluate the residual \f$ R_M(\mu)=T(\mu)-\hat T_M(\mu) \f$
     * where \f$ \hat T_M(\mu) \f$ is the DEIM approximation of size
     * \f$ M\f$.
     * \return a shared pointer on the Residual.
     */
    tensor_ptrtype residual( parameter_type const& mu, bool force_fem=false );


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
        if( (M < 0) || (M > M_M) )
            M = M_M;
        vectorN_type rhs (M);
        vectorN_type coeff (M);
        if ( M > 0 )
        {
            for ( int i=0; i<M; i++ )
                rhs(i) = evaluate( T, online ? M_indexR[i]:M_index[i], online );
            coeff = M_B.block(0,0,M,M).fullPivLu().solve( rhs );
        }

        return coeff;
    }


    //! evaluate V= V+a*vec
    void add( vector_ptrtype V, double const& a, vector_ptrtype vec ) { V->add( a, vec ); }
    //! evaluate M= M+a*mat
    void add( sparse_matrix_ptrtype M, double const& a, sparse_matrix_ptrtype mat )
        { M->addMatrix( a, mat ); }

    //! \return a shared pointer on a copy of \p V
    vector_ptrtype copyTensor( vector_ptrtype V )
        {
            //vector_ptrtype newV = V->clone();
            vector_ptrtype newV = backend()->newVector( V->mapPtr() );
            *newV = *V;
            return newV;
        }
    //! \return a shared pointer on a copy of \p M
    sparse_matrix_ptrtype copyTensor( sparse_matrix_ptrtype M )
        {
            sparse_matrix_ptrtype newM = backend()->newMatrix( M->mapColPtr(), M->mapRowPtr(), M->graph() );
            *newM=*M;
            return newM;
        }

    /**
     * Extract a submesh from the full dimension mesh.
     * The submesh only contained the elements related to the dofs selectioned in the greedy algorithm
     * Also build a interpolation space based on this submesh.
     * The interpolation mesh is then given to the online model
     */
    virtual void updateSubMesh()=0;

    //! Rebuild the correspondency map between the interpolation space Rh and the FE space Xh
    virtual void rebuildMap()=0;

    //! Perform an expansion of a RB solution using the projected basis stored in DEIM
    virtual element_type deimExpansion( vectorN_type const& urb )=0;


    virtual space_ptrtype newInterpolationSpace( mesh_ptrtype const& mesh )=0;


    //! Serialization method used to save the DB
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & __ar, const unsigned int __version );

protected :
    std::string M_prefix;
    std::string M_prefixModel;
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
    std::vector<int> M_n_block_rb;

    solutionsmap_type M_solutions;

    bool M_rebuild, M_nl_assembly, M_store_tensors, M_write_nl_solutions, M_last_solve_is_ok;
    std::string M_write_nl_directory;

    crb_ptrtype M_crb;
    bool M_crb_built,M_offline_step, M_restart,M_use_ser,M_ser_use_rb;
    int M_ser_frequency;

    std::vector<element_ptrtype> M_rb;
    block_rb_basis_type M_block_rb;
    space_ptrtype Rh;

    spaces_map_ptrtype M_map;
    parameter_type M_last_mu;
    vectorN_type M_last_beta;

}; // Class DEIMBase


template <typename ParameterSpaceType, typename SpaceType, typename TensorType>
DEIMBase<ParameterSpaceType,SpaceType,TensorType>::DEIMBase(  space_ptrtype Xh, parameterspace_ptrtype Dmu, sampling_ptrtype sampling, uuids::uuid const& uid, std::string const& modelname, std::string prefix, std::string const& dbfilename, std::string const& dbdirectory, worldcomm_ptr_t const& worldComm, std::string const& prefixModel ) :
    super ( ( boost::format( "%1%DEIM%2%" ) %(is_matrix ? "M":"") %prefix  ).str(),
            "deim",
            uid,
            worldComm ),
    M_prefix( prefix ),
    M_prefixModel(prefixModel),
    M_Dmu( Dmu ),
    M_trainset( sampling ),
    M_M(0),
    M_user_max( ioption(  prefixvm( M_prefix, "deim.dimension-max" ) ) ),
    M_n_rb( 0 ),
    M_tol( doption( prefixvm( M_prefix, "deim.greedy.rtol") ) ),
    M_Atol( doption( prefixvm( M_prefix, "deim.greedy.atol") ) ),
    M_max_value( -1 ),
    M_rebuild( boption( prefixvm( M_prefix, "deim.rebuild-database") ) ),
    M_nl_assembly(false),
    M_store_tensors( false ),
    M_write_nl_solutions( boption( prefixvm( M_prefix, "deim.elements.write") ) ),
    M_last_solve_is_ok( true ),
    M_write_nl_directory( soption(prefixvm( M_prefix, "deim.elements.directory") ) ),
    M_crb_built( false ),
    M_offline_step( true ),
    M_restart( true ),
    M_use_ser( ioption(_prefix=M_prefixModel,_name="ser.eim-frequency") || ioption(_prefix=M_prefixModel,_name="ser.rb-frequency") ),
    M_ser_use_rb( false ),
    M_map( new spaces_map_type ),
    M_last_mu( Dmu->element() )
{
    using Feel::cout;
    LOG(INFO) <<"DEIMBase constructor begin\n";

    if ( dbfilename.empty() )
    {
        this->setDBDirectory( modelname,uid );
        this->addDBSubDirectory( "deim"+M_prefix );
    }
    else
    {
        fs::path dbfilenamePath = dbfilename;
        this->setDBFilename( dbfilenamePath.filename().string() );
        if ( dbfilenamePath.is_relative() )
        {
            std::string mysubDir;
            if ( !dbfilenamePath.parent_path().empty() )
                mysubDir = dbfilenamePath.parent_path().string();

            fs::path dbdirectoryPath = dbdirectory;
            if ( !dbdirectoryPath.is_absolute() )
                dbdirectoryPath = fs::absolute( dbdirectoryPath );
            this->setDBDirectory( dbdirectoryPath.string() );

            dbfilenamePath = (dbdirectoryPath/dbfilenamePath).string();

            if ( !mysubDir.empty() )
                this->addDBSubDirectory( mysubDir );
        }
        else
            this->setDBDirectory( dbfilenamePath.parent_path().string() );
    }


    LOG(INFO) <<"DEIMBase constructor end\n";
} // DEIMBase Constructor


template <typename ParameterSpaceType, typename SpaceType, typename TensorType>
void
DEIMBase<ParameterSpaceType,SpaceType,TensorType>::run()
{
    using Feel::cout;
    LOG(INFO) <<"DEIM run function begin\n";

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
    cout << name() + " Offline sampling size = "<< M_trainset->size()<<std::endl;

    int sampling_size = M_trainset->size();
    if ( M_user_max>sampling_size )
    {
        cout << name()+" : Sampling size (="<< sampling_size
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
            addNewVector(mu);
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

    updateSubMesh();

    this->saveDB();

    cout <<"===========================================\n";
    cout << this->name() + " : Stopping greedy algorithm. Number of basis function : "<<M_M<<std::endl;

    LOG(INFO) <<"DEIM run function end\n";
    toc(this->name() + " : Offline Total Time");
} // run


template <typename ParameterSpaceType, typename SpaceType, typename TensorType>
void
DEIMBase<ParameterSpaceType,SpaceType,TensorType>::reAssembleFromDb()
{
    tic();
    LOG(INFO) << "DEIM reAssemblefromdb begin\n";
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
        if ( M_M>=M_user_max)
            break;
    }
    M_solutions.clear();
    rebuildMap();
    cout <<"===========================================\n";
    cout << this->name() + " : Stopping greedy algorithm. Number of basis function : "<<M_M<<std::endl;
    LOG(INFO) << "DEIM reAssemblefromdb end\n";
    toc( name() + " : Reassemble "+std::to_string(M_M)+" basis" );
} // reAssembleFromDb


template <typename ParameterSpaceType, typename SpaceType, typename TensorType>
void
DEIMBase<ParameterSpaceType,SpaceType,TensorType>::addNewVector( parameter_type const& mu )
{
    tic();
    LOG(INFO) << this->name() + " : addNewVector() start with "<<mu.toString();
    bool forced_fem = !boption(_prefix=M_prefixModel, _name="ser.use-rb-in-eim-basis-build" );
    tensor_ptrtype Phi = residual( mu, forced_fem );

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

    //this->saveDB();
    LOG(INFO) << this->name() + " : addNewVector() end";
    toc( this->name() +" : Add new vector" );
} // addNewvector


template <typename ParameterSpaceType, typename SpaceType, typename TensorType>
double
DEIMBase<ParameterSpaceType,SpaceType,TensorType>::evaluate( vector_ptrtype V, int const& index, bool seq )
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
} // evaluate for vectors


template <typename ParameterSpaceType, typename SpaceType, typename TensorType>
double
DEIMBase<ParameterSpaceType,SpaceType,TensorType>::evaluate( sparse_matrix_ptrtype M, std::pair<int,int> const& idx, bool seq )
{
    int i = idx.first;
    int j =idx.second;
    double value=0;

    if ( seq )
        return M->operator() ( i,j );
    M->close();

    /*int proc_number = M->mapRow().procOnGlobalCluster(i);

    if ( !M->closed() )
        M->close();
    if ( Environment::worldComm().globalRank()==proc_number )
    {
        // value = M->operator() ( i - M->mapRow().firstDofGlobalCluster(),
        //                         j - M->mapCol().firstDofGlobalCluster());
        auto searchGpDof = M->mapRow().searchGlobalProcessDof( i );
        CHECK( boost::get<0>( searchGpDof ) ) << "Did not find p_dof "<< i <<" when it should be here\n";
        auto ipdof = boost::get<1>( searchGpDof );

        searchGpDof = M->mapCol().searchGlobalProcessDof( j );
        CHECK( boost::get<0>( searchGpDof ) ) << "Did not find p_dof "<< j <<" when it should be here\n";
        auto jpdof = boost::get<1>( searchGpDof );

        value = M->operator()( ipdof, jpdof );
    }

     boost::mpi::broadcast( Environment::worldComm(), value, proc_number );*/
    int ncols=0;
    int ierr=0;
    const PetscScalar *petsc_row;
    const PetscInt    *petsc_cols;

    PetscInt m, n;
    ierr = MatGetOwnershipRange(toPETSc(M)->mat(),&m,&n);CHKERRABORT( M->comm(),ierr );
    int proc_number = M->mapRow().procOnGlobalCluster(i);

    if ( m<=i && i<n )
    {
        ierr = MatGetRow( toPETSc(M)->mat(), i, &ncols, &petsc_cols, &petsc_row );
        CHKERRABORT( M->comm(),ierr );

        for ( int k=0; k<ncols; k++ )
        {
            if ( petsc_cols[k]==j )
                value = static_cast<double>( petsc_row[k] );
        }
        if ( value!=0 && Environment::worldComm().globalRank()!=proc_number )
        {
            Feel::cout << "Entry found on proc "<< Environment::worldComm().globalRank()
                       <<", proc_number is = " << proc_number << ", (i,j)=("<< i<<","<<j<< ")\n";
        }

        ierr  = MatRestoreRow( toPETSc(M)->mat(), i, &ncols, &petsc_cols, &petsc_row );
        CHKERRABORT( M->comm(),ierr );
    }
    boost::mpi::broadcast( Environment::worldComm(), value, proc_number );


    return value;
} // evaluate for matrices


template <typename ParameterSpaceType, typename SpaceType, typename TensorType>
typename DEIMBase<ParameterSpaceType,SpaceType,TensorType>::vectormax_type
DEIMBase<ParameterSpaceType,SpaceType,TensorType>::vectorMaxAbs( vector_ptrtype V )
{
    auto newV = V->clone();
    *newV = *V;
    newV->abs();
    int index = 0;
    double max = newV->maxWithIndex( &index );

    double max_eval=math::abs(evaluate(V,index));
    DCHECK( max_eval==max ) << "evaluation of V(i)=" <<max_eval<<", maximum="<<max <<std::endl;
    return boost::make_tuple( max, index );
} // vectorMaxAbs for vectors


template <typename ParameterSpaceType, typename SpaceType, typename TensorType>
typename DEIMBase<ParameterSpaceType,SpaceType,TensorType>::matrixmax_type
DEIMBase<ParameterSpaceType,SpaceType,TensorType>::vectorMaxAbs( sparse_matrix_ptrtype M )
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
    {
        auto searchGpDof = V->map().searchGlobalProcessDof( i_row );
        CHECK( boost::get<0>( searchGpDof ) ) << "Did not find p_dof "<< i_row <<" when it should be here\n";
        auto gpdof = boost::get<1>( searchGpDof );
        i_col = idx[gpdof];
        //i_col = idx[i_row - V->map().firstDofGlobalCluster()];

    }
    boost::mpi::broadcast( Environment::worldComm(), i_col, proc_number );

    std::pair<int,int> index (i_row,i_col);
    double max_eval=math::abs(evaluate(M,index));
    CHECK( std::abs(max_eval-max)<1e-10 ) << "evaluation of M(i,j)=" <<max_eval<<", maximum="<<max <<std::endl;

    return boost::make_tuple( max, index );
} // vectorMaxAbs for matrices


template <typename ParameterSpaceType, typename SpaceType, typename TensorType>
typename DEIMBase<ParameterSpaceType,SpaceType,TensorType>::bestfit_type
DEIMBase<ParameterSpaceType,SpaceType,TensorType>::computeBestFit()
{
    tic();
    LOG(INFO) << this->name() + " : computeBestFit() start";
    double max=0;
    auto mu_max = M_trainset->max().template get<0>();
    std::vector<parameter_type> to_remove;
    for ( auto const& mu : *M_trainset )
    {
        tensor_ptrtype T = residual( mu );
        if ( M_last_solve_is_ok )
        {
            auto vec_max = vectorMaxAbs(T);
            auto norm = vec_max.template get<0>();

            if ( norm>max )
            {
                max = norm;
                mu_max = mu;
            }
        }
        else
            to_remove.push_back( mu );
    }
    for ( auto const& mu : to_remove )
    {
        Feel::cout << this->name() + " : solver did not converge for mu="<<mu.toString()<<", parameter removed from the sampling\n";
        auto it = std::find( M_trainset->begin(), M_trainset->end(), mu  );
        if ( it!=M_trainset->end() )
            M_trainset->erase( it );
    }
    LOG(INFO) << this->name() + " : computeBestFit() end";
    toc(this->name() + " : compute best fit");

    return boost::make_tuple( mu_max, max );
} //computeBestFit


template <typename ParameterSpaceType, typename SpaceType, typename TensorType>
typename DEIMBase<ParameterSpaceType,SpaceType,TensorType>::tensor_ptrtype
DEIMBase<ParameterSpaceType,SpaceType,TensorType>::residual( parameter_type const& mu, bool force_fem )
{
    LOG(INFO) << this->name() + " : residual() start with "<< mu.toString();
    tensor_ptrtype T;

    if( !M_store_tensors || M_use_ser || !M_solutions[mu.key()] )
    {
        T = this->assemble( mu, false, force_fem );
        T->close();
        if( M_store_tensors && M_last_solve_is_ok && !M_use_ser )
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
} //residual


template<class Archive, class SpaceType, typename RbType >
struct SerializeByBlock
{
    SerializeByBlock( Archive & ar, std::vector<int> const& n,  SpaceType const& Rh, RbType& rb  ) :
        m_ar( ar ),
        m_n( n ),
        m_Rh( Rh ),
        m_rb( rb )
        {}

    template <typename T>
    void operator()( T const& t ) const
    {
        if ( Archive::is_loading::value )
        {
            if ( fusion::at_c<T::value>( m_rb ).size()!=m_n[T::value] )
            {
                auto subRh = m_Rh->template functionSpace<T::value>();
                for ( int i=0; i<m_n[T::value]; i++ )
                    fusion::at_c<T::value>( m_rb ).push_back( subRh->elementPtr() );
                for ( int i=0; i<m_n[T::value]; i++ )
                    m_ar & BOOST_SERIALIZATION_NVP( unwrap_ptr(fusion::at_c<T::value>( m_rb )[i]) );
            }
        }
        else
            for ( int i=0; i<m_n[T::value]; i++ )
                m_ar & BOOST_SERIALIZATION_NVP( unwrap_ptr(fusion::at_c<T::value>( m_rb )[i]) );
    }


private :
    Archive& m_ar;
    std::vector<int> const& m_n;
    SpaceType const& m_Rh;
    RbType& m_rb;

};

template<class Archive, class SpaceType, typename RbType >
void
serializeByBlock( Archive & ar, std::vector<int> const& n,  SpaceType const& Rh, RbType& rb, typename std::enable_if<SpaceType::element_type::is_composite>::type* = nullptr  )
{
    SerializeByBlock<Archive,SpaceType,RbType> s( ar, n, Rh, rb );
    mpl::range_c< int, 0, SpaceType::element_type::nSpaces > range;
    fusion::for_each( range, s );
}
template<class Archive, class SpaceType, typename RbType >
void
serializeByBlock( Archive & ar, std::vector<int> const& n,  SpaceType const& Rh, RbType& rb, typename std::enable_if<!SpaceType::element_type::is_composite>::type* = nullptr  )
{}


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
    __ar & BOOST_SERIALIZATION_NVP( M_n_block_rb );
    DVLOG(2) << "M_n_block_rb saved/loaded\n";

    if ( Archive::is_loading::value )
    {
        if ( M_rb.size()!=M_n_rb  )
        {
            for(int i = 0; i < M_n_rb; i++ )
                M_rb.push_back( Rh->elementPtr() );
            for( int i = 0; i < M_n_rb; ++ i )
                __ar & BOOST_SERIALIZATION_NVP( unwrap_ptr(M_rb[i]) );
            DVLOG(2) << "M_rb saved/loaded\n";
        }
        for ( auto& mu : M_mus )
            mu.setParameterSpace( M_Dmu );
    }
    else
    {
        for( int i = 0; i < M_n_rb; ++ i )
            __ar & BOOST_SERIALIZATION_NVP( unwrap_ptr(M_rb[i]) );
    }

    if ( M_n_block_rb.size() )
        serializeByBlock<Archive,space_ptrtype,block_rb_basis_type>( __ar, M_n_block_rb, Rh, M_block_rb );

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
                                     _worldcomm= Environment::worldCommSeqPtr() );
            Rh = newInterpolationSpace(seqmesh);
        }

        boost::archive::binary_iarchive ia( ifs );
        ia >> *this;
        this->setIsLoaded( true );

        return true;
    }

    return false;
}




} // namespace Feel

#endif

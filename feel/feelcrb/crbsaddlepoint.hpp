/* -*- mode: c++ -*-

 This file is part of the Feel library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 2009-11-24

 Copyright (C) 2009-2012 Université Joseph Fourier (Grenoble I)

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
 * @file   crbsaddlepoint.hpp
 * @date   Mon Nov 24 11:25:55 2014
 *
 * @brief
 *
 *
 */

#ifndef __CRBSADDLEPOINT_H
#define __CRBSADDLEPOINT_H 1

// #include <boost/multi_array.hpp>
// #include <boost/tuple/tuple.hpp>
// #include "boost/tuple/tuple_io.hpp"
// #include <boost/format.hpp>
// #include <boost/foreach.hpp>
// #include <boost/bimap.hpp>
// #include <boost/bimap/support/lambda.hpp>
// #include <boost/archive/text_oarchive.hpp>
// #include <boost/archive/text_iarchive.hpp>
// #include <boost/math/special_functions/fpclassify.hpp>
// #include <fstream>

// #include <boost/serialization/vector.hpp>
// #include <boost/serialization/list.hpp>
// #include <boost/serialization/string.hpp>
// #include <boost/serialization/version.hpp>
// #include <boost/serialization/split_member.hpp>

// #include <vector>

// #include <Eigen/Core>
// #include <Eigen/LU>
// #include <Eigen/Dense>

// #include <feel/feelalg/solvereigen.hpp>

// #include <feel/feelcore/environment.hpp>
// #include <feel/feelcore/parameter.hpp>
// #include <feel/feelcrb/parameterspace.hpp>
// #include <feel/feelcrb/crbdb.hpp>
// #include <feel/feelcrb/crbscm.hpp>
// #include <feel/feelcore/serialization.hpp>
// #include <feel/feelfilters/exporter.hpp>

#include <feel/feel.hpp>
#include <feel/feelcrb/crb.hpp>

namespace Feel
{

/**
 * \class CRBSaddlePoint
 * \brief Certfified Reduced BAsis for Saddle Point Problems
 *
 * \author JB Wahl
 */


template<typename TruthModelType>
class CRBSaddlePoint :
        public CRB<TruthModelType>
{
    typedef CRB<TruthModelType> super_crb;

public:
    //@{ // Truth Model
    typedef TruthModelType model_type;
    typedef boost::shared_ptr<model_type> model_ptrtype;
    typedef typename model_type::mesh_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //@}

    //@{ /// Parameter Space
    typedef typename model_type::parameterspace_type parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef typename parameterspace_type::element_type parameter_type;
    typedef typename parameterspace_type::element_ptrtype parameter_ptrtype;
    typedef typename parameterspace_type::sampling_type sampling_type;
    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;
    //@}

    typedef boost::bimap< int, boost::tuple<double,double,double> > convergence_type;
    typedef double value_type;
    typedef typename convergence_type::value_type convergence;

    //@{ /// Function Space and Elements
    typedef typename model_type::space_type space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename model_type::functionspace_type functionspace_type;
    typedef typename model_type::functionspace_ptrtype functionspace_ptrtype;
    typedef typename model_type::element_type element_type;
    typedef typename model_type::element_ptrtype element_ptrtype;
    //@}

    //@{ Backend and Matrix
    typedef typename model_type::backend_type backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename model_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename model_type::vector_ptrtype vector_ptrtype;
    typedef typename model_type::beta_vector_type beta_vector_type;
    //@}

    //@{ /// Eigen Objects
    typedef Eigen::VectorXd vectorN_type;
    typedef Eigen::MatrixXd matrixN_type;
    typedef Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > map_dense_matrix_type;
    typedef Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic, 1> > map_dense_vector_type;
    typedef boost::tuple< std::vector<vectorN_type>,
                          std::vector<vectorN_type>,
                          std::vector<vectorN_type>,
                          std::vector<vectorN_type> > solutions_tuple;
    typedef boost::tuple< double,double,double,
                          std::vector< std::vector< double > >,
                          std::vector< std::vector< double > > > upper_bounds_tuple;
    typedef boost::tuple< double,double > matrix_info_tuple; //conditioning, determinant
    //@}

    //@{ /// Exporter
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    //@}

    //@{ /// Database
    typedef CRBElementsDB<model_type> crb_elements_db_type;
    typedef boost::shared_ptr<crb_elements_db_type> crb_elements_db_ptrtype;
    //@}


    /// Default constructor
    CRBSaddlePoint()
        :
        super_crb(),
        M_crbdb()
        {}


    /// constructor from command line options
    CRBSaddlePoint( std::string  name,
                    po::variables_map const& vm,
                    model_ptrtype const & model )
        :
        super_crb( name, vm , model ),
        M_crbdb( ( boost::format( "%1%" )
                   % ioption("crb.error-type") ).str(),
                 name,
                 ( boost::format( "%1%-%2%-%3%-saddlepoint" )
                   % name
                   % ioption("crb.output-index")
                   % ioption("crb.error-type") ).str(),
                 vm )
        {
            this->setTruthModel( model );
            if ( M_crbdb.loadDB() )
                LOG(INFO) << "Database " << M_crbdb.lookForDB()
                          << " available and loaded\n";

            // this will be in the offline step
            // (it's only when we enrich or create the database that we want to
            // have access to elements of the RB)
            /*this->M_elements_database.setMN( this->M_N );
            if( this->M_elements_database.loadDB() )
            {
                LOG(INFO) << "database for basis functions "
                          << this->M_elements_database.lookForDB() << " available and loaded\n";
                auto basis_functions = this->M_elements_database.wn();
                this->M_model->rBFunctionSpace()->setBasis( basis_functions );
            }
            else
                LOG( INFO ) <<"no database for basis functions loaded. Start from the begining";
             */
        }


    //! copy constructor
    CRBSaddlePoint( CRBSaddlePoint const & o )
        :
        super_crb( o ),
        M_crbdb( o )
        {}

    //! destructor
    ~CRBSaddlePoint()
        {}



    /// \name Accessors
    //@{
    WorldComm const& worldComm() const { return Environment::worldComm() ; }
    //@}


    /**
     * Returns the lower bound of the output
     *
     * \param mu \f$ \mu\f$ the parameter at which to evaluate the output
     * \param N the size of the reduced basis space to use
     * \param uN primal solution
     *
     *\return compute online the lower bound
     *\and also condition number of matrix A
     */
    boost::tuple<std::vector<double>,matrix_info_tuple>
    lb( size_type N, parameter_type const& mu,
        std::vector< vectorN_type >& uN, std::vector< vectorN_type >& uNdu ,
        std::vector<vectorN_type> & uNold, std::vector<vectorN_type> & uNduold,
        bool print_rb_matrix=false, int K=0 ) const;

    /**
     * Offline computation
     *
     * \return the convergence history (max error)
     */
    convergence_type offline();

    //@{ /// Database
    void saveDB();
    bool loadDB();
    //@}

private :
    void generateSampling();

    CRBDB M_crbdb;

    std::vector< std::vector<matrixN_type> > M_A11qm_pr;
    std::vector< std::vector<matrixN_type> > M_A12qm_pr;

    std::vector < std::vector<vectorN_type> > M_F1qm_pr;
    std::vector < std::vector<vectorN_type> > M_F2qm_pr;

    std::vector < std::vector<vectorN_type> > M_L1qm_pr;
    std::vector < std::vector<vectorN_type> > M_L2qm_pr;


    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const;

    template<class Archive>
    void load( Archive & ar, const unsigned int version ) ;

    BOOST_SERIALIZATION_SPLIT_MEMBER()

}; // class CRBSaddlePoint


template<typename TruthModelType>
typename CRBSaddlePoint<TruthModelType>::convergence_type
CRBSaddlePoint<TruthModelType>::offline()
{
    bool rebuild_database = boption( "crb.rebuild-database" );
    boost::timer ti;
    parameter_type mu ( this->M_Dmu );
    double delta_pr = 0;
    double delta_du = 0;
    size_type index;

    if ( Environment::isMasterRank() )
        std::cout << "Offline CRBSaddlePoint starts, this may take a while until Database is computed.." << std::endl;
    LOG(INFO) << "[CRBSaddlePoint::offline] Starting offline for output " << this->M_output_index << std::endl;
    LOG(INFO) << "[CRBSaddlePoint::offline] initialize underlying finite element model\n";

    if ( rebuild_database || this->M_N==0 )
    {
        ti.restart();
        // Generate sampling depending of the "crb.sampling-mode" option
        generateSampling();

        LOG(INFO) << " -- sampling init done in " << ti.elapsed() << "s";
        ti.restart();

        LOG( INFO )<<"[CRBTrilinear offline] M_error_type = "<<this->M_error_type;

        // empty sets
        this->M_WNmu->clear();
        if( this->M_error_type == CRB_NO_RESIDUAL )
            mu = this->M_Dmu->element();
        else
        {
            // start with M_C = { arg min mu, mu \in Xi }
            boost::tie( mu, index ) = this->M_Xi->min();
        }

        this->M_N=0;
        this->M_maxerror = 1e10;

        LOG(INFO) << "[CRBSaddlePoint::offline] allocate reduced basis data structures\n";



    } // if rebuild_database

    return this->M_rbconv;
} // offline()


template<typename TruthModelType>
typename boost::tuple<std::vector<double>,
                      typename CRBSaddlePoint<TruthModelType>::matrix_info_tuple >
CRBSaddlePoint<TruthModelType>::lb( size_type N, parameter_type const& mu,
                                  std::vector< vectorN_type >& uN,
                                  std::vector< vectorN_type >& uNdu,
                                  std::vector<vectorN_type> & uNold,
                                  std::vector<vectorN_type> & uNduold,
                                  bool print_rb_matrix, int K ) const
{
    std::vector<double>output_vector(1);
    auto matrix_info = boost::make_tuple( 0., 0. );

    return boost::make_tuple( output_vector, matrix_info);

} // lb()


template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::generateSampling()
{

    int proc_number = worldComm().globalRank();
    int total_proc = worldComm().globalSize();
    bool all_proc_same_sampling = boption( "crb.all-procs-have-same-sampling" );
    int sampling_size = ioption("crb.sampling-size");
    std::string sampling_mode = soption("crb.sampling-mode");
    std::string file_name = ( boost::format("M_Xi_%1%_"+sampling_mode+"-proc%2%on%3%")
                              %sampling_size
                              %proc_number
                              %total_proc ).str();

    if ( all_proc_same_sampling )
        file_name += "-all_proc_same_sampling";

    std::ifstream file ( file_name );

    if ( !file )
    {
        std::string supersamplingname =(boost::format("Dmu-%1%-generated-by-master-proc") %sampling_size ).str();

        if( sampling_mode == "log-random" )
            this->M_Xi->randomize( sampling_size , all_proc_same_sampling , supersamplingname );
        else if( sampling_mode == "log-equidistribute" )
            this->M_Xi->logEquidistribute( sampling_size , all_proc_same_sampling , supersamplingname );
        else if( sampling_mode == "equidistribute" )
            this->M_Xi->equidistribute( sampling_size , all_proc_same_sampling , supersamplingname );
        else
            throw std::logic_error( "[CRBTrilinear::offline] ERROR invalid option crb.sampling-mode, please select between log-random, log-equidistribute or equidistribute" );

        this->M_Xi->writeOnFile(file_name);
    }
    else
    {
        this->M_Xi->clear();
        this->M_Xi->readFromFile(file_name);
    }

    this->M_WNmu->setSuperSampling( this->M_Xi );
}


template<typename TruthModelType>
template<class Archive>
void
CRBSaddlePoint<TruthModelType>::save( Archive & ar, const unsigned int version ) const
{
    int proc_number = this->worldComm().globalRank();

    LOG(INFO) <<"[CRBSaddlepoint::save] version : "<<version<<std::endl;

    ar & BOOST_SERIALIZATION_NVP( this->M_output_index );
    ar & BOOST_SERIALIZATION_NVP( this->M_N );
    ar & BOOST_SERIALIZATION_NVP( this->M_rbconv );
    ar & BOOST_SERIALIZATION_NVP( this->M_error_type );
    ar & BOOST_SERIALIZATION_NVP( this->M_Xi );
    ar & BOOST_SERIALIZATION_NVP( this->M_WNmu );
    ar & BOOST_SERIALIZATION_NVP( this->M_Aqm_pr );
    ar & BOOST_SERIALIZATION_NVP( this->M_Fqm_pr );
    ar & BOOST_SERIALIZATION_NVP( this->M_Lqm_pr );

    ar & BOOST_SERIALIZATION_NVP( this->M_A11qm_pr );
    ar & BOOST_SERIALIZATION_NVP( this->M_A12qm_pr );
    ar & BOOST_SERIALIZATION_NVP( this->M_F1qm_pr );
    ar & BOOST_SERIALIZATION_NVP( this->M_F2qm_pr );
    ar & BOOST_SERIALIZATION_NVP( this->M_L1qm_pr );
    ar & BOOST_SERIALIZATION_NVP( this->M_L2qm_pr );

    ar & BOOST_SERIALIZATION_NVP( this->M_current_mu );
    ar & BOOST_SERIALIZATION_NVP( this->M_no_residual_index );

    ar & BOOST_SERIALIZATION_NVP( this->M_maxerror );
} //save( ... )


template<typename TruthModelType>
template<class Archive>
void
CRBSaddlePoint<TruthModelType>::load( Archive & ar, const unsigned int version )
{

    int proc_number = this->worldComm().globalRank();

    LOG(INFO) <<"[CRBSaddlePoint::load] version"<< version <<std::endl;

    ar & BOOST_SERIALIZATION_NVP( this->M_output_index );
    ar & BOOST_SERIALIZATION_NVP( this->M_N );
    ar & BOOST_SERIALIZATION_NVP( this->M_rbconv );
    ar & BOOST_SERIALIZATION_NVP( this->M_error_type );
    ar & BOOST_SERIALIZATION_NVP( this->M_Xi );
    ar & BOOST_SERIALIZATION_NVP( this->M_WNmu );
    ar & BOOST_SERIALIZATION_NVP( this->M_Aqm_pr );
    ar & BOOST_SERIALIZATION_NVP( this->M_Fqm_pr );
    ar & BOOST_SERIALIZATION_NVP( this->M_Lqm_pr );

    ar & BOOST_SERIALIZATION_NVP( this->M_A11qm_pr );
    ar & BOOST_SERIALIZATION_NVP( this->M_A12qm_pr );
    ar & BOOST_SERIALIZATION_NVP( this->M_F1qm_pr );
    ar & BOOST_SERIALIZATION_NVP( this->M_F2qm_pr );
    ar & BOOST_SERIALIZATION_NVP( this->M_L1qm_pr );
    ar & BOOST_SERIALIZATION_NVP( this->M_L2qm_pr );

    ar & BOOST_SERIALIZATION_NVP( this->M_current_mu );
    ar & BOOST_SERIALIZATION_NVP( this->M_no_residual_index );

    ar & BOOST_SERIALIZATION_NVP( this->M_maxerror );
} // load( ... )

template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::saveDB()
{
    super_crb::saveDB();

    fs::ofstream ofs( M_crbdb.dbLocalPath() / M_crbdb.dbFilename() );

    if ( ofs )
    {
        boost::archive::text_oarchive oa( ofs );
        // write class instance to archive
        oa << *this;
        // archive and stream closed when destructors are called
    }
} //saveDB()

template<typename TruthModelType>
bool
CRBSaddlePoint<TruthModelType>::loadDB()
{

    if ( this->rebuildDB() )
        return false;

    super_crb::loadDB();

    fs::path db = M_crbdb.lookForDB();

    if ( db.empty() )
        return false;

    if ( !fs::exists( db ) )
        return false;

    //std::cout << "Loading " << db << "...\n";
    fs::ifstream ifs( db );

    if ( ifs )
    {
        boost::archive::text_iarchive ia( ifs );
        // write class instance to archive
        ia >> *this;
        //std::cout << "Loading " << db << " done...\n";
        // archive and stream closed when destructors are called
        return true;
    }

    return false;
} // loadDB()

} // namespace Feel


namespace boost
{
namespace serialization
{
template< typename T>
struct version< Feel::CRBSaddlePoint<T> >
{
    // at the moment the version of the CRBTrilinear DB is 0. if any changes is done
    // to the format it is mandatory to increase the version number below
    // and use the new version number of identify the new entries in the DB
    typedef mpl::int_<0> type;
    typedef mpl::integral_c_tag tag;
    static const unsigned int value = version::type::value;
};
template<typename T> const unsigned int version<Feel::CRBSaddlePoint<T> >::value;

} // namespace serialization
} // namespace boost



#endif // CRBSADDLEPOINT_H

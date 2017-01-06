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
po::options_description crbSaddlePointOptions( std::string const& prefix="" );

template<typename TruthModelType>
class CRBSaddlePoint :
        public CRB<TruthModelType>
{
    typedef CRB<TruthModelType> super;

public:
    //@{ // Truth Model
    typedef TruthModelType model_type;
    typedef boost::shared_ptr<model_type> truth_model_ptrtype;
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




    CRBSaddlePoint( std::string const& name = "defaultname_crb",
                    WorldComm const& worldComm = Environment::worldComm() ) :
        super( name, worldComm )
        {}

    CRBSaddlePoint( std::string const& name, truth_model_ptrtype const & model ) :
        super( name, model )
        {}


    //@{ /// Database
    void saveDB();
    bool loadDB();
    //@}

private :
    void generateSampling();

    CRBDB M_crbdb;


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
void
CRBSaddlePoint<TruthModelType>::generateSampling()
{


}


template<typename TruthModelType>
template<class Archive>
void
CRBSaddlePoint<TruthModelType>::save( Archive & ar, const unsigned int version ) const
{
    int proc_number = this->worldComm().globalRank();

    LOG(INFO) <<"[CRBSaddlepoint::save] version : "<<version<<std::endl;


} //save( ... )


template<typename TruthModelType>
template<class Archive>
void
CRBSaddlePoint<TruthModelType>::load( Archive & ar, const unsigned int version )
{

    int proc_number = this->worldComm().globalRank();

    LOG(INFO) <<"[CRBSaddlePoint::load] version"<< version <<std::endl;


} // load( ... )

template<typename TruthModelType>
void
CRBSaddlePoint<TruthModelType>::saveDB()
{

} //saveDB()

template<typename TruthModelType>
bool
CRBSaddlePoint<TruthModelType>::loadDB()
{

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

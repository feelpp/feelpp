/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-06-04

  Copyright (C) 2014-2020 Feel++ Consortium

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
#if !defined(FEELPP_HOLO3_HBF_HPP)
#define FEELPP_HOLO3_HBF_HPP

//#include <fftw3.h>
#include <Eigen/Core>

#include <boost/bimap/bimap.hpp>
#include <boost/bimap/list_of.hpp>
#include <boost/bimap/unordered_set_of.hpp>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeltiming/tic.hpp>
#include <feel/feeldiscr/meshstructured.hpp>


/*
namespace boost
{
    template<class Archive, class T>
    inline void serialize(
        Archive & ar, 
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & t, 
        const unsigned int file_version
    ) 
    {
        ar & boost::serialization::make_array(t.data(), t.size()) ;
    }
}
*/
namespace boost {
    namespace serialization {

        template<class Archive, class T>
            void serialize(Archive & ar, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &t,
                    const unsigned int file_version)
            {
                split_free(ar, t, file_version); 
            }
        template<class Archive, class T>
            void save(Archive & ar, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &t,
                    const unsigned int file_version)
            {
                int x = t.rows();
                int y = t.cols();
                ar & x;
                ar & y;
                ar & boost::serialization::make_array(t.data(), t.size());
            }
        template<class Archive, class T>
            void load(Archive & ar, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &t,
                    const unsigned int file_version)
            {
                int x = t.rows();
                int y = t.cols();
                ar & x;
                ar & y;
                t.resize(x,y);
                ar & boost::serialization::make_array(t.data(), t.size());
            }
    } // namespace serialization
} // namespace boost


namespace Feel
{
template <typename T>
using holo3_image = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> ;

/**
 * Given a std::string path to a file, read the header and return the sizes of the HBF and the type of data
 */
std::tuple<int32_t, int32_t, std::string> readHBFHeaderAndSizes( std::string const& s );
/**
 * Given a std::ifstream, read the header and return the sizes of the HBF and the type of data
 */
std::tuple<int32_t, int32_t, std::string> readHBFHeaderAndSizes( std::istream& in );

/**
 * @brief read an HBF file containing an array of floats
 * @param s file name
 * @return eigen array of type T
 */
template <typename T> holo3_image<T> readHBF( std::string const& s );

/**
 * @brief read an HBF stream containing an array of floats
 * @param s hbf stream
 * @return eigen array of type T
 */
template <typename T> holo3_image<T> readHBF( std::istream& s );

/**
 * @brief write a HBF file containing an array of floats
 *
 * @param s file name
 * @param x eigen array of type T
 */
template <typename T> void writeHBF( std::string const& s, holo3_image<T> const& x, worldcomm_ptr_t wc = Environment::worldCommPtr() );

/**
 * @brief write a HBF stream containing an array of floats
 *
 * @param s stream to write the hbf to
 * @param x eigen array of type T
 */
template <typename T> void writeHBF( std::ostream& s, holo3_image<T> const& x );


holo3_image<float> cutHbf ( holo3_image<float> const& im, int n, int start );

namespace bimaps = boost::bimaps;

/**
 * @brief create bijection map between HBF indexing and Feel++ indexing
 */
class Hbf2Feelpp
{
public:
    typedef meta::Pch<Mesh<Hypercube<2>>,1>::ptrtype q1_space_ptrtype;
    typedef meta::Pch<Mesh<Hypercube<2>>,1>::type space_type;
    typedef space_type::element_type q1_element_type;
    
    typedef boost::bimap<bimaps::set_of<std::pair<int,int>>, bimaps::set_of<int> > dof_table;
    typedef dof_table::value_type dof_relation;
    typedef dof_table::left_iterator hbf_dof_iterator;
    typedef dof_table::left_const_iterator hbf_dof_const_iterator;
    typedef dof_table::right_iterator feelpp_dof_iterator;
    typedef dof_table::right_const_iterator feelpp_dof_const_iterator;

    /**
     * build the correspondence data structure
     * \param nx indicates the number of nodes in x direction
     * \param ny indicates the number of nodes in y direction
     * the number of cells is nx-1 and ny-1 in the x and y direction respectively
     */
    Hbf2Feelpp( int nx, int ny, q1_space_ptrtype Yh);

    /**
     * build a Q1 element from a nodal image
     */
    q1_element_type operator()( holo3_image<float> const& x );
    

    /**
     * from a Q1 element build a nodal image
     */
    template<typename ElementType>
    holo3_image<float>
    operator()( ElementType const& u )
        {
            holo3_image<float> x( M_rows, M_cols );
            for( auto const& dof : M_relation.left )
            {
                DCHECK( dof.first.first < M_rows ) << "invalid row index " << dof.first.first;
                DCHECK( dof.first.second < M_cols ) << "invalid col index " << dof.first.second;
                x( dof.first.first, dof.first.second ) = u(dof.second);
            }
            return x;
        }
private:
    int M_rows;
    int M_cols;
    q1_space_ptrtype M_Xh;
    dof_table M_relation;
  
};

class Hbf2FeelppStruc
{
public:

    //typedef meta::Pch<Mesh<Hypercube<2>>,1>::ptrtype q1_space_ptrtype;
    //typedef meta::Pch<Mesh<Hypercube<2>>,1>::type space_type;
    typedef meta::Pch<MeshStructured<Hypercube<2>>,1>::ptrtype q1_space_ptrtype;
    typedef meta::Pch<MeshStructured<Hypercube<2>>,1>::type space_type;
    typedef space_type::element_type q1_element_type;

    typedef boost::bimap<bimaps::unordered_set_of<std::pair<int,int>>, bimaps::unordered_set_of<int>, bimaps::list_of_relation > dof_table;
    typedef dof_table::value_type dof_relation;

    /**
     * build the correspondence data structure
     * \param nx indicates the number of nodes in x direction
     * \param ny indicates the number of nodes in y direction
     * the number of cells is nx-1 and ny-1 in the x and y direction respectively
     */
    Hbf2FeelppStruc( int nx, int ny, q1_space_ptrtype Yh );

    /**
     * build a Q1 element from a nodal image
     */
    q1_element_type operator()( holo3_image<float> const& x );
    q1_element_type operator()( holo3_image<float> const& x, q1_element_type& u);

    /**
     * build a Q0 element from a cell image
     */
    void cellToQ1( holo3_image<float> const& x, q1_element_type& u );
    
    /**
     * from a Q1 element build a nodal image
     */
    holo3_image<float>  operator()( q1_element_type const& u );

private:
    int M_rows;
    int M_cols;
    q1_space_ptrtype M_Xh;
    dof_table M_relation, M_relation_q0;
  
};


/**
 * @brief create bijection map between fine and coarse hbf based indexing grid
 */
class HbfFineToCoarse
{
public:
    typedef meta::Pch<Mesh<Hypercube<2>>,1>::ptrtype q1_space_ptrtype;
    typedef meta::Pch<Mesh<Hypercube<2>>,1>::type space_type;
    typedef space_type::element_type q1_element_type;

    typedef boost::bimap<bimaps::set_of<std::pair<int,int>>, bimaps::set_of<std::pair<int,int>> > dof_table;
    typedef dof_table::value_type dof_relation;
    typedef dof_table::left_iterator hbf_dof_iterator;
    typedef dof_table::left_const_iterator hbf_dof_const_iterator;
    typedef dof_table::right_iterator feelpp_dof_iterator;
    typedef dof_table::right_const_iterator feelpp_dof_const_iterator;

    /**
     * 
     */
    HbfFineToCoarse(int nx, int ny,int m );
                                                                                  
    /**
     * @brief [brief description]
     * @details [long description]
     * 
     * @param x [description]
     * @return [description]
     */
    holo3_image<float>
    operator()( holo3_image<float> const& u, std::string way );
   
    holo3_image<float>
    integ( holo3_image<float> const& u, std::string way,double pixel );
 
    /*
     holo3_image<float>
     changePict( holo3_image<float> const& u, std::string way )
     { 
     int itilde=-1;
     int jtilde=-1;
     if ( way == std::string("F2G") )
     {
     holo3_image<float> x(M_rows/SizeElem,M_cols/SizeElem);  
     for( int i = 0; i < M_rows; ++i )
     {
     if(i%SizeElem==0)
     for( int j = 0; j < M_cols; ++j )
     {
     if(j%SizeElem==0)
     {
     itilde=i/SizeElem;
     jtilde=j/SizeElem;
     x(itilde,jtilde)=u(i,j);
     }
     }
     }   
     return x;
     }
     else printf("No known way");
     };
     */


private:
    int M_rows;
    int M_cols;
    int SizeElem;
    q1_space_ptrtype M_Xh;
    dof_table M_relation;
};

class TransImage
{
    public :
    int T(holo3_image<float> im, std::pair<double,double> c);
    int T(holo3_image<float> im, std::pair<double,double> c, int L);
    
    private :
    double dx =8.9e-3;
    double dy =8.9e-3;
};

 
/**
 * @brief create bijection map between fine and coarse elem based indexing grid
 */
class ElemFineToCoarse
{
public:
    typedef meta::Pch<Mesh<Hypercube<2>>,1>::ptrtype q1_space_ptrtype;
    typedef meta::Pch<Mesh<Hypercube<2>>,1>::type space_type;
    typedef space_type::element_type q1_element_type;

    typedef boost::bimap<bimaps::set_of<std::pair<int,int>>, bimaps::set_of<int> > dof_table;
    typedef dof_table::value_type dof_relation;
    typedef dof_table::left_iterator hbf_dof_iterator;
    typedef dof_table::left_const_iterator hbf_dof_const_iterator;
    typedef dof_table::right_iterator feelpp_dof_iterator;
    typedef dof_table::right_const_iterator feelpp_dof_const_iterator;

    /**
     * 
     */
    ElemFineToCoarse(int nx, int ny,int N );
                                                                                  
    /**
     * @brief [brief description]
     * @details [long description]
     * 
     * @param x [description]
     * @return [description]
     */
    std::pair<double,double>
    operator()( std::pair<int,int> c);

    std::pair<int,int>
    operator()( std::pair<double,double> c, int num);

    std::pair<int,int>
    operator()( int num);
     
private:
    int M_rows;
    int M_cols;
    int SizeElem;
    q1_space_ptrtype M_Xh;
    dof_table M_relation;
};



} // Feel


#endif /* FEELPP_HOLO3_HBF_HPP */

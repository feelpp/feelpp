
#ifndef _EIGEN_MAP_HPP
#define _EIGEN_MAP_HPP 1

#include <Eigen/Core>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace ublas = boost::numeric::ublas;

namespace Feel {

template< int Rows, typename T >
auto
eigenMap( ublas::vector<T> const& v, int nrows = Rows)
{
    typedef Eigen::Matrix<T, Rows, 1, Eigen::ColMajor> EigenVectorType;
    return Eigen::Map< EigenVectorType >(
            const_cast<T*>( v.data().begin() ), nrows, 1
            );
}

template< int Rows, int Cols = Eigen::Dynamic, typename T >
auto
eigenMap( ublas::matrix<T, ublas::column_major> const& m, int nRows = Rows, int cols = Cols )
{
    typedef Eigen::Matrix<T, Rows, Cols, Eigen::ColMajor> EigenMatrixType;
    int nCols = cols;
    if( nCols == Eigen::Dynamic )
        nCols = m.size2();
    return Eigen::Map< EigenMatrixType >( 
            const_cast<T*>( m.data().begin() ), nRows, nCols
            );
}

} // namespace Feel

#endif // _EIGEN_MAP_HPP

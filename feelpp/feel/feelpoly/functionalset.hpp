/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2005-2006 EPFL
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/**
   \file functionalset.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-10-11
 */
#ifndef __FunctionalSet_H
#define __FunctionalSet_H 1

// clang-format off
#include <feel/feelcore/warnoff.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <Eigen/Core>
#include <feel/feelcore/warnon.hpp>
// clang-format on


#include <feel/feelpoly/functional.hpp>

namespace Feel
{
/**
 * \class FunctionalSet
 * \brief Set of functionals
 *
 * \ingroup Polynomial
 * @author Christophe Prud'homme
 * @see
 */
template<typename Space>
class FunctionalSet
{
public:


    /** @name Typedefs
     */
    //@{

    typedef Space space_type;
    typedef typename space_type::value_type value_type;


    typedef FunctionalSet<Space> functionalset_type;
    typedef functionalset_type self_type;
    typedef Functional<Space> functional_type;


    typedef typename space_type::matrix_type matrix_type;
    /** Canonical ordered dual matrix with a static or runtime expansion column axis. */
    using CoefficientStorage = ::Feel::detail::PolynomialCoefficientStorage<value_type,
        ::Feel::detail::polynomialCoefficientExtent<Space, orderSpecification<Space>()>()>;
    using eigen_matrix_type = Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using eigen_matrix_map_type = Eigen::Map<eigen_matrix_type>;
    using eigen_const_matrix_map_type = decltype( std::declval<CoefficientStorage const&>().matrix() );

    typedef std::vector<functional_type> fset_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /** @brief Construct an empty ordered dual set. */
    FunctionalSet()
        :
        M_space(),
        M_fset(),
        M_mat()
    {}

    /** @brief Construct an empty ordered dual set in a primal space. */
    FunctionalSet( space_type const& s )
        :
        M_space( s ),
        M_fset(),
        M_mat()
    {
    }
    /** @brief Assemble an ordered set from owned functionals. */
    FunctionalSet( space_type const& s, std::vector<functional_type> const& fset )
        :
        M_space( s ),
        M_fset( fset ),
        M_mat()
    {
        //std::cout << "FunctionalSet: " << fset[0].coeff() <<  "\n";
        this->setFunctionalSet( fset );
    }
    /** @brief Copy the primal space, mathematical functionals and assembled coefficients. */
    FunctionalSet( FunctionalSet const & fset )
        :
        M_space( fset.M_space ),
        M_fset( fset.M_fset ),
        M_mat( fset.M_mat )
    {}

    /** @brief Release the owned ordered dual set. */
    ~FunctionalSet()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    /** @brief Copy all ordered dual value state. */
    self_type& operator=( self_type const& fset )
    {
        if ( this != &fset )
        {
            M_space = fset.M_space;
            M_fset = fset.M_fset;
            M_mat = fset.M_mat;
        }

        return *this;
    }

    /**
     * \return the i-th functional
     */
    functional_type const& operator()( uint16_type i ) const
    {
        return M_fset[i];
    }

    /**
     * \return the value of the functional set applied to a polynomial
     */
    matrix_type operator()( space_type const& p ) const
    {
        auto dualView = M_mat.matrix();
        auto primalView = p.coefficients();
        decltype(auto) dual = space_type::polyset_type::toMatrix( dualView );
        decltype(auto) primal = space_type::polyset_type::toMatrix( primalView );
        if ( dual.cols() != primal.cols() )
            throw std::invalid_argument( "dual and primal expansion dimensions are incompatible" );
        return ::Feel::detail::materializePolynomialMatrix<matrix_type>( dual * primal.transpose() );
    }
    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the function space from which the functionals take
     * their values
     */
    space_type const& functionSpace() const
    {
        return M_space;
    }

    /**
     * This works only if the basis is orthonormal
     *
     * \return the representation of the functional set using basis
     * of the function space.
     */
    matrix_type rep() const
    {
        return M_mat.template legacy<matrix_type>();
    }

    /** \return basis-dependent matrix representation of the ordered dual set. */
    matrix_type dualMatrix() const
    {
        return rep();
    }

    /** @return A zero-copy view of the canonical ordered-dual coefficients. */
    [[nodiscard]] auto eigenDualMatrix() const&
    {
        return M_mat.matrix();
    }

    /** @brief Reject an ordered-dual view of a temporary owner. */
    auto eigenDualMatrix() const&& = delete;

    /** \return number of mathematical functionals in the ordered set. */
    [[nodiscard]] std::size_t size() const noexcept
    {
        return M_fset.size();
    }


    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set the function space
     */
    void setFunctionSpace( space_type const& __space )
    {
        M_space = __space;
    }

    /**
     * set the Functional set
     */
    void setFunctionalSet( std::vector<functional_type> const& fset )
    {
        constexpr Eigen::Index rowsPerFunctional = space_type::is_scalar ? 1 : space_type::nComponents;
        Eigen::Index const columns = fset.empty() ? M_mat.cols() : fset.front().rieszRepresentation().cols();
        auto const rows = ::Feel::detail::matrixElementCount( rowsPerFunctional, fset.size(), sizeof( value_type ) );
        eigen_matrix_type assembled( rows, columns );
        for ( std::size_t i = 0; i < fset.size(); ++i )
        {
            auto const functional = fset[i].rieszRepresentation();
            if ( functional.rows() != rowsPerFunctional || functional.cols() != columns )
                throw std::invalid_argument( "inconsistent functional shape in ordered dual set" );
            assembled.middleRows( i * rowsPerFunctional, rowsPerFunctional ) = functional;
        }
        M_mat.assign( assembled );
        M_fset = fset;
    }

    //@}

    /** @name  Methods
     */
    //@{


    //@}



protected:

private:
    /** Owned primal space from which the functionals take values. */
    space_type M_space;
    /** Mathematical functionals in their prescribed order. */
    fset_type M_fset;
    /** Sole assembled ordered-dual representation. */
    CoefficientStorage M_mat;
};
} // Feel
#endif /* __FunctionalSet_H */

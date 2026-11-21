/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2005-2006 EPFL
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/**
   \file functional.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-10-06
 */
#ifndef __Functional_H
#define __Functional_H 1

#include <boost/operators.hpp>
#include <Eigen/Core>
#include <feel/feelpoly/polynomialcoefficients.hpp>

// clang-format off
#include <feel/feelcore/warnoff.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <feel/feelcore/warnon.hpp>
// clang-format on

#include <feel/feelcore/feel.hpp>


namespace Feel
{
namespace ublas = boost::numeric::ublas;

/**
 * \class Functional
 * \brief represents a linear functional
 *
 * A functional is defined by a polynomial set and a set of
 * coefficients
 *
 * \ingroup Polynomial
 * @author Christophe Prud'homme
 * @see
*/
template<typename Space>
class Functional
    :
public boost::addable<Functional<Space> >
{
    typedef boost::addable<Functional<Space> > super;
public:


    /** @name Typedefs
     */
    //@{

    typedef typename Space::value_type value_type;
    typedef Functional<Space> self_type;
    typedef Space space_type;
    typedef Space polynomialset_type;
    typedef typename space_type::polynomial_type polynomial_type;
    typedef typename space_type::basis_type basis_type;
    typedef typename space_type::matrix_type matrix_type;

    static inline const uint16_type nComponents = space_type::nComponents;

    // representation type for the functionals
    typedef ublas::matrix<value_type> rep_type;
    /** Canonical Eigen owner retains the static expansion column axis. */
    using CoefficientStorage = ::Feel::detail::PolynomialCoefficientStorage<value_type,
        ::Feel::detail::polynomialCoefficientExtent<Space, orderSpecification<Space>()>()>;
    using eigen_representation_type = Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using eigen_representation_map_type = decltype( std::declval<CoefficientStorage const&>().matrix() );

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /** @brief Construct an empty functional in the default primal space. */
    Functional()
        :
        super(),
        M_p(),
        M_coeff()
    {}

    /** @brief Initialize the functional from the primal coefficient representation. */
    Functional( space_type const& P )
        :
        super(),
        M_p( P ),
        M_coeff( ::Feel::detail::polynomialCoefficientStorageFrom<CoefficientStorage>( M_p ) )

    {

    }

    /** @brief Own an explicit legacy Riesz representation in a primal space. */
    Functional( space_type const& P,
                matrix_type const& coeff )
        :
        super(),
        M_p( P ),
        M_coeff( coeff )
    {
    }

    /** @brief Copy primal-space and coefficient ownership. */
    Functional( Functional const & __f )
        :
        M_p( __f.M_p ),
        M_coeff( __f.M_coeff )
    {
    }

    /** @brief Release the owned primal space and coefficients. */
    ~Functional() = default;

    //@}

    /** @name Operator overloads
     */
    //@{

    /** @brief Copy the functional as a value, including its primal space. */
    self_type& operator=( self_type const& __f )
    {
        if ( this != &__f )
        {
            M_p = __f.M_p;
            M_coeff = __f.M_coeff;
        }

        return *this;
    }

    /**
     * add to another functional
     * it generates automatically operator+ thanks to addable
     */
    self_type& operator+=( const self_type& __f )
    {
        if ( M_coeff.rows() != __f.M_coeff.rows() || M_coeff.cols() != __f.M_coeff.cols() )
            throw std::invalid_argument( "functional addition requires matching coefficient shapes" );
        M_coeff.assign( M_coeff.matrix() + __f.M_coeff.matrix() );
        return *this;
    }

    /**
     * apply the functional to a polynomial
     *
     *
     * \param p polynomial
     * \return matrix resulting from the application of the functional to the polynomial
     */
    matrix_type operator()( polynomial_type const& p ) const
    {
        if ( p.coefficients().cols() != M_coeff.cols() )
            throw std::invalid_argument( "functional and polynomial expansion dimensions are incompatible" );
        return ::Feel::detail::materializePolynomialMatrix<matrix_type>( p.coefficients() * M_coeff.matrix().transpose() );
    }

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the dimension of the polynomial space
     */
    uint16_type size() const
    {
        return ::Feel::detail::checkedLocalIndexValue( M_coeff.cols() );
    }

    /**
     * \return the coefficient of the functional in the basis
     * associated with the polynomial space
     */
    rep_type coeff() const
    {
        return M_coeff.template legacy<rep_type>();
    }

    /**
     * \return zero-copy Eigen view of this functional's Riesz representation
     * in the current primal basis.
     *
     * The mathematical object remains `Functional`; this view is only a dense
     * representation used by assembly/factorization kernels.
     */
    [[nodiscard]] auto rieszRepresentation() const&
    {
        return M_coeff.matrix();
    }

    /** @brief Reject a borrowed Riesz view of a temporary functional. */
    auto rieszRepresentation() const&& = delete;

    //@}

    /** @name  Mutators
     */
    //@{

    /** @brief Replace the Riesz representation with owned legacy coefficients. */
    void setCoefficient( matrix_type const& __coeff )
    {
        M_coeff.assign( __coeff );
    }

    /** @brief Own native coefficients without retaining an Eigen expression. */
    template<typename Derived>
    void setCoefficient( Eigen::MatrixBase<Derived> const& coefficients )
    { M_coeff.assign( coefficients ); }


    //@}

    /** @name  Methods
     */
    //@{


    //@}



protected:

private:
    /** Owned primal polynomial space defining the Riesz representation. */
    space_type M_p;
    /** Sole Riesz coefficient owner. */
    CoefficientStorage M_coeff;
};

}
#endif /* __Functional_H */

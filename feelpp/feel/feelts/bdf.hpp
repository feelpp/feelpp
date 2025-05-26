/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2006-12-30

   Copyright (C) 2006-2008 Université Joseph Fourier (Grenoble)
   Copyright (C) 2011-2016 Feel++ Consortium

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 3.0 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file bdf.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-12-30
*/
#ifndef FEELPP_TS_BDF_H
#define FEELPP_TS_BDF_H

#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>



#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <feel/feelcore/parameter.hpp>

#include <feel/feelalg/glas.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelts/tsbase.hpp>

namespace Feel
{
namespace ublas = boost::numeric::ublas;

enum BDFTimeScheme { BDF_ORDER_ONE=1, BDF_ORDER_TWO, BDF_ORDER_THREE, BDF_ORDER_FOUR, BDF_MAX_ORDER = 4 };

/**
 * \class Bdf
 * \ingroup SpaceTime
 * \brief Backward Differentiation Formula (BDF) time discretization
 *
 * This class implements the Backward Differentiation Formula (BDF) method for solving time-dependent problems using a multistep, implicit time-stepping scheme.
 *
 * A differential equation of the form
 *
 * \f$ M \frac{du}{dt} = A u + f \f$
 *
 * is discretized in time using a polynomial \f$ p(t) \f$ of order \f$ n \f$ that interpolates the past state vectors \f$ (t_i, u_i) \f$ for \f$ i = k-n+1,\dots,k+1 \f$.
 *
 * The first derivative is approximated as:
 *
 * \f$ p'(t_{k+1}) = \frac{1}{\Delta t} \left( \alpha_0 u_{k+1} - \sum_{i=1}^n \alpha_i u_{k+1-i} \right) \f$
 *
 * Leading to the discrete equation:
 *
 * \f$ \frac{\alpha_0}{\Delta t} M u_{k+1} = A u_{k+1} + f_{k+1} + M \bar{p} \f$
 *
 * where:
 *
 * \f$ \bar{p} = \frac{1}{\Delta t} \sum_{i=1}^n \alpha_i u_{k+1-i} \f$
 *
 * The class provides access to:
 * - The BDF coefficients \f$ \alpha_i \f$
 * - The extrapolation coefficients \f$ \beta_i \f$ used to predict \f$ u_{k+1} \f$
 *
 * \f$ u_{k+1} \approx \sum_{i=0}^{n-1} \beta_i u_{k-i} \f$
 *
 * Additionally, the class supports the approximation of the second-order time derivative using:
 *
 * \f$ \frac{d^2 u}{dt^2}(t_{k+1}) \approx \frac{1}{\Delta t^2} \sum_{i=0}^{n+1} \alpha_i^{(2)} u_{k+1-i} \f$
 *
 * This is particularly useful for solving hyperbolic problems such as wave equations.
 *
 * The class manages the history of time states and time stamps, and supports both extrapolation and variational RHS assembly through `poly()`, `polyDeriv()`, and `polySecondDeriv()`.
 */
template<typename SpaceType>
class Bdf : public TSBase
{
    friend class boost::serialization::access;
    typedef TSBase super;
public:
    typedef Bdf<SpaceType> bdf_type;
    typedef std::shared_ptr<bdf_type> bdf_ptrtype;
    typedef SpaceType space_type;
    typedef std::shared_ptr<space_type>  space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename space_type::element_ptrtype element_ptrtype;
    typedef typename space_type::return_type return_type;
    typedef typename element_type::value_type value_type;
    typedef std::vector< element_ptrtype > unknowns_type;
    typedef typename node<value_type>::type node_type;

    typedef typename super::time_iterator time_iterator;

    /**
     * Constructor
     *
     * @param space approximation space
     * @param name name of the BDF
     */
    Bdf( space_ptrtype const& space, std::string const& name, std::string const& prefix="", po::variables_map const& vm =  Environment::vm() );

    Bdf( space_ptrtype const& __space, std::string const& name, std::string const& prefix, po::variables_map const& vm, int order,
         double ti, double tf, double dt, bool steady, bool reverse, bool restart, std::string const& restart_path, bool restart_at_last_save,
         bool save, int freq, bool rank_proc_in_files_name, std::string const& format, int n_consecutive_save );

    //! copy operator
    Bdf( Bdf const& b );

    ~Bdf() override;

    //! @return the FunctionSpace associated to BDF
    space_ptrtype const& functionSpace() const { return M_space; }

    //! apply a remesh by given the new functionspace and the interpolation matrix between previous and new function space
    void applyRemesh( space_ptrtype const& space, std::shared_ptr<MatrixSparse<double>> const& matInterp )
        {
            //! change space/fields and interpolate
            M_space = space;

            // not yet init : do nothing
            if ( M_history.empty() )
                return;

            // interpolate unknown fields
            unknowns_type oldUnknowns = M_history;
            for ( int k=0;k< M_history.size();++k )
            {
                M_history[k] = M_space->elementPtr();
                matInterp->multVector( unwrap_ptr(oldUnknowns[k]), unwrap_ptr(M_history[k]) );
            }

            //! recompute poly and polyDeriv
            M_poly.reset();
            M_polyDeriv.reset();
            this->computePolyAndPolyDeriv();

            // TODO : change repository where fields are saved
        }

    //! return a deep copy of the bdf object
    bdf_ptrtype deepCopy() const
    {
        auto b = bdf_ptrtype( new bdf_type( *this ) );

        for ( auto it = b->M_history.begin(), en = b->M_history.end(); it != en; ++ it )
        {
            *it = element_ptrtype( new element_type( M_space ) );
        }

        return b;
    }


    //! return the curent order used at current time time
    int timeOrder() const
    {
        return M_order_cur;
    }

    void setTimeOrder( int order )
    {
        M_order_cur = order;
        this->computePolyAndPolyDeriv();
    }

    //! return the order in time
    int bdfOrder() const
    {
        return M_order;
    }

    void setOrder( int order )
    {
        M_order = order;
        M_order_cur = M_order; // require when restart
    }

    //!return the prefix
    std::string bdfPrefix() const
    {
        return this->M_prefix;
    }

    //! return the number of iterations between order change
    int numberOfIterationsBetweenOrderChange() const
    {
        return M_iterations_between_order_change;
    }

    //! return the number of iterations since last order change
    int numberOfIterationsSinceOrderChange() const
    {
        return this->iteration()-M_last_iteration_since_order_change;
    }

    //! return number of consecutive save
    int numberOfConsecutiveSave() const { return M_numberOfConsecutiveSave; }

    //! set number of consecutive save
    void setNumberOfConsecutiveSave( int n )
    {
        if ( n > 0 )
            M_numberOfConsecutiveSave = n;
    }

    //! return a vector of the times prior to timeInitial() (included)
    std::map<int,double> priorTimes() const override
    {
        std::map<int,double> prior;
        for( int i = 0; i < this->M_order+2; ++i )
            prior[i]=timeInitial()-i*timeStep();
        return prior;
    }

    /**
       Initialize all the entries of the unknown vector to be derived with the
       vector u0 (duplicated)
    */
    void initialize( element_type const& u0 );

    /**
     * Initialize all the entries of the unknown vector to be derived with the
     * vector u0 (duplicated) and set the time step to dt
     */
    void initialize( std::vector<element_type> const& u0 );

    /**
       Initialize all the entries of the unknown vector to be derived with a
       set of vectors uv0
    */
    void initialize( unknowns_type const& uv0 );

    /**
       start the bdf
    */
    double start();

    /**
       start the bdf with a given initial value
       @param u0 initial value of the state vector
    */
    double start( element_type const& u0 );

    /**
       start the bdf with a given initial value
       @param u0 history of initial values of the state vector
    */
    double start( std::vector<element_type> const& u0 );

    /**
       start the bdf with a given initial value
       @param uv0 history of initial values of the state vector
    */
    double start( unknowns_type const& uv0 );

    /**
       restart the bdf
    */
    double restart();

    /**
       Update the vectors of the previous time steps by shifting on the right
       the old values.
       @param u_curr current (new) value of the state vector
    */
    template<typename container_type>
    void shiftRight( typename space_type::template Element<value_type, container_type> const& u_curr );

    /**
     * Move to next time step (solution not shifted).
     */
    double next() const override
    {
        double tcur = super::next();

        if ( M_strategyHighOrderStart == 1 )
        {
            if ( ( ( M_iteration - M_last_iteration_since_order_change ) == M_iterations_between_order_change ) &&
                 M_order_cur < M_order )
            {
                M_last_iteration_since_order_change = M_iteration;
                ++M_order_cur;
            }
        }

        return tcur;
    }

    /**
     * Move to next time step with solution shift.
     */
    template<typename container_type>
    double
    next( typename space_type::template Element<value_type, container_type> const& u_curr )
    {
        this->shiftRight( u_curr );

        double tcur = this->next();

        // do here because M_order_cur can change in the call of next()
        this->computePolyAndPolyDeriv();

        return tcur;
    }

    //! Return the first derivative at time t_{k}
    element_type const& firstDerivative() const
    {
        if (!M_firstDeriv)
            M_firstDeriv = M_space->elementPtr();
        M_firstDeriv->zero();
        for (int i = 0; i <= this->timeOrder(); ++i)
            M_firstDeriv->add(this->polyDerivCoefficient(i), *M_history[i]);
        return *M_firstDeriv;
    }

    //! Return the first derivative at time t_{k+1}
    element_type const& firstDerivative( element_type const& u ) const
    {
        if (!M_firstDeriv)
            M_firstDeriv = M_space->elementPtr();
        M_firstDeriv->zero();
        M_firstDeriv->add(this->polyDerivCoefficient(0), u);
        for (int i = 0; i < this->timeOrder(); ++i)
            M_firstDeriv->add(this->polyDerivCoefficient(i+1), *M_history[i]);
        return *M_firstDeriv;
    }


    //! Return the second derivative at time t_{k+1}
    element_type const& secondDerivative() const
    {
        if (!M_secondDeriv)
            M_secondDeriv = M_space->elementPtr();
        M_secondDeriv->zero();
        for (int i = 0; i < M_alpha2[this->timeOrder() - 1].size(); ++i)
            M_secondDeriv->add(this->polySecondDerivCoefficient(i), *M_history[i]);
        return *M_secondDeriv;
    }

    element_type const& extrapolation() const { return this->poly(); }

    /**
     * Return \f$ \alpha_i \f$
     */
    double polyCoefficient( int i ) const
    {
        CHECK( i >=0 && i < BDF_MAX_ORDER ) <<  "[BDF] invalid index " << i;
        return M_beta[this->timeOrder()-1][i];
    }

    /**
     * Return \f$ \frac{\alpha_i}{\Delta t} \f$
     */
    double polyDerivCoefficient( int i ) const
    {
        CHECK( i >=0 && i <= BDF_MAX_ORDER ) << "[BDF] invalid index " << i;
        return M_alpha[this->timeOrder()-1][i]/math::abs( this->timeStep() );
        //return M_alpha[this->timeOrder()-1][i]/this->timeStep();
    }

    /**
     * Return \f$ \frac{\alpha_i^{(2)}}{\Delta t^2} \f$
     */
    double polySecondDerivCoefficient(int i) const
    {
        int order = this->timeOrder();
        CHECK(order >= 2 && order <= BDF_MAX_ORDER);
        CHECK(i >= 0 && i < M_alpha2[order - 1].size());

        return M_alpha2[order - 1][i] / (this->timeStep() * this->timeStep());
    }

    //! Returns the right hand side \f$ \bar{p} \f$ of the time derivative formula
    element_type const& polyDeriv() const;

    //! Returns the right hand side \f$ \bar{p} \f$ of the time derivative formula
    element_ptrtype const& polyDerivPtr() const { return M_polyDeriv; }

    //! Returns the right hand side \f$ \bar{p}^{(2)} \f$ of the second order time derivative formula
    element_type const& polySecondDeriv() const
    {
        return *M_polySecondDeriv;
    }

    //! Returns the right hand side \f$ \bar{p}^{(2)} \f$ of the second order time derivative formula
    element_ptrtype const& polySecondDerivPtr() const { return M_polySecondDeriv; }

    //! Compute the polynomial extrapolation approximation of order n-1 of
    //! u^{n+1} defined by the n stored state vectors
    element_type const& poly() const;

    //! Compute the polynomial extrapolation approximation of order n-1 of
    //! u^{n+1} defined by the n stored state vectors
    element_ptrtype const& polyPtr() const { return M_poly; }

    //! Return a vector with the last n state vectors
    unknowns_type const& history() const { return M_history; }
    FEELPP_DEPRECATED unknowns_type const& unknowns() const { return M_history; }

    //! Return a vector with the last n state vectors
    unknowns_type& history() { return M_history; }
    FEELPP_DEPRECATED unknowns_type& unknowns() { return M_history; }

    //! Return the previous element at previous time i-1
    element_type const& history( int i ) const
    {
        CHECK( i >= 0 && i < M_history.size() ) << "[BDF] invalid index " << i;
        return *M_history[i];
    }
    FEELPP_DEPRECATED element_type& unknown( int i );

    //! Return the previous element at previous time i-1
    element_ptrtype historyPtr( int i ) const
    {
        CHECK( i >= 0 && i < M_history.size() ) << "[BDF] invalid index " << i;
        return M_history[i];
    }
    FEELPP_DEPRECATED element_ptrtype unknownPtr( int i );


    //! update field \u with derivative at previous time indexed by \i (i.e. curent_time - i - 1)
    void updateDerivative( element_type & u, int i = 0 ) const;

    element_type const& prior() const { return *M_history[0]; }

    element_type& prior() { return *M_history[0]; }

    template<typename container_type>
    FEELPP_DEPRECATED void setUnknown( int i,  typename space_type::template Element<value_type, container_type> const& e )
    {
        *M_history[i] = e;
    }
    template<typename container_type>
    void setHistory( int i,  typename space_type::template Element<value_type, container_type> const& e )
    {
        *M_history[i] = e;
    }

    /**
     * \brief Set history and update polynomial and derivatives
     *
     * Accepts n+1 unknowns (from u_{k+1}, u_k, ..., u_{k+1-n}) and fills M_history.
     * Automatically recomputes poly(), polyDeriv(), and polySecondDeriv().
     *
     * Usage:
     *   bdf->setHistory(u0, u1, u2);
     */
    template<typename... Elements>
    void setHistory(Elements const&... elems)
    {
        static_assert(sizeof...(elems) <= BDF_MAX_ORDER + 2, "Too many unknowns for BDF");

        std::array<element_type const*, sizeof...(elems)> args = { &elems... };
        M_history.resize( args.size() );
        for ( uint8_type __i = 0; __i < M_history.size(); ++__i )
        {
            M_history[__i] = M_space->elementPtr();
            *M_history[__i] = *args[__i];
            M_history[__i]->printMatlab(fmt::format("u{}", __i));
        }

        this->computePolyAndPolyDeriv();
    }
    void setHistory(std::vector<element_type> const& vec)
    {
        M_history.resize(vec.size());
        for (std::size_t i = 0; i < vec.size(); ++i)
        {
            M_history[i] = M_space->elementPtr();
            *M_history[i] = vec[i];
        }
    }
    void showMe( std::ostream& __out = std::cout ) const;

    //! Load current unknown in a file (hdf5, binary, ...)
    void loadCurrent();

    void print() const override
    {
        LOG(INFO) << "============================================================\n";
        LOG(INFO) << "BDF Information\n";
        LOG(INFO) << "   time step : " << this->timeStep() << "\n";
        LOG(INFO) << "time initial : " << this->timeInitial() << "\n";
        LOG(INFO) << "  time final : " << this->timeFinal() << "\n";
        LOG(INFO) << "  time order : " << this->timeOrder() << "\n";
    }


private:
    void init();

    //! Save current unknown in a file (hdf5, binary, ...)
    void saveCurrent();

    //! save/load Bdf metadata
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
        DVLOG(2) << "[BDF::serialize] saving/loading archive\n";
        ar & boost::serialization::base_object<TSBase>( *this );
    }

    //! compute BDF coefficients
    void computeCoefficients();

    //! compute extrapolation field and rhs part of bdf scheme
    void computePolyAndPolyDeriv();

private:

    //! bdf order
    int M_order;

    /**
     * strategy to start with high order scheme :
     * 0 start with order given
     * 1 start with order 1 and increase step by step up to order given
     */
    int M_strategyHighOrderStart;

    //! bdf order used at the current time
    mutable int M_order_cur;

    mutable int M_last_iteration_since_order_change;
    int M_iterations_between_order_change;

    //! space
    space_ptrtype M_space;

    //! Last n state vectors
    unknowns_type M_history;

    //! Coefficients \f$ \alpha_i \f$ of the time bdf discretization
    std::vector<ublas::vector<double>> M_alpha;
    //! Coefficients \f$ \alpha_i^{(2)} \f$ for second time derivative
    std::vector<ublas::vector<double>> M_alpha2;

    //! Coefficients \f$ \beta_i \f$ of the extrapolation
    std::vector<ublas::vector<double> > M_beta;

    //! extrapolation field and rhs part of bdf scheme
    mutable element_ptrtype M_poly, M_polyDeriv, M_firstDeriv;
    //! Storage of the second derivative vector
    mutable element_ptrtype M_polySecondDeriv, M_secondDeriv;

    int M_numberOfConsecutiveSave;
};

template <typename SpaceType>
Bdf<SpaceType>::Bdf( space_ptrtype const& __space,
                     std::string const& name,
                     std::string const& prefix,
                     po::variables_map const& vm )
    :
    super( name, prefix, __space->worldComm(), vm ),
    M_order( ioption(_prefix=prefix,_name="bdf.order",_vm=vm) ),
    M_strategyHighOrderStart( ioption(_prefix=prefix,_name="bdf.strategy-high-order-start",_vm=vm) ),
    M_order_cur( M_order ),
    M_iterations_between_order_change( ioption(_prefix=prefix,_name="bdf.iterations-between-order-change",_vm=vm) ),
    M_space( __space ),
    M_alpha( BDF_MAX_ORDER ),
    M_beta( BDF_MAX_ORDER ),
    M_numberOfConsecutiveSave( M_order+2 )
{
    computeCoefficients();

    CHECK( this->numberOfConsecutiveSave() >= this->bdfOrder() ) << "numberOfConsecutiveSave is too small, should be >= bdfOrder";
    M_history.resize( std::max(this->bdfOrder()+2, this->numberOfConsecutiveSave()) );
    for ( uint8_type __i = 0; __i < M_history.size(); ++__i )
    {
        M_history[__i] = element_ptrtype( new element_type( M_space ) );
        M_history[__i]->zero();
    }

    this->computePolyAndPolyDeriv();
}

template <typename SpaceType>
Bdf<SpaceType>::Bdf( space_ptrtype const& __space, std::string const& name, std::string const& prefix, po::variables_map const& vm, int order,
                     double ti, double tf, double dt, bool steady, bool reverse, bool restart, std::string const& restart_path, bool restart_at_last_save,
                     bool save, int freq, bool rank_proc_in_files_name, std::string const& format, int n_consecutive_save )
    :
    super( name, prefix, __space->worldComm(), vm, ti, tf, dt, steady, reverse, restart, restart_path, restart_at_last_save,
           save, freq, rank_proc_in_files_name, format ),
    M_order( order ),
    M_strategyHighOrderStart( ioption(_prefix=prefix,_name="bdf.strategy-high-order-start",_vm=vm) ),
    M_order_cur( M_order ),
    M_iterations_between_order_change( ioption(_prefix=prefix,_name="bdf.iterations-between-order-change",_vm=vm) ),
    M_space( __space ),
    M_alpha( BDF_MAX_ORDER ),
    M_beta( BDF_MAX_ORDER ),
    M_numberOfConsecutiveSave( n_consecutive_save )
{
    computeCoefficients();

    CHECK( this->numberOfConsecutiveSave() >= this->bdfOrder() ) << "numberOfConsecutiveSave is too small, should be >= bdfOrder";
    M_history.resize( std::max(this->bdfOrder()+2, this->numberOfConsecutiveSave()) );
    for ( uint8_type __i = 0; __i < M_history.size(); ++__i )
    {
        M_history[__i] = element_ptrtype( new element_type( M_space ) );
        M_history[__i]->zero();
    }

    this->computePolyAndPolyDeriv();
}


//! copy operator
template <typename SpaceType>
Bdf<SpaceType>::Bdf( Bdf const& b )
    :
        super( b ),
        M_order( b.M_order ),
        M_strategyHighOrderStart( b.M_strategyHighOrderStart ),
        M_order_cur( b.M_order_cur ),
        M_last_iteration_since_order_change( b.M_last_iteration_since_order_change ),
        M_iterations_between_order_change( b.M_iterations_between_order_change ),
        M_space( b.M_space ),
        M_history( b.M_history ),
        M_alpha( b.M_alpha ),
        M_beta( b.M_beta ),
        M_numberOfConsecutiveSave( b.M_numberOfConsecutiveSave )
{}

template <typename SpaceType>
void
Bdf<SpaceType>::computeCoefficients()
{
    M_alpha.resize(BDF_MAX_ORDER);
    M_beta.resize(BDF_MAX_ORDER);
    M_alpha2.resize(BDF_MAX_ORDER);

    for (int i = 0; i < BDF_MAX_ORDER; ++i)
    {
        M_alpha[i].clear();
        M_beta[i].clear();
        M_alpha2[i].clear();
    }

    // BDF1
    M_alpha[0].resize(2);
    M_alpha[0][0] = 1.0;
    M_alpha[0][1] = -1.0;

    M_beta[0].resize(1);
    M_beta[0][0] = 1.0;

    // BDF2
    M_alpha[1].resize(3);
    M_alpha[1][0] = 3.0 / 2.0;
    M_alpha[1][1] = -2.0;
    M_alpha[1][2] = 0.5;

    M_beta[1].resize(2);
    M_beta[1][0] = 2.0;
    M_beta[1][1] = -1.0;

    M_alpha2[1].resize(4);
    M_alpha2[1][0] = 2.0;
    M_alpha2[1][1] = -5.0;
    M_alpha2[1][2] = 4.0;
    M_alpha2[1][3] = -1.0;

    // BDF3
    M_alpha[2].resize(4);
    M_alpha[2][0] = 11.0 / 6.0;
    M_alpha[2][1] = -3.0;
    M_alpha[2][2] = 1.5;
    M_alpha[2][3] = -1.0 / 3.0;

    M_beta[2].resize(3);
    M_beta[2][0] = 3.0;
    M_beta[2][1] = -3.0;
    M_beta[2][2] = 1.0;

    M_alpha2[2].resize(5);
    M_alpha2[2][0] = 35.0 / 12.0;
    M_alpha2[2][1] = -104.0 / 12.0;
    M_alpha2[2][2] = 114.0 / 12.0;
    M_alpha2[2][3] = -56.0 / 12.0;
    M_alpha2[2][4] = 11.0 / 12.0;

    // BDF4
    M_alpha[3].resize(5);
    M_alpha[3][0] = 25.0 / 12.0;
    M_alpha[3][1] = -4.0;
    M_alpha[3][2] = 3.0;
    M_alpha[3][3] = -4.0 / 3.0;
    M_alpha[3][4] = 1.0 / 4.0;

    M_beta[3].resize(4);
    M_beta[3][0] = 4.0;
    M_beta[3][1] = -6.0;
    M_beta[3][2] = 4.0;
    M_beta[3][3] = -1.0;

    M_alpha2[3].resize(6);
    M_alpha2[3][0] = 45.0 / 12.0;
    M_alpha2[3][1] = -154.0 / 12.0;
    M_alpha2[3][2] = 214.0 / 12.0;
    M_alpha2[3][3] = -156.0 / 12.0;
    M_alpha2[3][4] = 61.0 / 12.0;
    M_alpha2[3][5] = -10.0 / 12.0;
}

template <typename SpaceType>
void
Bdf<SpaceType>::init()
{

    CHECK( this->numberOfConsecutiveSave() >= this->bdfOrder() ) << "numberOfConsecutiveSave is too small, should be >= bdfOrder";
    int sizeUnknowns = std::max(this->bdfOrder()+2, this->numberOfConsecutiveSave());
    if ( M_history.size() != sizeUnknowns )
    {
        M_history.resize( sizeUnknowns );
        for ( uint8_type __i = 0; __i < M_history.size(); ++__i )
        {
            if ( !M_history[__i] )
            {
                M_history[__i] = M_space->elementPtr();
                M_history[__i]->zero();
            }
        }
    }

    if ( this->path().empty() )
    {
        this->setPathSave( (boost::format("%3%bdf_o_%1%_dt_%2%")
                            %this->bdfOrder()
                            %this->timeStep()
                            %this->bdfPrefix()  ).str() );
    }
    // in super::init() M_iteration is set back to 0
    super::init();

    if ( !this->isRestart() )
    {
        M_last_iteration_since_order_change = 1;
        switch ( M_strategyHighOrderStart )
        {
        default :
        case 0 : M_order_cur = M_order; break;
        case 1 : M_order_cur = 1; break;
        }
    }
    else
    {
        M_last_iteration_since_order_change = 1;

        switch ( M_strategyHighOrderStart )
        {
        default :
        case 0 :
        {
            M_order_cur = M_order;
        }
        break;
        case 1 :
        {
            M_order_cur = 1;
            for ( int i = 2; i<=M_iteration; ++i )
            {
                if ( ( ( i - M_last_iteration_since_order_change ) == M_iterations_between_order_change ) &&
                     M_order_cur < M_order )
                {
                    M_last_iteration_since_order_change = i;
                    ++M_order_cur;
                }
                if ( M_order_cur == M_order )
                    break;
            }
        }
        break;
        }

        fs::path dirPath = ( this->restartPath().empty() )? this->path() : this->restartPath()/this->path();

        const int niteration = this->iterationNumber();
        int pmax = std::min( M_numberOfConsecutiveSave/*M_order*/, M_iteration+1 );

        for ( int p = 0; p < pmax; ++p )
        {
            int iteration = (M_iteration-p);
            // Load files beginning at last iteration.
            if( this->isReverseLoad() )
                iteration = niteration - iteration;
            CHECK( iteration >= 0 )
                << "BDF init: negative iteration: "<< iteration
                << "( niter:"<< niteration << ", iter:"<< M_iteration << ")";

            if ( fileFormat() == "hdf5")
            {
#ifdef FEELPP_HAS_HDF5
                fs::path fname =  dirPath / (boost::format("%1%-%2%.h5") % M_name % iteration).str();
                VLOG(1) << "BDF HDF5 load solution iteration " << iteration
                        << " time " << M_time
                        << " from " << fname.string();
                M_history[p]->loadHDF5( fname.string() );
#else
                CHECK( false ) << "hdf5 not detected";
#endif
            }
            else if ( this->fileFormat() == "binary")
            {
                // create and open a character archive for output
                std::ostringstream ostr;
                if( M_rankProcInNameOfFiles )
                    ostr << M_name << "-" << iteration <<"-proc"<<this->worldComm().globalRank()<<"on"<<this->worldComm().globalSize();
                else
                    ostr << M_name << "-" << iteration;
                DVLOG(2) << "[Bdf::init()] load file: " << ostr.str() << "\n";

                std::ifstream ifs;
                ifs.open( dirPath/ostr.str() );

                // load data from archive
                boost::archive::binary_iarchive ia( ifs );
                ia >> *M_history[p];
            }
        }
    }
}


template <typename SpaceType>
Bdf<SpaceType>::~Bdf()
{}


template <typename SpaceType>
void
Bdf<SpaceType>::initialize( element_type const& u0 )
{
    M_time_values_map.clear();
    std::ostringstream ostr;

    if( M_rankProcInNameOfFiles )
        ostr << M_name << "-" << 0<<"-proc"<<this->worldComm().globalRank()<<"on"<<this->worldComm().globalSize();
    else
        ostr << M_name << "-" << 0;
    //M_time_values_map.insert( std::make_pair( 0, boost::make_tuple( 0, ostr.str() ) ) );
    //M_time_values_map.push_back( 0 );
    M_time_values_map.push_back( M_Ti );
    std::for_each(M_history.begin(), M_history.end(), 
                  [u0](auto& element) { *element = u0; });
    this->computePolyAndPolyDeriv();
    this->saveCurrent();
}
template <typename SpaceType>
void
Bdf<SpaceType>::initialize( std::vector<element_type> const& u0 )
{
    M_time_values_map.clear();
    std::ostringstream ostr;

    if( M_rankProcInNameOfFiles )
        ostr << M_name << "-" << 0<<"-proc"<<this->worldComm().globalRank()<<"on"<<this->worldComm().globalSize();
    else
        ostr << M_name << "-" << 0;
    //M_time_values_map.insert( std::make_pair( 0, boost::make_tuple( 0, ostr.str() ) ) );
    //M_time_values_map.push_back( 0 );
    setHistory( u0 );
    this->computePolyAndPolyDeriv();
    this->saveCurrent();
}
template <typename SpaceType>
void
Bdf<SpaceType>::initialize( unknowns_type const& uv0 )
{
    M_time_values_map.clear();
    std::ostringstream ostr;

    if( M_rankProcInNameOfFiles )
        ostr << M_name << "-" << 0<<"-proc"<<this->worldComm().globalRank()<<"on"<<this->worldComm().globalSize();
    else
        ostr << M_name << "-" << 0;
    //M_time_values_map.insert( std::make_pair( 0, boost::make_tuple( 0, ostr.str() ) ) );
    //M_time_values_map.push_back( 0);
    M_time_values_map.push_back( M_Ti );

    if ( uv0.size() == 1 )
    {
        std::for_each( M_history.begin(), M_history.end(), 
                       [value = *uv0[0]]( auto& element ) { *element = value; } );
    }
    else if ( uv0.size() > 1 )
    {
        std::copy( uv0.begin(), uv0.end(), M_history.begin() );
    }

    this->computePolyAndPolyDeriv();
    this->saveCurrent();
}

template <typename SpaceType>
double
Bdf<SpaceType>::start()
{
    if ( this->isRestart() )
        return this->restart();

    this->init();
    this->initialize( unknowns_type(0) );
    double ti = super::start();
    return ti;
}

template <typename SpaceType>
double
Bdf<SpaceType>::start( element_type const& u0 )
{
    if ( this->isRestart() )
        return this->restart();

    this->init();
    this->initialize( u0 );
    double ti = super::start();
    return ti;
}

template <typename SpaceType>
double
Bdf<SpaceType>::start( std::vector<element_type> const& u0 )
{
    if ( this->isRestart() )
        return this->restart();

    this->init();
    this->initialize( u0 );
    double ti = super::start();
    return ti;
}

template <typename SpaceType>
double
Bdf<SpaceType>::start( unknowns_type const& uv0 )
{
    if ( this->isRestart() )
        return this->restart();

    this->init();
    this->initialize( uv0 );
    double ti = super::start();
    return ti;
}

template <typename SpaceType>
double
Bdf<SpaceType>::restart()
{
    this->init();
    this->computePolyAndPolyDeriv();

    double ti = super::restart();

    return ti;
}

template <typename SpaceType>
typename Bdf<SpaceType>::element_type&
Bdf<SpaceType>::unknown( int i )
{
    DVLOG(2) << "[Bdf::unknown] id: " << i << " l2norm = " << M_history[i]->l2Norm() << "\n";
    return *M_history[i];
}

template <typename SpaceType>
typename Bdf<SpaceType>::element_ptrtype
Bdf<SpaceType>::unknownPtr( int i )
{
    DVLOG(2) << "[Bdf::unknown] id: " << i << " l2norm = " << M_history[i]->l2Norm() << "\n";
    return M_history[i];
}


template <typename SpaceType>
void
Bdf<SpaceType>::saveCurrent()
{
    if (!this->saveInFile()) return;

    bool doSave=false;
    for ( uint8_type i = 0; i < M_numberOfConsecutiveSave/*this->timeOrder()*/ && !doSave; ++i )
        {
            int iterTranslate = M_iteration + M_numberOfConsecutiveSave/*this->timeOrder()*/-(i+1);
            if (iterTranslate % this->saveFreq()==0) doSave=true;
        }

    if (!doSave) return;

    TSBaseMetadata bdfsaver( *this );
    bdfsaver.save();

    {
        int iteration = M_iteration;

        if ( this->fileFormat() == "hdf5")
        {
#ifdef FEELPP_HAS_HDF5
            M_history[0]->saveHDF5( (M_path_save / (boost::format("%1%-%2%.h5")%M_name %iteration).str() ).string() );
#else
            CHECK( false ) << "hdf5 not detected";
#endif
        }
        else if ( this->fileFormat() == "binary")
        {
            std::ostringstream ostr;

            if( M_rankProcInNameOfFiles )
                ostr << M_name << "-" << iteration<<"-proc"<<this->worldComm().globalRank()<<"on"<<this->worldComm().globalSize();
            else
                ostr << M_name << "-" << iteration;
            // load data from archive
            std::ofstream ofs( M_path_save / ostr.str() );
            boost::archive::binary_oarchive oa( ofs );
            oa << *M_history[0];
        }

    }
}

//! Load current unknown in a file (hdf5, binary, ...)
//! The filename depends on the bdf name and the time iteration.
//!
//! Notes:
//!     - that if bdf is set to `setReverse(true)`, the iteration 0 will be loaded.
//!     - You might desire to read a bdf backward. Then you have to pass `setReverseLoad(true)`
//!       to load existing unknowns from last iteration! It is require after a restart if time
//!       loop sense changed.
template <typename SpaceType>
void
Bdf<SpaceType>::loadCurrent()
{
    //TSBaseMetadata bdfsaver( *this );
    //bdfsaver.save();

    {
        const int niteration = this->iterationNumber();
        int iteration = M_iteration;
        // Load files beginning at last iteration.
        if( this->isReverseLoad() )
            iteration = niteration - iteration + 1;
        LOG(INFO) << "BDF iteration: "<< iteration << "( niter:"<< niteration << ", iter:"<< M_iteration << ")";
        CHECK( iteration >= 0 )
            << "BDF loadCurrent: negative iteration: "<< iteration
            << "( niter:"<< niteration << ", iter:"<< M_iteration << ")";

        if ( this->fileFormat() == "hdf5")
        {
#ifdef FEELPP_HAS_HDF5
            fs::path fname =  M_path_save / (boost::format("%1%-%2%.h5")%M_name %iteration).str();
            LOG(INFO) << "BDF HDF5 load solution iteration " << iteration
                      << " time " << M_time
                      << " from " << fname.string();
            if ( fs::exists( fname ) )
                M_history[0]->loadHDF5( fname.string() );
            else
                throw std::invalid_argument( fname.string() + " not found" );
#else
            CHECK( false ) << "hdf5 not detected";
#endif
        }
        else if ( this->fileFormat() == "binary")
        {

            std::ostringstream ostr;

            if( M_rankProcInNameOfFiles )
                ostr << M_name << "-" << iteration<<"-proc"<<this->worldComm().globalRank()<<"on"<<this->worldComm().globalSize();
            else
                ostr << M_name << "-" << iteration;

            std::ifstream ifs( M_path_save / ostr.str() );

            // load data from archive
            boost::archive::binary_iarchive ia( ifs );
            ia >> *M_history[0];
        }
    }
}

template <typename SpaceType>
template <typename container_type>
void Bdf<SpaceType>::shiftRight( typename space_type::template Element<value_type, container_type> const& new_unk )
{
    DVLOG( 2 ) << "shiftRight: inserting time " << this->time() << "s\n";
    super::shiftRight();

    // Shift all previously stored BDF data
    auto it = std::next( M_history.rbegin() );
    std::for_each( M_history.rbegin(), std::prev( M_history.rend() ), 
                   [&it]( auto& element ) { *element = *(*it); ++it; } );

    // u(t^{n}) coefficient is in M_history[0]
    *M_history[0] = new_unk;

    // Log the l2 norm for each unknown
    int i = 0;
    for ( const auto& t : M_history )
    {
        DVLOG( 2 ) << "[Bdf::shiftright] id: " << i << " l2norm = " << t->l2Norm() << "\n";
        ++i;
    }

    // Save newly stored BDF data
    this->saveCurrent();
}

template <typename SpaceType>
typename Bdf<SpaceType>::element_type const&
Bdf<SpaceType>::polyDeriv() const
{
    return *M_polyDeriv;
}

template <typename SpaceType>
typename Bdf<SpaceType>::element_type const&
Bdf<SpaceType>::poly() const
{
    return *M_poly;
}

template <typename SpaceType>
void
Bdf<SpaceType>::computePolyAndPolyDeriv()
{
    if ( !M_poly )
        M_poly = M_space->elementPtr();
    if ( !M_polyDeriv )
        M_polyDeriv = M_space->elementPtr();

    M_poly->zero();
    for ( int i = 0; i < this->timeOrder(); ++i )
        M_poly->add(  this->polyCoefficient( i ),  *M_history[ i ] );

    M_polyDeriv->zero();
    for (int i = 1; i <= this->timeOrder(); ++i)
        M_polyDeriv->add(-this->polyDerivCoefficient(i), *M_history[i-1]);  // known part only

    // Compute second derivative polynomial
    if (!M_polySecondDeriv)
        M_polySecondDeriv = M_space->elementPtr();
    M_polySecondDeriv->zero();

    if (this->timeOrder() >= 2)
    {
        for (uint8_type i = 1; i < M_alpha2[this->timeOrder() - 1].size(); ++i)
            M_polySecondDeriv->add(-this->polySecondDerivCoefficient(i), *M_history[i-1]);
    }
}

template <typename SpaceType>
void
Bdf<SpaceType>::updateDerivative( element_type & u, int i ) const
{
    u = firstDerivative( u );
}

template <typename ... Ts>
auto bdf( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && space = args.get(_space);
    po::variables_map const& vm = args.get_else(_vm,Environment::vm());
    std::string const& prefix = args.get_else(_prefix,"");
    std::string const& name = args.get_else(_name,"bdf");
    int order = args.get_else_invocable( _order, [&prefix,&vm](){ return ioption(_prefix=prefix,_name="bdf.order",_vm=vm); } );
    double initial_time = args.get_else_invocable( _initial_time, [&prefix,&vm](){ return doption(_prefix=prefix,_name="bdf.time-initial",_vm=vm); } );
    double final_time = args.get_else_invocable( _final_time, [&prefix,&vm](){ return doption(_prefix=prefix,_name="bdf.time-final",_vm=vm); } );
    double time_step = args.get_else_invocable( _time_step, [&prefix,&vm](){ return doption(_prefix=prefix,_name="bdf.time-step",_vm=vm); } );
    int strategy = args.get_else_invocable( _strategy, [&prefix,&vm](){ return ioption(_prefix=prefix,_name="bdf.strategy",_vm=vm); } );
    bool steady = args.get_else_invocable( _steady, [&prefix,&vm](){ return boption(_prefix=prefix,_name="bdf.steady",_vm=vm); } );
    bool reverse = args.get_else_invocable( _reverse, [&prefix,&vm](){ return boption(_prefix=prefix,_name="bdf.reverse",_vm=vm); } );
    bool restart = args.get_else_invocable( _restart, [&prefix,&vm](){ return boption(_prefix=prefix,_name="bdf.restart",_vm=vm); } );
    std::string const& restart_path = args.get_else_invocable( _restart_path, [&prefix,&vm](){ return soption(_prefix=prefix,_name="bdf.restart.path",_vm=vm); } );
    bool restart_at_last_save = args.get_else_invocable( _restart_at_last_save, [&prefix,&vm](){ return boption(_prefix=prefix,_name="bdf.restart.at-last-save",_vm=vm); } );
    bool save = args.get_else_invocable( _save, [&prefix,&vm](){ return boption(_prefix=prefix,_name="bdf.save",_vm=vm); } );
    int freq = args.get_else_invocable( _freq, [&prefix,&vm](){ return ioption(_prefix=prefix,_name="bdf.save.freq",_vm=vm); } );
    std::string const& format = args.get_else_invocable( _format,[&prefix,&vm](){ return soption(_prefix=prefix,_name="bdf.file-format",_vm=vm); } );
    bool rank_proc_in_files_name = args.get_else_invocable( _rank_proc_in_files_name, [&prefix,&vm](){ return boption(_prefix=prefix,_name="bdf.rank-proc-in-files-name",_vm=vm); } );
    int n_consecutive_save = args.get_else( _n_consecutive_save, order );

    using _space_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(space)>>>;
    auto thebdf = std::shared_ptr<Bdf<_space_type> >( new Bdf<_space_type>( space,name,prefix,vm,order,
                                                                            initial_time, final_time, time_step, steady, reverse,
                                                                            restart, restart_path, restart_at_last_save,
                                                                            save, freq, rank_proc_in_files_name, format, n_consecutive_save
                                                                            ) );

    return thebdf;
}


}
#endif

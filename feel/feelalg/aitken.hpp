/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Goncalo Pena <goncalo.pena@epfl.ch>
       Date: 15-07-2008

  Copyright (C) 2007-2008 Universite Joseph Fourier (Grenoble I)

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


#ifndef __AitkenExtrapolation
#define __AitkenExtrapolation 1

#include <string>
#include <vector>
#include <map>

#include <boost/assign/list_of.hpp>
#include <boost/assert.hpp>
#include <boost/timer.hpp>
#include <boost/tuple/tuple.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelalg/vector.hpp>


namespace Feel
{
/**
 * \class Aitken
 * \brief Aitken relaxation method for fixed point iterations
 *
 *\code
 * auto Xh = space_type::New(mesh);
 * auto residual = Xh->element();
 * auto u_old->element();
 * auto u_new->element();
 * AitkenType relaxmethod = (AitkenType)this->vm()["relaxmethod"].as<int>();
 * Aitken<space_type> aitken( Xh, relaxmethod,init_theta, tol );
 * //where init_theta is the initial value of relaxation parameter
 * // if relaxmethod=0(use the AITKEN_STANDARD method)
 * // else if relaxmethod=1(use the AITKEN_METHOD_1 method)
 * // relaxmethod=2(use the FIXED_RELAXATION_METHOD method)
 * in this last case the relaxation parameter theta remains fixed
 * during the itérations(theta = init_theta) and the particular case (init_theta=1)
 * corresponds to the case without relaxation
 * // initialize aitken
 * aitken.initialize( residual, u_new );
 * aitken.restart();
 * while(!aitken.isFinished())
 * {
 *
 *     u_old = u_new;
 *     commpute u_new;
 *     residual = u_new-u_old;
 *     u_new = aitken.apply(residual, u_new);
 *     aitken.printInfo();
 *     ++aitken;
 *
 *   }
 *
 *\endcode
 *
 * \author Goncalo Pena
 * \author Christophe Prud'homme
 * \author Vincent Chabannes
 */
template< typename fs_type >
class Aitken
{

public:

    typedef Aitken<fs_type> self_type;

    typedef fs_type functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;

    typedef typename functionspace_type::element_type element_type;

    typedef typename functionspace_type::template Element<typename functionspace_type::value_type,
            typename VectorUblas<typename functionspace_type::value_type>::range::type > element_range_type;

    /**
     * convergence_iteration_type:
     *  - int: iteration number
     *  - double: residual 2-norm
     *  - double: computing time for the iteration
     */
    typedef std::map<std::string, double> convergence_iteration_type;
    typedef std::map<int, convergence_iteration_type> convergence_type;

    /**
     * Constructor
     */
    Aitken( functionspace_ptrtype _Xh,
            AitkenType _aitkenType = AITKEN_STANDARD,
            double _failsafeParameter = 1.0,
            double _tol=1.0e-6 )
        :
        Xh( _Xh ),
        failsafeParameter( _failsafeParameter ),
        previousParameter( _failsafeParameter ),
        maxParameter( _failsafeParameter ),
        previousResidual( Xh, "previous residual" ),
        previousElement( Xh, "previous element" ),
        currentResidual( Xh, "current residual" ),
        currentElement( Xh, "current element" ),
        M_cptIteration( 1 ),
        M_aitkenType( _aitkenType ),
        M_tolerance( _tol ),
        M_residualConvergence( 1. ),
        M_hasConverged( false )
    {
    }

    /**
     * copy constructor
     */
    Aitken( Aitken const& tc )
        :
        Xh( tc.Xh ),
        failsafeParameter( tc.failsafeParameter ),
        previousParameter( tc.previousParameter ),
        maxParameter( tc.maxParameter ),
        previousResidual( tc.previousResidual ),
        previousElement( tc.previousElement ),
        currentResidual( tc.currentResidual ),
        currentElement( tc.currentElement ),
        M_cptIteration( tc.M_cptIteration ),
        M_aitkenType( tc.M_aitkenType ),
        M_tolerance( tc.M_tolerance ),
        M_residualConvergence( tc.M_residualConvergence ),
        M_hasConverged( tc.M_hasConverged )
    {
    }

    /**
     * destructor
     */
    ~Aitken() {}

    /**
     * initiliaze the aitken algorithm
     */
    BOOST_PARAMETER_MEMBER_FUNCTION(
        ( void ),
        initialize,
        tag,
        ( required
          ( residual, */*(element_type const&)*/ )
          ( currentElt,*/*(element_type const& )*/ ) ) )
    {
        initializeimpl( residual,currentElt );
    }

    /**
     * Compute theta and do a relaxation step : u^{n+1} = theta*u^{n+1} + (1-theta)*u^{n}
     */
    BOOST_PARAMETER_MEMBER_FUNCTION(
        ( element_type ),
        apply,
        tag,
        ( required
          ( residual, * )
          ( currentElt,* )
        ) //required
        ( optional
          ( forceRelaxation, ( bool ), false  )
        )//optional
    )
    {
        element_type newElt( Xh );
        applyimpl( newElt,residual,currentElt,forceRelaxation );
        return newElt;
    }

    /**
     * Compute theta and do a relaxation step : u^{n+1} = theta*u^{n+1} + (1-theta)*u^{n}
     */
    template< typename eltType >
    element_type operator()( element_type const& residual, eltType const& elem,bool _forceRelax=false )
    {
        element_type newElt( Xh );
        applyimpl( newElt,residual,elem,_forceRelax );
        return newElt;
    }

    /**
     * Compute theta and do a relaxation step : u^{n+1} = theta*u^{n+1} + (1-theta)*u^{n}
     */
    BOOST_PARAMETER_MEMBER_FUNCTION(
        ( void ),
        apply2,
        tag,
        ( required
          ( newElt, * )
          ( residual, * )
          ( currentElt,* )
        ) //required
        ( optional
          ( forceRelaxation, ( bool ), false  )
        ) //optional
    )
    {
        applyimpl( newElt,residual,currentElt,forceRelaxation );
    }

    //! set Aitken method type
    void setType( AitkenType t )
    {
        M_aitkenType = t;
    }

    //! \return Aitken method type
    AitkenType type( ) const
    {
        return M_aitkenType;
    }

    /**
     * shift current step to previous step. After the call, we are ready for the next step.
     */
    void shiftRight();

    /**
     * shift current step to previous step. After the call, we are ready for the next step.
     */
    self_type & operator++();

    /**
     * reset the previous parameter
     */
    void restart();

    /**
     * get theta
     */
    double theta()
    {
        return previousParameter;
    }

    /**
     * get number of iterations
     */

    uint nIterations()
    {
        return  M_cptIteration;
    }

    bool isFinished()
    {
        return M_hasConverged;
    }

    double residualNorm()
    {
        return M_residualConvergence;
    }

    void printInfo();

    /**
     * save converegence history
     * \param fname name of the file to save the convergence history
     */
    void saveConvergenceHistory( std::string const& fname ) const;

    void forceConvergence( bool b )
    {
        M_hasConverged=b;
    }

private:
    /**
     * initiliaze the aitken algorithm
     */
    void initializeimpl( element_type const& residual, element_type const& elem );

    /**
     * initiliaze the aitken algorithm
     */
    void initializeimpl( element_type const& residual, element_range_type const& elem );

    /**
     * compute a residual norm for convergence
     */
    void computeResidualNorm();

    /**
     * Set the current element
     */
    void setElement( element_type const& residual, element_type const& elem );

    /**
     * Set the current element
     */
    void setElement( element_type const& residual, element_range_type const& elem );

    /**
     * Compute theta and do a relaxation step : u^{n+1} = theta*u^{n+1} + (1-theta)*u^{n}
     */
    void applyimpl( element_type & new_elem,element_type const& residual, element_type const& elem,bool _forceRelax=false );

    /**
     * Compute theta and do a relaxation step : u^{n+1} = theta*u^{n+1} + (1-theta)*u^{n}
     */
    void applyimpl( element_range_type new_elem,element_type const& residual, element_range_type const& elem,bool _forceRelax=false );

    /**
     * Compute Aitken parameter
     */
    void calculateParameter();

    /**
     * Compute Aitken parameter
     */
    void calculateParameter( mpl::int_<AITKEN_STANDARD> /**/ );

    /**
     * Compute Aitken parameter
     */
    void calculateParameter( mpl::int_<AITKEN_METHOD_1> /**/ );

    /**
     * Do a relaxation step : u^{n+1} = theta*u^{n+1} + (1-theta)*u^{n}
     */
    //void relaxationStep( element_type& new_elem );
    template< typename eltType >
    void relaxationStep( eltType& new_elem );

    //! \return the convergence history
    convergence_type const& convergenceHistory() const
    {
        return M_convergence;
    }


private:

    /**
     * function space
     */
    functionspace_ptrtype Xh;

    double failsafeParameter, previousParameter, maxParameter;

    element_type previousResidual, previousElement, currentResidual, currentElement;

    uint M_cptIteration;
    AitkenType M_aitkenType;
    double M_tolerance;
    double M_residualConvergence;
    bool M_hasConverged;
    convergence_type M_convergence;


    mutable boost::timer M_timer;
};


//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::initializeimpl( element_type const& residual, element_type const& elem )
{
    previousResidual = residual;
    previousElement = elem;
}

//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::initializeimpl( element_type const& residual, element_range_type const& elem )
{
    previousResidual = residual;
    previousElement.zero();
    previousElement.add( 1.,elem );
    /*previousElement = vf::project(previousElement.functionSpace(),
      elements(previousElement.mesh()),
      vf::idv(elem) );*/
}

//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::computeResidualNorm()
{
    auto oldEltL2Norm = previousElement.l2Norm();

    if ( oldEltL2Norm > 1e-8 )
        M_residualConvergence = currentResidual.l2Norm()/oldEltL2Norm;

    else
        M_residualConvergence = currentResidual.l2Norm();

    if ( M_residualConvergence <  M_tolerance ) M_hasConverged=true;

}

//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::setElement( element_type const& residual, element_type const& elem )
{
    currentResidual = residual;
    currentElement = elem;
}

//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::setElement( element_type const& residual, element_range_type const& elem )
{
    currentResidual = residual;
    currentElement.zero();
    currentElement.add( 1.,elem );
    /*currentElement = vf::project(currentElement.functionSpace(),
      elements(currentElement.mesh()),
      vf::idv(elem) );*/
}

//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::applyimpl( element_type & new_elem,element_type const& residual, element_type const& elem,bool _forceRelax )
{

    setElement( residual,elem );

    if ( M_cptIteration>=2 )
    {
        computeResidualNorm();

        if ( !M_hasConverged || _forceRelax )
        {
            calculateParameter();
        }
    }

    if ( !M_hasConverged || _forceRelax )
        relaxationStep( new_elem );

}

//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::applyimpl( element_range_type new_elem,element_type const& residual, element_range_type const& elem,bool _forceRelax )
{

    setElement( residual,elem );

    if ( M_cptIteration>=2 )
    {
        computeResidualNorm();

        if ( !M_hasConverged  || _forceRelax )
        {
            calculateParameter();
        }
    }

    if ( !M_hasConverged  || _forceRelax )
        relaxationStep( new_elem );

}

//-----------------------------------------------------------------------------------------//
#if 0
template< typename fs_type >
void
Aitken<fs_type>::relaxationStep( element_type& new_elem )
{
    new_elem = currentResidual;
    //new_elem.scale( -previousParameter );
    new_elem.scale( -( 1-previousParameter ) );

    new_elem += currentElement;
}
#else
template< typename fs_type >
template< typename eltType >
void
Aitken<fs_type>::relaxationStep( eltType& new_elem )
{
    new_elem = currentResidual;
    //new_elem.scale( -previousParameter );
    new_elem.scale( -( 1-previousParameter ) );

    //new_elem += currentElement;
    new_elem.add( 1.,currentElement );
}

#endif

//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::shiftRight()
{
    // store convergence history
    double rel_residual = 1;

    if ( M_cptIteration > 1 )
        rel_residual = currentResidual.l2Norm()/M_convergence.find( 1 )->second.find( "absolute_residual" )->second;

    convergence_iteration_type cit =
        boost::assign::map_list_of
        ( "absolute_residual", currentResidual.l2Norm() )
        ( "relative_residual", rel_residual )
        ( "time", M_timer.elapsed() )
        ( "relaxation_parameter",  previousParameter );
    M_convergence.insert( std::make_pair( M_cptIteration, cit ) );
    M_timer.restart();

    previousResidual = currentResidual;
    previousElement = currentElement;


    ++M_cptIteration;
}

//-----------------------------------------------------------------------------------------//

template< typename fs_type >
typename Aitken<fs_type>::self_type &
Aitken<fs_type>::operator++()
{
    shiftRight();

    return *this;
}

//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::restart()
{
    if ( !M_convergence.empty() )
    {
        std::string entry = "relaxation_parameter";
        auto it = std::max_element( M_convergence.begin(), M_convergence.end(),
                                    [entry] ( std::pair<int, convergence_iteration_type> const& x,
                                              std::pair<int, convergence_iteration_type> const& y )
        {
            return x.second.find( entry )->second < y.second.find( entry )->second;
        } );

        if ( it != M_convergence.end() )
            maxParameter = std::max( maxParameter,it->second.find( entry )->second );

        previousParameter = maxParameter;
    }

    M_cptIteration=1;
    M_hasConverged=false;
    M_convergence.clear();
    M_timer.restart();
}

//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::printInfo()
{
    std::cout << "[Aitken] iteration : "<< M_cptIteration
              <<" theta=" << previousParameter
              <<" residualNorm : " << M_residualConvergence
              << "\n";

    Feel::Log() << "[Aitken] iteration : "<< M_cptIteration
                <<" theta=" << previousParameter
                <<" residualNorm : " << M_residualConvergence
                << "\n";
}

template< typename fs_type >
void
Aitken<fs_type>::saveConvergenceHistory( std::string const& fname ) const
{
    std::ofstream ofs( fname.c_str() );

    if ( ofs )
    {
        ofs.setf( std::ios::scientific );
        BOOST_FOREACH( auto cit, M_convergence )
        {
            // iteration
            ofs << std::setw( 6 ) << cit.first << " ";
            BOOST_FOREACH( auto it, cit.second )
            {
                ofs << std::setw( 20 ) << std::setprecision( 16 ) << it.second << " ";
            }
            ofs << "\n";
        }
    }

    else
    {
        std::cerr << "[Aitken] convergence history filename " << fname << " could not be opened"
                  << std::endl;
    }
}
//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::calculateParameter()
{
    if ( M_aitkenType==AITKEN_STANDARD )
        calculateParameter( mpl::int_<AITKEN_STANDARD>() );

    if ( M_aitkenType==AITKEN_METHOD_1 )
        calculateParameter( mpl::int_<AITKEN_METHOD_1>() );

    if ( M_aitkenType==FIXED_RELAXATION_METHOD )
        previousParameter = failsafeParameter;
}

//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::calculateParameter( mpl::int_<AITKEN_STANDARD> /**/ )
{

    element_type aux( Xh, "aux" );

    aux = currentResidual;
    aux -= previousResidual;

    double scalar = inner_product( aux, aux );

    aux.scale( 1.0/scalar );

    element_type aux2( Xh, "aux2" );

    aux2 = currentElement;
    aux2 -= previousElement;

    scalar = inner_product( aux2, aux );

    if ( scalar > 1 )
        scalar = /*previousParameter*/failsafeParameter;

    if ( scalar < 0 )
        scalar = /*previousParameter*/failsafeParameter;

    previousParameter = 1-scalar;
}

//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::calculateParameter( mpl::int_<AITKEN_METHOD_1> /**/ )
{

    element_type aux( Xh, "aux" );

    aux = currentResidual;
    aux -= previousResidual;

    double scalar = inner_product( aux, aux );

    aux.scale( 1.0/scalar );

    /*element_type aux2( Xh, "aux2");

      aux2 = currentElement;
      aux2 -= previousElement;*/

    scalar = inner_product( previousResidual , aux );
    scalar = -previousParameter*scalar;

#if 1
    if ( scalar > 1 )
        scalar = /*previousParameter;*/failsafeParameter;

    if ( scalar < 0 )
        scalar = /*previousParameter;*/failsafeParameter;
#endif
    previousParameter = scalar;

}


//-----------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------//

template <typename fs_type >
void
operator++( boost::shared_ptr<Aitken<fs_type> > & aitk )
{
    aitk->shiftRight();
}

//-----------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------//

template<typename SpaceType>
boost::shared_ptr<Aitken<SpaceType> >
aitkenNew( boost::shared_ptr<SpaceType> const& _space,
           AitkenType _type,
           double _init_theta,
           double _tol )
{

    boost::shared_ptr<Aitken<SpaceType> > Aitk( new Aitken<SpaceType>( _space,_type,_init_theta,_tol ) );

    return Aitk;
}



template<typename Args>
struct compute_aitken_return
{
    typedef typename boost::remove_reference<typename parameter::binding<Args, tag::space>::type>::type::element_type space_type;

    typedef Aitken<space_type> type;
    typedef boost::shared_ptr<type> ptrtype;
};


BOOST_PARAMETER_FUNCTION(
    ( typename compute_aitken_return<Args>::type ),
    aitken,
    tag,
    ( required
      ( space,    *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ) )
    )//required
    ( optional
      ( type, ( AitkenType ), AITKEN_STANDARD  )
      ( initial_theta, *( boost::is_arithmetic<mpl::_> ), 1.0  )
      ( tolerance, *( boost::is_arithmetic<mpl::_> ), 1.0e-6  )
    )//optional
)
{
    return *aitkenNew( space,type,initial_theta,tolerance );
}

BOOST_PARAMETER_FUNCTION(
    ( typename compute_aitken_return<Args>::ptrtype ),
    aitkenPtr,
    tag,
    ( required
      ( space,    *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ) )
    )//required
    ( optional
      ( type, ( AitkenType ), AITKEN_STANDARD  )
      ( initial_theta, *( boost::is_arithmetic<mpl::_> ), 1.0  )
      ( tolerance, *( boost::is_arithmetic<mpl::_> ), 1.0e-6  )
    )//optional
)
{
    return aitkenNew( space,type,initial_theta,tolerance );
}




} // End namespace Feel

#endif

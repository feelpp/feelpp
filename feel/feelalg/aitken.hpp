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
            double _tol=1.0e-6,
            double _minParam = 1e-4,
            size_type _maxit=1000 )
        :
        M_Xh( _Xh ),
        M_failsafeParameter( _failsafeParameter ),
        M_previousParameter( _failsafeParameter ),
        M_maxParameter( _failsafeParameter ),
        M_minParameter( _minParam ),
        M_previousResidual( M_Xh, "previous residual" ),
        M_previousElement( M_Xh, "previous element" ),
        M_currentResidual( M_Xh, "current residual" ),
        M_currentElement( M_Xh, "current element" ),
        M_cptIteration( 1 ),
        M_aitkenType( _aitkenType ),
        M_tolerance( _tol ),
        M_residualConvergence( 1. ),
        M_hasConverged( false ),
        M_maxit( _maxit )
    {
    }

    /**
     * copy constructor
     */
    Aitken( Aitken const& tc )
        :
        M_Xh( tc.M_Xh ),
        M_failsafeParameter( tc.M_failsafeParameter ),
        M_previousParameter( tc.M_previousParameter ),
        M_maxParameter( tc.M_maxParameter ),
        M_minParameter( tc.M_minParameter ),
        M_previousResidual( tc.M_previousResidual ),
        M_previousElement( tc.M_previousElement ),
        M_currentResidual( tc.M_currentResidual ),
        M_currentElement( tc.M_currentElement ),
        M_cptIteration( tc.M_cptIteration ),
        M_aitkenType( tc.M_aitkenType ),
        M_tolerance( tc.M_tolerance ),
        M_residualConvergence( tc.M_residualConvergence ),
        M_hasConverged( tc.M_hasConverged ),
        M_maxit( tc.M_maxit )
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
        element_type newElt( M_Xh );
        applyimpl( newElt,residual,currentElt,forceRelaxation );
        return newElt;
    }

    /**
     * Compute theta and do a relaxation step : u^{n+1} = theta*u^{n+1} + (1-theta)*u^{n}
     */
    template< typename eltType >
    element_type operator()( element_type const& residual, eltType const& elem,bool _forceRelax=false )
    {
        element_type newElt( M_Xh );
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
    double theta() const
    {
        return M_previousParameter;
    }

    /**
     * set theta
     */
    void setTheta(double v)
    {
        M_previousParameter=v;
    }

    /**
     * get number of iterations
     */

    size_type nIterations() const
    {
        return  M_cptIteration;
    }

    bool isFinished() const
    {
        return M_hasConverged || ( this->nIterations() > this->maxit() );
    }

    bool hasConverged() const
    {
        return M_hasConverged;
    }

    double residualNorm() const
    {
        return M_residualConvergence;
    }

    size_type maxit() const
    {
        return M_maxit;
    }
    /**
     * elements access
     */
    element_type const& previousResidual() const { return M_previousResidual; }
    element_type const& previousElement() const { return M_previousElement; }
    element_type const& currentResidual() const { return M_currentResidual; }
    element_type const& currentElement() const { return M_currentElement; }


    void printInfo() const;

    /**
     * save converegence history
     * \param fname name of the file to save the convergence history
     */
    void saveConvergenceHistory( std::string const& fname ) const;

    void forceConvergence( bool b )
    {
        M_hasConverged=b;
    }


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
     * Compute Aitken parameter
     */
    void calculateParameter();

    /**
     * Do a relaxation step : u^{n+1} = theta*u^{n+1} + (1-theta)*u^{n}
     */
    template< typename eltType >
    void relaxationStep( eltType& new_elem );

    /**
     * Do a relaxation step : u^{n+1} = theta*u^{n+1} + (1-theta)*u^{n}
     */
    void relaxationStep();

    //! \return the convergence history
    convergence_type const& convergenceHistory() const
    {
        return M_convergence;
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
    void calculateParameter( mpl::int_<AITKEN_STANDARD> /**/ );

    /**
     * Compute Aitken parameter
     */
    void calculateParameter( mpl::int_<AITKEN_METHOD_1> /**/ );


private:

    /**
     * function space
     */
    functionspace_ptrtype M_Xh;

    double M_failsafeParameter, M_previousParameter, M_maxParameter;
    double M_minParameter;

    element_type M_previousResidual, M_previousElement, M_currentResidual, M_currentElement;

    size_type M_cptIteration;
    AitkenType M_aitkenType;
    double M_tolerance;
    double M_residualConvergence;
    bool M_hasConverged;
    size_type M_maxit;
    convergence_type M_convergence;


    mutable boost::timer M_timer;
};


//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::initializeimpl( element_type const& residual, element_type const& elem )
{
    M_previousResidual = residual;
    M_previousElement = elem;
}

//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::initializeimpl( element_type const& residual, element_range_type const& elem )
{
    M_previousResidual = residual;
    M_previousElement.zero();
    M_previousElement.add( 1.,elem );
    /*M_previousElement = vf::project(M_previousElement.functionSpace(),
      elements(M_previousElement.mesh()),
      vf::idv(elem) );*/
}

//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::computeResidualNorm()
{
    auto oldEltL2Norm = M_previousElement.l2Norm();

    if ( oldEltL2Norm > 1e-8 )
        M_residualConvergence = M_currentResidual.l2Norm()/oldEltL2Norm;

    else
        M_residualConvergence = M_currentResidual.l2Norm();

    if ( M_residualConvergence <  M_tolerance ) M_hasConverged=true;

}

//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::setElement( element_type const& residual, element_type const& elem )
{
    M_currentResidual = residual;
    M_currentElement = elem;
}

//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::setElement( element_type const& residual, element_range_type const& elem )
{
    M_currentResidual = residual;
    M_currentElement.zero();
    M_currentElement.add( 1.,elem );
    /*M_currentElement = vf::project(M_currentElement.functionSpace(),
      elements(M_currentElement.mesh()),
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
template< typename fs_type >
template< typename eltType >
void
Aitken<fs_type>::relaxationStep( eltType& new_elem )
{
#if 0
    new_elem = M_currentResidual;
    //new_elem.scale( -M_previousParameter );
    new_elem.scale( -( 1-M_previousParameter ) );

    //new_elem += M_currentElement;
    new_elem.add( 1.,M_currentElement );
#else
    new_elem = M_previousElement;
    new_elem.add(M_previousParameter,M_currentResidual);

    M_currentElement = new_elem;
#endif
}

//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::relaxationStep()
{
    M_currentElement = M_previousElement;
    M_currentElement.add(M_previousParameter,M_currentResidual);
}


//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::shiftRight()
{
    // store convergence history
    double rel_residual = 1;

    if ( M_cptIteration > 1 )
        rel_residual = M_currentResidual.l2Norm()/M_convergence.find( 1 )->second.find( "absolute_residual" )->second;

    convergence_iteration_type cit =
        boost::assign::map_list_of
        ( "absolute_residual", M_currentResidual.l2Norm() )
        ( "relative_residual", rel_residual )
        ( "time", M_timer.elapsed() )
        ( "relaxation_parameter",  M_previousParameter );
    M_convergence.insert( std::make_pair( M_cptIteration, cit ) );
    M_timer.restart();

    M_previousResidual = M_currentResidual;
    M_previousElement = M_currentElement;


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
            M_maxParameter = std::max( M_maxParameter,it->second.find( entry )->second );

        M_previousParameter = M_maxParameter;
    }

    M_cptIteration=1;
    M_hasConverged=false;
    M_convergence.clear();
    M_timer.restart();
}

//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::printInfo() const
{
    std::cout << "[Aitken] iteration : "<< M_cptIteration
              <<" theta=" << M_previousParameter
              <<" residualNorm : " << M_residualConvergence
              << "\n";

    LOG(INFO) << "[Aitken] iteration : "<< M_cptIteration
                <<" theta=" << M_previousParameter
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
        M_previousParameter = M_failsafeParameter;
}

//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::calculateParameter( mpl::int_<AITKEN_STANDARD> /**/ )
{

    element_type aux( M_Xh, "aux" );

    aux = M_currentResidual;
    aux -= M_previousResidual;

    double scalar = inner_product( aux, aux );

    aux.scale( 1.0/scalar );

    element_type aux2( M_Xh, "aux2" );

    aux2 = M_currentElement;
    aux2 -= M_previousElement;

    scalar = inner_product( aux2, aux );

    if ( scalar > 1 )
        scalar = /*M_previousParameter*/M_failsafeParameter;

    if ( scalar < 0 )
        scalar = /*M_previousParameter*/M_failsafeParameter;

    M_previousParameter = 1-scalar;
}

//-----------------------------------------------------------------------------------------//

template< typename fs_type >
void
Aitken<fs_type>::calculateParameter( mpl::int_<AITKEN_METHOD_1> /**/ )
{

    element_type aux( M_Xh, "aux" );

    aux = M_currentResidual;
    aux -= M_previousResidual;

    double scalar = inner_product( aux, aux );

    aux.scale( 1.0/scalar );

    /*element_type aux2( M_Xh, "aux2");

      aux2 = M_currentElement;
      aux2 -= M_previousElement;*/

    scalar = inner_product( M_previousResidual , aux );
    scalar = -M_previousParameter*scalar;

#if 1
    if ( scalar > 1 )
        scalar = /*M_previousParameter;*/M_failsafeParameter;

    if ( scalar < /*0*/M_minParameter )
        scalar = /*M_previousParameter;*/M_failsafeParameter;
#endif
    M_previousParameter = scalar;

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
aitkenImpl( boost::shared_ptr<SpaceType> const& _space,
            std::string _type,
            double _init_theta,
            double _tol,
            double _minParam,
            size_type _maxit )
{

    AitkenType myAitkenType=AITKEN_METHOD_1;
    if ( _type == "standard" )
        myAitkenType = AITKEN_STANDARD;
    else if ( _type == "method1" )
        myAitkenType = AITKEN_METHOD_1;
    else if ( _type == "fixed-relaxation" )
        myAitkenType = FIXED_RELAXATION_METHOD;
    else
        CHECK( false ) << "invalid aitken type " << _type;

    boost::shared_ptr<Aitken<SpaceType> > Aitk( new Aitken<SpaceType>( _space,myAitkenType,_init_theta,_tol,_minParam,_maxit ) );
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
      ( prefix,*,"" )
      ( type, *( boost::is_convertible<mpl::_,std::string> ), soption(_prefix=prefix,_name="aitken.type") /* AITKEN_STANDARD*/ )
      ( initial_theta, *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="aitken.initial_theta") /*1.0*/  )
      ( min_theta, *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="aitken.min_theta") /*1e-4*/  )
      ( tolerance, *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="aitken.tol") /*1.0e-6*/  )
      ( maxit,      ( size_type ), ioption(_prefix=prefix,_name="aitken.maxit") )
    )//optional
)
{
    return *aitkenImpl( space,type,initial_theta,tolerance,min_theta,maxit );
}

BOOST_PARAMETER_FUNCTION(
    ( typename compute_aitken_return<Args>::ptrtype ),
    aitkenPtr,
    tag,
    ( required
      ( space,    *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ) )
    )//required
    ( optional
      ( prefix,*,"" )
      ( type, *( boost::is_convertible<mpl::_,std::string> ), soption(_prefix=prefix,_name="aitken.type") /* AITKEN_STANDARD*/ )
      ( initial_theta, *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="aitken.initial_theta") /*1.0*/  )
      ( min_theta, *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="aitken.min_theta") /*1e-4*/  )
      ( tolerance, *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="aitken.tol") /*1.0e-6*/  )
      ( maxit,      ( size_type ), ioption(_prefix=prefix,_name="aitken.maxit") )
    )//optional
)
{
    return aitkenImpl( space,type,initial_theta,tolerance,min_theta,maxit );
}




} // End namespace Feel

#endif

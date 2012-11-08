/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-07-15

  Copyright (C) 2010-2011 Université Joseph Fourier (Grenoble I)

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
   \file simget.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-07-15
 */
#ifndef __Simget_H
#define __Simget_H 1

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/about.hpp>
#include <feel/feelcore/factory.hpp>
#include <feel/feelcore/singleton.hpp>

namespace Feel
{
namespace ptree = boost::property_tree;

class Application;

/**
 * \class Simget
 * \brief Simulation Object
 *
 * A Simget is an object that provides two flavors of \c run() member function
 * - \c run() without any argument which simulates a blackbox \f$ F \f$ whitout
 *   any outputs or inputs
 * - <tt> run( double* X, int P, double* Y, int N ) </tt> which simulates a
 *    blackbox with input/output relationship \f$ Y = F(X) \f$ with \f$ Y \in
 *    \mathbb{R}^N\f$ and \f$ X \in \mathbb{R}^P\f$.
 *
 * @author Christophe Prud'homme
 * @see
 */
class Simget
{
public:


    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * constructor with a \c variables_map
     */
    Simget();

    //! destructor
    virtual ~Simget() {}

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    Simget& operator=( Simget const & o )
    {
        if ( this != &o )
        {
        }

        return *this;
    }
    //@}

    /** @name Accessors
     */
    //@{

    //! return the name of the simget
    virtual std::string name() const
    {
        return M_about.appName();
    }

    //! \return the mpi communicator
    mpi::communicator comm() const
    {
        return M_comm;
    }

    //! \return the \c variables_map
    po::variables_map const& vm() const
    {
        return M_vm;
    }

    //! \return the \c AboutData object
    AboutData const& about() const
    {
        return M_about;
    }

    //! return the mesh size
    double meshSize() const
    {
        return M_meshSize;
    }
    //! return the mesh size
    double meshSizeInit() const
    {
        return M_meshSizeInit;
    }

    //! return the refinement level
    int level() const
    {
        return M_level;
    }

    //! return the statistics associated to the simget after calling run
    ptree::ptree const& stats() const
    {
        return M_stats;
    }
    //! return the statistics associated to the simget after calling run
    ptree::ptree& stats()
    {
        return M_stats;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    //! set the mesh size
    void setMeshSize( double h )
    {
        M_meshSize= h;
    }

    //! set the initial mesh size
    void setMeshSizeInit( double h )
    {
        M_meshSizeInit = h;
    }

    //! set the refinment level if applicable
    void setLevel( int level )
    {
        M_level= level;
    }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * simply execute the simget
     */
    virtual void run() = 0;

    /**
     * models the input/output relation \f$ Y=F(X) \f$
     */
    virtual void run( const double* X, unsigned long P, double* Y, unsigned long N ) { run(); };

    /**
     * print statistics from simget
     */
    void print( std::ostream& out, std::vector<ptree::ptree> & stats );

    //@}

protected:

    friend class Application;

    /**
     * change repository.
     */
    Simget& changeRepository( boost::format fmt );

protected:
    double M_meshSize;
    double M_meshSizeInit;
    int M_level;
    ptree::ptree M_stats;
private:
    mpi::communicator M_comm;
    po::variables_map M_vm;
    AboutData M_about;


};


/**
 * Simget factory
 */
typedef Singleton< Factory< Simget, std::string > > SimgetFactory;

#if 0
template<typename SimgetType>
Simget*
createSimget( AboutData cponst& about )
{
    return new SimgetType( about );
}

#define REGISTER_SIMGET_IN_FACTORY( simget, simgetname, about )         \
    SimgetFactory::instance().registerProduct( simgetname, &createSimget<simget>( about )

#if defined( FEELPP_HAS_OCT_H )
#define OCTNAME(name,dim, order) BOOST_PP_CAT(BOOST_PP_CAT(BOOST_PP_CAT(name,dim),_),order)

#define REGISTER_SIMGET( dim, order )                                   \
    DEFUN_DLD (OCTNAME(residualestimator_,dim,order), args, nargout, "Residual Estimator for the Laplacian") \
    {                                                                   \
        int nargin = args.length ();                                    \
        if (nargin != 1)                                                \
            print_usage ();                                             \
        else                                                            \
        {                                                               \
            NDArray A = args(0).array_value ();                         \
            dim_vector dims = A.dims();                                 \
            dim_vector ydims(1);                                        \
            ydims(0)=4;                                                 \
                                                                        \
            NDArray Y(ydims);                                           \
            {                                                           \
                static bool is_init = false;                            \
                if ( !is_init )                                         \
                {                                                       \
                    if (!Feel::Environment::initialized() )             \
                        new Feel::Environment();                        \
                    is_init = true;                                     \
                }                                                       \
                boost::shared_ptr<ResidualEstimator<dim,order> > OCTNAME(app_,dim,order)( new ResidualEstimator<dim,order>( makeAbout() ) ); \
                std::vector<double> x( dims(0) );                       \
                for( int j = 0; j < dims(0); ++j )                      \
                    {                                                   \
                        x[j] = A(j);                                    \
                        std::cout << "x["<< j << "]=" << x[j] << "\n";  \
                    }                                                   \
                std::vector<double> y( 4 );                             \
                OCTNAME(app_,dim,order)->run( x.data(), dims(1), y.data(), 4 ); \
                Y(0)=y[0];                                              \
                Y(1)=y[1];                                              \
                Y(2)=y[2];                                              \
                Y(3)=y[3];                                              \
            }                                                           \
            if (! error_state)                                          \
                return octave_value (Y);                                \
        }                                                               \
        return octave_value_list ();                                    \
    }



#endif /* FEELPP_HAS_OCT_H */
#endif // 0
}
#endif /* __Simget_H */

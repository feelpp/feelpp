/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel++ library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2006-07-04

  Copyright (C) 2006 EPFL
  Copyright (C) 2007-2010 Université Joseph Fourier (Grenoble 1)


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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file bench1.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2006-07-04
 */
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

//#include <boost/test/unit_test.hpp>
//using boost::unit_test::test_suite;

#include <boost/program_options.hpp>
#include <boost/lambda/bind.hpp>

#include <feel/feelconfig.h>


#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>


#include <feel/feelalg/backend.hpp>
#include <feel/feelvf/vf.hpp>

#if defined(FEELPP_HAS_GOOGLE_PROFILER_H)
#include <google/profiler.h>
#endif // FEELPP_HAS_GOOGLE_PROFILER_H




namespace Feel
{

using namespace Feel::vf;


/*!
 * \class Bench1
 * \brief Benchmark for assembly performance in 1D, 2D and 3D
 *
 * The benchmark is called as follows:
 * \code
 * bench1 --dim=1 --hsize=0.1
 * bench1 --dim=2 --hsize=0.1
 * bench1 --dim=3 --hsize=0.1
 * \endcode
 *
 * For a fixed \p hsize the bench is run in 1D, 2D or 3D (by default
 * 1D)
 */
class Bench1
    :
public Application
{
public:


    /** @name Typedefs
     */
    //@{

    typedef Application super;

    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_ptrtype vector_ptrtype;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    FEELPP_DONT_INLINE Bench1( int argc,
                               char** argv,
                               AboutData const& ad,
                               po::options_description const& od );

    FEELPP_DONT_INLINE ~Bench1();

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    void run();

    //@}



protected:

private:

    template<typename FSType> FEELPP_DONT_INLINE void A( boost::shared_ptr<FSType> const& Xh, mpl::int_<1> );
    template<typename FSType> FEELPP_DONT_INLINE void A( boost::shared_ptr<FSType> const& Xh, mpl::int_<2> );
    template<typename FSType> FEELPP_DONT_INLINE void A( boost::shared_ptr<FSType> const& Xh, mpl::int_<3> );
    template<typename FSType> FEELPP_DONT_INLINE void R( boost::shared_ptr<FSType> const& Xh );
    template<typename FSType> FEELPP_DONT_INLINE void D( boost::shared_ptr<FSType> const& Xh );
    template<typename FSType> FEELPP_DONT_INLINE void DR( boost::shared_ptr<FSType> const& Xh );

    template<typename FSType> FEELPP_DONT_INLINE void ADR( boost::shared_ptr<FSType> const& Xh, mpl::int_<1> );
    template<typename FSType> FEELPP_DONT_INLINE void ADR( boost::shared_ptr<FSType> const& Xh, mpl::int_<2> );
    template<typename FSType> FEELPP_DONT_INLINE void ADR( boost::shared_ptr<FSType> const& Xh, mpl::int_<3> );
    /**
     * 1D performance test
     */
    FEELPP_DONT_INLINE void run1d();

    /**
     * 2D performance test
     */
    FEELPP_DONT_INLINE void run2d();

    /**
     * 3D performance test
     */
    FEELPP_DONT_INLINE void run3d();

    /**
     * dimension independant code
     */
    template<typename MeshType, int Order>
    FEELPP_DONT_INLINE void bench1( boost::shared_ptr<MeshType> & mesh );

private:

    backend_ptrtype M_backend;
    double meshSize;
};


}

/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-08-10

  Copyright (C) 2009-2011 Université Joseph Fourier (Grenoble I)

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
   \file opusscm.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-08-10
 */
#include <boost/tuple/tuple_io.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>

#include <opusdata.hpp>
#include <opusmodelbase.hpp>
#include <opusmodelfactory.hpp>
#include <opusmodelrb.hpp>

#include <feel/feelcrb/crb.hpp>
#include <feel/feelcrb/crbmodel.hpp>



namespace Feel
{
template<int OrderU, int OrderP, int OrderT> class OpusModelRB;

po::options_description makeEadsSCMOptions();

/**
 * \class EadsSCMApp
 * \brief Eads SCM application
 *
 * This class implements the Opus application, getting the command
 * line arguments and running the actual code.
 *
 * @author Christophe Prud'homme
 */
class EadsSCMApp   : public Application
{
    typedef Application super;
public:

    typedef CRBModel<OpusModelRB<2,1,2> > opusmodel_type;
    //typedef CRBModel<Heat1D> opusmodel_type;
    typedef boost::shared_ptr<opusmodel_type> opusmodel_ptrtype;
    typedef CRBSCM<opusmodel_type> scm_type;
    typedef boost::shared_ptr<scm_type> scm_ptrtype;

    EadsSCMApp( AboutData const& ad, po::options_description const& od );
    EadsSCMApp( int argc, char** argv, AboutData const& ad, po::options_description const& od );

    void init();
    void run( std::ofstream& os, scm_type::parameter_type const& mu, int K );
    void run();

private:

    opusmodel_ptrtype M_opusmodel;
    scm_ptrtype M_scm;
}; // Opus

} // Feel


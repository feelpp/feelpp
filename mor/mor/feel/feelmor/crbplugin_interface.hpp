//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 09 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
#ifndef FEELPP_CRBPLUGIN_INTERFACE_HPP
#define FEELPP_CRBPLUGIN_INTERFACE_HPP 1


#include <string>
#include <feel/feelconfig.h>
#include <feel/feelcore/info.hpp>
#include <feel/feelcore/singleton.hpp>
#include <feel/feelmor/parameterspace.hpp>
#include <feel/feelmor/crbdata.hpp>
#include <feel/feelmor/crbenums.hpp>
#include <feel/feelmesh/meshbase.hpp>
#include <feel/feeldiscr/doftablebase.hpp>
#include <feel/feelmor/crbmodelbase.hpp>
#include <feel/feelalg/vector.hpp>

namespace Feel {

//!
//! Interface for CRB plugins
//!
class CRBPluginAPI
{
public:

    //!
    //! @return name of the plugin
    //!
    virtual std::string const& name() const = 0;

    //!
    //! load database
    //!
    virtual void loadDB( std::string const&, crb::load l = crb::load::rb ) = 0;

    //!
    //! load database from db id.
    //! we already know the \p name() of the model, we only need the \p id and the \p root repository
    //!
    virtual void loadDBFromId( std::string const& id, crb::load l = crb::load::rb, std::string const& root = Environment::rootRepository() ) = 0;

    //!
    //! load last created or modified database
    //! @param last defines whether to load last created or modified DB
    //! @param load defines the type of data data loaded from the DB
    //! @param root root repository of the DB
    //!
    virtual void loadDBLast ( crb::last last = crb::last::modified,
                              crb::load load = crb::load::rb,
                              std::string const& root = Environment::rootRepository() ) = 0;

    //!
    //! return true if some data from the DB is loaded, false otherwise
    //!
    virtual bool isDBLoaded() const = 0;

    //!
    //! return true if  the Reduced basis DB is loaded
    //!
    virtual bool isReducedBasisModelDBLoaded() const = 0;

    //!
    //! return true if the Finite Element DB is loaded
    //!
    virtual bool isFiniteElementModelDBLoaded() const = 0;

    //!
    //! return true if all the DB data is loaded
    //!
    virtual bool isAllLoaded() const = 0;

    //!
    //! @return the parameter space
    //!
    virtual std::shared_ptr<ParameterSpaceX> parameterSpace() const = 0;

    //!
    //! @return the crb model
    //!
    virtual std::shared_ptr<CRBModelBase> crbmodel() const = 0;

    //!
    //! @return the meshes
    //!
    virtual std::vector<std::shared_ptr<MeshBase<>>> meshes() const = 0;

    //!
    //! @return the doftable with datamap
    //!
    virtual std::pair<std::vector<std::shared_ptr<DofTableBase<>>>,std::shared_ptr<DataMap<>>> doftables() const = 0;

    //!
    //! @return an element of the fe space
    //!
    virtual std::shared_ptr<Vector<double>> feElement() const = 0;

    //!
    //! @return a list of sub-element which compose the element of the fe space
    //!
    virtual std::vector<std::shared_ptr<Vector<double>> > feSubElements( std::shared_ptr<Vector<double>> u ) const = 0;

    //!
    //! @return reduced basis functions primal
    //!
    virtual std::vector<std::shared_ptr<Vector<double>>> reducedBasisFunctionsPrimal() const = 0;

    //!
    //! @return reduced basis functions dual
    //!
    virtual std::vector<std::shared_ptr<Vector<double>>> reducedBasisFunctionsDual() const = 0;

    //!
    //! run the crb online code
    //! @param mu parameter at which CRB is evaluated
    //! @param time collection of timers (prediction, error bounds)
    //! @param eps max tolerance for the evaluation
    //! @param N number of basis functions
    //! @param print_rb_matrix print the reduced matrix
    //!
    virtual CRBResults run( ParameterSpaceX::Element const& mu,
                            vectorN_type & time, double eps , int N, bool print_rb_matrix ) const = 0;

    virtual CRBResults run( ParameterSpaceX::Element const& mu,
                            double eps , int N, bool print_rb_matrix ) const
        {
            vectorN_type t;
            return this->run( mu, t, eps, N, print_rb_matrix );
        }

    //!
    //! run plugin over a parameter sampling
    //!
    virtual std::vector<CRBResults> run( std::vector<ParameterSpaceX::Element> const& S,
                                         double eps , int N, bool print_rb_matrix ) const = 0;

    //!
    //! expand the rb field to fe field
    //!
    virtual void expansion( vectorN_type const& uRB, Vector<double> & uFE, int N=-1 ) const = 0;

    //!
    //! initialize the exporter
    //!
    virtual void initExporter()  = 0;

    //!
    //! write to disk the results of the run
    //! @param name the name of the field to exporter
    //! @param results the results to be exported
    //!
    virtual void exportField( std::string const & name, CRBResults const& results ) = 0;

    //!
    //! save exporter
    //!
    virtual void saveExporter() const = 0;

    //!
    //! virtual destructor
    //!
    virtual ~CRBPluginAPI() {}

protected:

    //!
    //! set the name of the plugin
    //!
    virtual void setName( std::string const& ) = 0;
};

//!
//! command line options for plugins
//!
po::options_description crbPluginOptions( std::string const& prefix = "" );

using crbpluginapi_create_t = std::shared_ptr<CRBPluginAPI> ();
using crbpluginapi_create_ft = boost::function<crbpluginapi_create_t>;

//!
//! create the plugin \p name from plugin library \p pluginlibname located in \p dirname
//! @param n name of the plugin
//! @param pluginlibname name of the plugin library
//! @param dirname location of the plugin library
//!
std::shared_ptr<CRBPluginAPI> factoryCRBPlugin( std::string const& n,
                                                  std::string const& pluginlibname = "",
                                                  std::string const& dirname = Info::libdir().string()
                                                  );

}

#endif

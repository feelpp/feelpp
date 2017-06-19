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
#include <feel/feelcore/singleton.hpp>
#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelcrb/crbdata.hpp>
#if defined(FEELPP_HAS_VTK)
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#endif

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
    virtual void loadDB( std::string ) = 0;
    
    //!
    //! @return the parameter space
    //!
    virtual boost::shared_ptr<ParameterSpaceX> parameterSpace() const = 0;

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
    
#if defined(FEELPP_HAS_VTK)
    //!
    //! exporter to VTK data structure 
    //!
    virtual vtkSmartPointer<vtkUnstructuredGrid> exporterVTK() const = 0;
#endif // FEELPP_HAS_VTK

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

using crbpluginapi_create_t = boost::shared_ptr<CRBPluginAPI> ();
using crbpluginapi_create_ft = boost::function<crbpluginapi_create_t>;
boost::shared_ptr<CRBPluginAPI> factoryCRBPlugin( std::string const& dirname, std::string const& n );
//!
//! 
//!
boost::function<crbpluginapi_create_t> makeCRBPluginCreator( std::string const& dirname, std::string const& pluginname );

#define FEELPP_CRB_PLUGIN( classname, strname )                         \
class FEELPP_EXPORT BOOST_PP_CAT( classname, Plugin ) : public CRBPlugin<classname>   \
{                                                                       \
public:                                                                 \
    using this_t = BOOST_PP_CAT(classname,Plugin);                      \
    BOOST_PP_CAT(classname,Plugin)()                                    \
        :                                                               \
        CRBPlugin<classname>( BOOST_PP_STRINGIZE( strname ) )           \
        {}                                                              \
                                                                        \
    /* Factory method */                                                \
    static boost::shared_ptr<this_t> create()                           \
        {                                                               \
            return boost::shared_ptr<this_t>( new this_t() );           \
        }                                                               \
};                                                                      \
                                                                        \
                                                                        \
BOOST_DLL_ALIAS( Feel::BOOST_PP_CAT(classname,Plugin)::create, create_crbplugin )


#define FEELPP_CRB_PLUGIN_TEMPLATE( classname, classtemplate, strname ) \
class FEELPP_EXPORT BOOST_PP_CAT( classname, Plugin ) : public CRBPlugin<classtemplate> \
{                                                                       \
public:                                                                 \
    using this_t = BOOST_PP_CAT(classname,Plugin);                      \
    BOOST_PP_CAT(classname,Plugin)()                                    \
        :                                                               \
        CRBPlugin<classtemplate>( BOOST_PP_STRINGIZE( strname ) )       \
        {}                                                              \
                                                                        \
    /* Factory method */                                                \
    static boost::shared_ptr<this_t> create()                           \
        {                                                               \
            return boost::shared_ptr<this_t>( new this_t() );           \
        }                                                               \
};                                                                      \
                                                                        \
                                                                        \
BOOST_DLL_ALIAS( Feel::BOOST_PP_CAT(classname,Plugin)::create, create_crbplugin )


}

#endif

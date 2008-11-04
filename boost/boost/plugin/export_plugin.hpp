// Copyright Vladimir Prus 2004.
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_EXPORT_PLUGIN_VP_2004_08_25
#define BOOST_EXPORT_PLUGIN_VP_2004_08_25

#include <string>
#include <map>

#include <boost/config.hpp>
#include <boost/any.hpp>

#include <boost/plugin/config.hpp>
#include <boost/plugin/concrete_factory.hpp>

///////////////////////////////////////////////////////////////////////////////
#define BOOST_PLUGIN_EXPORT(BaseType, ActualType, name)                       \
    extern "C" BOOST_PLUGIN_EXPORT_API std::map<std::string, boost::any>&     \
        BOOST_PLUGIN_API boost_exported_plugins_list();                       \
    namespace {                                                               \
        struct boost_plugin_exporter1 {                                       \
            boost_plugin_exporter1()                                          \
            {                                                                 \
                static boost::plugin::concrete_factory<BaseType, ActualType> cf;\
                boost::plugin::abstract_factory<BaseType>* w = &cf;           \
                boost_exported_plugins_list().insert(std::make_pair(name, w));\
            }                                                                 \
        } boost_plugin_exporter_instance1;                                    \
    }                                                                         \
    /**/

///////////////////////////////////////////////////////////////////////////////
#define BOOST_PLUGIN_EXPORT_LIST()                                            \
    extern "C" BOOST_PLUGIN_EXPORT_API std::map<std::string, boost::any>&     \
        BOOST_PLUGIN_API boost_exported_plugins_list()                        \
    {                                                                         \
        static std::map<std::string, boost::any> r;                           \
        return r;                                                             \
    }                                                                         \
    /**/

#endif


// Copyright Vladimir Prus 2004.
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_PLUGIN_FACTORY_VP_2004_08_25
#define BOOST_PLUGIN_FACTORY_VP_2004_08_25

#include <utility>
#include <stdexcept>
#include <string>
#include <utility>

#include <boost/config.hpp>
#include <boost/any.hpp>
#include <boost/function.hpp>

#include <boost/plugin/virtual_constructors.hpp>
#include <boost/plugin/abstract_factory.hpp>
#include <boost/plugin/dll.hpp>
#include <boost/plugin/export_plugin.hpp>

namespace boost { namespace plugin {

    ///////////////////////////////////////////////////////////////////////////
    namespace detail {

        template<typename BasePlugin>
        std::pair<abstract_factory<BasePlugin> *, dll_handle>
        get_abstract_factory(dll const& d, std::string const &class_name)
        {
            typedef typename remove_pointer<get_plugins_list_type>::type PointedType;
            typedef boost::function<void (get_plugins_list_type)> DeleterType;

            std::pair<get_plugins_list_type, DeleterType> f =
                d.get<get_plugins_list_type, DeleterType>("boost_exported_plugins_list");

            exported_plugins_type& e = (*f.first)();
            typename exported_plugins_type::iterator it = e.find(class_name);
            if (it != e.end()) {
#if 0
                std::cout << "found object " << class_name << "\n"
                          << "typeid(BasePlugin) = " << typeid( BasePlugin ).name() << "\n"
                          << "typeid(it) = " << it->second.type().name() << "\n"
                          << "abstract_factory<BasePlugin> " << typeid( abstract_factory<BasePlugin> ).name() << "\n";
#endif

#if BOOST_VERSION == 103301
                    abstract_factory<BasePlugin>** xw =
                     boost::unsafe_any_cast<abstract_factory<BasePlugin> *>(&(*it).second);
#else
                    abstract_factory<BasePlugin>** xw =
                     boost::any_cast<abstract_factory<BasePlugin> *>(&(*it).second);
#endif // BOOST_VERSION

                if (!xw) {
                    throw std::logic_error("Boost.Plugin: Can't cast to the right factor type\n");
                }

                abstract_factory<BasePlugin> *w = *xw;
                return make_pair(w, boost::shared_ptr<PointedType>(f.first, f.second));
            }
            else {
                BOOST_PLUGIN_OSSTREAM str;
                str << "Boost.Plugin: Class '" << class_name
                    << "' was not found in the shared library '"
                    << d.get_name() << "'.";

                throw std::logic_error(BOOST_PLUGIN_OSSTREAM_GETSTRING(str));
            }
        }

        ///////////////////////////////////////////////////////////////////////
        struct plugin_factory_item_base
        {
            void create(int******);

        protected:
            dll m_dll;
        };

        ///////////////////////////////////////////////////////////////////////
        template<typename BasePlugin, typename Base, typename Parameters>
        struct plugin_factory_item;

        template<typename BasePlugin, typename Base>
        struct plugin_factory_item<BasePlugin, Base, boost::mpl::list<> >
        :   public Base
        {
            BasePlugin* create(std::string const& name)
            {
                std::pair<abstract_factory<BasePlugin> *, dll_handle> r =
                    get_abstract_factory<BasePlugin>(this->m_dll, name);
                return r.first->create(r.second);
            }
        };

        template<typename BasePlugin, typename Base, typename A1>
        struct plugin_factory_item<BasePlugin, Base, boost::mpl::list<A1> >
        :   public Base
        {
            using Base::create;
            BasePlugin* create(std::string const& name, A1 a1)
            {
                std::pair<abstract_factory<BasePlugin> *, dll_handle> r =
                    get_abstract_factory<BasePlugin>(this->m_dll, name);
                return r.first->create(r.second, a1);
            }
        };

        template<typename BasePlugin, typename Base, typename A1, typename A2>
        struct plugin_factory_item<BasePlugin, Base, boost::mpl::list<A1, A2> >
        :   public Base
        {
            using Base::create;
            BasePlugin* create(std::string const& name, A1 a1, A2 a2)
            {
                std::pair<abstract_factory<BasePlugin> *, dll_handle> r =
                    get_abstract_factory<BasePlugin>(this->m_dll, name);
                return r.first->create(r.second, a1, a2);
            }
        };
    }

///////////////////////////////////////////////////////////////////////////////
//
//  Bring in the remaining plugin_factory_item definitions for parameter
//  counts greater 2
//
///////////////////////////////////////////////////////////////////////////////
#include <boost/plugin/detail/plugin_factory_impl.hpp>

    ///////////////////////////////////////////////////////////////////////////
    template<class BasePlugin>
    struct plugin_factory
    :   public boost::mpl::inherit_linearly <
            typename virtual_constructors<BasePlugin>::type,
            detail::plugin_factory_item<BasePlugin,
                boost::mpl::placeholders::_, boost::mpl::placeholders::_>,
            detail::plugin_factory_item_base
        >::type
    {
        plugin_factory(dll const& d)
        {
            this->m_dll = d;
        }
    };

///////////////////////////////////////////////////////////////////////////////
}}  // namespace boost::plugin

#endif

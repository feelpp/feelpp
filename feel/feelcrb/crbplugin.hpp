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
//! @date 10 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!

#ifndef FEELPP_CRBPLUGIN_HPP
#define FEELPP_CRBPLUGIN_HPP 1

#include <boost/dll/alias.hpp> // for BOOST_DLL_ALIAS
#include <boost/algorithm/string.hpp>
#include <feel/options.hpp>
#include <feel/feelcrb/crbplugin_interface.hpp>
#include <feel/feelcrb/crbenums.hpp>
#include <feel/feelcrb/modelcrbbase.hpp>
#include <feel/feelcrb/crb_trilinear.hpp>


namespace Feel {

namespace crbplugin_details
{

struct DofTablesComposite
{
    template <typename T>
    void operator()( boost::shared_ptr<T> const& x ) const
        {
            doftables.push_back( x->dof() );
        }

    mutable std::vector<boost::shared_ptr<DofTableBase>> doftables;
};

template <typename SpaceType>
std::vector<boost::shared_ptr<DofTableBase>>
doftables( boost::shared_ptr<SpaceType> const& space, typename std::enable_if< !SpaceType::is_composite >::type* = nullptr )
{
    std::vector<boost::shared_ptr<DofTableBase>> dt;
    dt.push_back( space->dof() );
    return dt;
}
template <typename SpaceType>
std::vector<boost::shared_ptr<DofTableBase>>
doftables( boost::shared_ptr<SpaceType> const& space, typename std::enable_if< SpaceType::is_composite >::type* = nullptr )
{
    Feel::crbplugin_details::DofTablesComposite dtc;
    boost::fusion::for_each( space->functionSpaces(),dtc );
    return dtc.doftables;
}


template <typename ElementType>
struct SubElementsComposite
{
    SubElementsComposite( boost::shared_ptr<ElementType> _uFE ) : uFE( _uFE ) {}

    template <typename T>
    void operator()( T const& t ) const
        {
            subelements.push_back( uFE->template elementPtr<T::value>() );
        }

    boost::shared_ptr<ElementType> uFE;
    mutable std::vector<boost::shared_ptr<Vector<typename ElementType::value_type>> > subelements;
};

template <typename ElementType>
std::vector<boost::shared_ptr<Vector<typename ElementType::value_type>> >
subelements( boost::shared_ptr<ElementType> uFE, typename std::enable_if< !ElementType::is_composite >::type* = nullptr )
{
    std::vector<boost::shared_ptr<Vector<typename ElementType::value_type>> > res;
    res.push_back( uFE );
    return res;
}
template <typename ElementType>
std::vector<boost::shared_ptr<Vector<typename ElementType::value_type>> >
subelements( boost::shared_ptr<ElementType> uFE, typename std::enable_if< ElementType::is_composite >::type* = nullptr )
{
    mpl::range_c<int,0,ElementType::functionspace_type::nSpaces> keySpaces;
    Feel::crbplugin_details::SubElementsComposite<ElementType> sec(uFE);
    boost::fusion::for_each( keySpaces, sec );
    return sec.subelements;
}
}


//!
//! Generic Plugin for CRB applications
//!
template<typename ModelT, template <class T> class CRBModelT = CRBModel, template <class T> class CRBT = CRB>
class CRBPlugin : public CRBPluginAPI
{
public:

    using model_t = ModelT;
    using crbmodel_type = CRBModelT<model_t>;
    using crb_type = CRBT<crbmodel_type> ;
    using crb_ptrtype = boost::shared_ptr<crb_type>;
    using mesh_t = typename model_t::mesh_type;
    using exporter_ptr_t = boost::shared_ptr<Exporter<mesh_t> >;
    
    CRBPlugin( std::string const& name )
        :
        M_name( name ),
        M_crb( boost::make_shared<crb_type>(crb::stage::online) )
        {
        }

    std::string const& name() const override
        {
            return M_name;
        }
    void loadDB( std::string const& filename, crb::load l ) override
        {
            if ( !fs::exists( filename ) )
                throw std::invalid_argument("file does not exist");

            if ( ( l == crb::load::all ) ||  (l == crb::load::fe ) )
                M_crb->setLoadBasisFromDB( true );
            else
                M_crb->setLoadBasisFromDB( false );
            M_crb->loadJson( filename );

            LOG(INFO) << "Loaded DB CRBPlugin " << filename;
        }

    void loadDBFromId( std::string const& id, crb::load l = crb::load::rb, std::string const& root = Environment::rootRepository() ) override
        {
            if ( !fs::exists( root ) )
                throw std::invalid_argument(std::string("root repository ") + root + " does not exist");
            fs::path crbdb = fs::path(root) / "crbdb";
            if ( !fs::exists( crbdb ) )
                throw std::invalid_argument(std::string("crbdb repository ") + crbdb.string() + " does not exist");

            fs::path dbbasedir = crbdb / fs::path(name()) ;
            if ( !fs::exists( dbbasedir ) && !fs::is_directory(dbbasedir) )
                throw std::invalid_argument(std::string("db directory ") + dbbasedir.string() + " does not exist");
            // either id provides the full directory or part of it
            // try first full path
            typedef std::vector<fs::path> vec;
            vec d;

            std::copy(fs::directory_iterator(dbbasedir), fs::directory_iterator(), std::back_inserter(d));
            std::sort(d.begin(), d.end());

            //std::cout << "dbbasedir=" << dbbasedir.string()<< " id=" << id << std::endl;
            for( auto& dbdir : d )
            {
                //std::cout << "dbdir = " << dbdir.string() << std::endl;
                
                if ( boost::ends_with( dbdir.string(), id ) )
                {
                    fs::path dbfilename = dbdir / fs::path(name() + ".crb.json");
                    if (!fs::exists(dbfilename))
                        continue;
                    loadDB( dbfilename.string(), l );
                    return;
                }
            }
            throw std::invalid_argument(std::string("Database for ") + name() + " with id " + id + " not found");
        }

    boost::shared_ptr<ParameterSpaceX> parameterSpace() const override
        {
            DCHECK( M_crb ) << "DB not loaded";
            return M_crb->Dmu();
        }

    boost::shared_ptr<CRBModelBase> crbmodel() const override
        {
            DCHECK( M_crb ) << "DB not loaded";
            return M_crb->model();
        }
    std::vector<boost::shared_ptr<MeshBase>> meshes() const override
        {
            DCHECK( M_crb ) << "DB not loaded";
            std::vector<boost::shared_ptr<MeshBase>> m;
            m.push_back( M_crb->model()->rBFunctionSpace()->functionSpace()->mesh() );
            // TODO composite case with several meshes
            return m;
        }


    std::pair<std::vector<boost::shared_ptr<DofTableBase>>,boost::shared_ptr<DataMap>> doftables() const override
        {
            DCHECK( M_crb ) << "DB not loaded";
            if ( M_crb->model() && M_crb->model()->rBFunctionSpace() && M_crb->model()->rBFunctionSpace()->functionSpace() )
                return std::make_pair( Feel::crbplugin_details::doftables( M_crb->model()->rBFunctionSpace()->functionSpace() ),
                                       M_crb->model()->rBFunctionSpace()->functionSpace()->dof() );
            else
                return std::make_pair( std::vector<boost::shared_ptr<DofTableBase>>(), boost::shared_ptr<DataMap>() );
        }

    boost::shared_ptr<Vector<double>> feElement() const override
        {
            DCHECK( M_crb ) << "DB not loaded";
            if ( M_crb->model() && M_crb->model()->rBFunctionSpace() && M_crb->model()->rBFunctionSpace()->functionSpace() )
                return M_crb->model()->rBFunctionSpace()->functionSpace()->elementPtr();
            else
                return boost::shared_ptr<Vector<double>>();
        }

    std::vector<boost::shared_ptr<Vector<double>> > feSubElements( boost::shared_ptr<Vector<double>> u ) const override
        {
            DCHECK( M_crb ) << "DB not loaded";
            auto uFE = boost::dynamic_pointer_cast< typename crbmodel_type::space_type::element_type >( u );
            CHECK( uFE ) << "dynamic_pointer_cast fails : wrong type of element u";
            return Feel::crbplugin_details::subelements( uFE );
        }

    std::vector<boost::shared_ptr<Vector<double>>> reducedBasisFunctionsPrimal() const override
        {
            DCHECK( M_crb ) << "DB not loaded";
            auto const& rbPrimal = M_crb->model()->rBFunctionSpace()->primalRB();
            int nBasis = rbPrimal.size();
            std::vector<boost::shared_ptr<Vector<double>>> res( nBasis );
            for ( int k=0;k<nBasis;++k )
            {
                auto u = M_crb->model()->rBFunctionSpace()->functionSpace()->elementPtr();
                *u = rbPrimal[k];
                res[k] = u;
            }
            return res;
        }
    std::vector<boost::shared_ptr<Vector<double>>> reducedBasisFunctionsDual() const override
        {
            DCHECK( M_crb ) << "DB not loaded";
            auto const& rbDual = M_crb->model()->rBFunctionSpace()->dualRB();
            int nBasis = rbDual.size();
            std::vector<boost::shared_ptr<Vector<double>>> res( nBasis );
            for ( int k=0;k<nBasis;++k )
            {
                auto u = M_crb->model()->rBFunctionSpace()->functionSpace()->elementPtr();
                *u = rbDual[k];
                res[k] = u;
            }
            return res;
        }

    CRBResults run( ParameterSpaceX::Element const& mu,
                    vectorN_type & time, double eps , int N, bool print_rb_matrix ) const override
        {
            DCHECK( M_crb ) << "DB not loaded";
            return M_crb->run( mu, time, eps, N, print_rb_matrix );
        }

    std::vector<CRBResults> run( std::vector<ParameterSpaceX::Element> const& S,
                                 double eps , int N, bool print_rb_matrix ) const override
        {
            
            return M_crb->run( S, eps, N, print_rb_matrix );
        }
    
    void expansion( vectorN_type const& uRB, Vector<double> & uFE,  int N ) const override
        {
            DCHECK( M_crb ) << "DB not loaded";
            auto uRBforExpansion = M_crb->model()->rBFunctionSpace()->element();
            uRBforExpansion.container() = uRB;
            M_crb->model()->rBFunctionSpace()->expansion( uRBforExpansion, uFE, N );
        }



    void initExporter() override
        {
            fieldExporter = exporter( _mesh=M_crb->model()->rBFunctionSpace()->mesh() );
        }

    void exportField( std::string const& name, CRBResults const& res ) override
        {
            auto const& solution = res.coefficients();
            auto uN = M_crb->model()->rBFunctionSpace()->element();
            //auto const& uN = solutions.template get<0>();
            uN.container() = solution;//.get<0>().back();
            auto sol = M_crb->model()->rBFunctionSpace()->expansion( uN );
            fieldExporter->add( name, sol );
        }

    void saveExporter() const override
        {
            fieldExporter->save();
        }

protected:
    void setName( std::string const& name ) override
        {
            M_name = name;
        }
protected:
    std::string M_name;
    crb_ptrtype M_crb;
    exporter_ptr_t fieldExporter;
};



#define FEELPP_CRB_PLUGIN( classname, strname )                         \
    class FEELPP_EXPORT BOOST_PP_CAT( classname, Plugin ) : public CRBPlugin<classname> \
{                                                                       \
public:                                                                 \
    using this_t = BOOST_PP_CAT(classname,Plugin);                      \
    BOOST_PP_CAT(classname,Plugin)()                                    \
        :                                                               \
        CRBPlugin<classname>( strname )                                 \
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

#define FEELPP_CRBTRILINEAR_PLUGIN( classname, strname )                         \
    class FEELPP_EXPORT BOOST_PP_CAT( classname, Plugin ) :             \
        public CRBPlugin<classname,CRBModelTrilinear,CRBTrilinear>      \
{                                                                       \
public:                                                                 \
    using this_t = BOOST_PP_CAT(classname,Plugin);                      \
    BOOST_PP_CAT(classname,Plugin)()                                    \
        :                                                               \
        CRBPlugin<classname,CRBModelTrilinear,CRBTrilinear>( strname )  \
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

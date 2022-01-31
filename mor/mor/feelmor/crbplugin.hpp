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
#include <feel/feelmor/crbplugin_interface.hpp>
#include <feel/feelmor/crbenums.hpp>
#include <feel/feelmor/modelcrbbase.hpp>
#include <feel/feelmor/crb_trilinear.hpp>
#include <feel/feelmor/crbsaddlepoint.hpp>
#include <feel/feelmor/crbmodelsaddlepoint.hpp>
#include <feel/feelmor/biotsavartrb.hpp>


namespace Feel {

namespace crbplugin_details
{

struct DofTablesComposite
{
    DofTablesComposite( std::vector<std::shared_ptr<DofTableBase<>>> & doftables )
        :
        M_doftables( doftables )
        {}

    template <typename T>
    void operator()( std::shared_ptr<T> const& x ) const
        {
            M_doftables.push_back( x->dof() );
        }
private :
    std::vector<std::shared_ptr<DofTableBase<>>> & M_doftables;
};

template <typename SpaceType>
std::vector<std::shared_ptr<DofTableBase<>>>
doftables( std::shared_ptr<SpaceType> const& space, typename std::enable_if< !SpaceType::is_composite >::type* = nullptr )
{
    std::vector<std::shared_ptr<DofTableBase<>>> dt;
    dt.push_back( space->dof() );
    return dt;
}
template <typename SpaceType>
std::vector<std::shared_ptr<DofTableBase<>>>
doftables( std::shared_ptr<SpaceType> const& space, typename std::enable_if< SpaceType::is_composite >::type* = nullptr )
{
    std::vector<std::shared_ptr<DofTableBase<>>> dt;
    boost::fusion::for_each( space->functionSpaces(), Feel::crbplugin_details::DofTablesComposite( dt ) );
    return dt;
}


template <typename ElementType>
struct SubElementsComposite
{
    explicit SubElementsComposite( std::shared_ptr<ElementType> const& uFE,
                                   std::vector<std::shared_ptr<Vector<typename ElementType::value_type>> > & subelements )
        :
        M_uFE( uFE ), M_subelements( subelements ) {}

    template <typename T>
    void operator()( T const& t ) const
        {
            M_subelements.push_back( M_uFE->template elementPtr<T::value>() );
        }
private :
    std::shared_ptr<ElementType> M_uFE;
    std::vector<std::shared_ptr<Vector<typename ElementType::value_type>> > & M_subelements;
};

template <typename ElementType>
std::vector<std::shared_ptr<Vector<typename ElementType::value_type>> >
subelements( std::shared_ptr<ElementType> uFE, typename std::enable_if< !ElementType::is_composite >::type* = nullptr )
{
    std::vector<std::shared_ptr<Vector<typename ElementType::value_type>> > res;
    res.push_back( uFE );
    return res;
}
template <typename ElementType>
std::vector<std::shared_ptr<Vector<typename ElementType::value_type>> >
subelements( std::shared_ptr<ElementType> uFE, typename std::enable_if< ElementType::is_composite >::type* = nullptr )
{
    std::vector<std::shared_ptr<Vector<typename ElementType::value_type>> > subelements;
    mpl::range_c<int,0,ElementType::functionspace_type::nSpaces> keySpaces;
    boost::fusion::for_each( keySpaces, Feel::crbplugin_details::SubElementsComposite<ElementType>(uFE, subelements ) );
    return subelements;
}
}


//!
//! Generic Plugin for CRB applications
//!
template<typename ModelT, template <class T> class CRBModelT = CRBModel, template <class T> class AlgoT = CRB, template <class T> class AlgoBaseT = CRB>
class CRBPlugin : public CRBPluginAPI
{
public:

    using model_t = ModelT;
    using crbmodel_type = CRBModelT<model_t>;
    using crb_type = AlgoBaseT<crbmodel_type> ;
    using crb_ptrtype = std::shared_ptr<crb_type>;
    using method_t = AlgoT<crbmodel_type>;
    using mesh_t = typename model_t::mesh_type;
    using exporter_ptr_t = std::shared_ptr<Exporter<mesh_t> >;
    
    CRBPlugin( std::string const& name )
        :
        M_name( name ),
        M_load( crb::load::none )
        {
            M_crb = method_t::New(name, crb::stage::online);
        }

    std::string const& name() const override
        {
            return M_name;
        }
    void loadDB( std::string const& filename, crb::load l ) override
        {
            M_load = l;
            M_crb->loadDB( filename, l );
        }

    void loadDBFromId( std::string const& id, crb::load l = crb::load::rb, std::string const& root = Environment::rootRepository() ) override
        {
            M_load = l;
            M_crb->loadDBFromId( id, l, root );
        }
    
    void loadDBLast( crb::last last = crb::last::modified, crb::load l = crb::load::rb, std::string const& root = Environment::rootRepository() ) override
        {
            M_load = l;
            M_crb->loadDBLast( last, l, root );
        }

    bool isDBLoaded() const override { return M_load != crb::load::none; }
    
    bool isReducedBasisModelDBLoaded() const override { return (M_load == crb::load::rb) || (M_load == crb::load::all); }

    bool isFiniteElementModelDBLoaded() const override { return (M_load == crb::load::fe) || (M_load == crb::load::all); }

    bool isAllLoaded() const override { return M_load == crb::load::all; }
    
    std::shared_ptr<ParameterSpaceX> parameterSpace() const override
        {
            DCHECK( M_crb ) << "DB not loaded";
            return M_crb->Dmu();
        }

    std::shared_ptr<CRBModelBase> crbmodel() const override
        {
            DCHECK( M_crb ) << "DB not loaded";
            return M_crb->model();
        }
    std::vector<std::shared_ptr<MeshBase<>>> meshes() const override
        {
            DCHECK( M_crb ) << "DB not loaded";
            std::vector<std::shared_ptr<MeshBase<>>> m;
            m.push_back( M_crb->model()->rBFunctionSpace()->functionSpace()->mesh() );
            // TODO composite case with several meshes
            return m;
        }


    std::pair<std::vector<std::shared_ptr<DofTableBase<>>>,std::shared_ptr<DataMap<>>> doftables() const override
        {
            DCHECK( M_crb ) << "DB not loaded";
            if ( M_crb->model() && M_crb->model()->rBFunctionSpace() && M_crb->model()->rBFunctionSpace()->functionSpace() )
                return std::make_pair( Feel::crbplugin_details::doftables( M_crb->model()->rBFunctionSpace()->functionSpace() ),
                                       M_crb->model()->rBFunctionSpace()->functionSpace()->dof() );
            else
                return std::make_pair( std::vector<std::shared_ptr<DofTableBase<>>>(), std::shared_ptr<DataMap<>>() );
        }

    std::shared_ptr<Vector<double>> feElement() const override
        {
            DCHECK( M_crb ) << "DB not loaded";
            if ( M_crb->model() && M_crb->model()->rBFunctionSpace() && M_crb->model()->rBFunctionSpace()->functionSpace() )
                return M_crb->model()->rBFunctionSpace()->functionSpace()->elementPtr();
            else
                return std::shared_ptr<Vector<double>>();
        }

    std::vector<std::shared_ptr<Vector<double>> > feSubElements( std::shared_ptr<Vector<double>> u ) const override
        {
            DCHECK( M_crb ) << "DB not loaded";
            auto uFE = std::dynamic_pointer_cast< typename crbmodel_type::space_type::element_type >( u );
            CHECK( uFE ) << "dynamic_pointer_cast fails : wrong type of element u";
            return Feel::crbplugin_details::subelements( uFE );
        }

    std::vector<std::shared_ptr<Vector<double>>> reducedBasisFunctionsPrimal() const override
        {
            DCHECK( M_crb ) << "DB not loaded";
            auto const& rbPrimal = M_crb->model()->rBFunctionSpace()->primalRB();
            int nBasis = rbPrimal.size();
            std::vector<std::shared_ptr<Vector<double>>> res( nBasis );
            for ( int k=0;k<nBasis;++k )
            {
                res[k] = rbPrimal[k];
            }
            return res;
        }
    std::vector<std::shared_ptr<Vector<double>>> reducedBasisFunctionsDual() const override
        {
            DCHECK( M_crb ) << "DB not loaded";
            auto const& rbDual = M_crb->model()->rBFunctionSpace()->dualRB();
            int nBasis = rbDual.size();
            std::vector<std::shared_ptr<Vector<double>>> res( nBasis );
            for ( int k=0;k<nBasis;++k )
            {
                res[k] = rbDual[k];;
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
    crb::load M_load;
};



#define FEELPP_CRB_PLUGIN( classname, strname )                         \
    class FEELPP_EXPORT BOOST_PP_CAT( classname, Plugin ) : public CRBPlugin<classname> \
{                                                                       \
public:                                                                 \
    using this_t = BOOST_PP_CAT(classname,Plugin);                      \
    BOOST_PP_CAT(classname,Plugin)()                                    \
        :                                                               \
        CRBPlugin<classname>( BOOST_PP_STRINGIZE( strname ) )           \
        {}                                                              \
                                                                        \
    /* Factory method */                                                \
    static std::shared_ptr<this_t> create()                           \
        {                                                               \
            return std::shared_ptr<this_t>( new this_t() );           \
        }                                                               \
};                                                                      \
                                                                        \
                                                                        \
BOOST_DLL_ALIAS( Feel::BOOST_PP_CAT(classname,Plugin)::create, BOOST_PP_CAT(create_crbplugin_,strname) )

#define FEELPP_CRBTRILINEAR_PLUGIN( classname, strname )                \
    class FEELPP_EXPORT BOOST_PP_CAT( classname, Plugin ) :             \
        public CRBPlugin<classname,CRBModelTrilinear,CRBTrilinear>      \
{                                                                       \
public:                                                                 \
    using this_t = BOOST_PP_CAT(classname,Plugin);                      \
    BOOST_PP_CAT(classname,Plugin)()                                    \
        :                                                               \
        CRBPlugin<classname,CRBModelTrilinear,CRBTrilinear>( BOOST_PP_STRINGIZE( strname ) ) \
        {}                                                              \
                                                                        \
    /* Factory method */                                                \
    static std::shared_ptr<this_t> create()                           \
        {                                                               \
            return std::shared_ptr<this_t>( new this_t() );           \
        }                                                               \
};                                                                      \
                                                                        \
                                                                        \
BOOST_DLL_ALIAS( Feel::BOOST_PP_CAT(classname,Plugin)::create, BOOST_PP_CAT(create_crbplugin_,strname) )


#define FEELPP_CRBSADDLEPOINT_PLUGIN( classname, strname )              \
    class FEELPP_EXPORT BOOST_PP_CAT( classname, Plugin ) :             \
        public CRBPlugin<classname,CRBModelSaddlePoint,CRBSaddlePoint>  \
{                                                                       \
public:                                                                 \
    using this_t = BOOST_PP_CAT(classname,Plugin);                      \
    BOOST_PP_CAT(classname,Plugin)()                                    \
        :                                                               \
        CRBPlugin<classname,CRBModelSaddlePoint,CRBSaddlePoint>( BOOST_PP_STRINGIZE( strname ) ) \
        {}                                                              \
                                                                        \
    /* Factory method */                                                \
    static std::shared_ptr<this_t> create()                           \
        {                                                               \
            return std::shared_ptr<this_t>( new this_t() );           \
        }                                                               \
};                                                                      \
                                                                        \
                                                                        \
BOOST_DLL_ALIAS( Feel::BOOST_PP_CAT(classname,Plugin)::create, BOOST_PP_CAT(create_crbplugin_,strname) )


#define FEELPP_BIOTSAVART_PLUGIN( model, classname, strname )            \
    class FEELPP_EXPORT BOOST_PP_CAT( classname, Plugin ) :             \
        public CRBPlugin<model,classname,BiotSavartRB,CRBPluginBase>                 \
{                                                                       \
public:                                                                 \
    using this_t = BOOST_PP_CAT(classname,Plugin);                      \
    BOOST_PP_CAT(classname,Plugin)()                                    \
        :                                                               \
        CRBPlugin<model,classname,BiotSavartRB,CRBPluginBase>( BOOST_PP_STRINGIZE( strname ) ) \
        {}                                                              \
                                                                        \
    /* Factory method */                                                \
    static std::shared_ptr<this_t> create()                           \
        {                                                               \
            return std::shared_ptr<this_t>( new this_t() );           \
        }                                                               \
};                                                                      \
                                                                        \
                                                                        \
BOOST_DLL_ALIAS( Feel::BOOST_PP_CAT(classname,Plugin)::create, BOOST_PP_CAT(create_crbplugin_,strname) )



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
    static std::shared_ptr<this_t> create()                           \
        {                                                               \
            return std::shared_ptr<this_t>( new this_t() );           \
        }                                                               \
};                                                                      \
                                                                        \
                                                                        \
BOOST_DLL_ALIAS( Feel::BOOST_PP_CAT(classname,Plugin)::create, BOOST_PP_CAT(create_crbplugin_,strname) )


}
#endif

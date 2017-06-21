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

#include <feel/options.hpp>
#include <feel/feelcrb/crbplugin_interface.hpp>
#include <feel/feelcrb/modelcrbbase.hpp>


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
template<typename ModelT>
class CRBPlugin : public CRBPluginAPI
{
public:

    using model_t = ModelT;
    typedef Feel::CRBModel< model_t > crbmodel_type;
    typedef Feel::CRB<crbmodel_type> crb_type;
    using mesh_t = typename model_t::mesh_type;
    using exporter_ptr_t = boost::shared_ptr<Exporter<mesh_t> >;

    CRBPlugin( std::string const& name )
        :
        M_name( name )
        {
        }

    std::string const& name() const override
        {
            return M_name;
        }
    void loadDB( std::string filename ) override
        {
            if ( !fs::exists( filename ) )
                throw std::invalid_argument("file does not exist");
            std::cout << "Loading " << filename << std::endl;
            boost::shared_ptr<model_t> model( new model_t );
            model->loadJson( filename, "crbmodel" );
            std::cout << "loaded model\n";
            boost::shared_ptr<crbmodel_type> crbmodel( new crbmodel_type( model, false ) );
            crbmodel->loadJson( filename, "crbmodel" );
            std::cout << "loaded crbmodel\n";

            crb.reset( new crb_type );
            crb->setTruthModel( crbmodel );
            crb->loadJson( filename );
            std::cout << "Loaded " << filename << std::endl;
        }

    boost::shared_ptr<ParameterSpaceX> parameterSpace() const override
        {
            DCHECK( crb ) << "DB not loaded";
            return crb->Dmu();
        }

    boost::shared_ptr<CRBModelBase> crbmodel() const override
        {
            DCHECK( crb ) << "DB not loaded";
            return crb->model();
        }
    std::vector<boost::shared_ptr<MeshBase>> meshes() const override
        {
            DCHECK( crb ) << "DB not loaded";
            std::vector<boost::shared_ptr<MeshBase>> m;
            m.push_back( crb->model()->rBFunctionSpace()->functionSpace()->mesh() );
            // TODO composite case with several meshes
            return m;
        }

    std::pair<std::vector<boost::shared_ptr<DofTableBase>>,boost::shared_ptr<DataMap>> doftables() const override
        {
            DCHECK( crb ) << "DB not loaded";
            if ( crb->model() && crb->model()->rBFunctionSpace() && crb->model()->rBFunctionSpace()->functionSpace() )
                return std::make_pair( Feel::crbplugin_details::doftables( crb->model()->rBFunctionSpace()->functionSpace() ),
                                       crb->model()->rBFunctionSpace()->functionSpace()->dof() );
            else
                return std::make_pair( std::vector<boost::shared_ptr<DofTableBase>>(), boost::shared_ptr<DataMap>() );
        }

    boost::shared_ptr<Vector<double>> feElement() const override
        {
            DCHECK( crb ) << "DB not loaded";
            if ( crb->model() && crb->model()->rBFunctionSpace() && crb->model()->rBFunctionSpace()->functionSpace() )
                return crb->model()->rBFunctionSpace()->functionSpace()->elementPtr();
            else
                return boost::shared_ptr<Vector<double>>();
        }

    std::vector<boost::shared_ptr<Vector<double>> > feSubElements( boost::shared_ptr<Vector<double>> u ) const override
        {
            DCHECK( crb ) << "DB not loaded";
            auto uFE = boost::dynamic_pointer_cast< typename crbmodel_type::space_type::element_type >( u );
            CHECK( uFE ) << "dynamic_pointer_cast fails : wrong type of element u";
            return Feel::crbplugin_details::subelements( uFE );
        }

    std::vector<boost::shared_ptr<Vector<double>>> reducedBasisFunctionsPrimal() const override
        {
            DCHECK( crb ) << "DB not loaded";
            auto const& rbPrimal = crb->model()->rBFunctionSpace()->primalRB();
            int nBasis = rbPrimal.size();
            std::vector<boost::shared_ptr<Vector<double>>> res( nBasis );
            for ( int k=0;k<nBasis;++k )
            {
                auto u = crb->model()->rBFunctionSpace()->functionSpace()->elementPtr();
                *u = rbPrimal[k];
                res[k] = u;
            }
            return res;
        }
    std::vector<boost::shared_ptr<Vector<double>>> reducedBasisFunctionsDual() const override
        {
            DCHECK( crb ) << "DB not loaded";
            auto const& rbDual = crb->model()->rBFunctionSpace()->dualRB();
            int nBasis = rbDual.size();
            std::vector<boost::shared_ptr<Vector<double>>> res( nBasis );
            for ( int k=0;k<nBasis;++k )
            {
                auto u = crb->model()->rBFunctionSpace()->functionSpace()->elementPtr();
                *u = rbDual[k];
                res[k] = u;
            }
            return res;
        }

    CRBResults run( ParameterSpaceX::Element const& mu, 
                    vectorN_type & time, double eps , int N, bool print_rb_matrix ) const override
        {
            DCHECK( crb ) << "DB not loaded";
            return crb->run( mu, time, eps, N, print_rb_matrix );
        }

    void expansion( vectorN_type const& uRB, Vector<double> & uFE,  int N ) const override
        {
            DCHECK( crb ) << "DB not loaded";
            auto uRBforExpansion = crb->model()->rBFunctionSpace()->element();
            uRBforExpansion.container() = uRB;
            crb->model()->rBFunctionSpace()->expansion( uRBforExpansion, uFE, N );
        }


    void initExporter() override
        {
            fieldExporter = exporter( _mesh=crb->model()->rBFunctionSpace()->mesh() );
        }

    void exportField( std::string const& name, CRBResults const& res ) override
        {
            auto const& solution = res.coefficients();
            auto uN = crb->model()->rBFunctionSpace()->element();
            //auto const& uN = solutions.template get<0>();
            uN.container() = solution;//.get<0>().back();
            auto sol = crb->model()->rBFunctionSpace()->expansion( uN );
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
    boost::shared_ptr<crb_type> crb;
    exporter_ptr_t fieldExporter;
};



}
#endif

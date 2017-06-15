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

    boost::shared_ptr<CRBModelBase> crbmodel() const
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

    CRBResults run( ParameterSpaceX::Element const& mu, 
                    vectorN_type & time, double eps , int N, bool print_rb_matrix ) const override
        {
            DCHECK( crb ) << "DB not loaded";
            return crb->run( mu, time, eps, N, print_rb_matrix );
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

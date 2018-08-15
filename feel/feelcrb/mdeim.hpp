/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

   This file is part of the Feel library

   Author(s): JB WAHL <wahl.jb@gmail.com>
   Date: 2018-03-27

   Copyright (C) 2018 Feel++ Consortium

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
   \file mdeim.hpp
   \author JB Wahl
   \date 2018-03-27
 */
#ifndef _FEELPP_MDEIM_HPP
#define _FEELPP_MDEIM_HPP 1

#include <feel/feelcrb/deimmodel.hpp>

namespace Feel
{

template <typename ModelType>
class MDEIM :
        public DEIMModel<ModelType,
                         typename Backend<typename ModelType::value_type>::sparse_matrix_type>
{
    typedef DEIMModel<ModelType, typename Backend<typename ModelType::value_type>::sparse_matrix_type> super_type;

public :
    typedef ModelType model_type;
    typedef std::shared_ptr<model_type> model_ptrtype;
    typedef typename super_type::parameter_type parameter_type;
    typedef typename super_type::parameterspace_ptrtype parameterspace_ptrtype;
    typedef typename super_type::sampling_ptrtype sampling_ptrtype;
    typedef typename super_type::tensor_ptrtype sparse_matrix_ptrtype;
    typedef typename super_type::element_type element_type;
    typedef typename super_type::mesh_type mesh_type;
    typedef typename super_type::space_type space_type;


    MDEIM() :
        super_type()
        {}

    MDEIM( model_ptrtype model, sampling_ptrtype sampling, std::string prefix,
           std::string const& dbfilename, std::string const& dbdirectory, int tag ) :
        super_type( model, sampling, prefix, dbfilename, dbdirectory, tag )
        {
            this->M_store_tensors = boption( prefixvm( this->M_prefix, "deim.store-matrices") );
            this->init();
        }

    ~MDEIM()
        {}

private :
    sparse_matrix_ptrtype modelAssemble( parameter_type const& mu, bool online=false ) override
        {
            if ( online )
                return this->M_online_model->assembleForMDEIM( mu, this->M_tag );
            return this->M_model->assembleForMDEIM( mu, this->M_tag );
        }

    sparse_matrix_ptrtype modelAssemble( parameter_type const& mu, element_type const& u, bool online=false ) override
        {
            if ( online )
                return this->M_online_model->assembleForMDEIMnl(mu,u,this->M_tag);
            return this->M_model->assembleForMDEIMnl(mu,u,this->M_tag);
        }

    void updateSubMesh() override
        {
            auto Xh = this->M_model->functionSpace();
            auto mesh = Xh->mesh();

            std::vector<int> my_index;
            for ( auto index : this->M_index )
            {
                my_index.push_back( index.first );
                if ( index.first!=index.second )
                    my_index.push_back( index.second );
            }
            this->updateEltsId( my_index );

            // create new submesh with the new elements and reread it in sequential
            auto submesh = createSubmesh( mesh, idelements(mesh,this->M_elts_ids.begin(), this->M_elts_ids.end()) );
            saveGMSHMesh( _mesh=submesh, _filename=this->name(true)+"-submesh.msh" );
            auto seqmesh = loadMesh( _mesh=new mesh_type,
                                     _filename=this->name(true)+"-submesh.msh",
                                     _worldcomm= Environment::worldCommSeqPtr() );
            Rh = space_type::New( seqmesh,
                                  _worldscomm=makeWorldsComm(space_type::nSpaces,Environment::worldCommSeqPtr()) );

            this->M_map->init( Rh, Xh );

            // on each proc : store the new indexR corresponding to the reduced space
            this->M_indexR.clear();
            for ( int i=0; i<this->M_index.size(); i++ )
            {
                auto i1 = this->M_index[i].first;
                auto i2 = this->M_index[i].second;

                int ir1 = this->M_map->clusterToSequential( i1 );
                int ir2 = this->M_map->clusterToSequential( i2 );
                CHECK( ir1>=0 )<< "No matching reduced index in the map for Xh index "<< i1 <<std::endl;
                CHECK( ir2>=0 )<< "No matching reduced index in the map for Xh index "<< i2 <<std::endl;

                this->M_indexR.push_back( std::make_pair(ir1, ir2) );
            }
            this->M_online_model->setFunctionSpaces( Rh );

        }

private :
    using super_type::Rh;
}; // class MDEIM




namespace detail
{
template <typename Args>
struct compute_mdeim_return
{
    typedef typename boost::remove_reference<typename boost::remove_pointer<typename parameter::binding<Args, tag::model>::type>::type>::type::element_type model1_type;
    typedef typename boost::remove_const<typename boost::remove_pointer<model1_type>::type>::type model_type;

    typedef MDEIM<model_type> type;
    typedef std::shared_ptr<type> ptrtype;
};
} // namespace detail


BOOST_PARAMETER_FUNCTION(
                         ( typename Feel::detail::compute_mdeim_return<Args>::ptrtype ), // 1. return type
                         mdeim,                        // 2. name of the function template
                         tag,                                        // 3. namespace of tag types
                         ( required
                           ( in_out(model),          * )
                           ) // required
                         ( optional
                           ( sampling, *, nullptr )
                           ( prefix, *( boost::is_convertible<mpl::_,std::string> ), "" )
                           ( filename, *( boost::is_convertible<mpl::_,std::string> ), "" )
                           ( directory, *( boost::is_convertible<mpl::_,std::string> ), "" )
                           ( tag, *( boost::is_convertible<mpl::_,int> ), 0 )
                           ) // optionnal
                         )
{
    typedef typename Feel::detail::compute_mdeim_return<Args>::type mdeim_type;
    return std::make_shared<mdeim_type>( model, sampling, prefix, filename, directory, tag );
}




} //namespace Feel


#endif

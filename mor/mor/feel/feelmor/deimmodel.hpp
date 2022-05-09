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
   \file deimmodel.hpp
   \author JB Wahl
   \date 2018-03-27
 */

#ifndef _FEELPP_DEIMMODEL_HPP
#define _FEELPP_DEIMMODEL_HPP 1

#include <feel/feelmor/deimbase.hpp>

namespace Feel
{

/**
 * Class inherited from DEIMBase
 * It contains the model
 * This class is still generic and support either matrices or vectors.
 */
template <typename ModelType, typename TensorType>
class DEIMModel :
        public DEIMBase<typename ModelType::parameterspace_type, typename ModelType::space_type,
                        TensorType>
{
    typedef DEIMModel<ModelType,TensorType> self_type;
    typedef DEIMBase<typename ModelType::parameterspace_type, typename ModelType::space_type, TensorType> super_type;

public :
    typedef ModelType model_type;
    typedef std::shared_ptr<model_type> model_ptrtype;
    typedef typename model_type::parameterspace_type parameterspace_type;
    typedef typename super_type::parameterspace_ptrtype parameterspace_ptrtype;
    typedef typename super_type::parameter_type parameter_type;
    typedef typename super_type::sampling_ptrtype sampling_ptrtype;
    typedef typename super_type::tensor_ptrtype tensor_ptrtype;
    typedef typename super_type::element_type element_type;
    typedef typename super_type::element_type element_ptrtype;
    typedef typename super_type::mesh_type mesh_type;
    typedef typename super_type::space_type space_type;
    typedef typename super_type::space_ptrtype space_ptrtype;
    typedef typename super_type::rbspace_ptrtype rbspace_ptrtype;
    typedef typename super_type::rangespace_type rangespace_type;
    typedef typename super_type::mesh_ptrtype mesh_ptrtype;

    static const bool by_block = model_type::by_block;
    static const int n_block = space_type::nSpaces;

    //! Default Constructor
    DEIMModel() : super_type()
        {}

    //! Constructor from model
    DEIMModel( model_ptrtype model, sampling_ptrtype sampling, std::string prefix,
               std::string const& dbfilename, std::string const& dbdirectory, int tag, po::variables_map vm );

    //! Destructor
    virtual ~DEIMModel() {}

    /**
     * Initialize the class: Check if the model is linear or not
     */
    void init();

    tensor_ptrtype assemble( parameter_type const& mu, bool online=false, bool force_fem=false ) override;
    tensor_ptrtype assemble( parameter_type const& mu, element_type const& u, bool online=false ) override
        {
            CHECK(this->M_nl_assembly) << "You called nl coefficient for DEIM but you implemented assembleForDEIM(mu)\n";
            return modelAssemble(mu,u,online);
        }

    void updateRb( rbspace_ptrtype const& XN, std::vector<std::vector<int>> const& subN ) override
        {
            M_subN = subN;
            return updateRb( XN, mpl::bool_<by_block>() );
        }

    //! \return the online model
    model_ptrtype & onlineModel() { return M_online_model; }
    //! \return the model
    model_ptrtype & model() { return M_model; }

    int tag() override { return M_tag; }

    int subN( int const& n_space ) const
        { return subN( n_space, WNmuSize() ); }
    int subN( int const& n_space, int const& N ) const
        {
            CHECK( M_subN.size()>0 )<<"No subN with this deim\n";
            CHECK( n_space<n_block ) <<"Invalid space number\n";
            CHECK( N>=0 )<<"Invalid size N="<<N<<std::endl;
            CHECK( N<M_subN[n_space].size() )<<"Invalid size N="<<N
                                             <<", size="<<M_subN[n_space].size()<<std::endl;
            return M_subN[n_space][N];
        }

    int dimension() const
        { return dimension( WNmuSize() ); }
    int dimension( int N ) const
        {
            CHECK( M_subN.size()>0 )<<"No subN with this deim\n";
            int dim = 0;
            for ( int n=0; n<n_block; n++ )
                dim += subN(n,N);
            return dim;
        }

    int WNmuSize() const
        {
            CHECK( M_subN.size()>0 )<<"No subN with this deim\n";
            int size = M_subN[0].size();
            for ( int i=0; i<M_subN.size(); i++ )
                CHECK( size==M_subN[i].size() ) <<"Space #"<<i<<" has different size ("
                                                <<size<<" versus "<< M_subN[i].size()<<")\n";
            return size;
        }

protected :
    //! UpdateRb without block structure
    void updateRb( rbspace_ptrtype const& XN, mpl::false_ );
    //! UpdateRb with block structure
    void updateRb( rbspace_ptrtype const& XN, mpl::true_ );

    element_type deimExpansion( vectorN_type const& urb ) override
        { return deimExpansion( urb, mpl::bool_<by_block>() ); }
    //! \return the expansion of \p urb
    element_type deimExpansion( vectorN_type const& urb, mpl::false_  );
    //! \return the expansion of \p urb using block strcture
    element_type deimExpansion( vectorN_type const& urb, mpl::true_  );

    void rebuildMap() override
        { this->M_map->init( Rh, M_model->functionSpace() ); }


    //! call the assemble function of the model with parameter \p mu
    virtual tensor_ptrtype modelAssemble( parameter_type const& mu, bool online=false )=0;
    //! call the assemble function of the model with parameter \p mu and FE solution \p u
    virtual tensor_ptrtype modelAssemble( parameter_type const& mu, element_type const& u, bool online=false )=0;

    //! Update the list of element containing the dofs in the \p index_list
    void updateEltsId( std::vector<int> const& index_list )
        {
            this->M_elts_ids.clear();
            return updateEltsId( index_list, mpl::bool_<space_type::is_composite>() );
        }
    //! Update the list of element containing the dofs in the \p index_list for non-composite spaces
    void updateEltsId( std::vector<int> const& index_list, mpl::false_ );
    //! Update the list of element containing the dofs in the \p index_list for composite spaces
    void updateEltsId( std::vector<int> const& index_list, mpl::true_ )
        {
            UdpateEltsIdComposite builder( this->M_model->functionSpace(), index_list );
            rangespace_type range;
            this->M_elts_ids = boost::fusion::fold(range, std::set<int>(), builder);
        }

    virtual space_ptrtype newInterpolationSpace( mesh_ptrtype const& mesh ) override
        {
            return space_type::New( _mesh=mesh,
                                    _range=M_online_model->functionspaceMeshSupport( mesh ) );
        }


protected :
    model_ptrtype M_model, M_online_model;
    int M_tag;
    std::vector<std::vector<int>> M_subN;
    using super_type::M_write_nl_solutions;
    using super_type::M_write_nl_directory;
    using super_type::M_rb;
    using super_type::Rh;

private :
    struct UpdateRbByBlock
    {
        UpdateRbByBlock( self_type* deim,  rbspace_ptrtype const& XN ):
            m_deim( deim ),
            m_XN( XN )
        {}

        template <typename T>
        void operator()( T const& t ) const
        {
            auto subXN = m_XN->template rbFunctionSpace<T::value>();
            auto subRh = m_deim->Rh->template functionSpace<T::value>();
            auto wn = subXN->primalRB();

            m_deim->template subRb<T::value>().resize( wn.size() );
            for( int i=0; i<wn.size(); i++ )
            {
                m_deim->template subRb<T::value>()[i] = subRh->elementPtr();
                m_deim->M_map->template project<T::value>( *(m_deim->template subRb<T::value>()[i]), *wn[i] );
            }

            m_deim->M_n_block_rb[T::value] = wn.size();
        }

    private :
        self_type* m_deim;
        rbspace_ptrtype m_XN;
    };

    struct ExpansionByBlock
    {
        ExpansionByBlock( self_type* deim, vectorN_type const& urb ) :
            m_deim( deim ),
            m_urb( urb ),
            N(0),
            m_start(0),
            U( m_deim->onlineModel()->functionSpace() )
        {
            for ( int i=0; i<m_deim->WNmuSize(); i++ )
                if ( m_deim->dimension(i)==m_urb.size() )
                    N=i;
            CHECK(N!=0) <<"Error trying to split urb\n";
        }

        template <typename T>
        void operator()( T const& t ) const
        {
            auto WN = m_deim->template subRb<T::value>();
            int Nwn = m_deim->subN(T::value,N);
            CHECK( m_start+Nwn<=m_urb.size() ) << "invalide expansion size, N<n="<<Nwn<<", size="
                                               <<m_urb.size() << ", space="<<T::value<<std::endl;

            auto  coeff = m_urb.segment( m_start, Nwn );
            auto u = U.template element<T::value>();
            u = Feel::expansion( WN, coeff, Nwn ).container();
            m_start += Nwn;
        }

        element_type& field() { return U; }

    private :
        self_type* m_deim;
        vectorN_type m_urb;
        int N;
        mutable int m_start;
        mutable element_type U;
    };

    struct UdpateEltsIdComposite
    {
        using result_type = std::set<int>;

        UdpateEltsIdComposite( space_ptrtype const& Xh, std::vector<int> const& index ) :
            m_Xh( Xh ),
            m_index( index )
        {}

        template <typename T>
        result_type operator()( result_type const& s, T const& t ) const
        {
            auto subXh = m_Xh->template functionSpace<T::value>();
            auto mesh = m_Xh->mesh();
            std::set<int> r1;
            std::set<int> r2;
            for ( auto index : m_index )
            {
                auto searchGpDof = m_Xh->dof()->searchGlobalProcessDof( index );
                if (  boost::get<0>( searchGpDof ) )
                {
                    size_type gpdof = boost::get<1>( searchGpDof );
                    int space = m_Xh->dof()->databaseIndexFromContainerId( gpdof );
                    if ( space==T::value )
                    {
                        gpdof = m_Xh->dof()->containerIdToDofId( space, gpdof );
                        CHECK( gpdof!=invalid_v<size_type> ) <<"Dof not found\n";
                        for ( auto const& dof : subXh->dof()->globalDof( gpdof ) )
                        {
                            size_type eltId = dof.second.elementId();
                            if ( mesh->element( eltId ).isGhostCell() )
                                continue;
                            r1.insert( eltId );
                        }
                    }
                }
            }
            std::set_union(std::begin(s), std::end(s),
                           std::begin(r1), std::end(r1),
                           std::inserter(r2, std::begin(r2)));
            return r2;
        }

    private:
        space_ptrtype m_Xh;
        std::vector<int> m_index;
    };


}; // class DEIMModel


template <typename ModelType, typename TensorType>
DEIMModel<ModelType,TensorType>::DEIMModel( model_ptrtype model, sampling_ptrtype sampling, std::string prefix, std::string const& dbfilename, std::string const& dbdirectory, int tag, po::variables_map vm ) :
    super_type( model->functionSpace(), model->parameterSpace(),
                sampling, model->uuid(), model->modelName(),
                prefix, dbfilename, dbdirectory, Environment::worldCommPtr(), model->prefix(), vm ),
    M_model( model ),
    M_tag( tag )
{
    this->M_online_model = model_ptrtype( new model_type( model->modelName() ) );
//    this->M_online_model = std::make_shared<model_type>( model->modelName(), model->uuid(), model->worldCommPtr(), model->prefix());
    //this->M_online_model->setModelOnlineDeim( prefixvm(this->M_prefix,"deim-online"), this->M_vm );
    this->M_online_model->setOnlineModel();

    if ( !this->M_rebuild )
    {
        if ( this->loadDB() )
            cout<< this->name() + " : Database loaded with " << this->M_M << " basis functions\n";
        else
        {
            cout << this->name() + " : No Database loaded : start greedy algorithm from beginning\n";
            this->M_rebuild=true;
        }
    }
    else
        cout << this->name() + " : option deim.rebuild-database=true : start greedy algorithm from beginning\n";
    if ( Rh )
        M_online_model->setFunctionSpaces( Rh );
    this->M_online_model->initOnlineModel(M_model);
} //DEIMModel


template <typename ModelType, typename TensorType>
void
DEIMModel<ModelType,TensorType>::init()
{
    if ( this->M_rebuild )
    {
        auto mu = this->M_trainset->max().template get<0>();
        auto T = this->assemble(mu);
        if (!T)
        {
            this->M_nl_assembly=true;
            auto u = M_model->functionSpace()->element();
            auto Tnl = this->assemble(mu,u);
            CHECK( Tnl ) << "You want to use DEIM but you did not implement assmbleForDEIM functions\n";
        }
    }
} // init


template <typename ModelType, typename TensorType>
typename DEIMModel<ModelType,TensorType>::tensor_ptrtype
DEIMModel<ModelType,TensorType>::assemble( parameter_type const& mu, bool online, bool force_fem )
{
    if ( this->M_nl_assembly )
    {
        if ( online )
            Feel::cout << this->name() + " : WARNING : Call of online nl assembly with no solution u\n";
        //CHECK(!online) << "Call of online nl assembly with no solution u\n";

        auto u = M_model->functionSpace()->element();
        bool need_solve = true;

        if ( this->M_ser_use_rb && this->M_crb && !online && !force_fem )
        {
            std::vector<vectorN_type> uN, uNdu, uNold, uNduold;
            auto o = this->M_crb->lb( this->M_crb->dimension(), mu, uN, uNdu , uNold, uNduold );
            int size = uN.size();

            u = this->M_crb->expansion( uN[size-1], this->M_crb->WNmuSize(), false );
        }
        else
        {
            if ( this->M_ser_use_rb && !this->M_crb )
                Feel::cout <<this->name() + " WARNING : Suppose to use crb expansion with no crb class ! u will be computed using model->solve\n";

            if ( M_write_nl_solutions )
            {
                need_solve = !u.load( _path=M_write_nl_directory,
                                      _suffix=std::to_string(mu.key()), _type="hdf5" );
                if ( need_solve )
                    LOG(INFO) << this->name() + " : Unable to load nl solution in directotry "
                              << M_write_nl_directory << ", for parameter : " << mu.toString()
                              <<" / " << mu.key()<< ". Solve function will be called.";
                else
                    LOG(INFO) << this->name() + " : NL solution loaded in directotry "
                              << M_write_nl_directory << ", for parameter : " << mu.toString()
                              <<" / " << mu.key();
            }

            if ( need_solve )
            {
                LOG(INFO) << this->name() + " : calling solve function for parameter " << mu.toString()
                          <<" / " << mu.key();
                auto o = M_model->safeSolve(mu);
                u = o.first;
                this->M_last_solve_is_ok = o.second;

                if ( M_write_nl_solutions && o.second )
                {
                    LOG(INFO) << this->name() + " : Writing solution on disk in directory "
                              << M_write_nl_directory << ", for parameter : " << mu.toString()
                              <<" / " << mu.key();
                    u.save( _path=M_write_nl_directory,
                            _suffix=std::to_string(mu.key()), _type="hdf5" );
                }
            }

        }

        return modelAssemble(mu,u);
    }
    return modelAssemble(mu,online);
} // assemble


template <typename ModelType, typename TensorType>
void
DEIMModel<ModelType,TensorType>::updateRb( rbspace_ptrtype const& XN, mpl::false_ )
{
    auto wn = XN->primalRB();
    M_rb.resize( wn.size() );

    for ( int i=0; i<wn.size(); i++ )
    {
        M_rb[i] = Rh->elementPtr();
        this->M_map->project( *M_rb[i], *wn[i] );
    }
    this->M_n_rb = M_rb.size();
} // updateRb no block

template <typename ModelType, typename TensorType>
void
DEIMModel<ModelType,TensorType>::updateRb( rbspace_ptrtype const& XN, mpl::true_ )
{
    BOOST_STATIC_ASSERT( space_type::is_composite );//<< "Use of CBR block but space is not composite\n";
    if ( !this->M_n_block_rb.size() )
        this->M_n_block_rb.resize( space_type::nSpaces );

    rangespace_type range;
    auto builder = UpdateRbByBlock( this, XN );
    boost::fusion::for_each( range, builder );
} // updateRb with block


template <typename ModelType, typename TensorType>
typename DEIMModel<ModelType,TensorType>::element_type
DEIMModel<ModelType,TensorType>::deimExpansion( vectorN_type const& urb, mpl::false_  )
{
    int N = urb.size();
    FEELPP_ASSERT( N <= M_rb.size() )( N )( M_rb.size() ).error( "invalid expansion size ( N and M_rb ) ");

    return Feel::expansion( M_rb, urb, N );
} // deimExpansion no block

template <typename ModelType, typename TensorType>
typename DEIMModel<ModelType,TensorType>::element_type
DEIMModel<ModelType,TensorType>::deimExpansion( vectorN_type const& urb, mpl::true_  )
{
    BOOST_STATIC_ASSERT( space_type::is_composite );//<< "Use of CBR block but space is not composite\n";

    rangespace_type range;
    auto builder = ExpansionByBlock( this, urb );
    boost::fusion::for_each( range, builder );
    return builder.field();
} // deimExpansion with block


template <typename ModelType, typename TensorType>
void
DEIMModel<ModelType,TensorType>::updateEltsId( std::vector<int> const& index_list, mpl::false_ )
{
    auto Xh = this->M_model->functionSpace();
    auto mesh = Xh->mesh();

    for ( auto index : index_list )
    {
        auto searchGpDof = Xh->dof()->searchGlobalProcessDof( index );
        if ( boost::get<0>( searchGpDof ) )
        {
            size_type gpdof = boost::get<1>( searchGpDof );
            for ( auto const& dof : Xh->dof()->globalDof( gpdof ) )
            {
                size_type eltId = dof.second.elementId();
                if ( mesh->element( eltId ).isGhostCell() )
                    continue;
                this->M_elts_ids.insert( eltId );
            }
        }
    }
}



} // namespace Feel

#endif

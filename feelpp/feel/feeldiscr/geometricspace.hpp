/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
   Date: 2016-09-09

   Copyright (C) 2016 Feel++ Consortium

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 3.0 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file geometricspace.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2016-09-09
*/
#ifndef FEELPP_GEOMETRIC_SPACE_HPP
#define FEELPP_GEOMETRIC_SPACE_HPP 1

namespace Feel
{

struct GeometricSpaceBase : public CommObject
{
    using super = CommObject;
    GeometricSpaceBase() : super( Environment::worldCommPtr() ) {}
    explicit GeometricSpaceBase( worldcomm_ptr_t const& w ) : super( w ) {}
    GeometricSpaceBase( GeometricSpaceBase const& ) = default;
    GeometricSpaceBase( GeometricSpaceBase && ) = default;
    GeometricSpaceBase& operator=( GeometricSpaceBase const& ) = default;
    GeometricSpaceBase& operator=( GeometricSpaceBase && ) = default;
    ~GeometricSpaceBase() override = default;

};
struct ContextGeometricBase {};

template <typename MeshType>
class GeometricSpace :
        public GeometricSpaceBase,
        public std::enable_shared_from_this<GeometricSpace<MeshType> >
{
    typedef GeometricSpace<MeshType> self_type;
public :

    using super = GeometricSpaceBase;
    typedef self_type geometricspace_type;
    typedef std::shared_ptr<geometricspace_type> geometricspace_ptrtype;

    typedef MeshType mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    typedef typename mesh_type::element_type element_type;
    typedef typename element_type::gm_type gm_type;
    static const size_type gmc_context_v = vm::POINT|vm::DYNAMIC;
    typedef typename gm_type::template Context<element_type> gmc_type;
    typedef std::shared_ptr<gmc_type> gmc_ptrtype;


    explicit GeometricSpace( worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() )
        :
        super( worldComm )
        {}

    explicit GeometricSpace( mesh_ptrtype const& mesh )
        :
        super( mesh->worldCommPtr() ),
        M_mesh( mesh )
        {}

    mesh_ptrtype const& mesh() const { return M_mesh; }

    void setMesh( mesh_ptrtype const& mesh )
        {
            M_mesh = mesh;
        }

    class ContextGeometric : public ContextGeometricBase
    {
    public :
        using gmc_type = typename GeometricSpace<MeshType>::gmc_type;
        typedef std::shared_ptr<gmc_type> geometric_mapping_context_ptrtype;

        ContextGeometric() = default;

        ContextGeometric( geometricspace_ptrtype const& Xh )
            :
            M_Xh( Xh )
            {}
        ContextGeometric( geometricspace_ptrtype const& Xh, geometric_mapping_context_ptrtype const& gmc )
            :
            M_Xh( Xh ),
            M_gmc( gmc )
            {}

        geometric_mapping_context_ptrtype const& gmContext() const
        {
            return M_gmc;
        }

        void setMeshGeoContext( mesh_ptrtype const& meshGeoContext )
        {
            M_meshGeoContext = meshGeoContext;
        }

        template <size_type CTX>
        void updateGmcContext( size_type dynctx = 0 )
            {
                if ( M_gmc )
                    M_gmc->template updateContext<gmc_context_v|CTX>( dynctx );
            }

    private :
        friend class boost::serialization::access;

        template<class Archive>
        void save( Archive & ar, const unsigned int version ) const
            {
                CHECK( M_gmc ) << "no gmc defined";
                auto const& meshEltCtx = M_gmc->element();
                ar & BOOST_SERIALIZATION_NVP( meshEltCtx );
                ar & boost::serialization::make_nvp( "gmContext", *M_gmc );
            }
        template<class Archive>
        void load( Archive & ar, const unsigned int version )
            {
                element_type meshEltCtx;
                ar & BOOST_SERIALIZATION_NVP( meshEltCtx );

                if ( M_meshGeoContext )
                {
                    // std::cout << "geospace ContextGeo with minimal mesh\n";
                    if ( !M_meshGeoContext->hasElement( meshEltCtx.id() ) )
                    {
                        M_meshGeoContext->addElement( meshEltCtx, false );
                    }
                    auto const& meshEltCtxRegister = M_meshGeoContext->element( meshEltCtx.id() );
                    M_gmc = M_meshGeoContext->gm()->template context<gmc_context_v>( meshEltCtxRegister, typename gmc_type::precompute_ptrtype{} );
                    //M_gmc.reset( new gmc_type( M_meshGeoContext->gm(),meshEltCtxRegister ) );
                    ar & boost::serialization::make_nvp( "gmContext", *M_gmc );
                }
                else if ( M_Xh && M_Xh->mesh() )
                {
                    // std::cout << "geospace ContextGeo with full mesh\n";
                    CHECK ( M_Xh->mesh()->hasElement( meshEltCtx.id()/*, meshEltCtx.processId()*/ ) ) << "fails because mesh doesnt have the element reloaded for gmc";
                    auto const& meshEltCtxRegister = M_Xh->mesh()->element( meshEltCtx.id()/*, meshEltCtx.processId()*/ );
                    //M_gmc.reset( new gmc_type( M_Xh->mesh()->gm(),meshEltCtxRegister ) );
                    M_gmc = M_Xh->mesh()->gm()->template context<gmc_context_v>( meshEltCtxRegister, typename gmc_type::precompute_ptrtype{} );
                    ar & boost::serialization::make_nvp( "gmContext", *M_gmc );

                }
            }
        BOOST_SERIALIZATION_SPLIT_MEMBER()

    private :
        geometricspace_ptrtype M_Xh;
        geometric_mapping_context_ptrtype M_gmc;
        mesh_ptrtype M_meshGeoContext;
    };

    class Context
        :
        public std::map<int,std::shared_ptr<ContextGeometric>>
    {
    public:
        typedef std::map<int,std::shared_ptr<ContextGeometric>> super_type;
        typedef typename super_type::iterator iterator;

        typedef geometricspace_type functionspace_type;

        typedef typename mesh_type::node_type node_type;
        typedef typename matrix_node<typename node_type::value_type>::type matrix_node_type;


        Context( geometricspace_ptrtype const& Xh )
            :
            M_Xh( Xh ),
            M_ctxHaveBeenMpiBroadcasted( false )
            {}

        geometricspace_ptrtype functionSpace() const
            {
                return M_Xh;
            }
        geometricspace_type* ptrFunctionSpace() const
            {
                return M_Xh.get();
            }

        int nPoints() const
        {
            return M_t.size();
        }

        bool ctxHaveBeenMpiBroadcasted() const { return M_ctxHaveBeenMpiBroadcasted; }

        void removeCtx()
        {
            this->clear();
            M_t.clear();
            M_t_proc.clear();
            M_meshGeoContext.reset();
        }


        std::pair<iterator, bool>
        add( node_type const& t )
        {
            int ptIdInCtx = M_t.size();
            M_t.push_back( t );

            std::pair<iterator, bool> ret = std::make_pair(this->end(),false);

            //rank of the current processor
            rank_type currentProcRank = M_Xh->mesh()->worldComm().globalRank();
            //total number of processors
            rank_type nprocs = M_Xh->mesh()->worldComm().globalSize();

            // localise t in space, find geometrical element in which t belongs
            matrix_node_type m( mesh_type::nRealDim, 1 );
            for(int i = 0; i < mesh_type::nRealDim; ++i )
                m(i,0) = t(i);
            auto loc =  M_Xh->mesh()->tool_localization();
            loc->setExtrapolation( false );
            auto analysis = loc->run_analysis( m, invalid_v<typename mesh_type::size_type> );
            auto found_points = analysis.template get<0>();
            bool found = found_points[0];

            std::vector<int> found_pt( nprocs, 0 );
            if( found )
                found_pt[currentProcRank]=1;
            std::vector<int> global_found_pt( nprocs, 0 );
            if ( nprocs > 1 )
                mpi::all_reduce( M_Xh->mesh()->comm(), found_pt.data(), found_pt.size(), global_found_pt.data(), std::plus<int>() );
            else
                global_found_pt[ 0 ] = found_pt[ 0 ];
            // only one proc has the point
            bool findOneProcess = false;
            for ( rank_type p=0;p<nprocs;++p )
            {
                if ( findOneProcess )
                    global_found_pt[p] = 0;
                if ( global_found_pt[p] == 1 )
                    findOneProcess=true;
            }
            if ( global_found_pt[ currentProcRank ] == 1 ) //we are on the proc that have the searched point
            {
                auto it = loc->result_analysis_begin();
                auto en = loc->result_analysis_end();
                DCHECK( boost::next(it) == en ) << "Logic problem in finding one point in the mesh\n";
                auto eid = it->first;
                auto xref = boost::get<1>( *(it->second.begin()) );
                DVLOG(2) << "found point " << t << " in element " << eid << " on proc "<<currentProcRank
                         << "  - reference coordinates " << xref;

                M_eltToUpdate[eid].push_back( std::make_tuple(ptIdInCtx, xref ) );
            }
#if 0
            if( found ) //we are on the proc that have the searched point
            {

                found_pt[currentProcRank]=1;

                auto it = loc->result_analysis_begin();
                auto en = loc->result_analysis_end();
                DCHECK( boost::next(it) == en ) << "Logic problem in finding one point in the mesh\n";
                auto eid = it->first;
                auto xref = boost::get<1>( *(it->second.begin()) );
                DVLOG(2) << "found point " << t << " in element " << eid << " on proc "<<currentProcRank
                         << "  - reference coordinates " << xref;

                typename mesh_type::gm_type::precompute_type::matrix_node_t_type p(mesh_type::nDim,1);

                ublas::column( p, 0 ) = xref;
                // compute for each basis function in reference element its
                // value at \hat{t} in reference element
                auto gmpc = M_Xh->mesh()->gm()->preCompute( M_Xh->mesh()->gm(), p );
                DVLOG(2) << "build precompute data structure for geometric mapping\n";
                // build geometric mapping
                auto gmc = M_Xh->mesh()->gm()->template context<gmc_context_v>( M_Xh->mesh()->element( eid ), gmpc );
                DVLOG(2) << "build geometric mapping context\n";

                int number = M_t.size()-1;
                ret = this->insert( std::make_pair( number , std::make_shared<ContextGeometric>( M_Xh, gmc ) ) );

                if ( nprocs > 1 )
                    mpi::all_reduce( M_Xh->mesh()->comm(), found_pt.data(), found_pt.size(), global_found_pt.data(), std::plus<int>() );
                else
                    global_found_pt[ 0 ] = found_pt[ 0 ];

            }//if( found )
            else
            {
                if ( nprocs > 1 )
                    mpi::all_reduce( M_Xh->mesh()->comm(), found_pt.data(), found_pt.size(), global_found_pt.data(), std::plus<int>() );

            }//not found case
#endif
            //verify that the point is on a proc
            bool found_on_a_proc = false;
            for (int i = 0 ; i < global_found_pt.size(); ++i )
            {
                if ( global_found_pt[i] != 0 )
                {
                    DVLOG(2) << "processor " << i << " has the point " << t << "\n";
                    found_on_a_proc = true;
                    M_t_proc.push_back(i);
                    break;
                }
            }
            CHECK( found_on_a_proc ) << "the point " << t << " was not found ! \n";
            
            //this->updateForUse();// TODO!!!!
            
            // update geo context for each process + define context mesh
            if ( false )
            {
                int mynumber = M_t.size()-1;
                this->syncCtx( mynumber );
            }

            return ret;
        }


        void updateForUse()
        {
            //M_eltToUpdate[eid].push_back( std::make_tuple(ptIdInCtx, xref ) );
            for ( auto const& [eid,data] : M_eltToUpdate )
            {
                typename mesh_type::gm_type::precompute_type::matrix_node_t_type xrefs(mesh_type::nDim,data.size());
                int k=0;
                std::vector<index_type> ptIds(data.size());
                for ( auto const& [ptIdInCtx,xref] : data )
                {
                    ublas::column( xrefs, k ) = xref;
                    ptIds[k] = ptIdInCtx;
                    ++k;
                }
                // compute for each basis function in reference element its
                // value at \hat{t} in reference element
                auto gmpc = M_Xh->mesh()->gm()->preCompute( M_Xh->mesh()->gm(), xrefs );
                DVLOG(2) << "build precompute data structure for geometric mapping\n";
                // build geometric mapping
                auto gmc = M_Xh->mesh()->gm()->template context<gmc_context_v>( M_Xh->mesh()->element( eid ), gmpc );
                DVLOG(2) << "build geometric mapping context\n";

                int number = this->size();//M_t.size()-1;
                /*ret =*/ this->insert( std::make_pair( number , std::make_shared<ContextGeometric>( M_Xh, gmc ) ) );

                M_ctxIdToPointIds.emplace( number, std::move(ptIds) );
            }

            M_eltToUpdate.clear();

        }


        // TODO : remove this and move into inherits
        std::map<int,std::vector<index_type>> const& ctxIdToPointIds() const { return M_ctxIdToPointIds; }

        template <size_type CTX>
        void updateGmcContext( size_type dynctx = 0 )
        {
            for ( auto& [ptId,geoCtx] : *this )
            {
                if ( geoCtx )
                    geoCtx->template updateGmcContext<CTX>( dynctx );
            }
        }

    private :
        void syncCtx( int ptId )
        {
            rank_type procId = M_t_proc[ptId];
            rank_type myrank = M_Xh->worldComm().rank();

            if ( !M_meshGeoContext )
                M_meshGeoContext = std::make_shared<mesh_type>( M_Xh->worldCommPtr() );

            std::shared_ptr<ContextGeometric> geoCtxReload;
            if ( myrank == procId )
            {
                // update additional mesh used in rb context
                CHECK( this->find( ptId ) != this->end() ) << "point id is not saved on this process";
                geoCtxReload = this->operator[]( ptId );
                CHECK( geoCtxReload ) << "geoCtxReload not init";
                geoCtxReload->setMeshGeoContext( M_meshGeoContext );
                auto const& modelMeshEltCtx = geoCtxReload->gmContext()->element();
                if ( !M_meshGeoContext->hasElement( modelMeshEltCtx.id(), modelMeshEltCtx.processId() ) )
                {
                    element_type meshEltCtx = modelMeshEltCtx;
                    M_meshGeoContext->addElement( meshEltCtx, false );
                }
            }
            else
            {
                // create a new geo context
                geoCtxReload.reset( new ContextGeometric( M_Xh ) );
                geoCtxReload->setMeshGeoContext( M_meshGeoContext );
            }
            // recv geo context from process which have localized the point
            mpi::broadcast( M_Xh->worldComm().globalComm(), *geoCtxReload, procId );

            // save geo context with process which doesnt have localized the point
            if ( myrank != procId )
                this->operator[]( ptId ) = geoCtxReload;

            //std::cout << "["<<M_Xh->worldComm().rank()<<"] M_meshGeoContext->numElements() : " << M_meshGeoContext->numElements() << "\n";

        }

        friend class boost::serialization::access;

        template<class Archive>
        void save( Archive & ar, const unsigned int version ) const
            {
                ar & BOOST_SERIALIZATION_NVP( M_t );
                ar & BOOST_SERIALIZATION_NVP( M_t_proc );

                std::vector<int> geoCtxKeys;
                for ( auto const& geoCtxPair : *this )
                    geoCtxKeys.push_back( geoCtxPair.first );
                ar & BOOST_SERIALIZATION_NVP( geoCtxKeys );

                for ( int geoCtxKey : geoCtxKeys )
                {
                    std::string geoctxNameInSerialization = (boost::format("geoSpaceContext_%1%")%geoCtxKey).str();
                    // std::cout << "geospace::Context save name : " << geoctxNameInSerialization << "\n";
                    ar & boost::serialization::make_nvp( geoctxNameInSerialization.c_str(), *(this->find( geoCtxKey )->second) );
                }
            }
        template<class Archive>
        void load( Archive & ar, const unsigned int version )
            {
                ar & BOOST_SERIALIZATION_NVP( M_t );
                ar & BOOST_SERIALIZATION_NVP( M_t_proc );

                std::vector<int> geoCtxKeys;
                ar & BOOST_SERIALIZATION_NVP( geoCtxKeys );

                if ( !M_meshGeoContext )
                    M_meshGeoContext = std::make_shared<mesh_type>( M_Xh->worldCommPtr() );

                for ( int geoCtxKey : geoCtxKeys )
                {
                    std::shared_ptr<ContextGeometric> geoCtxReload( new ContextGeometric( M_Xh ) );
                    geoCtxReload->setMeshGeoContext( M_meshGeoContext );
                    std::string geoctxNameInSerialization = (boost::format("geoSpaceContext_%1%")%geoCtxKey).str();
                    // std::cout << "geospace::Context load name : " << geoctxNameInSerialization << "\n";
                    ar & boost::serialization::make_nvp( geoctxNameInSerialization.c_str(), *geoCtxReload );
                    /*ret =*/ this->insert( std::make_pair( geoCtxKey , geoCtxReload ) );
                }
            }
        BOOST_SERIALIZATION_SPLIT_MEMBER()

    private :
        geometricspace_ptrtype M_Xh;
        std::vector<node_type> M_t;
        std::vector<int> M_t_proc;//point number i is on proc M_t_proc[i]
        mesh_ptrtype M_meshGeoContext;
        bool M_ctxHaveBeenMpiBroadcasted;

        //! internal container used by updateForUse
        std::map<size_type, std::vector<std::tuple<size_type, node_type> > > M_eltToUpdate;

        std::map<int,std::vector<index_type>> M_ctxIdToPointIds;

    };

    Context context() /*const*/ { return Context( this->shared_from_this() ); }

    std::shared_ptr<ContextGeometric> contextBasis( std::pair<int, std::shared_ptr<ContextGeometric>> const& p, Context const& c ) const { return p.second; }

private :
    mesh_ptrtype M_mesh;

};

#if 0
template <typename MeshType, typename ... CTX>
auto selectGeomapContext( std::shared_ptr<MeshType> const& mesh, CTX const& ... ctx )
{
    using geoelement_type = typename MeshType::element_type;
    using gmc_type = typename geoelement_type::gm_type::template Context<geoelement_type>;
    std::shared_ptr<gmc_type> res;
    hana::for_each( hana::make_tuple( ctx... ), [&mesh,&res]( auto const& e )
                    {
                        if constexpr ( std::is_same_v<std::decay_t<decltype(e->gmContext())>, std::shared_ptr<gmc_type> > )
                        {
                            if ( e->gmContext()->element().mesh()->isSameMesh( mesh ) )
                                res = e->gmContext();
                        }
                    } );
    return res;
}
#endif

} //namespace Feel

#endif // FEELPP_GEOMETRIC_SPACE_HPP

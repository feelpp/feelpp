/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

#ifndef FEELPP_FEELDISCR_DETAIL_FUNCTIONSPACELEGACYCOMPOSITE_HPP
#define FEELPP_FEELDISCR_DETAIL_FUNCTIONSPACELEGACYCOMPOSITE_HPP 1

// Internal FunctionSpace legacy-composite support fragment.
// This file is included from functionspace.hpp while namespace Feel::detail is open.

template<bool IsComposite>
struct LegacyCompositeFunctionSpacePolicy
{
    static constexpr bool legacy_composite_enabled = FEELPP_ENABLE_LEGACY_COMPOSITE_FUNCTIONSPACE;
    static constexpr bool is_legacy_composite = IsComposite;
    static constexpr bool uses_internal_composite = IsComposite;
    static constexpr bool is_product_backed_composite = false;
};

template<typename SpaceType, bool IsComposite>
class LegacyCompositeFunctionSpaceStorage
{};

template<typename SpaceType>
class LegacyCompositeFunctionSpaceStorage<SpaceType, true>
{
public:
    using functionspace_vector_type = typename SpaceType::functionspace_vector_type;

    functionspace_vector_type const& functionSpaces() const noexcept
    {
        return M_functionspaces;
    }

    functionspace_vector_type& functionSpaces() noexcept
    {
        return M_functionspaces;
    }

    void setFunctionSpaces( functionspace_vector_type const& functionspaces )
    {
        M_functionspaces = functionspaces;
    }

    template<typename ... FSpaceList>
    void setFunctionSpacesFromList( FSpaceList const&... fspacelist )
    {
        M_functionspaces = boost::fusion::make_vector( fspacelist... );
    }

    template<int I>
    typename mpl::at_c<functionspace_vector_type,I>::type
    functionSpace()
    {
        return boost::fusion::at_c<I>( M_functionspaces );
    }

    template<int I>
    typename mpl::at_c<functionspace_vector_type,I>::type const&
    functionSpace() const
    {
        return boost::fusion::at_c<I>( M_functionspaces );
    }

private:
    functionspace_vector_type M_functionspaces;
};

template<typename SpaceType>
struct InitializeSpace
{
    typedef typename SpaceType::functionspace_vector_type functionspace_vector_type;
    typedef typename SpaceType::mesh_ptrtype MeshPtrType;
    typedef typename SpaceType::mesh_support_vector_type mesh_support_vector_type;
    using globaldof_type = Dof<typename SpaceType::mesh_type::size_type>;
    InitializeSpace( functionspace_vector_type & functionspaces,
                     MeshPtrType const& mesh,
                     mesh_support_vector_type const& meshSupport,
                     std::vector<globaldof_type> const& dofindices,
                     worldscomm_ptr_t const & worldsComm,
                     std::vector<DofTableExtendedType> extendedDofTable )
        :
        M_functionspaces( functionspaces ),
        M_cursor( 0 ),
        M_worldsComm( worldsComm ),
        M_mesh( mesh ),
        M_meshSupport( meshSupport ),
        M_dofindices( dofindices ),
        M_extendedDofTable( extendedDofTable )
    {}
    template <typename T>
    void operator()( T const& t ) const
        {
            if constexpr ( is_shared_ptr<MeshPtrType>() )
            {
                typedef typename fusion::result_of::at_c<functionspace_vector_type,T::value>::type _subspace_ptrtype;
                typedef typename boost::remove_reference<_subspace_ptrtype>::type subspace_ptrtype;
                typedef typename subspace_ptrtype::element_type subspace_type;

                auto & subSpace = boost::fusion::at_c<T::value>( M_functionspaces );
                auto subMeshSupport = typename subspace_type::mesh_support_vector_type( boost::fusion::at_c<T::value>( M_meshSupport ) );
                subSpace = subspace_ptrtype( new subspace_type( M_mesh, subMeshSupport, M_dofindices,
                                                                makeWorldsComm( 1,M_worldsComm[M_cursor] ),
                                                                std::vector<DofTableExtendedType>( 1,M_extendedDofTable[M_cursor] ) ) );
                FEELPP_ASSERT( subSpace ).error( "invalid function space" );

                ++M_cursor;// warning M_cursor < nb color
            }
            else
            {
                typedef typename fusion::result_of::at_c<functionspace_vector_type,T::value>::type _subspace_ptrtype;
                typedef typename boost::remove_reference<_subspace_ptrtype>::type subspace_ptrtype;
                typedef typename subspace_ptrtype::element_type subspace_type;

                auto & subSpace = boost::fusion::at_c<T::value>( M_functionspaces );
                // look for T::mesh_ptrtype in MeshPtrType
                //auto m = *fusion::find<typename subspace_type::mesh_ptrtype>(M_mesh);
                auto m = boost::fusion::at_c<T::value>( M_mesh );
                auto subMeshSupport = typename subspace_type::mesh_support_vector_type( boost::fusion::at_c<T::value>( M_meshSupport ) );
                subSpace = subspace_ptrtype( new subspace_type( m, subMeshSupport, M_dofindices,
                                                                makeWorldsComm( 1,M_worldsComm[M_cursor] ),
                                                                std::vector<DofTableExtendedType>( 1,M_extendedDofTable[M_cursor] ) ) );
                FEELPP_ASSERT( subSpace ).error( "invalid function space" );

                ++M_cursor;// warning M_cursor < nb color
            }
        }
    functionspace_vector_type & M_functionspaces;
    mutable uint16_type M_cursor;
    worldscomm_ptr_t M_worldsComm;
    MeshPtrType M_mesh;
    mesh_support_vector_type const& M_meshSupport;
    std::vector<globaldof_type> const& M_dofindices;
    std::vector<DofTableExtendedType> M_extendedDofTable;
};
template<typename DofType>
struct updateDataMapProcessStandard
{
    typedef std::shared_ptr<DofType> result_type;

    updateDataMapProcessStandard( worldcomm_ptr_t const& worldComm,
                                  uint16_type nSpaces )
        :
        M_worldComm( worldComm ),
        M_cursor( 0 ),
        M_lastCursor( nSpaces-1 )
    {}

    template <typename T>
    result_type operator()( result_type const& r, std::shared_ptr<T> & x ) const
    {
        M_subdm.push_back( x->mapPtr() );
        if ( M_cursor == M_lastCursor )
        {
            result_type dm = std::make_shared<DofType>( M_subdm,M_worldComm );
            return dm;
        }
        ++M_cursor;
        return r;
    }

    worldcomm_ptr_t M_worldComm;
    mutable uint16_type M_cursor;
    uint16_type M_lastCursor;
    mutable std::vector<datamap_ptrtype<>> M_subdm;
};





struct NbDof
{
    typedef size_type result_type;
    NbDof( size_type start = 0, size_type size = invalid_v<size_type> )
        :
        M_cursor( start ),
        M_finish( size )
    {}
    template<typename Sig>
    struct result;

    template<typename T, typename S>
#if BOOST_VERSION < 104200
    struct result<NbDof( T,S )>
#else
    struct result<NbDof( S,T )>
#endif
:
    boost::remove_reference<S>
    {};
    template <typename T>
    size_type
    operator()( T const& x, size_type s ) const
    {
        size_type ret = s;

        if ( !x )
            return ret;

        if ( M_cursor < M_finish )
            ret += x->nDof();

        ++M_cursor;
        return ret;
    }

    template <typename T>
    size_type
    operator()( size_type s, T const& x ) const
    {
        return this->operator()( x, s );
    }
private:
    mutable size_type M_cursor;
    size_type M_finish;
};

#if 0
struct NLocalDof
{
    NLocalDof( size_type start = 0, size_type size = invalid_v<size_type> )
        :
        M_cursor( start ),
        M_finish( size )
    {}
    template<typename Sig>
    struct result;

    template<typename T, typename S>
#if BOOST_VERSION < 104200
    struct result<NLocalDof( T,S )>
#else
    struct result<NLocalDof( S,T )>
#endif
:
    boost::remove_reference<S>
    {};
    template <typename T>
    size_type
    operator()( T const& x, size_type s ) const
    {
        size_type ret = s;

        if ( M_cursor < M_finish )
            ret += x->nLocalDof();

        ++M_cursor;
        return ret;
    }
    template <typename T>
    size_type
    operator()( size_type s, T const& x ) const
    {
        return this->operator()( x, s );
    }
private:
    mutable size_type M_cursor;
    size_type M_finish;
};
#else // MPI
template< typename IsWithGhostType>
struct NLocalDof
{

    NLocalDof( worldscomm_ptr_t const & worldsComm = Environment::worldsComm(1),
               bool useOffSubSpace = false,
               size_type start = 0, size_type size = invalid_v<size_type> )
        :
        M_cursor( start ),
        M_finish( size ),
        M_worldsComm( worldsComm ),
        M_useOffSubSpace( useOffSubSpace )
    {}
    template<typename Sig>
    struct result;

    template<typename T, typename S>
#if BOOST_VERSION < 104200
    struct result<NLocalDof( T,S )>
#else
    struct result<NLocalDof( S,T )>
#endif
:
    boost::remove_reference<S>
    {};

    template <typename T>
    size_type
    nLocalDof( T const& x, mpl::bool_<true> /**/ ) const
    {
        return x->nLocalDofWithGhost();
    }

    template <typename T>
    size_type
    nLocalDof( T const& x, mpl::bool_<false> /**/ ) const
    {
        return x->nLocalDofWithoutGhost();
    }

    template <typename T>
    size_type
    operator()( T const& x, size_type s ) const
    {
        size_type ret = s;

        if ( M_cursor < M_finish )
        {
            if ( M_useOffSubSpace )
            {
                ret += nLocalDof( x, mpl::bool_<IsWithGhostType::value>() );
            }

            else
            {
                if ( M_worldsComm[M_cursor]->isActive() )
                    ret += nLocalDof( x, mpl::bool_<IsWithGhostType::value>() );
            }
        }

        ++M_cursor;
        return ret;
    }

    template <typename T>
    size_type
    operator()( size_type s, T const& x ) const
    {
        return this->operator()( x, s );
    }
private:
    mutable size_type M_cursor;
    size_type M_finish;
    worldscomm_ptr_t M_worldsComm;
    bool M_useOffSubSpace;
};
#endif // end MPI


template< typename IsWithGhostType>
struct NLocalDofOnProc
{

    NLocalDofOnProc( const int proc,
                     worldscomm_ptr_t const & worldsComm = Environment::worldsComm(1),
                     bool useOffSubSpace = false,
                     size_type start = 0, size_type size = invalid_v<size_type> )
        :
        M_proc(proc),
        M_cursor( start ),
        M_finish( size ),
        M_worldsComm( worldsComm ),
        M_useOffSubSpace( useOffSubSpace )
    {}

    template<typename Sig>
    struct result;

    template<typename T, typename S>
#if BOOST_VERSION < 104200
    struct result<NLocalDofOnProc( T,S )>
#else
    struct result<NLocalDofOnProc( S,T )>
#endif
:
    boost::remove_reference<S>
    {};

    template <typename T>
    size_type
    nLocalDof( T const& x, mpl::bool_<true> /**/ ) const
    {
        return x->nLocalDofWithGhostOnProc(M_proc);
    }

    template <typename T>
    size_type
    nLocalDof( T const& x, mpl::bool_<false> /**/ ) const
    {
        return x->nLocalDofWithoutGhostOnProc(M_proc);
    }

    template <typename T>
    size_type
    operator()( T const& x, size_type s ) const
    {
        size_type ret = s;

        if ( M_cursor < M_finish )
        {
            if ( M_useOffSubSpace )
            {
                ret += nLocalDof( x, mpl::bool_<IsWithGhostType::value>() );
            }

            else
            {
                if ( M_worldsComm[M_cursor]->isActive() )
                    ret += nLocalDof( x, mpl::bool_<IsWithGhostType::value>() );
            }
        }

        ++M_cursor;
        return ret;
    }

    template <typename T>
    size_type
    operator()( size_type s, T const& x ) const
    {
        return this->operator()( x, s );
    }
private:
    int M_proc;
    mutable size_type M_cursor;
    size_type M_finish;
    worldscomm_ptr_t M_worldsComm;
    bool M_useOffSubSpace;
}; // NLocalDofOnProc


template<int i,typename SpaceCompositeType>
struct InitializeContainersOff
{
    explicit InitializeContainersOff( std::shared_ptr<SpaceCompositeType> const& _space )
        :
        M_cursor( 0 ),
        M_space( _space )
    {}
    template <typename T>
    void operator()( std::shared_ptr<T> & x ) const
    {
        if ( M_cursor==i && !x )
            x = std::shared_ptr<T>( new T( M_space->template functionSpace<i>()->dof() ) );

        ++M_cursor;// warning M_cursor < nb color
    }
    mutable uint16_type M_cursor;
    std::shared_ptr<SpaceCompositeType> M_space;
};


template<int i,typename SpaceCompositeType>
struct SendContainersOn
{
    SendContainersOn( std::shared_ptr<SpaceCompositeType> const& _space,
                      std::vector<double> const& _dataToSend )
        :
        M_cursor( 0 ),
        M_space( _space ),
        M_dataToSend( _dataToSend )
    {}
    template <typename T>
    void operator()( std::shared_ptr<T> & x ) const
    {
        if ( M_cursor!=i )
        {
            int locRank=M_space->worldComm().localRank();
            int globRank=M_space->worldComm().localColorToGlobalRank( M_cursor,locRank );
            int tag = 0;
            //std::cout << "\n I am proc " << M_space->worldComm().globalRank()
            //          << " I send to proc " << globRank << std::endl;
            M_space->worldComm().globalComm().send( globRank,tag,M_dataToSend );
        }

        ++M_cursor;// warning M_cursor < nb color
    }
    mutable uint16_type M_cursor;
    std::shared_ptr<SpaceCompositeType> M_space;
    std::vector<double> M_dataToSend;
};


template<int i,typename SpaceCompositeType>
struct RecvContainersOff
{
    explicit RecvContainersOff( std::shared_ptr<SpaceCompositeType> const& _space )
        :
        M_cursor( 0 ),
        M_space( _space )
    {}
    template <typename T>
    void operator()( std::shared_ptr<T> & x ) const
    {
        if ( M_cursor==i )
        {
            std::vector<double> dataToRecv( M_space->template functionSpace<i>()->nLocalDof() );
            int locRank=M_space->worldComm().localRank();
            int globRank=M_space->worldComm().localColorToGlobalRank( i,locRank );
            int tag = 0;//locRank;
            //std::cout << "\n I am proc " << M_space->worldComm().globalRank()
            //          << " I recv to proc " << globRank << std::endl;
            M_space->worldComm().globalComm().recv( globRank,tag,dataToRecv );
            std::copy( dataToRecv.begin(), dataToRecv.end(), x->begin() );
        }

        ++M_cursor;// warning M_cursor < nb color
    }
    mutable uint16_type M_cursor;
    std::shared_ptr<SpaceCompositeType> M_space;
};



template< typename map_type >
struct searchIndicesBySpace
{
    searchIndicesBySpace()
    {}

    searchIndicesBySpace( map_type& /*u*/ )
    {}

    template<typename T>
    searchIndicesBySpace( T const& fspace, map_type& u )
    {
        u = getIndicesFromSpace( fspace,u );
    }

    template<typename Sig>
    struct result;

    template<typename T, typename M>
#if BOOST_VERSION < 104200
    struct result<searchIndicesBySpace( T,M )>
#else
    struct result<searchIndicesBySpace( M,T )>
#endif
:
    boost::remove_reference<M>
    {};

    template < typename T >
    map_type getIndicesFromSpace( T const& fspace, map_type t ) const
    {
        if ( fspace->mesh()->numElements() == 0 )
            return t;

        size_type nProc = fspace->dof()->nProcessors();

        //search for the biggest index already in t; this will give the shift for the dofs
        std::vector< size_type > max_per_space;

        for ( size_type j=0; j<t.size(); j++ )
        {
            size_type _end = t[j].size();

            if ( _end )
                max_per_space.push_back( t[j][_end-1] );
        }

        //from all max indices found, determine the biggest
        size_type max_index = 0;

        if ( t.size() )
            max_index = *max_element( max_per_space.begin(), max_per_space.end() ) + 1;

        //std::cout << "maximum index " << max_index << "\n";

        //loop in all processors
        for ( size_type i=0; i<nProc; i++ )
        {
            /*
              std::cout << "Processor " << i << " has dofs"
              << " from " << fspace->dof()->firstDof(i)
              << " to " << fspace->dof()->lastDof(i) << "\n";
            */

            size_type _first = fspace->dof()->firstDof( i );
            size_type _last  = fspace->dof()->lastDof( i );

            //the dofs numbering for the current space start at max_index+1
            for ( size_type j=_first; j<=_last; j++ )
                t[i].push_back( max_index + j );
        }

        return t;
    }
    template <typename T>
    map_type
#if BOOST_VERSION < 104200
    operator()( T const& fspace, map_type t ) const
#else
    operator()( map_type t, T const& fspace ) const
#endif
    {
        return getIndicesFromSpace( fspace,t );
    }
};

// get start for each proc ->( proc0 : 0 ), (proc1 : sumdofproc0 ), (proc2 : sumdofproc0+sumdofproc1 ) ....
struct computeStartOfFieldSplit
{
    typedef boost::tuple< uint16_type , size_type > result_type;

    template<typename T>
    result_type operator()( result_type const &  previousRes, T const& t )
    {
        auto cptSpaces = previousRes.get<0>();
        auto start = previousRes.get<1>();

        for (int proc=0;proc<t->dof()->worldComm().globalSize();++proc)
            {
                if (proc < t->dof()->worldComm().globalRank())
                    start+=t->dof()->nLocalDofWithoutGhost(proc);
            }
        return boost::make_tuple( ++cptSpaces, start );
    }
};

struct hasSubSpaceWithComponentsSplit
{
    typedef bool result_type;
    template<typename T>
    result_type operator()( result_type const &  previousRes, T const& t )
    {
        //return ( T::element_type::dof_type::is_product && T::element_type::dof_type::nComponents > 1 ) || previousRes;
        return t->map().hasIndexSplitWithComponents() || previousRes;
    }
};

// compute split
template<bool UseComponentsSplit>
struct computeNDofForEachSpace
{
    computeNDofForEachSpace(size_type startSplit)
        :
        M_indexSplit( new IndexSplit() ),
        M_startSplit(startSplit)
    {}

    std::shared_ptr<IndexSplit> const& indexSplit() const { return M_indexSplit; }

    typedef boost::tuple< uint16_type, size_type, IndexSplit > result_type;

    template<typename T>
    void operator()( T const& t ) const
    {
        this->operator()( t, mpl::bool_<UseComponentsSplit>() );
    }
    template<typename T>
    void operator()( T const& t, mpl::false_ ) const
    {
        M_indexSplit->addSplit( M_startSplit, t->map().indexSplit() );
    }
    template<typename T>
    void operator()( T const& t, mpl::true_ ) const
    {
        M_indexSplit->addSplit( M_startSplit, t->map().indexSplitWithComponents() );
    }

    mutable std::shared_ptr<IndexSplit> M_indexSplit;
    size_type M_startSplit;
};

struct rebuildDofPointsTool
{

    template <typename T>
    void operator()( std::shared_ptr<T> & x ) const
    {
        x->dof()->rebuildDofPoints( *x->mesh() );
    }
};

struct BasisName
{
    typedef std::string result_type;

    template<typename T>
    result_type operator()( result_type const & previousRes, T const& t )
    {
        std::ostringstream os;

        if ( previousRes.size() )
            os << previousRes << "_" << t->basis()->familyName();

        else
            os << t->basis()->familyName();

        return os.str();
    }
};

struct BasisOrder
{
    typedef std::vector<int> result_type;

    template<typename T>
    result_type operator()( result_type const & previousRes, T const& t )
    {
        std::vector<int> res( previousRes );
        res.push_back( t->nSubFunctionSpace() );
        return res;
    }
};

template<typename SpaceType>
struct LegacyCompositeFunctionSpaceDofOps
{
    using functionspace_vector_type = typename SpaceType::functionspace_vector_type;
    using dof_type = typename SpaceType::dof_type;
    using dof_ptrtype = typename SpaceType::dof_ptrtype;

    static size_type nDof( functionspace_vector_type const& functionspaces )
    {
        return fusion::accumulate( functionspaces, size_type( 0 ), NbDof() );
    }

    static size_type nLocalDofWithGhost( functionspace_vector_type const& functionspaces,
                                         worldscomm_ptr_t const& worldsComm )
    {
        return fusion::accumulate( functionspaces, size_type( 0 ), NLocalDof<mpl::bool_<true> >( worldsComm ) );
    }

    static size_type nLocalDofWithoutGhost( functionspace_vector_type const& functionspaces,
                                            worldscomm_ptr_t const& worldsComm )
    {
        return fusion::accumulate( functionspaces, size_type( 0 ), NLocalDof<mpl::bool_<false> >( worldsComm ) );
    }

    static size_type nLocalDofWithGhostOnProc( functionspace_vector_type const& functionspaces,
                                               int proc,
                                               worldscomm_ptr_t const& worldsComm )
    {
        return fusion::accumulate( functionspaces, size_type( 0 ), NLocalDofOnProc<mpl::bool_<true> >( proc, worldsComm ) );
    }

    static size_type nLocalDofWithoutGhostOnProc( functionspace_vector_type const& functionspaces,
                                                  int proc,
                                                  worldscomm_ptr_t const& worldsComm )
    {
        return fusion::accumulate( functionspaces, size_type( 0 ), NLocalDofOnProc<mpl::bool_<false> >( proc, worldsComm ) );
    }

    static size_type nDofStart( functionspace_vector_type const& functionspaces, size_type i )
    {
        return fusion::accumulate( functionspaces, size_type( 0 ), NbDof( 0, i ) );
    }

    static size_type nLocalDofWithGhostStart( functionspace_vector_type const& functionspaces,
                                              worldscomm_ptr_t const& worldsComm,
                                              size_type i )
    {
        return fusion::accumulate( functionspaces, size_type( 0 ), NLocalDof<mpl::bool_<true> >( worldsComm, true, 0, i ) );
    }

    static size_type nLocalDofWithoutGhostStart( functionspace_vector_type const& functionspaces,
                                                 worldscomm_ptr_t const& worldsComm,
                                                 size_type i )
    {
        return fusion::accumulate( functionspaces, size_type( 0 ), NLocalDof<mpl::bool_<false> >( worldsComm, true, 0, i ) );
    }

    static size_type nLocalDofWithGhostOnProcStart( functionspace_vector_type const& functionspaces,
                                                    int proc,
                                                    worldscomm_ptr_t const& worldsComm,
                                                    size_type i )
    {
        return fusion::accumulate( functionspaces, size_type( 0 ), NLocalDofOnProc<mpl::bool_<true> >( proc, worldsComm, true, 0, i ) );
    }

    static size_type nLocalDofWithoutGhostOnProcStart( functionspace_vector_type const& functionspaces,
                                                       int proc,
                                                       worldscomm_ptr_t const& worldsComm,
                                                       size_type i )
    {
        return fusion::accumulate( functionspaces, size_type( 0 ), NLocalDofOnProc<mpl::bool_<false> >( proc, worldsComm, true, 0, i ) );
    }

    static std::string basisName( functionspace_vector_type const& functionspaces )
    {
        return fusion::accumulate( functionspaces, std::string(), BasisName() );
    }

    static std::vector<int> basisOrder( functionspace_vector_type const& functionspaces )
    {
        return fusion::accumulate( functionspaces, std::vector<int>(), BasisOrder() );
    }

    static std::shared_ptr<IndexSplit> buildDofIndexSplit( functionspace_vector_type const& functionspaces )
    {
        auto startSplit = boost::fusion::fold( functionspaces, boost::make_tuple( 0, 0 ), computeStartOfFieldSplit() ).template get<1>();
        auto computeSplit = computeNDofForEachSpace<false>( startSplit );
        boost::fusion::for_each( functionspaces, computeSplit );
        return computeSplit.indexSplit();
    }

    static std::shared_ptr<IndexSplit> buildDofIndexSplitWithComponents( functionspace_vector_type const& functionspaces )
    {
        bool hasCompSplit = boost::fusion::fold( functionspaces, false, hasSubSpaceWithComponentsSplit() );
        if ( hasCompSplit )
        {
            auto startSplit = boost::fusion::fold( functionspaces, boost::make_tuple( 0, 0 ), computeStartOfFieldSplit() ).template get<1>();
            auto computeSplit = computeNDofForEachSpace<true>( startSplit );
            boost::fusion::for_each( functionspaces, computeSplit );
            return computeSplit.indexSplit();
        }
        return std::shared_ptr<IndexSplit>();
    }

    static void installDofIndexSplits( dof_ptrtype const& dof, functionspace_vector_type const& functionspaces )
    {
        CHECK( dof ) << "invalid legacy composite dof table";

        dof->setIndexSplit( buildDofIndexSplit( functionspaces ) );

        auto indexSplitWithComponents = buildDofIndexSplitWithComponents( functionspaces );
        if ( indexSplitWithComponents )
            dof->setIndexSplitWithComponents( indexSplitWithComponents );
    }

    static dof_ptrtype buildCompositeDof( functionspace_vector_type& functionspaces,
                                          dof_ptrtype const& dof,
                                          worldcomm_ptr_t const& worldComm,
                                          uint16_type nSubFunctionSpace )
    {
        auto dofInitTool = updateDataMapProcessStandard<dof_type>( worldComm, nSubFunctionSpace );
        return fusion::fold( functionspaces, dof, dofInitTool );
    }

    static void rebuildDofPoints( functionspace_vector_type& functionspaces )
    {
        fusion::for_each( functionspaces, rebuildDofPointsTool() );
    }
};

template<typename SpaceType>
struct LegacyCompositeFunctionSpaceElementOps
{
    using functionspace_vector_type = typename SpaceType::functionspace_vector_type;

    template<int I>
    static size_type firstSubSpaceActiveDofCount( SpaceType const& space, size_type fallback )
    {
        if constexpr ( SpaceType::is_composite )
            return space.template functionSpace<I>()->dof()->nLocalDofWithoutGhost();
        else
            return fallback;
    }

    template<int I, typename ElementType>
    static typename ElementType::template sub_element<I>::type&
    subElement( typename ElementType::element_vector_type& elements )
    {
        CHECK( fusion::at_c<I>( elements ).second ) << " has not element \n";
        return *( fusion::at_c<I>( elements ).second );
    }

    template<int I, typename ElementType>
    static typename ElementType::template sub_element<I>::type const&
    subElement( typename ElementType::element_vector_type const& elements )
    {
        CHECK( fusion::at_c<I>( elements ).second ) << " has not element \n";
        return *( fusion::at_c<I>( elements ).second );
    }

    template<int I, typename ElementType>
    static typename ElementType::template sub_element<I>::ptrtype&
    subElementPtr( typename ElementType::element_vector_type& elements )
    {
        CHECK( fusion::at_c<I>( elements ).second ) << " has not element \n";
        return fusion::at_c<I>( elements ).second;
    }

    template<int I, typename ElementType>
    static typename ElementType::template sub_element<I>::ptrtype const&
    subElementPtr( typename ElementType::element_vector_type const& elements )
    {
        CHECK( fusion::at_c<I>( elements ).second ) << " has not element \n";
        return fusion::at_c<I>( elements ).second;
    }

    template<int I, typename ElementType>
    static typename ElementType::template sub_element<I>::type
    buildSubElementView( ElementType& element,
                         boost::optional<typename ElementType::container_vector_type>& containersOffProcess,
                         std::string const& name,
                         bool updateOffViews )
    {
        size_type nbdof_start = element.functionSpace()->nLocalDofWithoutGhostStart( I );
        size_type nbdofWithGhost_start = element.functionSpace()->nLocalDofWithGhostStart( I );
        size_type startDofIndexGhost = nbdofWithGhost_start - nbdof_start;
        //if ( !Cont::is_shallow_array_adaptor_vector )
            //startDofIndexGhost += element.functionSpace()->dof()->nLocalDofWithoutGhost();

        typename mpl::at_c<functionspace_vector_type,I>::type space( element.functionSpace()->template functionSpace<I>() );
        DVLOG(2) << "Element <" << I << ">::start :  "<< nbdof_start << "\n";
        DVLOG(2) << "Element <" << I << ">::size :  "<<  space->nDof()<< "\n";
        DVLOG(2) << "Element <" << I << ">::local size :  "<<  space->nLocalDof()<< "\n";
        DVLOG(2) << "Element <" << -1 << ">::size :  "<<  element.size() << "\n";

        if ( element.functionSpace()->template functionSpace<I>()->worldComm().isActive() )
        {
            typename ElementType::ct_type ct( element,
                                              ublas::range( nbdof_start, nbdof_start + space->dof()->nLocalDofWithoutGhost() ),
                                              ublas::range( startDofIndexGhost, startDofIndexGhost + space->dof()->nLocalGhosts() ),
                                              element.functionSpace()->template functionSpace<I>()->dof() );

            if ( element.worldComm().globalSize() > 1 && updateOffViews && !element.functionSpace()->hasEntriesForAllSpaces() )
            {
                std::vector<double> dataToSend( ct.begin(), ct.end() );

                if ( !containersOffProcess )
                    containersOffProcess = typename ElementType::container_vector_type();

                fusion::for_each( *containersOffProcess, SendContainersOn<I,SpaceType>( element.functionSpace(), dataToSend ) );
            }

            DVLOG(2) << "Element <" << I << ">::range.size :  "<< ct.size() << "\n";
            DVLOG(2) << "Element <" << I << ">::range.start :  "<< ct.start() << "\n";
            return typename ElementType::template sub_element<I>::type( space, ct, name );
        }
        else
        {
            if ( !containersOffProcess )
                containersOffProcess = typename ElementType::container_vector_type();

            fusion::for_each( *containersOffProcess, InitializeContainersOff<I,SpaceType>( element.functionSpace() ) );

            if ( element.worldComm().globalSize() > 1 && updateOffViews && !element.functionSpace()->hasEntriesForAllSpaces() )
                fusion::for_each( *containersOffProcess, RecvContainersOff<I,SpaceType>( element.functionSpace() ) );

            typename ElementType::ct_type ct( *fusion::at_c<I>( *containersOffProcess ),
                                              ublas::range( 0, space->nLocalDof() ),
                                              ublas::range( 0, 0 ),
                                              element.functionSpace()->template functionSpace<I>()->dof() );

            DVLOG(2) << "Element <" << I << ">::range.size :  "<<  ct.size()<< "\n";
            DVLOG(2) << "Element <" << I << ">::range.start :  "<<  ct.start()<< "\n";

            return typename ElementType::template sub_element<I>::type( space, ct, name );
        }
    }
};

template<typename ElementType>
struct LegacyCompositeInitializeElement
{
    explicit LegacyCompositeInitializeElement( ElementType* element )
        :
        M_element( element )
    {}

    template<typename T>
    void operator()( T& x ) const
    {
        using key_type = typename T::first_type;
        using myelt_type = typename T::second_type::element_type;
        std::string name = (boost::format("%1%_%2%")%M_element->name() %key_type::value).str();

        if ( M_element->functionSpace() )
        {
            if ( !x.second || ( !x.second->functionSpace() ) )
            {
                auto e = M_element->template elementImpl<key_type::value>( name );
                auto sp = std::make_shared<myelt_type>( e );
                x = std::make_pair( key_type(), sp );
            }
        }
        else if ( !x.second )
        {
            x = std::make_pair( key_type(), nullptr );
        }
    }

    ElementType* M_element;
};

template<typename SpaceType>
struct createWorldsComm
{

    typedef typename SpaceType::mesh_ptrtype mesh_ptrtype;
    typedef typename SpaceType::meshes_list meshes_list;
    static inline const bool useMeshesList = !boost::is_base_of<MeshBase<>, meshes_list >::value;

    struct UpdateWorldsComm
    {
        UpdateWorldsComm( createWorldsComm<SpaceType> & cwc )
            :
            M_cwc( cwc )
            {}
        template<typename T>
        void operator()( T const& t) const
            {
                M_cwc.M_worldsComm.push_back( t->worldComm().shared_from_this() );
            }
        createWorldsComm<SpaceType> & M_cwc;
    };

    createWorldsComm( mesh_ptrtype const& mesh )
        {
            this->init<useMeshesList>( mesh );
        }
    template<bool _UseMeshesList >
    void init( mesh_ptrtype const& mesh, typename std::enable_if< !_UseMeshesList >::type* = nullptr )
        {
            M_worldsComm.resize( SpaceType::nSpaces, mesh->worldComm().shared_from_this() );
        }
    template<bool _UseMeshesList >
    void init( mesh_ptrtype const& mesh, typename std::enable_if< _UseMeshesList >::type* = nullptr )
        {
            boost::fusion::for_each( mesh, UpdateWorldsComm( *this ) );
        }
    worldscomm_ptr_t worldsComm()       { return M_worldsComm; }
    worldscomm_ptr_t worldsComm() const { return M_worldsComm; }

    worldscomm_ptr_t M_worldsComm;
};

template<typename SpaceType>
std::vector<DofTableExtendedType>
createInfoExtendedDofTable( DofTableExtendedType b )
{
    return std::vector<DofTableExtendedType>( SpaceType::nSpaces,b );
}
template<typename SpaceType>
std::vector<DofTableExtendedType>
createInfoExtendedDofTable( std::vector<DofTableExtendedType> const& b )
{
    CHECK( b.size() == SpaceType::nSpaces ) << "invalid extended doftable info vector size : " << b.size() << " should be : " << SpaceType::nSpaces;
    return b;
}

template<typename SpaceType>
struct createMeshSupport
{
    typedef typename SpaceType::mesh_support_vector_type mesh_support_vector_type;
    typedef typename fusion::result_of::at_c<mesh_support_vector_type,0>::type _mesh_support_ptrtype;
    typedef typename boost::remove_reference<_mesh_support_ptrtype>::type mesh_support_ptrtype;
    typedef typename mesh_support_ptrtype::element_type mesh_support_type;
    typedef typename SpaceType::mesh_ptrtype mesh_ptrtype;
    typedef typename mesh_support_type::range_elements_type range_elements_type;

    typedef typename SpaceType::meshes_list meshes_list;
    static inline const bool useMeshesList = !boost::is_base_of<MeshBase<>, meshes_list >::value;

    struct HasAllMeshSupportDefined
    {
        typedef bool result_type;
        template<typename T>
        result_type operator()( result_type const& r, T const& t) const
            {
                if ( !t )
                    return false;
                else
                    return r;
            }
    };

    struct UpdateMeshSupport
    {
        UpdateMeshSupport( createMeshSupport<SpaceType> & cms )
            :
            M_cms( cms )
            {}
        template<typename T>
        void operator()( T const& t) const
            {
                this->updateImpl<T,useMeshesList>( t );
            }
        template<typename T,bool _UseMeshesList >
            requires (!_UseMeshesList)
        void updateImpl( T const& t ) const
            {
                auto & meshSupport = boost::fusion::at_c<T::value>( M_cms.M_meshSupportVector );
                if ( meshSupport )
                    return;
                CHECK( M_cms.M_meshSupport0 ) << "no mesh support defined";
                meshSupport = M_cms.M_meshSupport0;
            }
        template<typename T,bool _UseMeshesList >
            requires _UseMeshesList
        void updateImpl( T const& t ) const
            {
                auto & meshSupport = boost::fusion::at_c<T::value>( M_cms.M_meshSupportVector );
                if ( meshSupport )
                    return;
                typedef typename fusion::result_of::at_c<mesh_support_vector_type,T::value>::type _submesh_support_ptrtype;
                typedef typename boost::remove_reference<_submesh_support_ptrtype>::type submesh_support_ptrtype;
                typedef typename submesh_support_ptrtype::element_type submesh_support_type;
                auto const& mesh = boost::fusion::at_c<T::value>( M_cms.M_mesh );
                meshSupport.reset( new submesh_support_type(mesh) );
            }

        createMeshSupport<SpaceType> & M_cms;
    };

    createMeshSupport( mesh_ptrtype const& mesh, mesh_support_vector_type const& meshSupport )
        :
        M_mesh( mesh ),
        M_meshSupportVector( meshSupport )
        {
            this->init<useMeshesList>(mesh);
        }
    template<typename RangeType>
        requires is_range_v<RangeType>
    createMeshSupport( mesh_ptrtype const& mesh, RangeType && rangeMeshElt )
        :
        M_mesh( mesh )
        {
            this->init2<useMeshesList>(mesh,std::forward<RangeType>(rangeMeshElt));
        }
    createMeshSupport( mesh_ptrtype const& mesh, mesh_support_ptrtype const& meshSupport )
        :
        M_mesh( mesh )
        {
            M_meshSupport0 = meshSupport;
            mpl::range_c<int,0,SpaceType::nSpaces> keySpaces;
            boost::fusion::for_each( keySpaces, UpdateMeshSupport( *this ) );
        }

    template<bool _UseMeshesList >
    void init( mesh_ptrtype const& mesh, typename std::enable_if< !_UseMeshesList >::type* = nullptr )
        {
            HasAllMeshSupportDefined hasMSFunctor;
            bool hasMS = boost::fusion::fold( M_meshSupportVector, true, hasMSFunctor );
            if ( !hasMS )
                M_meshSupport0.reset( new mesh_support_type(mesh) );

            mpl::range_c<int,0,SpaceType::nSpaces> keySpaces;
            boost::fusion::for_each( keySpaces, UpdateMeshSupport( *this ) );
        }
    template<bool _UseMeshesList >
    void init( mesh_ptrtype const& mesh, typename std::enable_if< _UseMeshesList >::type* = nullptr )
        {
            mpl::range_c<int,0,SpaceType::nSpaces> keySpaces;
            boost::fusion::for_each( keySpaces, UpdateMeshSupport( *this ) );
        }
    template<bool _UseMeshesList, typename RangeType >
    void init2( mesh_ptrtype const& mesh, RangeType && rangeMeshElt )
        {
            if constexpr ( _UseMeshesList )
            {
                CHECK( false ) << fmt::format( "MeshSupport not allowed in Mesh List" );
            }
            else
            {
                if ( std::forward<RangeType>( rangeMeshElt ).container() )
                {
                    M_meshSupport0.reset( new mesh_support_type(mesh,std::forward<RangeType>(rangeMeshElt) ) );
                }
                else
                {
                    M_meshSupport0.reset( new mesh_support_type(mesh) );
                }

                mpl::range_c<int,0,SpaceType::nSpaces> keySpaces;
                boost::fusion::for_each( keySpaces, UpdateMeshSupport( *this ) );
            }
        }

    mesh_ptrtype const& M_mesh;
    mesh_support_vector_type M_meshSupportVector;
    mesh_support_ptrtype M_meshSupport0;
};

template<typename SpaceType>
struct FunctionSpaceMeshSupport
{
    typedef typename SpaceType::mesh_support_vector_type mesh_support_vector_type;

    struct UpdateMeshSupport
    {
        UpdateMeshSupport( FunctionSpaceMeshSupport<SpaceType> & fsms )
            :
            M_fsms( fsms )
            {}

        template<typename T>
        void operator()( T const& t) const
            {
                this->updateImpl<T,SpaceType::is_composite>( t );
            }

        template<typename T,bool _IsComposite >
        void updateImpl( T const& t, typename std::enable_if< !_IsComposite >::type* = nullptr ) const
            {
                auto doftable = M_fsms.M_space.dof();
                if ( !doftable )
                    return;
                if ( doftable->hasMeshSupport() )
                {
                    auto & meshSupport = boost::fusion::at_c<T::value>( M_fsms.M_meshSupportVector );
                    meshSupport = doftable->meshSupport();
                }
            }
        template<typename T,bool _IsComposite >
        void updateImpl( T const& t, typename std::enable_if< _IsComposite >::type* = nullptr ) const
            {
                auto subspace = M_fsms.M_space.template functionSpace<T::value>();
                if ( !subspace )
                    return;
                auto doftable = subspace->dof();
                if ( !doftable )
                    return;
                if ( doftable->hasMeshSupport() )
                {
                    auto & meshSupport = boost::fusion::at_c<T::value>( M_fsms.M_meshSupportVector );
                    meshSupport = doftable->meshSupport();
                }
            }
        FunctionSpaceMeshSupport<SpaceType> & M_fsms;
    };

    FunctionSpaceMeshSupport( SpaceType const& space )
        :
        M_space( space )
        {
            mpl::range_c<int,0,SpaceType::nSpaces> keySpaces;
            boost::fusion::for_each( keySpaces, UpdateMeshSupport( *this ) );
        }

    SpaceType const& M_space;
    mesh_support_vector_type M_meshSupportVector;
};

#endif // FEELPP_FEELDISCR_DETAIL_FUNCTIONSPACELEGACYCOMPOSITE_HPP

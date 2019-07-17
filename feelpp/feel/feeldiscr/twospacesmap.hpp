#ifndef TWOSPACESMAP_HPP
#define TWOSPACESMAP_HPP

#include <feel/feeldiscr/functionspace.hpp>


namespace Feel
{

/**
 * This class provides a correspondancy map between two spaces and particularly
 * between a sequential space Xs defined on a submesh of the parallel space Xp.
 * For any dof in Xs the class provides the index which has the corresponding dof in parallel
 * and the index of this dof in Xp.
 *
 */
template <typename SpaceType>
class TwoSpacesMap
{
public :
    typedef SpaceType space_type;
    typedef std::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    template <int T>
    using subspace_type = typename space_type::template sub_functionspace<T>::type;
    template <int T>
    using subspace_ptrtype = std::shared_ptr<subspace_type<T>>;
    template<int T>
    using subelement_type = typename subspace_type<T>::element_type;


    static const bool is_composite = space_type::is_composite;
    static const int n_spaces = space_type::nSpaces;

    typedef typename mpl::range_c< int, 0, n_spaces > rangespace_type;


    //! Default Constructor
    TwoSpacesMap( int s=1 ) :
        M_shift( s )
    {}

    //! Constructor from spaces
    TwoSpacesMap( space_ptrtype _Xs, space_ptrtype _Xp, int s=1 ) :
        M_shift( s )
    {
        init( _Xs, _Xp );
    }


    //! initialize the map between \p _Xs and \p _Xp
    void init( space_ptrtype _Xs, space_ptrtype _Xp );

    //! build the map, depending if the spaces are composite or not
    void build() { build( mpl::bool_<is_composite>() ); }
    //! build the map, non-composite spaces
    void build( mpl::false_ );
    //! build the map, composite spaces
    void build( mpl::true_ );


    /**
     * \return the sequential index corresponding to the given global cluster dof \p p_dof
     * The sequential index is broadcasted to all proc
     */
    int clusterToSequential( size_type const& p_dof )
        {
            int proc_number = Xp->dof()->procOnGlobalCluster( p_dof );
            int s_dof = -1;
            if ( Environment::worldComm().globalRank()==proc_number )
            {

                auto searchGpDof = Xp->dof()->searchGlobalProcessDof( p_dof );
                CHECK( boost::get<0>( searchGpDof ) ) << "Did not find p_dof "<< p_dof <<" when it should be here\n";
                auto gpdof = boost::get<1>( searchGpDof );
                s_dof = parallelToSequential( gpdof );
            }
            boost::mpi::broadcast( Environment::worldComm(), s_dof, proc_number );

            return s_dof;
        }


    /**
     * \return the sequential index corresponding to the given global process dof \p p_dof
     * Return -1 if the \p p_dof is not found
     */
    int parallelToSequential( size_type const& p_dof )
        {
            auto it = M_p_to_s.find( p_dof );
            if ( it==M_p_to_s.end() )
                return -1;
            else
                return it->second;
        }


    /**
     * \return a pair <proc,p_dof> corresponding to the given \p s_dof.
     * proc is the number of the process where is located the global process p_dof
     */
    std::pair<int,int> sequentialToParallel( size_type const& s_dof )
        {
            auto it = M_s_to_p.find( s_dof );
            if ( it==M_s_to_p.end() )
                return std::make_pair<int,int>( -1, -1 );
            else
                return it->second;
        }


    //! project the parallel vector \p up on the sequential vector \p us
    void project( element_type& us, element_type const& up )
    {
        us.zero();
        for ( size_type i=0; i<us.size(); i++ )
        {
            auto e = sequentialToParallel( i );
            if ( e.first >= 0 )
                us(i) = up( e.second );
        }

        gather( us );
    }

    //! Project the parallel vector \p up from the N-th subspace on the sequential vector \p us
    template <int N>
    void project( subelement_type<N>& us, subelement_type<N> const& up )
    {
        us.zero();
        for( size_type i=0; i<us.size(); i++ )
        {
            auto it = M_s_to_p_comp[N].find( i );
            if ( it!=M_s_to_p_comp[N].end() )
            {
                us(i) = up( it->second.second );
            }
        }
        gather( us );
    }


private :
    //! Gather the contribution to the sequenetial from all process. Used for the projection
    template <typename ElementType>
    void gather( ElementType& us )
    {
        int world_size = Environment::worldComm().globalSize();
        auto Rh = us.functionSpace();
        auto ut = Rh->element();
        std::vector<ElementType> all_u( world_size, ut );

        mpi::all_gather( Environment::worldComm().globalComm(),
                         us,
                         all_u );

        for ( auto ut : all_u )
        {
            for ( int j=0; j<us.size(); j++ )
            {
                if ( ut(j)!=0 )
                    us(j) = ut(j);
                if ( ut(j)!=0 && us(j)!=0 )
                    CHECK( us(j)-ut(j)<1e-12 ) << "proc=" << Environment::worldComm().globalRank() <<", us(j)=" <<us(j)<<", ut(j)="<<ut(j)<<std::endl;
            }
        }
    }


    struct BuildForComposite
    {
        static const int n_spaces = space_type::nSpaces;

        BuildForComposite( space_ptrtype _Xs, space_ptrtype _Xp, int s,
                           std::vector<std::map<int,std::pair<int,int>>> & s_to_p_comp,
                           std::vector<std::map<int,int>> & p_to_s_comp,
                           std::map<int,std::pair<int,int>> & s_to_p,
                           std::map<int,int> & p_to_s ) :
            compositeXs( _Xs ), compositeXp( _Xp ),
            M_shift( s ), m_s_to_p_comp( s_to_p_comp ), m_p_to_s_comp( p_to_s_comp ), m_s_to_p( s_to_p ), m_p_to_s( p_to_s )
        {
            m_s_to_p_comp.resize( n_spaces );
            m_p_to_s_comp.resize( n_spaces );
        }

        template <typename T>
        void operator()( T const& t ) const
        {

            int q = T::value;
            auto Xs = compositeXs->template functionSpace<T::value>();
            auto Xp = compositeXp->template functionSpace<T::value>();

            m_s_to_p_comp[q].clear();
            m_p_to_s_comp[q].clear();

            std::map<std::set<int>,int> elts_map_s;
            for ( auto const& eltWrap : elements(Xs->mesh()) )
            {
                auto const& elt = unwrap_ref( eltWrap );
                std::set<int> pts_id;
                for ( int p=0; p<elt.nPoints(); p++ )
                    pts_id.insert( elt.point(p).id() );
                elts_map_s[pts_id] = elt.id();
            }
            // loop on the elements of the full mesh splited between all procs
            for ( auto const& eltWrap : elements(Xp->mesh()) )
            {
                auto const& elt = unwrap_ref( eltWrap );
                std::set<int> pts_id;
                for ( int p=0; p<elt.nPoints(); p++ )
                    pts_id.insert( elt.point(p).id() + M_shift );

                // check if this elements is in the submesh of the sequential space
                auto map_it = elts_map_s.find( pts_id );
                if ( map_it!=elts_map_s.end() ) // the element exists in the seq mesh
                {
                    int eid_s = map_it->second;

                    // loop on each ldof of the element : get the globaldof id associated
                    // and put it in the maps
                    for ( auto const& ldof : Xp->dof()->localDof(elt.id()) )
                    {
                        int gdof_s = Xs->dof()->localToGlobalId( eid_s, ldof.first.localDof() );
                        int gdof_p = ldof.second.index() ;

                        auto s_to_p_it = m_s_to_p_comp[q].find( gdof_s );
                        if ( s_to_p_it==m_s_to_p_comp[q].end() )
                        {
                            m_s_to_p_comp[q][gdof_s] = std::make_pair( Xp->worldComm().globalRank(),
                                                                       gdof_p );
                            m_p_to_s_comp[q][gdof_p] = gdof_s;

                            int vec_gdof_s = compositeXs->dof()->dofIdToContainerId( q, gdof_s );
                            int vec_gdof_p = compositeXp->dof()->dofIdToContainerId( q, gdof_p );
                            m_s_to_p[vec_gdof_s] = std::make_pair( Xp->worldComm().globalRank(),
                                                                   vec_gdof_p );
                            m_p_to_s[vec_gdof_p] = vec_gdof_s;
                        }
                        //#if !defined( NDEBUG )
                        else
                        {
                            int old_gdof_p = m_s_to_p_comp[q][gdof_s].second;
                            CHECK( old_gdof_p==gdof_p ) <<"Error on proc "<< Environment::worldComm().globalRank() <<", the dof_s has the corrspondency to gdof_p :" << old_gdof_p <<", and the new one is "<< gdof_p <<std::endl;
                        }
                        //#endif
                    }
                }
            }
        }

    private :
        space_ptrtype compositeXs, compositeXp;
        int M_shift;

        std::vector<std::map<int,std::pair<int,int>>> & m_s_to_p_comp;
        std::vector<std::map<int,int>> & m_p_to_s_comp;
        std::map<int,std::pair<int,int>> &  m_s_to_p;
        std::map<int,int> & m_p_to_s;

    };


private :
    space_ptrtype Xs, Xp;
    int M_shift;

    std::map<int,std::pair<int,int>>  M_s_to_p;
    std::map<int,int> M_p_to_s;

    std::vector<std::map<int,std::pair<int,int>>>  M_s_to_p_comp;
    std::vector<std::map<int,int>> M_p_to_s_comp;

}; // class TwoSpaceMap


template <typename SpaceType>
void
TwoSpacesMap<SpaceType>::init( space_ptrtype _Xs, space_ptrtype _Xp )
{
    Xs = _Xs;
    Xp = _Xp;

    if ( n_spaces>1 )
    {
        M_s_to_p_comp.resize(n_spaces);
        M_p_to_s_comp.resize(n_spaces);
    }

    build();
}


template <typename SpaceType>
void
TwoSpacesMap<SpaceType>::build( mpl::false_ )
{
    M_s_to_p.clear();
    M_p_to_s.clear();

    // create a map between points id and elements in Xs (supposed to be small)
    std::map<std::set<int>,int> elts_map_s;
    for ( auto const& eltWrap : elements(Xs->mesh()) )
    {
        auto const& elt = unwrap_ref( eltWrap );
        std::set<int> pts_id;
        for ( int p=0; p<elt.nPoints(); p++ )
            pts_id.insert( elt.point(p).id() );
        elts_map_s[pts_id] = elt.id();
    }

    // loop on the elements of the full mesh splited between all procs
    for ( auto const& eltWrap : elements(Xp->mesh()) )
    {
        auto const& elt = unwrap_ref( eltWrap );
        std::set<int> pts_id;
        for ( int p=0; p<elt.nPoints(); p++ )
            pts_id.insert( elt.point(p).id() + M_shift );

        // check if this elements is in the submesh of the sequential space
        auto map_it = elts_map_s.find( pts_id );
        if ( map_it!=elts_map_s.end() ) // the element exists in the seq mesh
        {
            int eid_s = map_it->second;

            // loop on each ldof of the element : get the globaldof id associated
            // and put it in the maps
            for ( auto const& ldof : Xp->dof()->localDof(elt.id()) )
            {
                int gdof_s = Xs->dof()->localToGlobalId( eid_s, ldof.first.localDof() );
                int gdof_p = ldof.second.index() ;

                auto s_to_p_it = M_s_to_p.find( gdof_s );
                if ( s_to_p_it==M_s_to_p.end() )
                {
                    M_s_to_p[gdof_s] = std::make_pair( Xp->worldComm().globalRank(),
                                                       gdof_p );
                    M_p_to_s[gdof_p] = gdof_s;
                }
#if !defined( NDEBUG )
                else
                {
                    int old_gdof_p = M_s_to_p[gdof_s].second;
                    CHECK( old_gdof_p==gdof_p ) <<"Error on proc "<< Environment::worldComm().globalRank() <<", the dof_s has the corrspondency to gdof_p :" << old_gdof_p <<", and the new one is "<< gdof_p <<std::endl;
                }
#endif
            }
        }
    }

} // build non-composite


template <typename SpaceType>
void
TwoSpacesMap<SpaceType>::build( mpl::true_ )
{
    rangespace_type range;
    boost::fusion::for_each( range, BuildForComposite( Xs, Xp, M_shift, M_s_to_p_comp, M_p_to_s_comp, M_s_to_p, M_p_to_s ) );
}




} //namespace Feel
#endif

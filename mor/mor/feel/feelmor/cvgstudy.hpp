#ifndef CVG_STUDY_H
#define CVG_STUDY_H 1

#include <feel/feelmor/crb.hpp>
#include <feel/feelmor/crbmodel.hpp>
#include <feel/feelmor/eim.hpp>
#include <feel/feelmor/deim.hpp>
#include <feel/feelmor/ser.hpp>

namespace Feel
{

template<typename ModelType,
         template < typename ReducedMethod > class RM=CRB,
         template < typename ModelInterface > class RmModel=CRBModel>
class CvgStudy
{
public :
    using model_type = ModelType;
    using crb_model_type = RmModel<model_type>;
    using crb_model_ptrtype = std::shared_ptr<crb_model_type>;
    using crb_type = RM<crb_model_type>;
    using crb_ptrtype = std::shared_ptr<crb_type>;
    using mesh_type = typename model_type::mesh_type;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;
    using sampling_type = typename crb_type::sampling_type;
    using sampling_ptrtype = std::shared_ptr<sampling_type>;
    using ser_type = SER<crb_type>;
    using ser_ptrtype = std::shared_ptr<ser_type>;

    typedef typename model_type::element_type element_type;
    typedef typename model_type::functionspace_type functionspace_type;
    static const int nb_spaces = functionspace_type::nSpaces;
    typedef typename mpl::if_< boost::is_same< mpl::int_<nb_spaces> , mpl::int_<2> > , fusion::vector< mpl::int_<0>, mpl::int_<1> >  ,
                               typename mpl::if_ < boost::is_same< mpl::int_<nb_spaces> , mpl::int_<3> > ,
                                                   fusion::vector < mpl::int_<0> , mpl::int_<1> , mpl::int_<2> >,
                                                   typename mpl::if_< boost::is_same< mpl::int_<nb_spaces> , mpl::int_<4> >,
                                                                      fusion::vector< mpl::int_<0>, mpl::int_<1>, mpl::int_<2>, mpl::int_<3> >,
                                                                      fusion::vector< mpl::int_<0>, mpl::int_<1>, mpl::int_<2>, mpl::int_<3>, mpl::int_<4> >
                                                                      >::type >::type >::type index_vector_type;

    static const bool is_composite = functionspace_type::is_composite;

    CvgStudy()
        :
        crb_model( std::make_shared<crb_model_type>(crb::stage::offline) ),
        crb( crb_type::New( crb_model->model()->modelName(), crb_model, crb::stage::offline) )
    {
        crb_model = std::make_shared<crb_model_type>(crb::stage::offline);
        crb = crb_type::New( crb_model->model()->modelName(),
                             crb_model, crb::stage::offline);
        if( crb_model->hasEim() && crb_model->useSER() )
        {
            ser = std::make_shared<ser_type>( crb, crb_model );
            ser->run();
        }
        else
            crb->offline();
    }

    void run( int const& s )
    {
        using Feel::cout;
        int N = crb->WNmuSize();
        int i=0;
        auto mesh = crb_model->model()->functionSpace()->mesh();
        int size = s;
        auto zero_mean = crb->hasZeroMean();

        sampling_ptrtype sampling( new sampling_type( crb_model->parameterSpace() ) );
        if ( size==-1 )
        {
            sampling = crb->wnmu();
            size = sampling->size();
        }
        else
        {
            sampling->clear();
            sampling->randomize( size, true );
        }

        std::vector<std::vector<double>> output_error( size );
        std::vector<std::vector<double>> l2_error( size );
        std::vector<std::vector<double>> timers( size );
        std::vector<std::vector<std::vector<double>>> l2_composite;
        std::vector<bool> converged ( size, true );
        std::map<std::string,std::vector<std::vector<double>>> timers_rb;

        boost::mpi::timer timer;
        if ( is_composite )
        {
            l2_composite.resize( nb_spaces );
            for ( int i=0; i<nb_spaces; i++ )
                l2_composite[i].resize( size );
        }


        M_e = exporter( _mesh=crb->model()->functionSpace()->mesh(), _name="cvg-study" );

        int j=1;
        cout << "Convergence Study Starts\n"
             << "==================================================\n";
        crb->setCheckCvg(true);
        for ( auto mu : *sampling )
        {
            cout << "["<<j<<"] mu="<<mu.toString() <<std::endl;
            j++;
            tic();
            auto u_fem = crb->femSolve(mu);//crb_model->model()->solve( mu );
            auto u_rb = u_fem.functionSpace()->element();

            // TODO fix in general case
            if ( zero_mean[1] )
            {
                auto p = u_fem.template element<1>();
                auto pv = p.functionSpace()->element();
                pv.on( elements(mesh), idv(p) );
                double mean_p = p.min();//mean( elements(mesh), idv(p) )(0,0);
                p.on( elements(mesh), idv(pv)-cst(mean_p) );
            }

            double output_fem = crb_model->output( crb->outputIndex(), mu, u_fem );
            double u_l2 = l2Norm(u_fem);

            if ( u_l2!=0 )
            {
                vectorN_type l2_c;
                if ( is_composite )
                {
                    l2_c = l2Composite(u_fem);
                }
                toc("FEM Solve");

                cout << std::setw(5) << "N"
                     << std::setw(25) << "outputFem"
                     << std::setw(25) << "outputRB"
                     << std::setw(25) << "outputError"
                     << std::setw(25) << "l2Error";
                if ( is_composite )
                    for ( int k=0; k<nb_spaces; k++ )
                        cout << std::setw(25) << "l2Error_"<<k;
                cout << std::endl;

                std::vector<vectorN_type> Un, Undu, Unold, Unduold;
                Un.resize(1);
                Undu.resize(1);
                Unold.resize(1);
                Unduold.resize(1);
                std::vector<double> output_vector;
                output_vector.resize( 1, 0. );
                for ( int n=1; n<=N; n++ )
                {
                    timer.restart();
                    crb->onlineSolve( n, mu, Un, Undu, Unold, Unduold, output_vector, 0, false, true  );
                    double time = timer.elapsed();
                    timers[i].push_back( time );
                    auto t_map = crb->timerMap();
                    for ( auto m : t_map )
                    {
                        std::string name = m.first;
                        if ( timers_rb[name].size()==0 )
                            timers_rb[name].resize(size);
                        timers_rb[name][i].push_back(m.second);
                    }

                    //auto output_vector=o.template get<0>();
                    double output_rb = output_vector[0];
                    output_error[i].push_back( math::abs((output_fem-output_rb)/output_fem) );

                    vectorN_type u = Un[0];
                    u_rb = crb->expansion( u, n );
                    // TODO fix in general case
                    if ( zero_mean[1] )
                    {
                        auto p = u_rb.template element<1>();
                        auto pv = p.functionSpace()->element();
                        pv.on( elements(mesh), idv(p) );
                        double mean_p = p.min();//mean( elements(mesh), idv(p) )(0,0);
                        p.on( elements(mesh), idv(pv)-cst(mean_p) );
                    }

                    if ( !crb->onlineHasConverged() )
                    {
                        Feel::cout << "Not converged !\n";
                        converged[i] = false;
                        break;
                    }


                    auto u_diff = u_fem-u_rb;
                    l2_error[i].push_back( l2Norm(u_diff)/u_l2 );


                    vectorN_type diff_composite;
                    if ( is_composite )
                    {
                        diff_composite = l2Composite(u_diff);
                        for ( int k=0; k<nb_spaces; k++ )
                            l2_composite[k][i].push_back( diff_composite(k)/l2_c(k) );
                    }

                    cout << std::setw(5) << n
                         << std::setw(25) << output_fem
                         << std::setw(25) << output_rb
                         << std::setw(25) << output_error[i][n-1]
                         << std::setw(25) << l2_error[i][n-1];
                    if ( is_composite )
                        for ( int k=0; k<nb_spaces; k++ )
                            cout << std::setw(25) << l2_composite[k][i][n-1];
                    cout<<std::endl;
                }
                exportSolutions( u_rb, u_fem, i );
                cout << "==================================================\n";
                i++;
            }
            else
            {
                l2_error.pop_back();
                output_error.pop_back();
                if ( is_composite )
                    for ( int k=0; k<nb_spaces; k++ )
                        l2_composite[k].pop_back();
            }
        }

        M_e->save();
        this->writeOnFile( l2_error, "l2", converged );
        this->writeOnFile( output_error, "output", converged );
        this->writeOnFile( timers, "timers_rb", converged, false );
        for ( auto m : timers_rb )
        {
            if ( m.second.size()!=0 )
                this->writeOnFile( m.second, "timer"+m.first, converged, false );
        }
        if ( is_composite )
            for ( int k=0; k<nb_spaces; k++ )
                this->writeOnFile( l2_composite[k], "l2_"+std::to_string(k), converged );
    }


private :
    void exportSolutions( element_type const& u_rb, element_type const& u_fem, int const& i )
        {
            exportSolutions( u_rb, u_fem, i, mpl::bool_< is_composite >() );
        }
    void exportSolutions( element_type const& u_rb, element_type const& u_fem, int const& i, mpl::bool_<false> )
        {
            M_e->add( "u_rb"+std::to_string(i), u_rb );
            M_e->add( "u_fem"+std::to_string(i),u_fem );
        }
    void exportSolutions( element_type const& u_rb, element_type const& u_fem, int const& i, mpl::bool_<true> )
        {
            ExportSolutionsComposite e( u_rb, u_fem, M_e, i );
            index_vector_type index_vector;
            fusion::for_each( index_vector, e );
        }

    void writeOnFile( std::vector<std::vector<double>>& data, std::string const& name,
                      std::vector<bool> const& converged, bool geo_mean=true )
    {
        for ( int i= converged.size()-1; i>=0; i-- )
        {
            if ( !converged[i] )
                data.erase( data.begin() + i );
        }

        int N = crb->WNmuSize();
        int size = data.size();

        if ( Environment::isMasterRank() )
        {
            std::ofstream error_dat ( name+"_error.dat", std::ios::trunc );
            std::ofstream stats_dat ( name+"_stats.dat", std::ios::trunc );

            for ( int i=0; i<size; i++ )
            {
                for ( int n=0; n<N; n++ )
                    error_dat << std::setw(25) << data[i][n];
                error_dat << std::endl;
            }
            error_dat.close();

            stats_dat << std::setw(5)<< "N" << std::setw(25)<< "min"
                      << std::setw(25)<< "max" << std::setw(25)<< "mean\n";

            for ( int n=0; n<N; n++ )
            {
                double data_min = data[0][n];
                double data_max = data[0][n];
                double data_mean = 0;
                if ( geo_mean )
                    data_mean=1;
                for ( int i=0; i<size; i++ )
                {
                    data_min = std::min( data_min, data[i][n] );
                    data_max = std::max( data_max, data[i][n] );
                    if ( geo_mean )
                        data_mean *= std::pow(data[i][n], 1./size);
                    else
                        data_mean += data[i][n]/(double)size;
                }
                stats_dat << std::setw(5)<< n+1 << std::setw(25)<< data_min
                          << std::setw(25)<< data_max << std::setw(25)<< data_mean <<std::endl;
            }
            stats_dat.close();

        }
    }

    double l2Norm( element_type const& u )
        {
            return l2Norm( u, mpl::bool_< is_composite >() );
        }
    double l2Norm( element_type const& u, mpl::bool_<false> )
        {
            return normL2( _range=u.functionSpace()->template rangeElements<0>(), _expr=idv(u) );
        }
    double l2Norm( element_type const& u, mpl::bool_<true>)
        {
            ComputeNormL2InCompositeCase compute_normL2_in_composite_case( u );
            index_vector_type index_vector;
            fusion::for_each( index_vector, compute_normL2_in_composite_case );
            return compute_normL2_in_composite_case.norm();
        }

    vectorN_type l2Composite( element_type const& u )
        {
            return l2Composite( u, mpl::bool_< is_composite >() );
        }
    vectorN_type l2Composite( element_type const& u, mpl::bool_<false> )
        {
            vectorN_type dummy;
            return dummy;
        }
    vectorN_type l2Composite( element_type const& u, mpl::bool_<true>)
        {
            ComputeNormL2InCompositeCase compute_normL2_in_composite_case( u );
            index_vector_type index_vector;
            fusion::for_each( index_vector, compute_normL2_in_composite_case );
            return compute_normL2_in_composite_case.norms();
        }

    struct ComputeNormL2InCompositeCase
    {
        explicit ComputeNormL2InCompositeCase( element_type const composite_u )
            :
            M_composite_u( composite_u )
            {}

        template< typename T >
        void
        operator()( const T& t ) const
            {
                int i = T::value;
                if( i == 0 )
                    M_vec.resize( 1 );
                else
                    M_vec.conservativeResize( i+1 );

                auto u = M_composite_u.template element< T::value >();
                double norm  = normL2(_range=u.functionSpace()->template rangeElements<0>(),_expr=( idv(u) ) );
                M_vec(i)= norm;
            }

        double norm()
            { return std::sqrt(M_vec.sum()); }

        vectorN_type const& norms()
            { return M_vec; }

        mutable vectorN_type M_vec;
        element_type M_composite_u;
    };

    struct ExportSolutionsComposite
    {
        ExportSolutionsComposite( element_type const& u_rb, element_type const& u_fem, std::shared_ptr<Exporter<mesh_type>> e, int const& i ) :
            m_urb( u_rb ),
            m_ufem( u_fem ),
            m_e( e ),
            m_i( i )
            {}

        template <typename T>
        void operator()( T const& t ) const
            {
                auto suburb = m_urb.template element<T::value>();
                auto subufem = m_ufem.template element<T::value>();
                m_e->add( "u_rb"+std::to_string(m_i)+"_"+std::to_string(T::value), suburb );
                m_e->add( "u_fem"+std::to_string(m_i)+"_"+std::to_string(T::value), subufem );
            }

    private:
        element_type m_urb,m_ufem;
        std::shared_ptr<Exporter<mesh_type>> m_e;
        int m_i;
    };

private :
    crb_model_ptrtype crb_model;
    crb_ptrtype crb;
    ser_ptrtype ser;
    std::shared_ptr<Exporter<mesh_type>> M_e;
};

} //namespace Feel
#endif

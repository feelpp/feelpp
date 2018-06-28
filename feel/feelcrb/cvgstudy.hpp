#ifndef CVG_STUDY_H
#define CVG_STUDY_H 1

#include <feel/feelcrb/crb.hpp>
#include <feel/feelcrb/crbmodel.hpp>
#include <feel/feelcrb/eim.hpp>
#include <feel/feelcrb/deim.hpp>
#include <feel/feelcrb/ser.hpp>

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
    using crb_model_ptrtype = boost::shared_ptr<crb_model_type>;
    using crb_type = RM<crb_model_type>;
    using crb_ptrtype = boost::shared_ptr<crb_type>;
    using mesh_type = typename model_type::mesh_type;
    using mesh_ptrtype = boost::shared_ptr<mesh_type>;
    using sampling_type = typename crb_type::sampling_type;
    using sampling_ptrtype = boost::shared_ptr<sampling_type>;
    using ser_type = SER<crb_type>;
    using ser_ptrtype = boost::shared_ptr<ser_type>;

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
    {
        crb_model = boost::make_shared<crb_model_type>(crb::stage::offline);
        crb = crb_type::New( crb_model->model()->modelName(),
                             crb_model, crb::stage::offline);
        if( crb_model->hasEim() && crb_model->useSER() )
        {
            ser = boost::make_shared<ser_type>( crb, crb_model );
            ser->run();
        }
        else
            crb->offline();
    }

    void run( int const& size )
    {
        using Feel::cout;
        int N = crb->dimension();
        int i=0;
        auto mesh = crb_model->model()->functionSpace()->mesh();
        std::vector<std::vector<double>> output_error( size );
        std::vector<std::vector<double>> l2_error( size );
        std::vector<std::vector<std::vector<double>>> l2_composite;

        if ( is_composite )
        {
            l2_composite.resize( nb_spaces );
            for ( int i=0; i<nb_spaces; i++ )
                l2_composite[i].resize( size );
        }

        sampling_ptrtype sampling( new sampling_type( crb_model->parameterSpace() ) );
        sampling->clear();
        sampling->randomize( size, true );

        cout << "Convergence Study Starts\n"
             << "==================================================\n";

        for ( auto mu : *sampling )
        {
            cout << "mu="<<mu.toString() <<std::endl;

            tic();
            auto u_fem = crb->offlineSolve(mu);//crb_model->model()->solve( mu );
            double output_fem = crb_model->output( crb->outputIndex(), mu, u_fem );
            double u_l2 = l2Norm(u_fem);
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


            for ( int n=1; n<=N; n++ )
            {
                std::vector<vectorN_type> Un, Undu, Unold, Unduold;
                auto o = crb->lb( n, mu, Un, Undu, Unold, Unduold );

                auto output_vector=o.template get<0>();
                double output_rb = output_vector[0];
                output_error[i].push_back( math::abs((output_fem-output_rb)/output_fem) );

                vectorN_type u = Un[0];
                auto u_rb = crb->expansion( u, n );
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
            cout << "==================================================\n";
            i++;
        }

        this->writeOnFile( l2_error, "l2");
        this->writeOnFile( output_error, "output");
        if ( is_composite )
            for ( int k=0; k<nb_spaces; k++ )
                this->writeOnFile( l2_composite[k], "l2_"+std::to_string(k) );
    }


private :
    void writeOnFile( std::vector<std::vector<double>> const& data, std::string const& name )
    {
        int N = crb->dimension();
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
                double data_mean = data[0][n];
                for ( int i=1; i<size; i++ )
                {
                    data_min = std::min( data_min, data[i][n] );
                    data_max = std::max( data_max, data[i][n] );
                    data_mean += data[i][n];
                }
                data_mean = data_mean/size;
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
            return norml2( u.functionSpace()->template rangeElements<0>(), idv(u) );
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
            ComputeNormL2InCompositeCase compute_normL2_in_composite_case( u );
            index_vector_type index_vector;
            fusion::for_each( index_vector, compute_normL2_in_composite_case );
            return compute_normL2_in_composite_case.norms();
        }

    struct ComputeNormL2InCompositeCase
    {
        ComputeNormL2InCompositeCase( element_type const composite_u )
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
                M_vec(i)= norm ;
            }

        double norm()
            { return M_vec.sum(); }

        vectorN_type const& norms()
            { return M_vec; }

        mutable vectorN_type M_vec;
        element_type M_composite_u;
    };


private :
    crb_model_ptrtype crb_model;
    crb_ptrtype crb;
    ser_ptrtype ser;
};

} //namespace Feel
#endif

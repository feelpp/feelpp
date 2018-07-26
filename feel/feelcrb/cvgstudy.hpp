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
    using crb_model_ptrtype = std::shared_ptr<crb_model_type>;
    using crb_type = RM<crb_model_type>;
    using crb_ptrtype = std::shared_ptr<crb_type>;
    using mesh_type = typename model_type::mesh_type;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;
    using sampling_type = typename crb_type::sampling_type;
    using sampling_ptrtype = std::shared_ptr<sampling_type>;
    using ser_type = SER<crb_type>;
    using ser_ptrtype = std::shared_ptr<ser_type>;

    CvgStudy()
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

    void run( int const& size )
    {
        using Feel::cout;
        int N = crb->dimension();
        int i=0;
        auto mesh = crb_model->model()->functionSpace()->mesh();
        std::vector<std::vector<double>> output_error( size );
        std::vector<std::vector<double>> l2_error( size );

        sampling_ptrtype sampling( new sampling_type( crb_model->parameterSpace() ) );
        sampling->clear();
        sampling->randomize( size, true );

        cout << "Convergence Study Starts\n"
             << "==================================================\n";

        for ( auto mu : *sampling )
        {
            cout << "mu="<<mu.toString() <<std::endl;

            tic();
            auto u_fem = crb_model->model()->solve( mu );
            double output_fem = crb_model->output( crb->outputIndex(), mu, u_fem );
            double u_l2 = normL2( elements(mesh), idv(u_fem) );
            toc("FEM Solve");

            cout << std::setw(5) << "N"
                 << std::setw(25) << "outputFem"
                 << std::setw(25) << "outputRB"
                 << std::setw(25) << "outputError"
                 << std::setw(25) << "l2Error\n";

            for ( int n=1; n<=N; n++ )
            {
                std::vector<vectorN_type> Un, Undu, Unold, Unduold;
                auto o = crb->lb( n, mu, Un, Undu, Unold, Unduold );

                auto output_vector=o.template get<0>();
                double output_rb = output_vector[0];
                output_error[i].push_back( math::abs((output_fem-output_rb)/output_fem) );

                vectorN_type u = Un[0];
                auto u_rb = crb->expansion( u, n );
                l2_error[i].push_back( normL2( elements(mesh), idv(u_rb)-idv(u_fem) )/u_l2 );

                cout << std::setw(5) << n
                     << std::setw(25) << output_fem
                     << std::setw(25) << output_rb
                     << std::setw(25) << output_error[i][n-1]
                     << std::setw(25) << l2_error[i][n-1] <<std::endl;
            }
            cout << "==================================================\n";
            i++;
        }

        this->writeOnFile( l2_error, "l2");
        this->writeOnFile( output_error, "output");
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



private :
    crb_model_ptrtype crb_model;
    crb_ptrtype crb;
    ser_ptrtype ser;
};

} //namespace Feel
#endif

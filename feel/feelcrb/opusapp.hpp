/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-06-18

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file opusapp.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-06-18
 */
#ifndef __OpusApp_H
#define __OpusApp_H 1

#include <feel/feel.hpp>

#include <boost/assign/std/vector.hpp>
#include <feel/feelcrb/crb.hpp>
#include <feel/feelcrb/eim.hpp>
#include <feel/feelcrb/crbmodel.hpp>
#include <boost/serialization/version.hpp>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>


namespace Feel
{
po::options_description opusapp_options( std::string const& prefix );
std::string _o( std::string const& prefix, std::string const& opt )
{
    std::string o = prefix;

    if ( !o.empty() )
        o += ".";

    return o + opt;
}


enum class SamplingMode
{
    RANDOM = 0, EQUIDISTRIBUTED = 1, LOGEQUIDISTRIBUTED = 2
};

#define prec 4
#define Pdim 7
#define fill ' '
#define dmanip std::scientific << std::setprecision( prec )
#define hdrmanip(N) std::setw(N) << std::setfill(fill) << std::right
#define tabmanip(N) std::setw(N) << std::setfill(fill) << std::right << dmanip


/**
 * \class OpusApp
 * \brief Certified Reduced Basis application
 *
 * @author Christophe Prud'homme
 */
    template<typename ModelType,
             template < typename ReducedMethod > class RM=CRB,
             template < typename ModelInterface > class Model=CRBModel>
class OpusApp   : public Application
{
    typedef Application super;
public:

    typedef double value_type;

    //! model type
    typedef ModelType model_type;
    typedef boost::shared_ptr<ModelType> model_ptrtype;

    //! function space type
    typedef typename model_type::functionspace_type functionspace_type;
    typedef typename model_type::functionspace_ptrtype functionspace_ptrtype;

    typedef typename model_type::element_type element_type;

    typedef Eigen::VectorXd vectorN_type;

#if 0
    //old
    typedef CRBModel<ModelType> crbmodel_type;
    typedef boost::shared_ptr<crbmodel_type> crbmodel_ptrtype;
    typedef CRB<crbmodel_type> crb_type;
    typedef boost::shared_ptr<crb_type> crb_ptrtype;
#endif

    typedef Model<ModelType> crbmodel_type;
    typedef boost::shared_ptr<crbmodel_type> crbmodel_ptrtype;
    typedef RM<crbmodel_type> crb_type;
    typedef boost::shared_ptr<crb_type> crb_ptrtype;

    typedef CRBModel<ModelType> crbmodelbilinear_type;

    typedef typename ModelType::parameter_type parameter_type;
    typedef std::vector< parameter_type > vector_parameter_type;
    OpusApp()
        :
        super(),
        M_mode( ( CRBModelMode )this->vm()[_o( this->about().appName(),"run.mode" )].template as<int>() )
        {
            this->init();
        }

    OpusApp( AboutData const& ad, po::options_description const& od )
        :
        super( ad, opusapp_options( ad.appName() ).add( od ).add( crbOptions() ).add( feel_options() ) ),
        M_mode( ( CRBModelMode )this->vm()[_o( this->about().appName(),"run.mode" )].template as<int>() )
        {
            this->init();
        }

    OpusApp( AboutData const& ad, po::options_description const& od, CRBModelMode mode )
        :
        super( ad, opusapp_options( ad.appName() ).add( od ).add( crbOptions() ).add( feel_options() ) ),
        M_mode( mode )
        {
            this->init();
        }

    OpusApp( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, opusapp_options( ad.appName() ).add( od ).add( crbOptions() ).add( feel_options() ) ),
        M_mode( ( CRBModelMode )this->vm()[_o( this->about().appName(),"run.mode" )].template as<int>() )
        {
            this->init();
        }
    OpusApp( int argc, char** argv, AboutData const& ad, po::options_description const& od, CRBModelMode mode )
        :
        super( argc, argv, ad, opusapp_options( ad.appName() ).add( od ).add( crbOptions() ).add( feel_options() ) ),
        M_mode( mode )
        {
            this->init();
        }
    void init()
        {
            try
            {
                M_current_path = fs::current_path();

                std::srand( static_cast<unsigned>( std::time( 0 ) ) );
                std::cout << this->about().appName() << std::endl;
                LOG(INFO) << "[OpusApp] constructor " << this->about().appName()  << "\n";

                // Check if user have given a name for result files repo
                // Note : this name is also use for database storage
                std::string results_repo_name;
                if( this->vm().count("crb.results-repo-name") )
                    results_repo_name = this->vm()["crb.results-repo-name"].template as<std::string>();
                else
                    results_repo_name = "default_repo";

                LOG(INFO) << "Name for results repo : " << results_repo_name << "\n";
                this->changeRepository( boost::format( "%1%/%2%/" )
                                        % this->about().appName()
                                        % results_repo_name
                                        );

                LOG(INFO) << "[OpusApp] ch repo" << "\n";
                this->setLogs();
                LOG(INFO) << "[OpusApp] set Logs" << "\n";
                LOG(INFO) << "[OpusApp] mode:" << ( int )M_mode << "\n";
                model = crbmodel_ptrtype( new crbmodel_type( this->vm(),M_mode ) );
                LOG(INFO) << "[OpusApp] get model done" << "\n";

                crb = crb_ptrtype( new crb_type( this->about().appName(),
                                                 this->vm() ,
                                                 model ) );

                LOG(INFO) << "[OpusApp] get crb done" << "\n";

                //VLOG(1) << "[OpusApp] get crb done" << "\n";
                //crb->setTruthModel( model );
                //VLOG(1) << "[OpusApp] constructor done" << "\n";
            }

            catch ( boost::bad_any_cast const& e )
            {
                std::cout << "[OpusApp] a bad any cast occured, probably a nonexistant or invalid  command line/ config options\n";
                std::cout << "[OpusApp] exception reason: " << e.what() << "\n";
            }

        }

    void setMode( std::string const& mode )
        {
            if ( mode == "pfem" ) M_mode = CRBModelMode::PFEM;

            if ( mode == "crb" ) M_mode = CRBModelMode::CRB;

            if ( mode == "scm" ) M_mode = CRBModelMode::SCM;

            if ( mode == "scm_online" ) M_mode = CRBModelMode::SCM_ONLINE;

            if ( mode == "crb_online" ) M_mode = CRBModelMode::CRB_ONLINE;
        }
    void setMode( CRBModelMode mode )
        {
            M_mode = mode;
        }

    void loadDB()
        {
            bool use_predefined = this->vm()["crb.use-predefined-WNmu"].template as<bool>();
            std::string file_name = ( boost::format("SamplingWNmu") ).str();
            int NlogEquidistributed = this->vm()["crb.use-logEquidistributed-WNmu"].template as<int>();
            int Nequidistributed = this->vm()["crb.use-equidistributed-WNmu"].template as<int>();
            typename crb_type::sampling_ptrtype Sampling( new typename crb_type::sampling_type( model->parameterSpace() ) );
            if( NlogEquidistributed+Nequidistributed > 0 )
            {
                if( NlogEquidistributed > 0 )
                    Sampling->logEquidistribute( NlogEquidistributed  );
                if( Nequidistributed > 0 )
                    Sampling->equidistribute( Nequidistributed  );
                Sampling->writeOnFile(file_name);
            }

            if ( M_mode == CRBModelMode::PFEM )
                return;

            if ( !crb->scm()->isDBLoaded() || crb->scm()->rebuildDB() )
            {
                if ( M_mode == CRBModelMode::SCM )
                {
                    if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
                        std::cout << "No SCM DB available, do scm offline computations first...\n";
                    if( crb->scm()->doScmForMassMatrix() )
                        crb->scm()->setScmForMassMatrix( true );

                    crb->scm()->offline();
                }
            }

            if ( !crb->isDBLoaded() || crb->rebuildDB() )
            {
                if ( M_mode == CRBModelMode::CRB )
                    //|| M_mode == CRBModelMode::SCM )
                {
                    if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
                        std::cout << "No CRB DB available, do crb offline computations...\n";
                    crb->offline();
                }

                else if ( M_mode != CRBModelMode::SCM )
                    throw std::logic_error( "CRB/SCM Database could not be loaded" );
            }

            if( crb->isDBLoaded() )
            {
                bool do_offline = false;
                int current_dimension = crb->dimension();
                int dimension_max = this->vm()["crb.dimension-max"].template as<int>();
                int sampling_size = 0;
                if( use_predefined )
                    sampling_size = Sampling->readFromFile(file_name);

                if( sampling_size > current_dimension )
                    do_offline = true;

                if( current_dimension < dimension_max && !use_predefined )
                    do_offline=true;

                if( do_offline )
                    crb->offline();
            }
        }

    FEELPP_DONT_INLINE
    void run()
        {
            int proc_number =  Environment::worldComm().globalRank();

            if ( this->vm().count( "help" ) )
            {
                std::cout << this->optionsDescription() << "\n";
                return;
            }

            this->loadDB();
            int run_sampling_size = this->vm()[_o( this->about().appName(),"run.sampling.size" )].template as<int>();
            SamplingMode run_sampling_type = ( SamplingMode )this->vm()[_o( this->about().appName(),"run.sampling.mode" )].template as<int>();
            int output_index = this->vm()["crb.output-index"].template as<int>();
            //int output_index = this->vm()[_o(this->about().appName(),"output.index")].template as<int>();

            typename crb_type::sampling_ptrtype Sampling( new typename crb_type::sampling_type( model->parameterSpace() ) );

            switch ( run_sampling_type )
            {
            default:
            case SamplingMode::RANDOM:
                Sampling->randomize( run_sampling_size  );
                break;

            case SamplingMode::EQUIDISTRIBUTED:
                Sampling->equidistribute( run_sampling_size  );
                break;

            case SamplingMode::LOGEQUIDISTRIBUTED:
                Sampling->logEquidistribute( run_sampling_size  );
                break;
            }


            std::map<CRBModelMode,std::vector<std::string> > hdrs;
            using namespace boost::assign;
            std::vector<std::string> pfemhdrs = boost::assign::list_of( "FEM Output" )( "FEM Time" );
            std::vector<std::string> crbhdrs = boost::assign::list_of( "FEM Output" )( "FEM Time" )( "RB Output" )( "Error Bounds" )( "CRB Time" )( "Rel. error" )( "Conditionning" )( "l2_error" )( "h1_error" );
            std::vector<std::string> scmhdrs = boost::assign::list_of( "Lb" )( "Lb Time" )( "Ub" )( "Ub Time" )( "FEM" )( "FEM Time" )( "Rel.(FEM-Lb)" );
            std::vector<std::string> crbonlinehdrs = boost::assign::list_of( "RB Output" )( "Error Bounds" )( "CRB Time" );
            std::vector<std::string> scmonlinehdrs = boost::assign::list_of( "Lb" )( "Lb Time" )( "Ub" )( "Ub Time" )( "Rel.(FEM-Lb)" );
            hdrs[CRBModelMode::PFEM] = pfemhdrs;
            hdrs[CRBModelMode::CRB] = crbhdrs;
            hdrs[CRBModelMode::SCM] = scmhdrs;
            hdrs[CRBModelMode::CRB_ONLINE] = crbonlinehdrs;
            hdrs[CRBModelMode::SCM_ONLINE] = scmonlinehdrs;
            std::ostringstream ostr;

            //if( boost::is_same<  crbmodel_type , crbmodelbilinear_type >::value )
            {
                if( crb->printErrorDuringOfflineStep() )
                    crb->printErrorsDuringRbConstruction();
                if ( crb->showMuSelection() )
                    crb->printMuSelection();
            }

            auto exporter = Exporter<typename crbmodel_type::mesh_type>::New( "ensight" );
            exporter->step( 0 )->setMesh( model->functionSpace()->mesh() );

            printParameterHdr( ostr, model->parameterSpace()->dimension(), hdrs[M_mode] );

            std::ofstream file_summary_of_simulations( ( boost::format( "summary_of_simulations_%d" ) %crb->dimension() ).str().c_str() ,std::ios::out | std::ios::app );

            int curpar = 0;
            if( crb->useWNmu() )
                Sampling = crb->wnmu();

            /* Example of use of the setElements (but can use write in the file SamplingForTest)
            vector_parameter_type V;
            parameter_type UserMu( model->parameterSpace() );
            double j=0.1;
            //for(int i=1; i<101; i++)  { UserMu(0)=j;  UserMu(1)=1; UserMu(2)=1.5; UserMu(3)=2; UserMu(4)=3; UserMu(5)=4; UserMu(6)=4.5; UserMu(7)=5; UserMu(8)=6; V.push_back(UserMu ); j+=0.1;}
            //for(int i=10; i<100; i+=10)    { UserMu(0)=i;  UserMu(1)=1; V.push_back(UserMu );}
            //for(int i=1e2; i<1e3; i+=1e2)  { UserMu(0)=i;  UserMu(1)=1; V.push_back(UserMu );}
            //for(int i=1e3; i<1e4; i+=1e3)  { UserMu(0)=i;  UserMu(1)=1; V.push_back(UserMu );}
            //for(int i=1e4; i<1e5; i+=1e4)  { UserMu(0)=i;  UserMu(1)=1; V.push_back(UserMu );}
            //for(int i=1e5; i<1e6; i+=1e5)  { UserMu(0)=i;  UserMu(1)=1; V.push_back(UserMu );}
            //UserMu(0)=1e6;  UserMu(1)=1; V.push_back(UserMu );
            //Sampling->setElements( V );
            */

            /**
             * note that in the file SamplingForTest we expect :
             * mu_0= [ value0 , value1 , ... ]
             * mu_1= [ value0 , value1 , ... ]
             **/
            if( this->vm()["crb.use-predefined-test-sampling"].template as<bool>() )
            {
                std::string file_name = ( boost::format("SamplingForTest") ).str();
                std::ifstream file ( file_name );
                if( file  )
                {
                    Sampling->clear();
                    Sampling->readFromFile( file_name ) ;
                }
                else
                    throw std::logic_error( "[OpusApp] file SamplingForTest was not found" );

            }



            //Statistics
            vectorN_type l2_error_vector( Sampling->size() );
            vectorN_type h1_error_vector( Sampling->size() );
            vectorN_type relative_error_vector( Sampling->size() );
            vectorN_type time_fem_vector ( Sampling->size() );
            vectorN_type time_crb_vector ( Sampling->size() );
            vectorN_type relative_estimated_error_vector;
            if( crb->errorType()!=2 )
                relative_estimated_error_vector.resize( Sampling->size() );


            BOOST_FOREACH( auto mu, *Sampling )
            {

                int size = mu.size();
                if( proc_number == Environment::worldComm().masterRank() )
                {
                    std::cout << "(" << curpar++ << "/" << Sampling->size() << ") mu = [ ";
                    for ( int i=0; i<size-1; i++ ) std::cout<< mu[i] <<" , ";
                    std::cout<< mu[size-1]<<" ]\n ";
                }

                std::ostringstream mu_str;
                //if too many parameters, it will crash
                int sizemax=8;
                if( size < sizemax )
                    sizemax=size;
                for ( int i=0; i<sizemax-1; i++ ) mu_str << std::scientific << std::setprecision( 5 ) << mu[i] <<",";
                mu_str << std::scientific << std::setprecision( 5 ) << mu[size-1];


                LOG(INFO) << "mu=" << mu << "\n";
                mu.check();

                if( this->vm()["crb.script-mode"].template as<bool>() )
                {
                    unsigned long N = mu.size() + 5;
                    unsigned long P = 2;
                    std::vector<double> X( N );
                    std::vector<double> Y( P );
                    for(int i=0; i<mu.size(); i++)
                        X[i] = mu[i];

                    int N_dim = this->vm()["crb.dimension"].template as<int>();
                    int N_dimMax = this->vm()["crb.dimension-max"].template as<int>();
                    int Nwn;
                    if( N_dim > 0 )
                        Nwn = N_dim;
                    else
                        Nwn = N_dimMax;

                    X[N-5] = output_index;
                    X[N-4] = Nwn;
                    X[N-3] = this->vm()["crb.online-tolerance"].template as<double>();
                    X[N-2] = this->vm()["crb.error-type"].template as<int>();
                    //X[N-1] = this->vm()["crb.compute-variance"].template as<int>();
                    X[N-1] = 0;
                    bool compute_variance = this->vm()["crb.compute-variance"].template as<bool>();
                    if ( compute_variance )
                        X[N-1] = 1;

                    this->run( X.data(), X.size(), Y.data(), Y.size() );

                    std::cout << "output = " << Y[0] << std::endl;

                    std::ofstream res(this->vm()["result-file"].template as<std::string>() );
                    res << "output="<< Y[0] << "\n";
                }
                else
                {
                    switch ( M_mode )
                        {
                        case  CRBModelMode::PFEM:
                            {
                                LOG(INFO) << "PFEM mode" << std::endl;
                                if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
                                    std::cout << "PFEM mode" << std::endl;
                                boost::mpi::timer ti;

                                auto u_fem = model->solve( mu );
                                std::ostringstream u_fem_str;
                                u_fem_str << "u_fem(" << mu_str.str() << ")";
                                u_fem.setName( u_fem_str.str()  );

                                LOG(INFO) << "compute output\n";
                                google::FlushLogFiles(google::GLOG_INFO);

                                exporter->step(0)->add( u_fem.name(), u_fem );
                                //model->solve( mu );
                                std::vector<double> o = boost::assign::list_of( model->output( output_index,mu ) )( ti.elapsed() );
                                if(proc_number == Environment::worldComm().masterRank() ) std::cout << "output=" << o[0] << "\n";
                                printEntry( ostr, mu, o );

                                std::ofstream res(this->vm()["result-file"].template as<std::string>() );
                                res << "output="<< o[0] << "\n";

                            }
                            break;

                        case  CRBModelMode::CRB:
                            {
                                LOG(INFO) << "CRB mode\n";
                                if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
                                    std::cout << "CRB mode\n";


                                boost::mpi::timer ti;
                                //in the case we don't do the offline step, we need the affine decomposition
                                model->computeAffineDecomposition();

                                ti.restart();
                                LOG(INFO) << "solve crb\n";
                                google::FlushLogFiles(google::GLOG_INFO);

                                //dimension of the RB (not necessarily the max)
                                int N =  this->vm()["crb.dimension"].template as<int>();

                                auto o = crb->run( mu,  this->vm()["crb.online-tolerance"].template as<double>() , N);
                                double time_crb = ti.elapsed();

                                //auto u_crb = crb->expansion( mu , N );
                                auto uN_0 = o.template get<4>();
                                element_type u_crb;
                                if( model->isSteady()) // Re-use uN given by lb in crb->run
                                    u_crb = crb->expansion( uN_0 , N ); // Re-use uN given by lb in crb->run
                                else
                                    u_crb = crb->expansion( mu , N );

                                std::ostringstream u_crb_str;
                                u_crb_str << "u_crb(" << mu_str.str() << ")";
                                u_crb.setName( u_crb_str.str()  );
                                LOG(INFO) << "export u_crb \n";
                                exporter->step(0)->add( u_crb.name(), u_crb );

                                double relative_error = -1;
                                double relative_estimated_error = -1;
                                double condition_number = o.template get<3>();
                                double l2_error = -1;
                                double h1_error = -1;
                                double output_fem = -1;
                                double time_fem = -1;

                                bool compute_fem = this->vm()["crb.compute-fem-during-online"].template as<bool>();
                                element_type u_fem;

                                if ( compute_fem )
                                {
                                    bool use_newton = this->vm()["crb.use-newton"].template as<bool>();
                                    ti.restart();
                                    LOG(INFO) << "solve u_fem\n";
                                    google::FlushLogFiles(google::GLOG_INFO);

                                    //auto u_fem = model->solveRB( mu );
                                    //auto u_fem = model->solveFemUsingOfflineEim( mu );

                                    //TODO : add use-newton condition
                                    if( boost::is_same<  crbmodel_type , crbmodelbilinear_type >::value && ! use_newton )
                                        u_fem = model->solveFemUsingOnlineEimPicard( mu );
                                    else
                                        u_fem = model->solve( mu );

                                    std::ostringstream u_fem_str;
                                    u_fem_str << "u_fem(" << mu_str.str() << ")";
                                    u_fem.setName( u_fem_str.str()  );

                                    LOG(INFO) << "compute output\n";
                                    google::FlushLogFiles(google::GLOG_INFO);

                                    LOG(INFO) << "export u_fem \n";
                                    exporter->step(0)->add( u_fem.name(), u_fem );

                                    std::vector<double> ofem = boost::assign::list_of( model->output( output_index,mu ) )( ti.elapsed() );


                                    relative_error = std::abs( ofem[0]-o.template get<0>() ) /ofem[0];
                                    relative_estimated_error = o.template get<1>() / ofem[0];

                                    //compute || u_fem - u_crb||_L2
                                    LOG(INFO) << "compute error \n";
                                    auto u_error = model->functionSpace()->element();
                                    std::ostringstream u_error_str;
                                    u_error = (( u_fem - u_crb ).pow(2)).sqrt()  ;
                                    u_error_str << "u_error(" << mu_str.str() << ")";
                                    u_error.setName( u_error_str.str()  );
                                    exporter->step(0)->add( u_error.name(), u_error );
                                    LOG(INFO) << "L2(fem)=" << l2Norm( u_fem )    << "\n";
                                    LOG(INFO) << "H1(fem)=" << h1Norm( u_fem )    << "\n";
                                    l2_error = l2Norm( u_error )/l2Norm( u_fem );
                                    h1_error = h1Norm( u_error )/h1Norm( u_fem );

                                    output_fem = ofem[0];
                                    time_fem = ofem[1];

                                }//compute-fem-during-online

                                if ( crb->errorType()==2 )
                                {
                                    std::vector<double> v = boost::assign::list_of( output_fem )( time_fem )( o.template get<0>() )( relative_estimated_error )( time_crb )( relative_error )( condition_number )( l2_error )( h1_error );
                                    if( proc_number == Environment::worldComm().masterRank() )
                                    {
                                        std::cout << "output=" << o.template get<0>() << " with " << o.template get<2>() << " basis functions\n";
                                        printEntry( file_summary_of_simulations, mu, v );
                                        printEntry( ostr, mu, v );
                                        //file_summary_of_simulations.close();

                                        if ( this->vm()["crb.compute-stat"].template as<bool>() && compute_fem )
                                        {
                                            relative_error_vector[curpar-1] = relative_error;
                                            l2_error_vector[curpar-1] = l2_error;
                                            h1_error_vector[curpar-1] = h1_error;
                                            time_fem_vector[curpar-1] = time_fem;
                                            time_crb_vector[curpar-1] = time_crb;
                                        }

                                        std::ofstream res(this->vm()["result-file"].template as<std::string>() );
                                        res << "output="<< o.template get<0>() << "\n";

                                    }

                                }//end of crb->errorType==2
                                else
                                {
                                    //if( ! boost::is_same<  crbmodel_type , crbmodelbilinear_type >::value )
                                    //    throw std::logic_error( "ERROR TYPE must be 2 when using CRBTrilinear (no error estimation)" );
                                    std::vector<double> v = boost::assign::list_of( output_fem )( time_fem )( o.template get<0>() )( relative_estimated_error )( ti.elapsed() ) ( relative_error )( condition_number )( l2_error )( h1_error ) ;
                                    if( proc_number == Environment::worldComm().masterRank() )
                                    {
                                        std::cout << "output=" << o.template get<0>() << " with " << o.template get<2>() << " basis functions  (relative error estimation on this output : " << relative_estimated_error<<") \n";
                                        //std::ofstream file_summary_of_simulations( ( boost::format( "summary_of_simulations_%d" ) % o.template get<2>() ).str().c_str() ,std::ios::out | std::ios::app );
                                        printEntry( file_summary_of_simulations, mu, v );
                                        printEntry( ostr, mu, v );
                                        //file_summary_of_simulations.close();

                                        if ( this->vm()["crb.compute-stat"].template as<bool>() && compute_fem )
                                        {
                                            relative_error_vector[curpar-1] = relative_error;
                                            l2_error_vector[curpar-1] = l2_error;
                                            h1_error_vector[curpar-1] = h1_error;
                                            time_fem_vector[curpar-1] = time_fem;
                                            time_crb_vector[curpar-1] = time_crb;
                                            relative_estimated_error_vector[curpar-1] = relative_estimated_error;
                                        }
                                        std::ofstream res(this->vm()["result-file"].template as<std::string>() );
                                        res << "output="<< o.template get<0>() << "\n";
                                    }//end of proc==master
                                }//end of else (errorType==2)
                                if (this->vm()["crb.cvg-study"].template as<bool>() && compute_fem )
                                {
                                    LOG(INFO) << "start convergence study...\n";
                                    std::map<int, boost::tuple<double,double,double,double> > conver;
                                    for( int N = 1; N < crb->dimension(); N++ )
                                    {
                                        LOG(INFO) << "N=" << N << "...\n";
                                        auto o = crb->run( mu,  this->vm()["crb.online-tolerance"].template as<double>() , N);
                                        auto u_crb = crb->expansion( mu , N );
                                        auto u_error = model->functionSpace()->element();
                                        u_error = u_fem - u_crb;
                                        double rel_err = std::abs( output_fem-o.template get<0>() ) /output_fem;
                                        double l2_error = l2Norm( u_error )/l2Norm( u_fem );
                                        double h1_error = h1Norm( u_error )/h1Norm( u_fem );
                                        double condition_number = o.template get<3>();
                                        conver[N]=boost::make_tuple( rel_err, l2_error, h1_error , condition_number );
                                        LOG(INFO) << "N=" << N << " " << rel_err << " " << l2_error << " " << h1_error << " " <<condition_number<<"\n";
                                        if ( proc_number == Environment::worldComm().masterRank() )
                                            std::cout << "N=" << N << " " << rel_err << " " << l2_error << " " << h1_error << " " <<condition_number<<std::endl;
                                        LOG(INFO) << "N=" << N << " done.\n";
                                    }
                                    if( proc_number == Environment::worldComm().masterRank() )
                                    {
                                        LOG(INFO) << "save in logfile\n";
                                        std::string mu_str;
                                        for ( int i=0; i<mu.size(); i++ )
                                            mu_str= mu_str + ( boost::format( "_%1%" ) %mu[i] ).str() ;
                                        std::string file_name = "convergence"+mu_str+".dat";
                                        std::ofstream conv( file_name );
                                        BOOST_FOREACH( auto en, conver )
                                            conv << en.first << "\t" << en.second.get<0>()  << "\t" << en.second.get<1>() << "\t" << en.second.get<2>() << "\t"<< en.second.get<3>() << "\n";
                                    }
                                }//end of cvg-study
                            }//case CRB
                            break;
                        case  CRBModelMode::CRB_ONLINE:
                            {
                                std::cout << "CRB Online mode\n";
                                boost::mpi::timer ti;
                                ti.restart();
                                auto o = crb->run( mu,  this->vm()["crb.online-tolerance"].template as<double>() );

                                if ( crb->errorType()==2 )
                                    {
                                        std::vector<double> v = boost::assign::list_of( o.template get<0>() )( ti.elapsed() );
                                        std::cout << "output=" << o.template get<0>() << " with " << o.template get<2>() << " basis functions\n";
                                        printEntry( ostr, mu, v );
                                    }

                                else
                                    {
                                        std::vector<double> v = boost::assign::list_of( o.template get<0>() )( o.template get<1>() )( ti.elapsed() );
                                        std::cout << "output=" << o.template get<0>() << " with " << o.template get<2>() << " basis functions  (error estimation on this output : " << o.template get<1>()<<") \n";
                                        printEntry( ostr, mu, v );
                                    }

                            }
                            break;

                        case  CRBModelMode::SCM:
                            {
                                std::cout << "SCM mode\n";
                                auto o = crb->scm()->run( mu, crb->scm()->KMax() );
                                printEntry( ostr, mu, o );
                            }
                            break;

                        case  CRBModelMode::SCM_ONLINE:
                            {
                                std::cout << "SCM online mode\n";
                                auto o = crb->scm()->run( mu, crb->scm()->KMax() );
                                printEntry( ostr, mu, o );
                            }
                            break;

                        }

                    std::cout << "------------------------------------------------------------\n";
                }
            }
            exporter->save();
            if( proc_number == Environment::worldComm().masterRank() ) std::cout << ostr.str() << "\n";

            bool compute_fem = this->vm()["crb.compute-fem-during-online"].template as<bool>();
            if ( this->vm()["crb.compute-stat"].template as<bool>() && compute_fem )
            {
                LOG( INFO ) << "compute statistics \n";
                Eigen::MatrixXf::Index index_max_l2;
                Eigen::MatrixXf::Index index_min_l2;
                Eigen::MatrixXf::Index index_max_h1;
                Eigen::MatrixXf::Index index_min_h1;
                Eigen::MatrixXf::Index index_max_time_fem;
                Eigen::MatrixXf::Index index_min_time_fem;
                Eigen::MatrixXf::Index index_max_time_crb;
                Eigen::MatrixXf::Index index_min_time_crb;
                Eigen::MatrixXf::Index index_max_estimated_error;
                Eigen::MatrixXf::Index index_min_estimated_error;
                Eigen::MatrixXf::Index index_max_output_error;
                Eigen::MatrixXf::Index index_min_output_error;

                double max_l2 = l2_error_vector.maxCoeff(&index_max_l2);
                double min_l2 = l2_error_vector.minCoeff(&index_min_l2);
                double mean_l2 = l2_error_vector.mean();
                double max_h1 = h1_error_vector.maxCoeff(&index_max_h1);
                double min_h1 = h1_error_vector.minCoeff(&index_min_h1);
                double mean_h1 = h1_error_vector.mean();
                double max_time_fem = time_fem_vector.maxCoeff(&index_max_time_fem);
                double min_time_fem = time_fem_vector.minCoeff(&index_min_time_fem);
                double mean_time_fem = time_fem_vector.mean();
                double max_time_crb = time_crb_vector.maxCoeff(&index_max_time_crb);
                double min_time_crb = time_crb_vector.minCoeff(&index_min_time_crb);
                double mean_time_crb = time_crb_vector.mean();
                double max_output_error = relative_error_vector.maxCoeff(&index_max_output_error);
                double min_output_error = relative_error_vector.minCoeff(&index_min_output_error);
                double mean_output_error = relative_error_vector.mean();
                double max_estimated_error = 0;
                double min_estimated_error = 0;
                double mean_estimated_error = 0;

                if( crb->errorType()!=2 )
                {
                    max_estimated_error = relative_estimated_error_vector.maxCoeff(&index_max_estimated_error);
                    min_estimated_error = relative_estimated_error_vector.minCoeff(&index_min_estimated_error);
                    mean_estimated_error = relative_estimated_error_vector.mean();
                }
                if( proc_number == Environment::worldComm().masterRank() )
                {
                    file_summary_of_simulations <<"\n\nStatistics\n";
                    file_summary_of_simulations <<"max of output error : "<<max_output_error<<" at the "<<index_max_output_error+1<<"^th simulation\n";
                    file_summary_of_simulations <<"min of output error : "<<min_output_error<<" at the "<<index_min_output_error+1<<"^th simulation\n";
                    file_summary_of_simulations <<"mean of output error : "<<mean_output_error<<"\n\n";
                    file_summary_of_simulations <<"max of estimated output error : "<<max_estimated_error<<" at the "<<index_max_estimated_error+1<<"^th simulation\n";
                    file_summary_of_simulations <<"min of estimated output error : "<<min_estimated_error<<" at the "<<index_min_estimated_error+1<<"^th simulation\n";
                    file_summary_of_simulations <<"mean of estimated output error : "<<mean_estimated_error<<"\n\n";
                    file_summary_of_simulations <<"max of L2 error : "<<max_l2<<" at the "<<index_max_l2+1<<"^th simulation\n";
                    file_summary_of_simulations <<"min of L2 error : "<<min_l2<<" at the "<<index_min_l2+1<<"^th simulation\n";
                    file_summary_of_simulations <<"mean of L2 error : "<<mean_l2<<"\n\n";
                    file_summary_of_simulations <<"max of H1 error : "<<max_h1<<" at the "<<index_max_h1+1<<"^th simulation\n";
                    file_summary_of_simulations <<"min of H1 error : "<<min_h1<<" at the "<<index_min_h1+1<<"^th simulation\n";
                    file_summary_of_simulations <<"mean of H1 error : "<<mean_h1<<"\n\n";
                    file_summary_of_simulations <<"max of time FEM : "<<max_time_fem<<" at the "<<index_max_time_fem+1<<"^th simulation\n";
                    file_summary_of_simulations <<"min of time FEM : "<<min_time_fem<<" at the "<<index_min_time_fem+1<<"^th simulation\n";
                    file_summary_of_simulations <<"mean of time FEM : "<<mean_time_fem<<"\n\n";
                    file_summary_of_simulations <<"max of time CRB : "<<max_time_crb<<" at the "<<index_max_time_crb+1<<"^th simulation\n";
                    file_summary_of_simulations <<"min of time CRB : "<<min_time_crb<<" at the "<<index_min_time_crb+1<<"^th simulation\n";
                    file_summary_of_simulations <<"mean of time CRB : "<<mean_time_crb<<"\n\n";
                }
            }//end of compute-stat
        }
    void run( const double * X, unsigned long N,
              double * Y, unsigned long P )
        {

            switch ( M_mode )
            {
            case  CRBModelMode::PFEM:
            {
                model->run( X, N, Y, P );
            }
            break;

            case  CRBModelMode::SCM:
            case  CRBModelMode::SCM_ONLINE:
            {
                crb->scm()->run( X, N, Y, P );
            }
            break;

            case  CRBModelMode::CRB:
            case  CRBModelMode::CRB_ONLINE:
            {
                crb->run( X, N, Y, P );
            }
            break;
            }

            fs::current_path( M_current_path );
        }

private:
    int printParameterHdr( std::ostream& os, int N, std::vector<std::string> outputhdrs )
        {
            for ( int i = 0; i < N; ++i )
            {
                std::ostringstream s;
                s << "mu" << i;
                os  << hdrmanip( prec+7 ) << s.str();
            }

            BOOST_FOREACH( auto output, outputhdrs )
            {
                os << hdrmanip( 15 ) << output;
            }
            os << "\n";

            return N*( prec+7 )+outputhdrs.size()*15;
        }
    void printEntry( std::ostream& os,
                     typename ModelType::parameter_type const& mu,
                     std::vector<double> const& outputs )
        {
            for ( int i = 0; i < mu.size(); ++i )
                os  << std::right <<std::setw( prec+7 ) << dmanip << mu[i];

            BOOST_FOREACH( auto o, outputs )
            {
                os << tabmanip( 15 ) << o;
            }
            os << "\n";
        }


    //double l2Norm( typename ModelType::parameter_type const& mu, int N )
    double l2Norm( element_type const& u )
    {
        static const bool is_composite = functionspace_type::is_composite;
        return l2Norm( u, mpl::bool_< is_composite >() );
    }
    double l2Norm( element_type const& u, mpl::bool_<false> )
    {
        auto mesh = model->functionSpace()->mesh();
        return math::sqrt( integrate( elements(mesh), (vf::idv(u))*(vf::idv(u)) ).evaluate()(0,0) );
    }
    double l2Norm( element_type const& u, mpl::bool_<true>)
    {
        auto mesh = model->functionSpace()->mesh();
        auto uT = u.template element<1>();
        return math::sqrt( integrate( elements(mesh), (vf::idv(uT))*(vf::idv(uT)) ).evaluate()(0,0) );
    }
    //double h1Norm( typename ModelType::parameter_type const& mu, int N )
    double h1Norm( element_type const& u )
    {
        static const bool is_composite = functionspace_type::is_composite;
        return h1Norm( u, mpl::bool_< is_composite >() );
    }
    double h1Norm( element_type const& u, mpl::bool_<false> )
    {
        auto mesh = model->functionSpace()->mesh();
        double l22 = integrate( elements(mesh), (vf::idv(u))*(vf::idv(u)) ).evaluate()(0,0);
        double semih12 = integrate( elements(mesh), (vf::gradv(u))*trans(vf::gradv(u)) ).evaluate()(0,0);
        return math::sqrt( l22+semih12 );
    }
    double h1Norm( element_type const& u, mpl::bool_<true>)
    {
        auto mesh = model->functionSpace()->mesh();
        auto u_femT = u.template element<1>();
        double l22 = integrate( elements(mesh), (vf::idv(u_femT))*(vf::idv(u_femT)) ).evaluate()(0,0);
        double semih12 = integrate( elements(mesh), (vf::gradv(u_femT))*trans(vf::gradv(u_femT))).evaluate()(0,0);
        return math::sqrt( l22+semih12 );
    }

private:
    CRBModelMode M_mode;
    crbmodel_ptrtype model;
    crb_ptrtype crb;

    fs::path M_current_path;
}; // OpusApp

} // Feel

#endif /* __OpusApp_H */


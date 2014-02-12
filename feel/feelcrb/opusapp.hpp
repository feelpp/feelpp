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
#include <boost/range/join.hpp>
#include <boost/regex.hpp>
#include <boost/assign/list_of.hpp>

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

    typedef typename crb_type::sampling_ptrtype sampling_ptrtype;

    OpusApp()
        :
        super(),
        M_mode( ( CRBModelMode )option(_name=_o( this->about().appName(),"run.mode" )).template as<int>() )
        {
            this->init();
        }

    OpusApp( AboutData const& ad, po::options_description const& od )
        :
        super( ad, opusapp_options( ad.appName() ).add( od ).add( crbOptions() ).add( feel_options() ).add( eimOptions() ).add( podOptions() ) ),
        M_mode( ( CRBModelMode )option(_name=_o( this->about().appName(),"run.mode" )).template as<int>() )
        {
            this->init();
        }

    OpusApp( AboutData const& ad, po::options_description const& od, CRBModelMode mode )
        :
        super( ad, opusapp_options( ad.appName() ).add( od ).add( crbOptions() ).add( feel_options() ).add( eimOptions() ).add( podOptions() ) ),
        M_mode( mode )
        {
            this->init();
        }

    OpusApp( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( ad, opusapp_options( ad.appName() ).add( od ).add( crbOptions() ).add( feel_options() ).add( eimOptions() ).add( podOptions() )  ),
        M_mode( ( CRBModelMode )option(_name=_o( this->about().appName(),"run.mode" )).template as<int>() )
        {
            this->init();
        }
    OpusApp( int argc, char** argv, AboutData const& ad, po::options_description const& od, CRBModelMode mode )
        :
        super( ad, opusapp_options( ad.appName() ).add( od ).add( crbOptions() ).add( feel_options() ).add( eimOptions() ).add( podOptions() ) ),
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
                if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
                    std::cout << this->about().appName() << std::endl;
                LOG(INFO) << "[OpusApp] constructor " << this->about().appName()  << "\n";

                // Check if user have given a name for result files repo
                // Note : this name is also use for database storage
                std::string results_repo_name;
                if( this->vm().count("crb.results-repo-name") )
                    results_repo_name = option(_name="crb.results-repo-name").template as<std::string>();
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
            int proc_number = Environment::worldComm().globalRank();
            int global_size = Environment::worldComm().globalSize();
            std::string pslogfile = ( boost::format("PsLogCrbOffline-%1%_%2%") %global_size %proc_number ).str();

            bool only_master=option(_name="crb.system-memory-evolution").template as<bool>();
            bool all_procs  =option(_name="crb.system-memory-evolution-on-all-procs").template as<bool>();
            bool only_one_proc=only_master * ( Environment::worldComm().globalRank()==Environment::worldComm().masterRank() );
            bool write_memory_evolution = all_procs || only_one_proc ;

            bool crb_use_predefined = option(_name="crb.use-predefined-WNmu").template as<bool>();
            std::string file_name;
            int NlogEquidistributed = option(_name="crb.use-logEquidistributed-WNmu").template as<int>();
            int Nequidistributed = option(_name="crb.use-equidistributed-WNmu").template as<int>();
            int NlogEquidistributedScm = option(_name="crb.scm.use-logEquidistributed-C").template as<int>();
            int NequidistributedScm = option(_name="crb.scm.use-equidistributed-C").template as<int>();
            typename crb_type::sampling_ptrtype Sampling( new typename crb_type::sampling_type( model->parameterSpace() ) );
            bool all_proc_have_same_sampling=true;

            if( crb_use_predefined )
            {
                file_name = ( boost::format("SamplingWNmu") ).str();
            }
            if( NlogEquidistributed+Nequidistributed > 0 )
            {
                file_name = ( boost::format("SamplingWNmu") ).str();
                if( NlogEquidistributed > 0 )
                    Sampling->logEquidistribute( NlogEquidistributed,all_proc_have_same_sampling );
                if( Nequidistributed > 0 )
                    Sampling->equidistribute( Nequidistributed,all_proc_have_same_sampling );
                Sampling->writeOnFile(file_name);
            }
            if( NlogEquidistributedScm+NequidistributedScm > 0 )
            {
                file_name = ( boost::format("SamplingC") ).str();
                if( NlogEquidistributedScm > 0 )
                    Sampling->logEquidistribute( NlogEquidistributedScm,all_proc_have_same_sampling );
                if( NequidistributedScm > 0 )
                    Sampling->equidistribute( NequidistributedScm,all_proc_have_same_sampling );
                Sampling->writeOnFile(file_name);
            }

            if ( M_mode == CRBModelMode::PFEM )
                return;

            if ( !crb->scm()->isDBLoaded() || crb->scm()->rebuildDB() )
            {
                if ( M_mode == CRBModelMode::SCM )
                {
                    if( proc_number == Environment::worldComm().masterRank() )
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
                    if( proc_number == Environment::worldComm().masterRank() )
                        std::cout << "No CRB DB available, do crb offline computations...\n";
                    crb->setOfflineStep( true );
                    crb->offline();
                    if( write_memory_evolution )
                        this->generateMemoryEvolution(pslogfile);
                }

                else if ( M_mode != CRBModelMode::SCM )
                    throw std::logic_error( "CRB/SCM Database could not be loaded" );
            }

            if( crb->isDBLoaded() )
            {
                int Nrestart = option(_name="crb.restart-from-N").template as<int>();
                bool do_offline = false;
                int current_dimension = crb->dimension();
                int dimension_max = option(_name="crb.dimension-max").template as<int>();
                int sampling_size = 0;
                if( crb_use_predefined )
                    sampling_size = Sampling->readFromFile(file_name);

                if( sampling_size > current_dimension )
                    do_offline = true;

                if( current_dimension < dimension_max && !crb_use_predefined )
                    do_offline=true;

                if( Nrestart > 1 )
                    do_offline=true;

                if( ! do_offline )
                {
                    crb->loadSCMDB();
                }

                if( do_offline )
                {
                    crb->setOfflineStep( true );
                    crb->offline();
                    if( write_memory_evolution )
                        this->generateMemoryEvolution(pslogfile);
                }
            }
        }


    FEELPP_DONT_INLINE
    void run()
        {

            bool export_solution = option(_name=_o( this->about().appName(),"export-solution" )).template as<bool>();
            int proc_number =  Environment::worldComm().globalRank();

            if ( this->vm().count( "help" ) )
            {
                std::cout << this->optionsDescription() << "\n";
                return;
            }

            this->loadDB();
            int run_sampling_size = option(_name=_o( this->about().appName(),"run.sampling.size" )).template as<int>();
            SamplingMode run_sampling_type = ( SamplingMode )option(_name=_o( this->about().appName(),"run.sampling.mode" )).template as<int>();
            int output_index = option(_name="crb.output-index").template as<int>();
            //int output_index = option(_name=_o(this->about().appName(),"output.index")).template as<int>();

            typename crb_type::sampling_ptrtype Sampling( new typename crb_type::sampling_type( model->parameterSpace() ) );

            int n_eval_computational_time = option(_name="eim.computational-time-neval").template as<int>();
            bool compute_fem = option(_name="crb.compute-fem-during-online").template as<bool>();
            bool compute_stat =  option(_name="crb.compute-stat").template as<bool>();

            bool use_predefined_sampling = option(_name="crb.use-predefined-test-sampling").template as<bool>();
            bool select_parameter_via_one_feel=option( _name="crb.select-parameter-via-one-feel").template as<bool>();
            bool sampling_is_already_generated=false;
            if( select_parameter_via_one_feel )
            {
                run_sampling_size=1;
                Sampling->clear();
                //in this case we want to visualize RB solution with parameters from one feel interface
                compute_fem=false;
                compute_stat=false;

                //parameters are given by a vector of double
                std::string string_parameters = option(_name="crb.user-parameters").template as< std::string >();
                std::vector< std::string > str;
                boost::split( str, string_parameters, boost::is_any_of(" "), boost::token_compress_on );
                parameter_type user_mu ( model->parameterSpace() );
                double user_parameter_size = str.size();
                double mu_size = user_mu.size();
                CHECK( user_parameter_size == mu_size )<<"[OpusApp] Error : parameters must have "<<mu_size<<" components and "<<user_parameter_size<<" have been given by the user \n";
                for(int i=0; i<mu_size; i++)
                {
                    double mu = boost::lexical_cast<double>( str[i] );
                    user_mu( i ) = mu;
                }
                Sampling->addElement( user_mu );
                sampling_is_already_generated=true;
            }


            std::string vary_only_parameter_components = option(_name="crb.vary-only-parameter-components").template as<std::string>();
            std::vector< std::string > str;
            boost::split( str, vary_only_parameter_components, boost::is_any_of(" "), boost::token_compress_on );
            int number_str=str.size();
            CHECK( number_str < 4 )<<"Error when using option crb.vary-only-parameter-components, at maximum we can vary 2 components of the parameter";
            int vary_mu_comp0=-1,vary_mu_comp1=-1;
            if( number_str > 1 )
            {
                Sampling->clear();
                compute_fem=false;
                compute_stat=false;
                export_solution=false;
                int size=-1;
                //here only one component vary
                if( number_str == 2 )
                {
                    vary_mu_comp0 = boost::lexical_cast<int>( str[0] );
                    size = boost::lexical_cast<int>( str[1] );
                    run_sampling_size=size;
                }
                else
                {
                    vary_mu_comp0 = boost::lexical_cast<int>( str[0] );
                    vary_mu_comp1 = boost::lexical_cast<int>( str[1] );
                    size = boost::lexical_cast<int>( str[2] );
                    run_sampling_size=size*size;
                }

                parameter_type user_mu ( model->parameterSpace() );
                double mu_size = user_mu.size();
                CHECK( vary_mu_comp0 < mu_size )<<"[OpusApp] error using crb.vary-only-parameter-components, the component "<<vary_mu_comp0<<" can't vary because parameter have a total of only "<<mu_size<<" components\n";
                if( number_str == 3 )
                {
                    CHECK( vary_mu_comp1 < mu_size )<<"[OpusApp] error using crb.vary-only-parameter-components, the component "<<vary_mu_comp1<<" can't vary because parameter have a total of only "<<mu_size<<" components\n";
                }

                std::vector< int > sampling_each_direction ( mu_size );
                for(int i=0; i<mu_size; i++)
                {
                    if( i == vary_mu_comp0 || i == vary_mu_comp1 )
                        sampling_each_direction[i]=size;
                    else
                        sampling_each_direction[i]=0;
                }
                Sampling->logEquidistributeProduct( sampling_each_direction );
                sampling_is_already_generated=true;
            }

            if( n_eval_computational_time > 0 )
            {
                compute_fem = false;
                auto eim_sc_vector = model->scalarContinuousEim();
                auto eim_sd_vector = model->scalarDiscontinuousEim();
                int size1 = eim_sc_vector.size();
                int size2 = eim_sd_vector.size();
                if( size1 + size2 == 0 )
                    throw std::logic_error( "[OpusApp] no eim object detected" );

                std::string appname = this->about().appName();
                for(int i=0; i<size1; i++)
                    eim_sc_vector[i]->computationalTimeStatistics(appname);
                for(int i=0; i<size2; i++)
                    eim_sd_vector[i]->computationalTimeStatistics(appname);

                run_sampling_size = 0;
            }
            n_eval_computational_time = option(_name="crb.computational-time-neval").template as<int>();
            if( n_eval_computational_time > 0 )
            {
                if( ! option(_name="crb.cvg-study").template as<bool>() )
                {
                    compute_fem = false;
                    run_sampling_size = 0;
                }
                std::string appname = this->about().appName();
                crb->computationalTimeStatistics( appname );
            }


            //here we can be interested by computing FEM and CRB solutions
            //so it is important that every proc has the same sampling (for FEM solution)
            bool all_proc_have_same_sampling=true;
            if( ! sampling_is_already_generated )
            {
                switch ( run_sampling_type )
                {
                default:
                case SamplingMode::RANDOM:
                    Sampling->randomize( run_sampling_size , all_proc_have_same_sampling );
                    break;

                case SamplingMode::EQUIDISTRIBUTED:
                    Sampling->equidistribute( run_sampling_size , all_proc_have_same_sampling );
                    break;

                case SamplingMode::LOGEQUIDISTRIBUTED:
                    Sampling->logEquidistribute( run_sampling_size , all_proc_have_same_sampling );
                    break;
                }
            }// ! select_parameter_via_one_feel


            std::map<CRBModelMode,std::vector<std::string> > hdrs;
            using namespace boost::assign;
            std::vector<std::string> pfemhdrs = boost::assign::list_of( "FEM Output" )( "FEM Time" );
            std::vector<std::string> crbhdrs = boost::assign::list_of( "FEM Output" )( "FEM Time" )( "RB Output" )( "Error Bounds" )( "CRB Time" )( "output error" )( "Conditionning" )( "l2_error" )( "h1_error" );
            std::vector<std::string> scmhdrs = boost::assign::list_of( "Lb" )( "Lb Time" )( "Ub" )( "Ub Time" )( "FEM" )( "FEM Time" )( "output error" );
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
                if ( crb->showMuSelection() && Environment::worldComm().globalRank()==Environment::worldComm().masterRank() )
                    crb->printMuSelection();
            }

            auto e = exporter( _mesh= model->functionSpace()->mesh()  );

            printParameterHdr( ostr, model->parameterSpace()->dimension(), hdrs[M_mode] );

            int crb_error_type = option(_name="crb.error-type").template as<int>();

            int dim=0;
            if( M_mode==CRBModelMode::CRB )
            {
                dim=crb->dimension();
                if( crb->useWNmu() )
                    Sampling = crb->wnmu();

                if( option(_name="crb.run-on-scm-parameters").template as<bool>() )
                {
                    Sampling = crb->scm()->c();
                    if( crb_error_type!=1 )
                        throw std::logic_error( "[OpusApp] The SCM has not been launched, you can't use the option crb.run-on-scm-parameters. Run the SCM ( option crb.error-type=1 ) or comment this option line." );
                }
            }
            if( M_mode==CRBModelMode::SCM )
            {
                dim=crb->scm()->KMax();
                if( option(_name="crb.scm.run-on-C").template as<bool>() )
                    Sampling = crb->scm()->c();
            }

            std::ofstream file_summary_of_simulations( ( boost::format( "summary_of_simulations_%d" ) %dim ).str().c_str() ,std::ios::out | std::ios::app );

            vectorN_type outputs_storage( run_sampling_size );
            vectorN_type mu0_storage( run_sampling_size );
            vectorN_type estimated_error_outputs_storage( run_sampling_size );

            int curpar = 0;

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

            // Script write current mu in cfg => need to write in SamplingForTest
            if( option(_name="crb.script-mode").template as<bool>() )
                {
                    // Sampling will be the parameter given by OT
                    if( proc_number == Environment::worldComm().masterRank() )
                            buildSamplingFromCfg();

                    Environment::worldComm().barrier();

                    compute_fem = false;
                    compute_stat = false;
                }

            /**
             * note that in the file SamplingForTest we expect :
             * mu_0= [ value0 , value1 , ... ]
             * mu_1= [ value0 , value1 , ... ]
             **/
            if(  use_predefined_sampling || option(_name="crb.script-mode").template as<bool>() )
            {
                std::string file_name = ( boost::format("SamplingForTest") ).str();
                std::ifstream file ( file_name );
                if( file )
                {
                    Sampling->clear();
                    Sampling->readFromFile( file_name ) ;
                }
                else
                {
                    VLOG(2) << "proc number : " << proc_number << "can't find file \n";
                    throw std::logic_error( "[OpusApp] file SamplingForTest was not found" );
                }
            }

            //Statistics
            vectorN_type l2_error_vector;
            vectorN_type h1_error_vector;
            vectorN_type relative_error_vector;
            vectorN_type time_fem_vector;
            vectorN_type time_crb_vector;
            vectorN_type relative_estimated_error_vector;

            vectorN_type scm_relative_error;

            bool solve_dual_problem = option(_name="crb.solve-dual-problem").template as<bool>();

            if (option(_name="crb.cvg-study").template as<bool>() && compute_fem )
            {
                //int Nmax = crb->dimension();
                int Nmax = option("crb.dimension-max").template as<int>();
                vector_sampling_for_primal_efficiency_under_1.resize(Nmax);
                for(int N=0; N<Nmax; N++)
                {
                    sampling_ptrtype sampling_primal ( new typename crb_type::sampling_type( model->parameterSpace() ) );
                    vector_sampling_for_primal_efficiency_under_1[N]=sampling_primal;
                }
                if( solve_dual_problem )
                {
                    vector_sampling_for_dual_efficiency_under_1.resize(Nmax);
                    for(int N=0; N<Nmax; N++)
                    {
                        sampling_ptrtype sampling_dual ( new typename crb_type::sampling_type( model->parameterSpace() ) );
                        vector_sampling_for_dual_efficiency_under_1[N]=sampling_dual;
                    }
                }
            }

            if( M_mode==CRBModelMode::CRB )
            {

                l2_error_vector.resize( Sampling->size() );
                h1_error_vector.resize( Sampling->size() );
                relative_error_vector.resize( Sampling->size() );
                time_fem_vector.resize( Sampling->size() );
                time_crb_vector.resize( Sampling->size() );

                if( crb->errorType()!=2 )
                    relative_estimated_error_vector.resize( Sampling->size() );

                crb->setOfflineStep( false );

                if (option(_name="eim.cvg-study").template as<bool>() )
                {
                    compute_fem=false;
                }

            }
            if( M_mode==CRBModelMode::SCM )
            {
                if (option(_name="crb.scm.cvg-study").template as<bool>() )
                    this->initializeConvergenceScmMap( Sampling->size() );

                scm_relative_error.resize( Sampling->size() );
            }

            int crb_dimension = option(_name="crb.dimension").template as<int>();
            int crb_dimension_max = option(_name="crb.dimension-max").template as<int>();
            double crb_online_tolerance = option(_name="crb.online-tolerance").template as<double>();
            bool crb_compute_variance  = option(_name="crb.compute-variance").template as<bool>();

            double output_fem = -1;


            //in the case we don't do the offline step, we need the affine decomposition
            model->computeAffineDecomposition();

            //compute beta coeff for reference parameters
            auto ref_mu = model->refParameter();
            double dt = model->timeStep();
            double ti = model->timeInitial();
            double tf = model->timeFinal();
            int K = ( tf - ti )/dt;
            std::vector< std::vector< std::vector< double > > > ref_betaAqm;
            for(int time_index=0; time_index<K; time_index++)
            {
                double time = time_index*dt;
                ref_betaAqm.push_back( model->computeBetaQm( ref_mu , time ).template get<1>() );
            }
            auto ref_betaMqm = model->computeBetaQm( ref_mu , tf ).template get<0>() ;

            int sampling_size = Sampling->size();

            BOOST_FOREACH( auto mu, *Sampling )
            {
                int size = mu.size();

                element_type u_crb; // expansion of reduced solution
                element_type u_crb_dual; // expansion of reduced solution ( dual )
#if !NDEBUG
                if( proc_number == Environment::worldComm().masterRank() )
                {
                    std::cout << "(" << curpar << "/" << Sampling->size() << ") mu = [ ";
                    for ( int i=0; i<size-1; i++ ) std::cout<< mu[i] <<" , ";
                    std::cout<< mu[size-1]<<" ]\n ";
                }
#endif
                curpar++;

                std::ostringstream mu_str;
                //if too many parameters, it will crash
                int sizemax=7;
                if( size < sizemax )
                    sizemax=size;
                for ( int i=0; i<sizemax-1; i++ ) mu_str << std::scientific << std::setprecision( 5 ) << mu[i] <<",";
                mu_str << std::scientific << std::setprecision( 5 ) << mu[size-1];

#if !NDEBUG
                LOG(INFO) << "mu=" << mu << "\n";
                mu.check();
#endif
                if( option(_name="crb.script-mode").template as<bool>() )
                {
                    unsigned long N = mu.size() + 5;
                    unsigned long P = 2;
                    std::vector<double> X( N );
                    std::vector<double> Y( P );
                    for(int i=0; i<mu.size(); i++)
                        X[i] = mu[i];

                    int N_dim = crb_dimension;
                    int N_dimMax = crb_dimension_max;
                    int Nwn;
                    if( N_dim > 0 )
                        Nwn = N_dim;
                    else
                        Nwn = N_dimMax;
                    X[N-5] = output_index;
                    X[N-4] = Nwn;
                    X[N-3] = crb_online_tolerance;
                    X[N-2] = crb_error_type;
                    //X[N-1] = option(_name="crb.compute-variance").template as<int>();
                    X[N-1] = 0;
                    bool compute_variance = crb_compute_variance;
                    if ( compute_variance )
                        X[N-1] = 1;

                    this->run( X.data(), X.size(), Y.data(), Y.size() );
                    //std::cout << "output = " << Y[0] << std::endl;

                    std::ofstream res(option(_name="result-file").template as<std::string>() );
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

                                model->computeAffineDecomposition();
                                auto u_fem =  model->solveFemUsingAffineDecompositionFixedPoint( mu );
                                std::ostringstream u_fem_str;
                                u_fem_str << "u_fem(" << mu_str.str() << ")";
                                u_fem.setName( u_fem_str.str()  );

                                LOG(INFO) << "compute output\n";
                                if( export_solution )
                                    e->add( u_fem.name(), u_fem );
                                //e->step(0)->add( u_fem.name(), u_fem );
                                //model->solve( mu );
                                std::vector<double> o = boost::assign::list_of( model->output( output_index,mu , u_fem, true) )( ti.elapsed() );
                                if(proc_number == Environment::worldComm().masterRank() ) std::cout << "output=" << o[0] << "\n";
                                printEntry( ostr, mu, o );

                                std::ofstream res(option(_name="result-file").template as<std::string>() );
                                res << "output="<< o[0] << "\n";

                            }
                            break;

                        case  CRBModelMode::CRB:
                            {
                                LOG(INFO) << "CRB mode\n";
                                if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
                                    std::cout << "CRB mode -- "<<curpar<<"/"<<sampling_size<<std::endl;


                                boost::mpi::timer ti;

                                ti.restart();
                                LOG(INFO) << "solve crb\n";
                                //google::FlushLogFiles(google::GLOG_INFO);

                                //dimension of the RB (not necessarily the max)
                                int N =  option(_name="crb.dimension").template as<int>();

                                bool print_rb_matrix = option(_name="crb.print-rb-matrix").template as<bool>();
                                double online_tol = option(_name="crb.online-tolerance").template as<double>();
                                auto o = crb->run( mu, online_tol , N, print_rb_matrix);
                                double time_crb = ti.elapsed();

                                auto WN = crb->wn();
                                auto WNdu = crb->wndu();
                                //auto u_crb = crb->expansion( mu , N );
                                auto solutions=o.template get<2>();
                                auto uN = solutions.template get<0>();
                                auto uNdu = solutions.template get<1>();

                                int size = uN.size();

                                //if( model->isSteady()) // Re-use uN given by lb in crb->run
                                u_crb = crb->expansion( uN[size-1] , N , WN ); // Re-use uN given by lb in crb->run
                                if( solve_dual_problem )
                                    u_crb_dual = crb->expansion( uNdu[0] , N , WNdu );
                                //else
                                //    u_crb = crb->expansion( mu , N , WN );

                                std::ostringstream u_crb_str;
                                u_crb_str << "u_crb(" << mu_str.str() << ")";
                                u_crb.setName( u_crb_str.str()  );
                                LOG(INFO) << "export u_crb \n";
                                if( export_solution )
                                {
                                    e->add( u_crb.name(), u_crb );
                                }

                                double relative_error = -1;
                                double relative_estimated_error = -1;
                                auto matrix_info = o.template get<3>();
                                double condition_number = matrix_info.template get<0>();
                                double l2_error = -1;
                                double h1_error = -1;
                                double l2_dual_error = -1;
                                double h1_dual_error = -1;
                                double time_fem = -1;

                                element_type u_fem ;
                                element_type u_dual_fem ;

                                auto all_upper_bounds = o.template get<6>();
                                double output_estimated_error = all_upper_bounds.template get<0>();
                                double solution_estimated_error = all_upper_bounds.template get<1>();
                                double dual_solution_estimated_error = all_upper_bounds.template get<2>();

                                double ocrb = o.template get<0>();
                                if( vary_mu_comp0 > -1 )
                                {
                                    double x = mu(vary_mu_comp0);
                                    double mu0 = mu(vary_mu_comp0);
                                    outputs_storage(curpar-1)=ocrb;
                                    mu0_storage(curpar-1)=mu0;
                                    estimated_error_outputs_storage(curpar-1)=output_estimated_error;
                                }

                                if ( compute_fem )
                                {
									bool use_newton = option(_name="crb.use-newton").template as<bool>();

                                    ti.restart();
                                    LOG(INFO) << "solve u_fem\n";

                                    //auto u_fem = model->solveRB( mu );
                                    //auto u_fem = model->solveFemUsingOfflineEim( mu );

                                    if( boost::is_same<  crbmodel_type , crbmodelbilinear_type >::value && ! use_newton )
                                        //use affine decomposition
                                        u_fem = model->solveFemUsingAffineDecompositionFixedPoint( mu );
                                    else
                                        u_fem = model->solve( mu );

                                    std::ostringstream u_fem_str;
                                    u_fem_str << "u_fem(" << mu_str.str() << ")";
                                    u_fem.setName( u_fem_str.str()  );

                                    if( export_solution )
                                    {
                                        LOG(INFO) << "export u_fem \n";
                                        e->add( u_fem.name(), u_fem );
                                    }
                                    std::vector<double> ofem = boost::assign::list_of( model->output( output_index,mu, u_fem ) )( ti.elapsed() );

                                    relative_error = std::abs( ofem[0]- ocrb) /ofem[0];
                                    relative_estimated_error = output_estimated_error / ofem[0];

                                    //compute || u_fem - u_crb||_L2
                                    LOG(INFO) << "compute error \n";
                                    auto u_error = model->functionSpace()->element();
                                    auto u_dual_error = model->functionSpace()->element();
                                    std::ostringstream u_error_str;
                                    u_error = (( u_fem - u_crb ).pow(2)).sqrt()  ;
                                    u_error_str << "u_error(" << mu_str.str() << ")";
                                    u_error.setName( u_error_str.str()  );
                                    if( export_solution )
                                        e->add( u_error.name(), u_error );
                                    LOG(INFO) << "L2(fem)=" << l2Norm( u_fem )    << "\n";
                                    LOG(INFO) << "H1(fem)=" << h1Norm( u_fem )    << "\n";
                                    l2_error = l2Norm( u_error )/l2Norm( u_fem );
                                    h1_error = h1Norm( u_error )/h1Norm( u_fem );

                                    output_fem = ofem[0];
                                    time_fem = ofem[1];

                                    if( boost::is_same<  crbmodel_type , crbmodelbilinear_type >::value && ! use_newton )
                                    {
                                        if( solve_dual_problem )
                                        {
                                            u_dual_fem =  model->solveFemDualUsingAffineDecompositionFixedPoint( mu );

                                            u_dual_error = model->functionSpace()->element();
                                            u_dual_error = (( u_dual_fem - u_crb_dual ).pow(2)).sqrt() ;
                                            l2_dual_error = l2Norm( u_dual_error )/l2Norm( u_dual_fem );
                                            h1_dual_error = h1Norm( u_dual_error )/h1Norm( u_dual_fem );
                                        }
                                    }

                                }//compute-fem-during-online

                                if ( crb->errorType()==2 )
                                {
                                    double ocrb = o.template get<0>();
                                    std::vector<double> v = boost::assign::list_of( output_fem )( time_fem )( ocrb )( relative_estimated_error )( time_crb )( relative_error )( condition_number )( l2_error )( h1_error );
                                    if( proc_number == Environment::worldComm().masterRank() )
                                    {
                                        std::cout << "output=" << ocrb << " with " << o.template get<1>() << " basis functions\n";
                                        printEntry( file_summary_of_simulations, mu, v );
                                        printEntry( ostr, mu, v );
                                        //file_summary_of_simulations.close();

                                        if ( option(_name="crb.compute-stat").template as<bool>() && compute_fem )
                                        {
                                            relative_error_vector[curpar-1] = relative_error;
                                            l2_error_vector[curpar-1] = l2_error;
                                            h1_error_vector[curpar-1] = h1_error;
                                            time_fem_vector[curpar-1] = time_fem;
                                            time_crb_vector[curpar-1] = time_crb;
                                        }

                                        std::ofstream res(option(_name="result-file").template as<std::string>() );
                                        res << "output="<< ocrb << "\n";

                                    }

                                }//end of crb->errorType==2
                                else
                                {
                                    //if( ! boost::is_same<  crbmodel_type , crbmodelbilinear_type >::value )
                                    //    throw std::logic_error( "ERROR TYPE must be 2 when using CRBTrilinear (no error estimation)" );
                                    double ocrb = o.template get<0>();
                                    std::vector<double> v = boost::assign::list_of( output_fem )( time_fem )( ocrb )( relative_estimated_error )( ti.elapsed() ) ( output_estimated_error )( condition_number )( l2_error )( h1_error ) ;
                                    if( proc_number == Environment::worldComm().masterRank() )
                                    {
                                        std::cout << "output=" << ocrb << " with " << o.template get<1>() << " basis functions  (error estimation on this output : " << output_estimated_error<<") \n";
                                        //std::ofstream file_summary_of_simulations( ( boost::format( "summary_of_simulations_%d" ) % o.template get<2>() ).str().c_str() ,std::ios::out | std::ios::app );
                                        printEntry( file_summary_of_simulations, mu, v );
                                        printEntry( ostr, mu, v );
                                        //file_summary_of_simulations.close();

                                        if ( option(_name="crb.compute-stat").template as<bool>() && compute_fem )
                                        {
                                            relative_error_vector[curpar-1] = relative_error;
                                            l2_error_vector[curpar-1] = l2_error;
                                            h1_error_vector[curpar-1] = h1_error;
                                            time_fem_vector[curpar-1] = time_fem;
                                            time_crb_vector[curpar-1] = time_crb;
                                            relative_estimated_error_vector[curpar-1] = relative_estimated_error;
                                        }
                                        std::ofstream res(option(_name="result-file").template as<std::string>() );
                                        res << "output="<< ocrb << "\n";
                                    }//end of proc==master
                                }//end of else (errorType==2)

                                if (option(_name="eim.cvg-study").template as<bool>() )
                                {
                                    bool check_name = false;
                                    std::string how_compute_unknown = option(_name=_o( this->about().appName(),"how-compute-unkown-for-eim" )).template as<std::string>();
                                    if( how_compute_unknown == "CRB-with-ad")
                                    {
                                        LOG( INFO ) << "convergence eim with CRB-with-ad ";
                                        check_name = true;
                                        this->studyEimConvergence( mu , u_crb , curpar );
                                    }
                                    if( how_compute_unknown == "FEM-with-ad")
                                    {
                                        LOG( INFO ) << "convergence eim with FEM-with-ad ";
                                        check_name = true;
                                        //fem computed via solveFemUsingAffineDecomposition use the affine decomposition
                                        auto fem_with_ad = model->solveFemUsingAffineDecompositionFixedPoint( mu );
                                        this->studyEimConvergence( mu , fem_with_ad , curpar );
                                    }
                                    if( how_compute_unknown == "FEM-without-ad")
                                    {
                                        LOG( INFO ) << "convergence eim with FEM-without-ad ";
                                        check_name = true;
                                        auto fem_without_ad = model->solve( mu );
                                        this->studyEimConvergence( mu , fem_without_ad , curpar );
                                    }
                                    if( ! check_name )
                                        throw std::logic_error( "OpusApp error with option how-compute-unknown-for-eim, please use CRB-with-ad, FEM-with-ad or FEM-without-ad" );
                                }

                                if (option(_name="crb.cvg-study").template as<bool>() && compute_fem )
                                {

                                    LOG(INFO) << "start convergence study...\n";
                                    std::map<int, boost::tuple<double,double,double,double,double,double,double> > conver;

                                    std::ofstream fileL2 ( "CrbConvergenceL2.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileH1 ( "CrbConvergenceH1.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileOutputError ( "CrbConvergenceOutputError.dat",std::ios::out | std::ios::app );
                                    std::ofstream fileOutputEstimatedError ( "CrbConvergenceOutputErrorEstimated.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileOutputErrorBoundEfficiency ( "CrbConvergenceOutputErrorBoundEfficiency.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileSolutionErrorBoundEfficiency ("CrbConvergencePrimalSolutionErrorBoundEfficiency.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileSolutionError ("CrbConvergencePrimalSolutionError.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileSolutionErrorEstimated ("CrbConvergencePrimalSolutionErrorEstimated.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileSolutionDualErrorBoundEfficiency ("CrbConvergenceDualSolutionErrorBoundEfficiency.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileSolutionDualError ("CrbConvergenceDualSolutionError.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileSolutionDualErrorEstimated ("CrbConvergenceDualSolutionErrorEstimated.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream filePrimalResidualNorm ("CrbConvergencePrimalResidualNorm.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileDualResidualNorm ( "CrbConvergenceDualResidualNorm.dat" ,std::ios::out | std::ios::app );


                                    int Nmax = option("crb.dimension-max").template as<int>();
                                    if( Environment::worldComm().isMasterRank() )
                                    {
                                            fileL2 << Nmax <<"\t";
                                            fileH1 << Nmax <<"\t";
                                            fileOutputError << Nmax <<"\t";
                                            fileOutputEstimatedError << Nmax << "\t";
                                            fileOutputErrorBoundEfficiency <<  Nmax << "\t";
                                            fileSolutionErrorBoundEfficiency << Nmax << "\t";
                                            fileSolutionError << Nmax << "\t";
                                            fileSolutionErrorEstimated << Nmax << "\t";
                                            fileSolutionDualErrorBoundEfficiency << Nmax << "\t" ;
                                            fileSolutionDualError << Nmax << "\t";
                                            fileSolutionDualErrorEstimated << Nmax << "\t";
                                            filePrimalResidualNorm << Nmax << "\t";
                                            fileDualResidualNorm << Nmax << "\t";

                                    }
                                    std::string str = "\t";
                                    for( int N = 1; N <= Nmax ; N++ )
                                    {
                                        auto o= crb->run( mu,  online_tol , N, print_rb_matrix);
                                        auto ocrb = o.template get<0>();
                                        auto solutions=o.template get<2>();
                                        auto u_crb = solutions.template get<0>();
                                        auto u_crb_du = solutions.template get<1>();
                                        int size = u_crb.size();
                                        auto uN = crb->expansion( u_crb[size-1], N, WN );

                                        element_type uNdu;

                                        auto u_error = u_fem - uN;
                                        auto u_dual_error = model->functionSpace()->element();
                                        if( solve_dual_problem )
                                        {
                                            uNdu = crb->expansion( u_crb_du[0], N, WNdu );
                                            u_dual_error = u_dual_fem - uNdu;
                                        }

                                        auto all_upper_bounds = o.template get<6>();
                                        output_estimated_error = all_upper_bounds.template get<0>();
                                        solution_estimated_error = all_upper_bounds.template get<1>();
                                        dual_solution_estimated_error = all_upper_bounds.template get<2>();

                                        //auto o = crb->run( mu,  option(_name="crb.online-tolerance").template as<double>() , N);
                                        double rel_err = std::abs( output_fem-ocrb ) /output_fem;

                                        double output_relative_estimated_error = output_estimated_error / output_fem;

                                        double primal_residual_norm = o.template get<4>();
                                        double dual_residual_norm = o.template get<5>();

                                        double solution_error=0;
                                        double dual_solution_error=0;
                                        double square_solution_error=0;
                                        double square_dual_solution_error=0;
                                        double ref_primal=0;
                                        double ref_dual=0;
                                        if( model->isSteady() )
                                        {
                                            //let ufem-ucrb = e
                                            //||| e |||_mu = sqrt( a( e , e ; mu ) ) = solution_error
                                            for(int q=0; q<model->Qa();q++)
                                            {
                                                for(int m=0; m<model->mMaxA(q); m++)
                                                {
                                                    solution_error +=  ref_betaAqm[0][q][m]*model->Aqm(q,m,u_error,u_error) ;
                                                    ref_primal +=  ref_betaAqm[0][q][m]*model->Aqm(q,m,u_fem,u_fem);
                                                }
                                            }

                                            if( solve_dual_problem )
                                            {
                                                for(int q=0; q<model->Qa();q++)
                                                {
                                                    for(int m=0; m<model->mMaxA(q); m++)
                                                    {
                                                        dual_solution_error += ref_betaAqm[0][q][m]*model->Aqm(q,m,u_dual_error,u_dual_error);
                                                        ref_dual += ref_betaAqm[0][q][m]*model->Aqm(q,m,u_dual_fem,u_dual_fem);
                                                    }
                                                }
                                                square_dual_solution_error = dual_solution_error;
                                                dual_solution_error = math::sqrt( dual_solution_error );
                                                ref_dual = math::sqrt( ref_dual );
                                            }
                                            square_solution_error = solution_error;
                                            solution_error = math::sqrt( solution_error );
                                            ref_primal = math::sqrt( ref_primal );
                                            //dual_solution_error = math::sqrt( model->scalarProduct( u_dual_error, u_dual_error ) );

                                        }
                                        else
                                        {
                                            double dt = model->timeStep();
                                            double ti = model->timeInitial();
                                            double tf = model->timeFinal();
                                            int K = ( tf - ti )/dt;

                                            for(int q=0; q<model->Qm();q++)
                                            {
                                                for(int m=0; m<model->mMaxM(q); m++)
                                                {
                                                    solution_error +=  ref_betaMqm[q][m]*model->Mqm(q,m,u_error,u_error);
                                                    ref_primal +=  ref_betaMqm[q][m]*model->Mqm(q,m,u_fem,u_fem);
                                                }
                                            }
                                            for(int time_index=0; time_index<K; time_index++)
                                            {
                                                double t=time_index*dt;
                                                for(int q=0; q<model->Qa();q++)
                                                {
                                                    for(int m=0; m<model->mMaxA(q); m++)
                                                    {
                                                        solution_error +=  ref_betaAqm[time_index][q][m]*model->Aqm(q,m,u_error,u_error) * dt;
                                                        ref_primal +=  ref_betaAqm[time_index][q][m]*model->Aqm(q,m,u_fem,u_fem) * dt;
                                                    }
                                                }
                                            }
                                            square_solution_error = solution_error;
                                            solution_error = math::sqrt( solution_error );
                                            ref_primal = math::sqrt( ref_primal );

                                            if( solve_dual_problem )
                                            {
                                                ti = model->timeFinal()+dt;
                                                tf = model->timeInitial()+dt;
                                                dt -= dt;

                                                for(int q=0; q<model->Qm();q++)
                                                {
                                                    for(int m=0; m<model->mMaxM(q); m++)
                                                    {
                                                        dual_solution_error +=  ref_betaMqm[q][m]*model->Mqm(q,m,u_dual_error,u_dual_error);
                                                        ref_dual +=  ref_betaMqm[q][m]*model->Mqm(q,m,u_dual_fem,u_dual_fem);
                                                    }
                                                }

                                                for(int time_index=0; time_index<K; time_index++)
                                                {
                                                    double t=time_index*dt;
                                                    for(int q=0; q<model->Qa();q++)
                                                    {
                                                        for(int m=0; m<model->mMaxA(q); m++)
                                                        {
                                                            dual_solution_error +=  ref_betaAqm[time_index][q][m]*model->Aqm(q,m,u_dual_error,u_dual_error) * dt;
                                                            ref_dual +=  ref_betaAqm[time_index][q][m]*model->Aqm(q,m,u_dual_fem,u_dual_fem) * dt;
                                                        }
                                                    }
                                                }
                                                square_dual_solution_error = dual_solution_error;
                                                dual_solution_error = math::sqrt( dual_solution_error );
                                                ref_dual = math::sqrt( ref_dual );

                                            }//if solve-dual

                                        }//transient case

                                        double l2_error = l2Norm( u_error )/l2Norm( u_fem );
                                        double h1_error = h1Norm( u_error )/h1Norm( u_fem );
                                        auto matrix_info = o.template get<3>();
                                        double condition_number = matrix_info.template get<0>();
                                        double output_error_bound_efficiency = output_relative_estimated_error / rel_err;

                                        double relative_primal_solution_error = solution_error / ref_primal ;
                                        double relative_primal_solution_estimated_error = solution_estimated_error / ref_primal;
                                        double relative_primal_solution_error_bound_efficiency = relative_primal_solution_estimated_error / relative_primal_solution_error;

                                        if( relative_primal_solution_error_bound_efficiency < 1 )
                                        {
                                            vector_sampling_for_primal_efficiency_under_1[N-1]->push_back( mu , 1);
                                        }

                                        double relative_dual_solution_error = 1;
                                        double relative_dual_solution_estimated_error = 1;
                                        double relative_dual_solution_error_bound_efficiency = 1;
                                        if( solve_dual_problem )
                                        {
                                            relative_dual_solution_error = dual_solution_error / ref_dual ;
                                            relative_dual_solution_estimated_error = dual_solution_estimated_error / ref_dual;
                                            relative_dual_solution_error_bound_efficiency = relative_dual_solution_estimated_error / relative_dual_solution_error;

                                            if( relative_dual_solution_error_bound_efficiency < 1 )
                                            {
                                                vector_sampling_for_dual_efficiency_under_1[N-1]->push_back( mu , 0);
                                            }

                                        }
                                        conver[N]=boost::make_tuple( rel_err, l2_error, h1_error , relative_estimated_error, condition_number , output_error_bound_efficiency , relative_primal_solution_error_bound_efficiency );

                                        //LOG(INFO) << "N=" << N << " " << rel_err << " " << l2_error << " " << h1_error << " " <<condition_number<<"\n";
                                        if ( proc_number == Environment::worldComm().masterRank() )
                                        {
                                            std::cout << "N=" << N << "Output =  "<< output_fem <<" OutputError = "<<rel_err <<" OutputErrorEstimated = "<<relative_estimated_error
                                                      <<"  L2Error = "<< l2_error << "  H1Error = " << h1_error <<std::endl;

                                            if( N == Nmax )
                                                str="\n";
                                            fileL2 << l2_error <<str;
                                            fileH1 << h1_error <<str;
                                            fileOutputError << rel_err <<str;
                                            fileOutputEstimatedError << output_relative_estimated_error << str;
                                            fileOutputErrorBoundEfficiency <<  output_error_bound_efficiency << str;
                                            fileSolutionErrorBoundEfficiency << relative_primal_solution_error_bound_efficiency << str;
                                            fileSolutionError << relative_primal_solution_error << str;
                                            fileSolutionErrorEstimated << relative_primal_solution_estimated_error << str;
                                            fileSolutionDualErrorBoundEfficiency << relative_dual_solution_error_bound_efficiency << str ;
                                            fileSolutionDualError << relative_dual_solution_error << str;
                                            fileSolutionDualErrorEstimated <<  relative_dual_solution_estimated_error << str;
                                            filePrimalResidualNorm << primal_residual_norm << str;
                                            fileDualResidualNorm <<  dual_residual_norm << str;
                                        }
                                        if( option(_name="crb.compute-matrix-information").template as<bool>() )
                                        {
                                            auto matrix_info = o.template get<3>();// conditioning of primal reduced matrix + determinant
                                            double conditioning = matrix_info.template get<0>();
                                            double determinant = matrix_info.template get<1>();
                                            LOG( INFO ) << " primal reduced matrix information ";
                                            LOG( INFO ) << std::setprecision(15)<<"mu : \n"<<mu;
                                            LOG( INFO ) << std::setprecision(15)<<"conditioning : "<<conditioning;
                                            LOG( INFO ) << std::setprecision(15)<<"determinant : "<<determinant;
                                        }
                                        if( relative_primal_solution_error_bound_efficiency < 1 )
                                        {
                                            LOG( INFO ) << "N : "<<N;
                                            LOG( INFO ) << std::setprecision(15)<<"efficiency of error estimation on primal solution is "<<relative_primal_solution_error_bound_efficiency<<" ( should be >= 1 )";
                                            LOG( INFO ) << std::setprecision(15)<<"mu : \n"<<mu;
                                            LOG( INFO ) << std::setprecision(15)<<"relative_primal_solution_estimated_error : "<<relative_primal_solution_estimated_error;
                                            LOG( INFO ) << std::setprecision(15)<<"relative_primal_solution_error : "<<relative_primal_solution_error;
                                            LOG( INFO ) << std::setprecision(15)<<"primal_solution_estimated_error : "<<solution_estimated_error;
                                            LOG( INFO ) << std::setprecision(15)<<"primal_solution_error : "<<solution_error;
                                            LOG( INFO ) << std::setprecision(15)<<"square error : "<<square_solution_error;
                                            //LOG( INFO ) << std::setprecision(15)<<"u_crb : \n"<<u_crb[size-1];
                                            LOG( INFO ) << std::setprecision(15)<<"primal solution norme  : "<<uN.l2Norm();
                                        }
                                        if( relative_dual_solution_error_bound_efficiency < 1 )
                                        {
                                            LOG( INFO ) <<std::setprecision(15)<< "efficiency of error estimation on dual solution is "<<relative_dual_solution_error_bound_efficiency<<" ( should be >= 1 )";
                                            LOG( INFO ) <<std::setprecision(15)<< "mu : \n"<<mu;
                                            LOG( INFO ) <<std::setprecision(15)<<"relative_dual_solution_estimated_error : "<<relative_dual_solution_estimated_error;
                                            LOG( INFO ) <<std::setprecision(15)<<"relative_dual_solution_error : "<<relative_dual_solution_error;
                                            LOG( INFO ) <<std::setprecision(15)<<"dual_solution_estimated_error : "<<dual_solution_estimated_error;
                                            LOG( INFO ) <<std::setprecision(15)<<"dual_solution_error : "<<dual_solution_error;
                                            LOG( INFO ) << std::setprecision(15)<<"square error : "<<square_dual_solution_error;
                                            //LOG( INFO ) << std::setprecision(15)<<"u_crb_du : \n"<<u_crb_du[0];
                                            LOG( INFO ) << std::setprecision(15)<<"dual solution norme  : "<<uNdu.l2Norm();
                                        }
                                        if( output_error_bound_efficiency < 1 )
                                        {
                                            LOG( INFO ) <<std::setprecision(15)<<"efficiency of error estimation on output is "<<output_error_bound_efficiency<<" ( should be >= 1 )";
                                            LOG( INFO ) <<std::setprecision(15)<< "mu : \n"<<mu;
                                            LOG( INFO ) <<std::setprecision(15)<< "output_relative_estimated_error : "<<output_relative_estimated_error;
                                            LOG( INFO ) <<std::setprecision(15)<< "output_relative_error : "<<rel_err;
                                            LOG( INFO ) <<std::setprecision(15)<< "output_estimated_error : "<<output_estimated_error;
                                            LOG( INFO ) <<std::setprecision(15)<< "output_error : "<< std::abs( output_fem-ocrb );
                                        }

                                        if ( option(_name="crb.check.residual").template as<bool>()  && solve_dual_problem  )
                                        {
                                            std::vector < std::vector<double> > primal_residual_coefficients = all_upper_bounds.template get<3>();
                                            std::vector < std::vector<double> > dual_residual_coefficients = all_upper_bounds.template get<4>();
                                            if( model->isSteady() )
                                            {
                                                crb->checkResidual( mu , primal_residual_coefficients, dual_residual_coefficients, uN, uNdu );
                                            }
                                            else
                                            {
                                                std::vector< element_type > uNelement;
                                                std::vector< element_type > uNelement_old;
                                                std::vector< element_type > uNelement_du;
                                                std::vector< element_type > uNelement_du_old;

                                                auto u_crb_old = solutions.template get<2>();
                                                auto u_crb_du_old = solutions.template get<3>();

                                                //size is the number of time step
                                                for(int t=0; t<size; t++)
                                                {
                                                    uNelement.push_back( crb->expansion( u_crb[t], N, WN ) );
                                                    uNelement_old.push_back( crb->expansion( u_crb_old[t], N, WN ) );
                                                    uNelement_du.push_back( crb->expansion( u_crb_du[t], N, WNdu ) );
                                                    uNelement_du_old.push_back( crb->expansion( u_crb_du_old[t], N, WNdu ) );
                                                }//loop over time step

                                                crb->compareResidualsForTransientProblems(N, mu ,
                                                                                          uNelement, uNelement_old, uNelement_du, uNelement_du_old,
                                                                                          primal_residual_coefficients, dual_residual_coefficients );
                                            }//transient case
                                        }//check residuals computations

                                    }//loop over basis functions ( N )
#if 0
                                        LOG(INFO) << "save in logfile\n";
                                        std::string mu_str;
                                        for ( int i=0; i<mu.size(); i++ )
                                            mu_str= mu_str + ( boost::format( "_%1%" ) %mu[i] ).str() ;
                                        std::string file_name = "convergence"+mu_str+".dat";
                                        std::ofstream conv( file_name );
                                        BOOST_FOREACH( auto en, conver )
                                            conv << en.first << "\t" << en.second.get<0>()  << "\t" << en.second.get<1>() << "\t" << en.second.get<2>() <<
                                            "\t"<< en.second.get<3>() << "\t"<< en.second.get<4>()<< "\t" <<en.second.get<5>()<< "\n";
#endif
                                }//end of cvg-study
                            }//case CRB
                            break;
                        case  CRBModelMode::CRB_ONLINE:
                            {
                                std::cout << "CRB Online mode\n";
                                boost::mpi::timer ti;
                                ti.restart();
                                auto o = crb->run( mu,  option(_name="crb.online-tolerance").template as<double>() );

                                if ( crb->errorType()==2 )
                                    {
                                        std::vector<double> v = boost::assign::list_of( o.template get<0>() )( ti.elapsed() );
                                        std::cout << "output=" << o.template get<0>() << " with " << o.template get<1>() << " basis functions\n";
                                        printEntry( ostr, mu, v );
                                    }

                                else
                                    {
                                        auto all_upper_bounds = o.template get<6>();
                                        double output_estimated_error = all_upper_bounds.template get<0>();
                                        double relative_estimated_error = output_estimated_error / output_fem;
                                        std::vector<double> v = boost::assign::list_of( o.template get<0>() )( output_estimated_error )( ti.elapsed() );
                                        std::cout << "output=" << o.template get<0>() << " with " << o.template get<1>() <<
                                            " basis functions  (relative error estimation on this output : " << relative_estimated_error<<") \n";
                                        printEntry( ostr, mu, v );
                                    }

                            }
                            break;

                        case  CRBModelMode::SCM:
                            {

                                if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
                                    std::cout << "SCM mode\n";
                                int kmax = crb->scm()->KMax();
                                auto o = crb->scm()->run( mu, kmax );
                                printEntry( file_summary_of_simulations, mu, o );
                                printEntry( ostr, mu, o );

                                scm_relative_error[curpar-1] = o[6];

                                if (option(_name="crb.scm.cvg-study").template as<bool>()  )
                                {
                                    LOG(INFO) << "start scm convergence study...\n";
                                    std::map<int, boost::tuple<double> > conver;
                                    for( int N = 1; N <= kmax; N++ )
                                    {
                                        auto o = crb->scm()->run( mu, N);
                                        double relative_error = o[6];
                                        conver[N]=boost::make_tuple( relative_error );
                                        if ( proc_number == Environment::worldComm().masterRank() )
                                            std::cout << "N=" << N << " " << relative_error <<std::endl;
                                        M_mapConvSCM["RelativeError"][N-1](curpar - 1) = relative_error;
                                    }
                                    if( proc_number == Environment::worldComm().masterRank() )
                                    {
                                        LOG(INFO) << "save in logfile\n";
                                        std::string mu_str;
                                        for ( int i=0; i<mu.size(); i++ )
                                            mu_str= mu_str + ( boost::format( "_%1%" ) %mu[i] ).str() ;
                                        std::string file_name = "convergence-scm-"+mu_str+".dat";
                                        std::ofstream conv( file_name );
                                        BOOST_FOREACH( auto en, conver )
                                            conv << en.first << "\t" << en.second.get<0>()  ;
                                    }
                                }//end of cvg-study

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

                    LOG( INFO ) << "------------------------------------------------------------";
                }
            }

            model->generateGeoFileForOutputPlot( outputs_storage , mu0_storage, estimated_error_outputs_storage );

            //model->computationalTimeEimStatistics();
            if( export_solution )
                e->save();

            if( proc_number == Environment::worldComm().masterRank() ) std::cout << ostr.str() << "\n";

            if (option(_name="eim.cvg-study").template as<bool>() && M_mode==CRBModelMode::CRB)
                this->doTheEimConvergenceStat( Sampling->size() );

            if (option(_name="crb.cvg-study").template as<bool>() && compute_fem && M_mode==CRBModelMode::CRB )
                this->doTheCrbConvergenceStat( Sampling->size() );

            if (option(_name="crb.scm.cvg-study").template as<bool>() && M_mode==CRBModelMode::SCM )
                this->doTheScmConvergenceStat( Sampling->size() );

            if ( compute_stat && compute_fem && M_mode==CRBModelMode::CRB )
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
            }//end of compute-stat CRB
            if( M_mode==CRBModelMode::SCM )
            {
                LOG( INFO ) << "compute statistics \n";
                Eigen::MatrixXf::Index index_max_error;
                Eigen::MatrixXf::Index index_min_error;
                double max_error = scm_relative_error.maxCoeff(&index_max_error);
                double min_error = scm_relative_error.minCoeff(&index_min_error);
                double mean_error = scm_relative_error.mean();
                if( proc_number == Environment::worldComm().masterRank() )
                {
                    file_summary_of_simulations <<"\n\nStatistics\n";
                    file_summary_of_simulations <<"max of relative error : "<<max_error<<" at the "<<index_max_error+1<<"^th simulation\n";
                    file_summary_of_simulations <<"min of relative error : "<<min_error<<" at the "<<index_min_error+1<<"^th simulation\n";
                    file_summary_of_simulations <<"mean of relative error : "<<mean_error<<"\n\n";
                }

            }//end of compute-stat SCM


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


    /*
     * read the file generated by PsLogger
     * and generate a file to plot memory evolution in function of number of elements
     * ( during offline step )
     */
    void generateMemoryEvolution(std::string pslogname)
    {
        std::ifstream file_read (pslogname);
        std::vector< double > memory;
        std::vector< double > percent;

        // number of elements in the RB
        int N = option("crb.dimension-max").template as<int>();
        //int N = crb->dimension();
        memory.resize( N );
        percent.resize( N );
        if( file_read )
        {
            double mem;
            double per;
            std::string str;
            //first line
            for(int i=0; i<5; i++)
                file_read >> str;
            //second line
            for(int i=0; i<3; i++)
                file_read >> str;
            //for each elements in the RB
            for(int i=0; i<N; i++)
            {
                file_read >> str ;
                file_read >> mem;
                file_read >> per;
                file_read >> str;
                memory[i] = mem;
                percent[i] = per;
            }
            file_read.close();
        }//if file
        else
        {
            throw std::logic_error( "[OpusApp::generateMemoryEvaluation] ERROR loading the file generated by ps-log (does not exist)" );
        }//no file

        int proc_number = Environment::worldComm().globalRank();
        int global_size = Environment::worldComm().globalSize();
        std::string file_name = ( boost::format("MemoryEvolution-%1%_%2%") %global_size %proc_number ).str();
        std::ofstream file_write;
        file_write.open( file_name,std::ios::out );
        file_write << "N \t Memory \t Percent \n";
        for(int i=0; i<N; i++)
            file_write << i << "\t"<< memory[i]<< "\t" << percent[i] << "\n";
        file_write.close();
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



    void initializeConvergenceScmMap( int sampling_size )
    {
        auto N = crb->scm()->KMax();
        M_mapConvSCM["RelativeError"] = std::vector<vectorN_type>(N);

        for(int j=0; j<N; j++)
        {
            M_mapConvSCM["RelativeError"][j].resize(sampling_size);
        }
    }

    void studyEimConvergence( typename ModelType::parameter_type const& mu , element_type & model_solution , int mu_number)
    {
        auto eim_sc_vector = model->scalarContinuousEim();
        auto eim_sd_vector = model->scalarDiscontinuousEim();

        for(int i=0; i<eim_sc_vector.size(); i++)
        {
            auto eim = eim_sc_vector[i];
            std::vector< std::string > all_file_name;
            all_file_name.push_back( eim->name()+"-EimConvergenceL2.dat");
            all_file_name.push_back( eim->name()+"-EimConvergenceL2estimated.dat");
            all_file_name.push_back( eim->name()+"-EimConvergenceL2ratio.dat");
            all_file_name.push_back( eim->name()+"-EimConvergenceLINF.dat");
            all_file_name.push_back( eim->name()+"-EimConvergenceLINFestimated.dat");
            all_file_name.push_back( eim->name()+"-EimConvergenceLINFratio.dat");
            eim->studyConvergence( mu , model_solution , all_file_name );
        }

        for(int i=0; i<eim_sd_vector.size(); i++)
        {
            std::vector< std::string > all_file_name;
            auto eim = eim_sd_vector[i];
            all_file_name.push_back( eim->name()+"-EimConvergenceL2.dat");
            all_file_name.push_back( eim->name()+"-EimConvergenceL2estimated.dat");
            all_file_name.push_back( eim->name()+"-EimConvergenceL2ratio.dat");
            all_file_name.push_back( eim->name()+"-EimConvergenceLINF.dat");
            all_file_name.push_back( eim->name()+"-EimConvergenceLINFestimated.dat");
            all_file_name.push_back( eim->name()+"-EimConvergenceLINFratio.dat");
            eim->studyConvergence( mu , model_solution , all_file_name );
        }

        //std::vector<double> l2error;
        //for(int j=0; j<l2error.size(); j++)
        //    M_mapConvEIM[eim->name()][j][mu_number-1] = l2error[j];

    }

    void doTheEimConvergenceStat( int sampling_size )
    {
        auto eim_sc_vector = model->scalarContinuousEim();
        auto eim_sd_vector = model->scalarDiscontinuousEim();

        for(int i=0; i<eim_sc_vector.size(); i++)
        {
            auto eim = eim_sc_vector[i];

            //load information contained in files
            std::vector< vectorN_type > L2, L2estimated, L2ratio, LINF, LINFestimated, LINFratio;
            std::string filename = eim->name()+"-EimConvergenceL2.dat";
            model->readConvergenceDataFromFile( L2, filename );
            filename = eim->name()+"-EimConvergenceL2estimated.dat";
            model->readConvergenceDataFromFile( L2estimated, filename );
            filename = eim->name()+"-EimConvergenceL2ratio.dat";
            model->readConvergenceDataFromFile( L2ratio, filename );
            filename = eim->name()+"-EimConvergenceLINF.dat";
            model->readConvergenceDataFromFile( LINF, filename );
            filename = eim->name()+"-EimConvergenceLINFestimated.dat";
            model->readConvergenceDataFromFile( LINFestimated, filename );
            filename = eim->name()+"-EimConvergenceLINFratio.dat";
            model->readConvergenceDataFromFile( LINFratio, filename );

            //write files containing statistics
            filename = "cvg-eim-"+eim->name()+"-L2.dat";
            model->writeConvergenceStatistics( L2 , filename);
            filename = "cvg-eim-"+eim->name()+"-L2estimated.dat";
            model->writeConvergenceStatistics( L2estimated , filename);
            filename = "cvg-eim-"+eim->name()+"-L2ratio.dat";
            model->writeConvergenceStatistics( L2ratio , filename);
            filename = "cvg-eim-"+eim->name()+"-LINF.dat";
            model->writeConvergenceStatistics( LINF , filename);
            filename = "cvg-eim-"+eim->name()+"-LINFestimated.dat";
            model->writeConvergenceStatistics( LINFestimated , filename);
            filename = "cvg-eim-"+eim->name()+"-LINFratio.dat";
            model->writeConvergenceStatistics( LINFratio , filename);
        }

        for(int i=0; i<eim_sd_vector.size(); i++)
        {
            auto eim = eim_sd_vector[i];

            //load information contained in files
            std::vector< vectorN_type > L2, L2estimated, L2ratio, LINF, LINFestimated, LINFratio;
            std::string filename = eim->name()+"-EimConvergenceL2.dat";
            model->readConvergenceDataFromFile( L2, filename );
            filename = eim->name()+"-EimConvergenceL2estimated.dat";
            model->readConvergenceDataFromFile( L2estimated, filename );
            filename = eim->name()+"-EimConvergenceL2ratio.dat";
            model->readConvergenceDataFromFile( L2ratio, filename );
            filename = eim->name()+"-EimConvergenceLINF.dat";
            model->readConvergenceDataFromFile( LINF, filename );
            filename = eim->name()+"-EimConvergenceLINFestimated.dat";
            model->readConvergenceDataFromFile( LINFestimated, filename );
            filename = eim->name()+"-EimConvergenceLINFratio.dat";
            model->readConvergenceDataFromFile( LINFratio, filename );

            //write files containing statistics
            filename = "cvg-eim-"+eim->name()+"-L2.dat";
            model->writeConvergenceStatistics( L2 , filename);
            filename = "cvg-eim-"+eim->name()+"-L2estimated.dat";
            model->writeConvergenceStatistics( L2estimated , filename);
            filename = "cvg-eim-"+eim->name()+"-L2ratio.dat";
            model->writeConvergenceStatistics( L2ratio , filename);
            filename = "cvg-eim-"+eim->name()+"-LINF.dat";
            model->writeConvergenceStatistics( LINF , filename);
            filename = "cvg-eim-"+eim->name()+"-LINFestimated.dat";
            model->writeConvergenceStatistics( LINFestimated , filename);
            filename = "cvg-eim-"+eim->name()+"-LINFratio.dat";
            model->writeConvergenceStatistics( LINFratio , filename);

        }
    }

    void doTheCrbConvergenceStat( int sampling_size )
    {
        int N = option("crb.dimension-max").template as<int>();

        std::vector< vectorN_type > L2, H1, OutputError, OutputErrorEstimated, OutputErrorBoundEfficiency;
        std::vector< vectorN_type > PrimalSolutionError, PrimalSolutionErrorEstimated, PrimalSolutionErrorBoundEfficiency;
        std::vector< vectorN_type > DualSolutionError, DualSolutionErrorEstimated, DualSolutionErrorBoundEfficiency;

        //load information contained in files
        std::string filename = "CrbConvergenceL2.dat";
        model->readConvergenceDataFromFile( L2, filename );
        filename = "CrbConvergenceH1.dat";
        model->readConvergenceDataFromFile( H1, filename );
        filename = "CrbConvergenceOutputError.dat";
        model->readConvergenceDataFromFile( OutputError, filename );
        filename = "CrbConvergenceOutputErrorEstimated.dat";
        model->readConvergenceDataFromFile( OutputErrorEstimated, filename );
        filename = "CrbConvergenceOutputErrorBoundEfficiency.dat";
        model->readConvergenceDataFromFile( OutputErrorBoundEfficiency, filename );
        filename = "CrbConvergencePrimalSolutionError.dat";
        model->readConvergenceDataFromFile( PrimalSolutionError , filename );
        filename = "CrbConvergencePrimalSolutionErrorEstimated.dat";
        model->readConvergenceDataFromFile( PrimalSolutionErrorEstimated , filename );
        filename = "CrbConvergencePrimalSolutionErrorBoundEfficiency.dat";
        model->readConvergenceDataFromFile( PrimalSolutionErrorBoundEfficiency , filename );
        filename = "CrbConvergenceDualSolutionError.dat";
        model->readConvergenceDataFromFile( DualSolutionError , filename );
        filename = "CrbConvergenceDualSolutionErrorEstimated.dat";
        model->readConvergenceDataFromFile( DualSolutionErrorEstimated , filename );
        filename = "CrbConvergenceDualSolutionErrorBoundEfficiency.dat";
        model->readConvergenceDataFromFile( DualSolutionErrorBoundEfficiency , filename );

        //write files containing statistics
        filename = "cvg-crb-L2.dat";
        model->writeConvergenceStatistics( L2 , filename);
        filename = "cvg-crb-H1.dat";
        model->writeConvergenceStatistics( H1 , filename);
        filename = "cvg-crb-OutputError.dat";
        model->writeConvergenceStatistics( OutputError , filename);
        filename = "cvg-crb-OutputErrorEstimated.dat";
        model->writeConvergenceStatistics( OutputErrorEstimated , filename);
        filename = "cvg-crb-OutputErrorBoundEfficiency.dat";
        model->writeConvergenceStatistics( OutputErrorBoundEfficiency , filename);
        filename = "cvg-crb-PrimalSolutionError.dat";
        model->writeConvergenceStatistics( PrimalSolutionError , filename);
        filename = "cvg-crb-PrimalSolutionErrorEstimated.dat";
        model->writeConvergenceStatistics( PrimalSolutionErrorEstimated , filename);
        filename = "cvg-crb-PrimalSolutionErrorBoundEfficiency.dat";
        model->writeConvergenceStatistics( PrimalSolutionErrorBoundEfficiency , filename);
        filename = "cvg-crb-DualSolutionError.dat";
        model->writeConvergenceStatistics( DualSolutionError , filename);
        filename = "cvg-crb-DualSolutionErrorEstimated.dat";
        model->writeConvergenceStatistics( DualSolutionErrorEstimated , filename);
        filename = "cvg-crb-DualSolutionErrorBoundEfficiency.dat";
        model->writeConvergenceStatistics( DualSolutionErrorBoundEfficiency , filename);

        //now we have error and estimated error statistics for all parameters (min max ect ... )
        //but it is interesting to have also the efficiency associated to parameter that generated the max/min
        //error on the output or on primal/dual solution
        //so that what we do here
        filename = "cvg-crb-OutputErrorBoundEfficiency2.dat";
        model->writeVectorsExtremumsRatio( OutputError, OutputErrorEstimated, filename );
        filename = "cvg-crb-PrimalSolutionErrorBoundEfficiency2.dat";
        model->writeVectorsExtremumsRatio( PrimalSolutionError, PrimalSolutionErrorEstimated, filename );
        filename = "cvg-crb-DualSolutionErrorBoundEfficiency2.dat";
        model->writeVectorsExtremumsRatio( DualSolutionError, DualSolutionErrorEstimated, filename );

        int Nmax = vector_sampling_for_primal_efficiency_under_1.size();
        for(int N=0; N<Nmax; N++)
        {
            std::string file_name_primal = (boost::format("Sampling-Primal-Problem-Bad-Efficiency-N=%1%") %(N+1) ).str();
            if( vector_sampling_for_primal_efficiency_under_1[N]->size() > 1 )
                vector_sampling_for_primal_efficiency_under_1[N]->writeOnFile(file_name_primal);
        }
        Nmax = vector_sampling_for_dual_efficiency_under_1.size();
        for(int N=0; N<Nmax; N++)
        {
            std::string file_name_dual = (boost::format("Sampling-Dual-Problem-Bad-Efficiency-N=%1%") %(N+1) ).str();
            if( vector_sampling_for_dual_efficiency_under_1[N]->size() > 1 )
                vector_sampling_for_dual_efficiency_under_1[N]->writeOnFile(file_name_dual);
        }

    }

    void doTheScmConvergenceStat( int sampling_size )
    {
        auto N = crb->scm()->KMax();
        std::list<std::string> list_error_type = boost::assign::list_of("RelativeError");
        BOOST_FOREACH( auto error_name, list_error_type)
        {
            std::ofstream conv;
            std::string file_name = "cvg-scm-"+ error_name +"-stats.dat";

            if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
            {
                conv.open(file_name, std::ios::app);
                conv << "Nb_basis" << "\t" << "Min" << "\t" << "Max" << "\t" << "Mean" << "\t" << "Variance" << "\n";
            }

            for(int j=0; j<N; j++)
            {
                double mean = M_mapConvSCM[error_name][j].mean();
                double variance = 0.0;
                for( int k=0; k < sampling_size; k++)
                    variance += (M_mapConvSCM[error_name][j](k) - mean)*(M_mapConvSCM[error_name][j](k) - mean)/sampling_size;

                if( Environment::worldComm().globalRank()  == Environment::worldComm().masterRank() )
                {
                    conv << j+1 << "\t"
                         << M_mapConvSCM[error_name][j].minCoeff() << "\t"
                         << M_mapConvSCM[error_name][j].maxCoeff() << "\t"
                         << mean << "\t" << variance << "\n";
                }
            }
            conv.close();
        }
    }


    // Script write current mu in cfg => need to write it in SamplingForTest
    void buildSamplingFromCfg()
    {
        // Size of mu
        int mu_size = model->parameterSpace()->dimension();

        // Clear SamplingForTest is exists, and open a new one
        fs::path input_file ("SamplingForTest");
        if( fs::exists(input_file) )
            std::remove( "SamplingForTest" );

        std::ofstream input( "SamplingForTest" );
        input << "mu= [ ";

        // Check if cfg file is readable
        std::ifstream cfg_file( option(_name="config-file").template as<std::string>() );
        if(!cfg_file)
            std::cout << "[Script-mode] Config file cannot be read" << std::endl;

        // OT writes values of mu in config file => read it and copy in SamplingForTest with specific syntax
        for(int i=1; i<=mu_size; i++)
            {
                // convert i into string
                std::ostringstream oss;
                oss << i;
                std::string is = oss.str();

                // Read cfg file, collect line with current mu_i
                std::string cfg_line_mu, tmp_content;
                std::ifstream cfg_file( option(_name="config-file").template as<std::string>() );
                while(cfg_file)
                    {
                        std::getline(cfg_file, tmp_content);
                        if(tmp_content.compare(0,2+is.size(),"mu"+is) == 0)
                            cfg_line_mu += tmp_content;
                    }

                //Regular expression : corresponds to one set in xml file (mu<i>=<value>)
                std::string expr_s = "mu"+is+"[[:space:]]*=[[:space:]]*([0-9]+(.?)[0-9]*(e(\\+|-)[0-9]+)?)[[:space:]]*";
                boost::regex expression( expr_s );

                //Match mu<i>=<value> in cfg file and copy to SamplingForTest
                boost::smatch what;
                auto is_match = boost::regex_match(cfg_line_mu, what, expression);
                //std::cout << "is match ?" << is_match << std::endl;
                if(is_match)
                    {
                        // what[0] is the complete string mu<i>=<value>
                        // what[1] is the submatch <value>
                        //std::cout << "what 1 = " << what[1] << std::endl;
                        if( i!=mu_size )
                            input << what[1] << " , ";
                        else
                            input << what[1] << " ]";
                    }
            }
    }

private:
    CRBModelMode M_mode;
    crbmodel_ptrtype model;
    crb_ptrtype crb;

    // For SCM convergence study
    std::map<std::string, std::vector<vectorN_type> > M_mapConvSCM;

    //vector of sampling to stock parameters for which the efficiency is under 1
    std::vector< sampling_ptrtype > vector_sampling_for_primal_efficiency_under_1;
    std::vector< sampling_ptrtype > vector_sampling_for_dual_efficiency_under_1;

    fs::path M_current_path;
}; // OpusApp

} // Feel

#endif /* __OpusApp_H */


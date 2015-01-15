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
   \author Cecile Daversin <daversin@math.unistra.fr>
   \author Stephane Veys
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
    RANDOM = 0, EQUIDISTRIBUTED = 1, LOGEQUIDISTRIBUTED = 2, READFROMCOMMANDLINE = 3
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
    typedef typename ModelType::mesh_type mesh_type;
    typedef std::vector< parameter_type > vector_parameter_type;

    typedef typename crb_type::sampling_ptrtype sampling_ptrtype;

    static const int nb_spaces = functionspace_type::nSpaces;
    typedef typename mpl::if_< boost::is_same< mpl::int_<nb_spaces> , mpl::int_<2> > , fusion::vector< mpl::int_<0>, mpl::int_<1> >  ,
                               typename mpl::if_ < boost::is_same< mpl::int_<nb_spaces> , mpl::int_<3> > ,
                                                   fusion::vector < mpl::int_<0> , mpl::int_<1> , mpl::int_<2> >,
                                                   typename mpl::if_< boost::is_same< mpl::int_<nb_spaces> , mpl::int_<4> >,
                                                                      fusion::vector< mpl::int_<0>, mpl::int_<1>, mpl::int_<2>, mpl::int_<3> >,
                                                                      fusion::vector< mpl::int_<0>, mpl::int_<1>, mpl::int_<2>, mpl::int_<3>, mpl::int_<4> >
                                                                      >::type >::type >::type index_vector_type;

    OpusApp()
        :
        super(),
        M_mode( ( CRBModelMode )ioption(_name=_o( this->about().appName(),"run.mode" )) )
        {
            this->init();
        }

    OpusApp( AboutData const& ad, po::options_description const& od )
        :
        super( ad, opusapp_options(ad.appName())
               .add( od )
               .add( crbOptions() )
               .add( eimOptions() )
               .add( podOptions() )),
        M_mode( ( CRBModelMode )ioption(_name=_o( this->about().appName(),"run.mode" )) )
        {
            this->init();
        }

    OpusApp( AboutData const& ad, po::options_description const& od, CRBModelMode mode )
        :
        super( ad, opusapp_options( ad.appName() )
               .add( od )
               .add( crbOptions() )
               .add( eimOptions() )
               .add( podOptions() )),
        M_mode( mode )
        {
            this->init();
        }

    OpusApp( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( ad, opusapp_options( ad.appName() )
               .add( od )
               .add( crbOptions() )
               .add( eimOptions() )
               .add( podOptions() )),
        M_mode( ( CRBModelMode )ioption(_name=_o( this->about().appName(),"run.mode" )) )
        {
            this->init();
        }
    OpusApp( int argc, char** argv, AboutData const& ad, po::options_description const& od, CRBModelMode mode )
        :
        super( ad, opusapp_options( ad.appName() )
               .add( od )
               .add( crbOptions() )
               .add( eimOptions() )
               .add( podOptions() )),
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
                    results_repo_name = soption(_name="crb.results-repo-name");
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

            bool only_master=boption(_name="crb.system-memory-evolution");
            bool all_procs  =boption(_name="crb.system-memory-evolution-on-all-procs");
            bool only_one_proc=only_master * ( Environment::worldComm().globalRank()==Environment::worldComm().masterRank() );
            bool write_memory_evolution = all_procs || only_one_proc ;

            bool crb_use_predefined = boption(_name="crb.use-predefined-WNmu");
            std::string file_name;
            int NlogEquidistributed = ioption(_name="crb.use-logEquidistributed-WNmu");
            int Nequidistributed = ioption(_name="crb.use-equidistributed-WNmu");
            int NlogEquidistributedScm = ioption(_name="crb.scm.use-logEquidistributed-C");
            int NequidistributedScm = ioption(_name="crb.scm.use-equidistributed-C");
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
                int Nrestart = ioption(_name="crb.restart-from-N");
                bool do_offline = false;
                int current_dimension = crb->dimension();
                int dimension_max = ioption(_name="crb.dimension-max");
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
            bool export_solution = boption(_name=_o( this->about().appName(),"export-solution" ));
            int proc_number =  Environment::worldComm().globalRank();

            bool load_elements_db= boption(_name="crb.load-elements-database");
            bool rebuild_db= boption(_name="crb.rebuild-database");

            int exportNameSize = ioption(_name="crb.export-name-max-size"); //paraview reads max 49 characters

            if ( this->vm().count( "help" ) )
            {
                std::cout << this->optionsDescription() << "\n";
                return;
            }

            //check options (does it make sens ?)
            bool option_checked=true;
            if( !load_elements_db && rebuild_db )
                option_checked=false;
            CHECK( option_checked )<<"options crb.load-elements-database : "<<load_elements_db<<" and crb.rebuild-database : "<<rebuild_db<<". If you don't want to load elements database maybe you want to apply RB approximation on a laptop wherease the RB was built on a super-computer ? If it's the case put crb.rebuild-database=false !! Else, you have to choose if you want to rebuild a RB database or to reload an existing one but not the elements database.\n";

            if( ! load_elements_db  )
            {
                M_mode = CRBModelMode::CRB_ONLINE;
                if( Environment::worldComm().isMasterRank() )
                {
                    std::cout<<"[OpusApp Information] You have choosen to reload an existing RB database without loading elments database. If the RB was built on an other computer make sure that database have been moved on in the right repositories.\n";
                }
            }

            this->loadDB();

            int run_sampling_size = ioption(_name=_o( this->about().appName(),"run.sampling.size" ));
            SamplingMode run_sampling_type = ( SamplingMode )ioption(_name=_o( this->about().appName(),"run.sampling.mode" ));
            int output_index = ioption(_name="crb.output-index");
            //int output_index = ioption(_name=_o(this->about().appName(),"output.index"));

            typename crb_type::sampling_ptrtype Sampling( new typename crb_type::sampling_type( model->parameterSpace() ) );

            int n_eval_computational_time = ioption(_name="eim.computational-time-neval");
            bool compute_fem = boption(_name="crb.compute-fem-during-online");
            bool compute_stat =  boption(_name="crb.compute-stat");

            bool use_predefined_sampling = boption(_name="crb.use-predefined-test-sampling");
            //bool select_parameter_via_one_feel=boption( _name="crb.select-parameter-via-one-feel");
            bool sampling_is_already_generated=false;


            if( ! load_elements_db )
            {
                compute_fem=false;
                export_solution=false;
                compute_stat=false;
            }
            bool select_parameter_via_one_feel=false;
            parameter_type user_mu_onefeel ( model->parameterSpace() );
            std::string string_parameters = soption(_name="crb.user-parameters");
            if( string_parameters == "" || string_parameters == " " )
            {
                select_parameter_via_one_feel=false;
            }
            else
            {
                select_parameter_via_one_feel=true;
            }
            if( select_parameter_via_one_feel )
            {
                run_sampling_size=1;
                Sampling->clear();
                //in this case we want to visualize RB solution with parameters from one feel interface
                compute_fem=false;
                compute_stat=false;

                //CHECK( load_elements_db )<<"[OpusApp] You have specified to not load elements database so it is impossible to visualize RB solution, use crb.load-elements-database=true !\n";

                //parameters are given by a vector of double
                std::vector< std::string > str;
                boost::split( str, string_parameters, boost::is_any_of(" "), boost::token_compress_on );
                double user_parameter_size = str.size();
                double mu_size = user_mu_onefeel.size();
                CHECK( user_parameter_size == mu_size )<<"[OpusApp] Error : parameters must have "<<mu_size<<" components and "<<user_parameter_size<<" have been given by the user \n";
                for(int i=0; i<mu_size; i++)
                {
                    double mu = boost::lexical_cast<double>( str[i] );
                    user_mu_onefeel( i ) = mu;
                }
                Sampling->addElement( user_mu_onefeel );
                sampling_is_already_generated=true;
            }

            std::string vary_only_parameter_components = soption(_name="crb.vary-only-parameter-components");
            std::vector< std::string > str;
            boost::split( str, vary_only_parameter_components, boost::is_any_of(" "), boost::token_compress_on );
            int number_str=str.size();
            CHECK( number_str < 5 )<<"Error when using option crb.vary-only-parameter-components, at maximum we can vary 2 components of the parameter";
            int vary_mu_comp0=-1,vary_mu_comp1=-1;
            int cutting_direction0=0;
            int cutting_direction1=0;
            bool vary_comp_time=false;
            if( vary_only_parameter_components!="" )
            {
                Sampling->clear();
                compute_fem=false;
                compute_stat=false;
                export_solution=false;
                int size=-1;
                if( number_str == 1 )
                {
                    //user want only to make time vary
                    CHECK( str[0] == "t" )<<"Error ! option crb.vary-only-parameter-components = "<<str[0]<<" but should be only 't' in this format";
                    vary_comp_time=true;
                }
                //here only one component vary
                if( number_str == 2 )
                {
                    if( str[0] == "t" )
                    {
                        //note that in this case we accept
                        //to have number of runs in this direction
                        //but we ignore it
                        vary_comp_time=true;
                    }
                    else
                    {
                        vary_mu_comp0 = boost::lexical_cast<int>( str[0] );
                    }
                    cutting_direction0 = boost::lexical_cast<int>( str[1] );
                }
                if( number_str == 3 )
                {
                    //in this configuration we have
                    //a parameter component and associated number of runs
                    //the time
                    if( str[0] == "t" )
                    {
                        vary_comp_time=true;
                        vary_mu_comp1 = boost::lexical_cast<int>( str[1] );
                        cutting_direction1 = boost::lexical_cast<int>( str[2] );
                    }
                    if( str[2] == "t" )
                    {
                        vary_comp_time=true;
                        vary_mu_comp1 = boost::lexical_cast<int>( str[0] );
                        cutting_direction1 =boost::lexical_cast<int>( str[1] );
                    }
                    if( str[0] != "t" && str[2] != "t" )
                    {
                       bool go=false;
                       CHECK( go )<<"A problem appears in the option crb.vary-only-parameter-components = "<<vary_only_parameter_components<<" No time 't' found !\n";
                    }
                }
                if( number_str==4 )
                {
                    if( str[0] == "t" )
                    {
                        vary_comp_time=true;
                        vary_mu_comp1 = boost::lexical_cast<int>( str[2] );
                        cutting_direction1 = boost::lexical_cast<int>( str[3] );
                    }
                    if( str[2] == "t" )
                    {
                        vary_comp_time=true;
                        vary_mu_comp1 = boost::lexical_cast<int>( str[0] );
                        cutting_direction1 = boost::lexical_cast<int>( str[1] );
                    }
                    if( str[0] != "t" && str[2] != "t" )
                    {
                        vary_mu_comp0 = boost::lexical_cast<int>( str[0] );
                        cutting_direction0 = boost::lexical_cast<int>( str[1] );
                        vary_mu_comp1 = boost::lexical_cast<int>( str[2] );
                        cutting_direction1 = boost::lexical_cast<int>( str[3] );
                    }
                }

                parameter_type user_mu ( model->parameterSpace() );
                double mu_size = user_mu.size();
                CHECK( vary_mu_comp0 < mu_size )<<"[OpusApp] error using crb.vary-only-parameter-components, the component "<<vary_mu_comp0<<" can't vary because parameter have a total of only "<<mu_size<<" components\n";
                if( number_str == 3 )
                {
                    CHECK( vary_mu_comp1 < mu_size )<<"[OpusApp] error using crb.vary-only-parameter-components, the component "<<vary_mu_comp1<<" can't vary because parameter have a total of only "<<mu_size<<" components\n";
                }
                //the sampling will b generated latter
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
            n_eval_computational_time = ioption(_name="crb.computational-time-neval");
            if( n_eval_computational_time > 0 )
            {
                if( ! boption(_name="crb.cvg-study") )
                {
                    compute_fem = false;
                    run_sampling_size = 0;
                }
                std::string appname = this->about().appName();
                //in the case we don't do the offline step, we need the affine decomposition
                model->computeAffineDecomposition();
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
                case SamplingMode::READFROMCOMMANDLINE:
                    std::vector<double> mu_list = option(_name=_o( this->about().appName(),"run.parameter" )).template as<std::vector<double>>();
                    parameter_type _mu = crb->Dmu()->element();
                    if ( crbmodel_type::ParameterSpaceDimension != mu_list.size() )
                        throw std::logic_error( "Parameter given by the command line option as note the expected size" );
                    Sampling->clear();
                    for ( int i=0 ; i<mu_list.size() ; i++ )
                        _mu(i)=mu_list[i];
                    Sampling->setElements( {_mu} );
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

                bool eim_mu_selection = boption(_name="eim.show-mu-selection");
                bool eim_t_selection = boption(_name="eim.show-t-selection");
                bool eim_offline_error = boption(_name="eim.show-offline-error");
                if( eim_mu_selection || eim_t_selection || eim_offline_error )
                {
                    auto eim_sc_vector = model->scalarContinuousEim();
                    auto eim_sd_vector = model->scalarDiscontinuousEim();

                    for(int i=0; i<eim_sc_vector.size(); i++)
                    {
                        auto eim = eim_sc_vector[i];
                        if( eim_mu_selection )
                            eim->printMuSelection();
                        if( eim_t_selection )
                            eim->printInterpolationPointsSelection();
                        if( eim_offline_error )
                            eim->printOfflineError();
                    }

                    for(int i=0; i<eim_sd_vector.size(); i++)
                    {
                        auto eim = eim_sd_vector[i];
                        if( eim_mu_selection )
                            eim->printMuSelection();
                        if( eim_t_selection )
                            eim->printInterpolationPointsSelection();
                        if( eim_offline_error )
                            eim->printOfflineError();
                    }
                }
            }

            auto e = exporter( _mesh= model->functionSpace()->mesh()  );

            printParameterHdr( ostr, model->parameterSpace()->dimension(), hdrs[M_mode] );

            int crb_error_type = ioption(_name="crb.error-type");

            int dim=0;
            if( M_mode==CRBModelMode::CRB )
            {
                dim=crb->dimension();
                if( crb->useWNmu() )
                    Sampling = crb->wnmu();

                if( boption(_name="crb.run-on-scm-parameters") )
                {
                    Sampling = crb->scm()->c();
                    if( crb_error_type!=1 )
                        throw std::logic_error( "[OpusApp] The SCM has not been launched, you can't use the option crb.run-on-scm-parameters. Run the SCM ( option crb.error-type=1 ) or comment this option line." );
                }
            }
            if( M_mode==CRBModelMode::SCM )
            {
                dim=crb->scm()->KMax();
                if( boption(_name="crb.scm.run-on-C") )
                    Sampling = crb->scm()->c();
            }

            std::ofstream file_summary_of_simulations( ( boost::format( "summary_of_simulations_%d" ) %dim ).str().c_str() ,std::ios::out | std::ios::app );

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
            if( boption(_name="crb.script-mode") )
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
            if(  use_predefined_sampling || boption(_name="crb.script-mode") )
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
            vectorN_type time_crb_vector_prediction;
            vectorN_type time_crb_vector_error_estimation;
            vectorN_type relative_estimated_error_vector;

            vectorN_type scm_relative_error;

            bool solve_dual_problem = boption(_name="crb.solve-dual-problem");

            if (boption(_name="crb.cvg-study") && compute_fem )
            {
                //int Nmax = crb->dimension();
                int Nmax = ioption("crb.dimension-max");
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
                time_crb_vector_prediction.resize( Sampling->size() );
                time_crb_vector_error_estimation.resize( Sampling->size() );

                if( crb->errorType()!=2 )
                    relative_estimated_error_vector.resize( Sampling->size() );

                crb->setOfflineStep( false );

                if (boption(_name="eim.cvg-study") )
                {
                    compute_fem=false;
                }

            }

            if( M_mode==CRBModelMode::SCM )
            {
                if (boption(_name="crb.scm.cvg-study") )
                    this->initializeConvergenceScmMap( Sampling->size() );

                scm_relative_error.resize( Sampling->size() );
            }

            int crb_dimension = ioption(_name="crb.dimension");
            int crb_dimension_max = ioption(_name="crb.dimension-max");
            double crb_online_tolerance = doption(_name="crb.online-tolerance");
            bool crb_compute_variance  = boption(_name="crb.compute-variance");

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
            std::vector< std::vector< std::vector< double > > > ref_betaLinearDecompositionAqm;
            for(int time_index=0; time_index<K; time_index++)
            {
                double time = time_index*dt;
                ref_betaAqm.push_back( model->computeBetaQm( ref_mu , time ).template get<1>() );
                ref_betaLinearDecompositionAqm.push_back( model->computeBetaLinearDecompositionA( ref_mu , time ) );
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
                if( boption(_name="crb.script-mode") )
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
                    //X[N-1] = ioption(_name="crb.compute-variance");
                    X[N-1] = 0;
                    bool compute_variance = crb_compute_variance;
                    if ( compute_variance )
                        X[N-1] = 1;

                    this->run( X.data(), X.size(), Y.data(), Y.size() );
                    //std::cout << "output = " << Y[0] << std::endl;

                    std::string resultFileName = soption(_name="result-file");
                    std::ofstream res(resultFileName);
                    res << "output="<< Y[0] << "\n";
                    res.close();
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
                                    {
                                        std::string exportName = u_fem.name().substr(0,exportNameSize) + "-" + std::to_string(curpar);
                                        e->add( exportName, u_fem );
                                    }
                                //model->solve( mu );
                                std::vector<double> o = boost::assign::list_of( model->output( output_index,mu , u_fem, true) )( ti.elapsed() );
                                if(proc_number == Environment::worldComm().masterRank() ) std::cout << "output=" << o[0] << "\n";
                                printEntry( ostr, mu, o );

                                std::ofstream res(soption(_name="result-file") );
                                res << "output="<< o[0] << "\n";

                            }
                            break;

                        case  CRBModelMode::CRB:
                            {
                                LOG(INFO) << "CRB mode\n";
                                if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
                                    std::cout << "CRB mode -- "<<curpar<<"/"<<sampling_size<<std::endl;


                                boost::mpi::timer ti;

                                LOG(INFO) << "solve crb\n";
                                //google::FlushLogFiles(google::GLOG_INFO);

                                //dimension of the RB (not necessarily the max)
                                int N =  ioption(_name="crb.dimension");

                                bool print_rb_matrix = boption(_name="crb.print-rb-matrix");
                                double online_tol = doption(_name="crb.online-tolerance");
                                vectorN_type time_crb;
                                ti.restart();

                                auto o = crb->run( mu, time_crb, online_tol , N, print_rb_matrix);
                                double time_crb_prediction=time_crb(0);
                                double time_crb_error_estimation=time_crb(1);

                                auto WN = crb->wn();
                                auto WNdu = crb->wndu();
                                //auto u_crb = crb->expansion( mu , N );
                                auto solutions=o.template get<2>();
                                auto uN = solutions.template get<0>();
                                auto uNdu = solutions.template get<1>();

                                int size = uN.size();

                                // Re-use uN given by lb in crb->run

                                u_crb = crb->expansion( uN[size-1] , N , WN );
                                if( solve_dual_problem )
                                    u_crb_dual = crb->expansion( uNdu[0] , N , WNdu );

                                std::ostringstream u_crb_str;
                                u_crb_str << "u_crb(" << mu_str.str() << ")";
                                u_crb.setName( u_crb_str.str()  );
                                LOG(INFO) << "export u_crb \n";
                                if( export_solution )
                                {
                                    if( select_parameter_via_one_feel )
                                    {
                                        model->adaptMesh( mu );
                                    }
                                    std::string exportName = u_crb.name().substr(0,exportNameSize) + "-" + std::to_string(curpar);
                                    e->add( exportName, u_crb );
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
                                auto vector_output_estimated_error = all_upper_bounds.template get<0>();
                                int last_time = vector_output_estimated_error.size()-1;
                                //output estimated error for last time
                                double output_estimated_error = vector_output_estimated_error[ last_time ];
                                double solution_estimated_error = all_upper_bounds.template get<1>();
                                double dual_solution_estimated_error = all_upper_bounds.template get<2>();

                                auto output_vector=o.template get<0>();
                                double output_vector_size=output_vector.size();
                                double ocrb = output_vector[output_vector_size-1];//output at last time
                                double time_fem_solve=-1;

                                if ( compute_fem )
                                {
									bool use_newton = boption(_name="crb.use-newton");

                                    LOG(INFO) << "solve u_fem\n";
                                    ti.restart();

                                    //auto u_fem = model->solveRB( mu );
                                    //auto u_fem = model->solveFemUsingOfflineEim( mu );

                                    if( boost::is_same<  crbmodel_type , crbmodelbilinear_type >::value )
                                    {
                                        if( boption(_name="crb.solve-fem-monolithic") )
                                        {
                                            u_fem = model->solveFemMonolithicFormulation( mu );
                                        }
                                        else
                                        {
                                            //use affine decomposition
                                            u_fem = model->solveFemUsingAffineDecompositionFixedPoint( mu );
                                        }
                                    }
                                    else
                                        u_fem = model->solve( mu );

                                    time_fem_solve=ti.elapsed();

                                    std::ostringstream u_fem_str;
                                    u_fem_str << "u_fem(" << mu_str.str() << ")";
                                    u_fem.setName( u_fem_str.str()  );

                                    if( export_solution )
                                    {
                                        LOG(INFO) << "export u_fem \n";
                                        std::string exportName = u_fem.name().substr(0,exportNameSize) + "-" + std::to_string(curpar);
                                        e->add( exportName, u_fem );
                                    }

                                    ti.restart();
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
                                        {
                                            std::string exportName = u_error.name().substr(0,exportNameSize) + "-" + std::to_string(curpar);
                                            e->add( exportName, u_error );
                                        }

                                    LOG(INFO) << "L2(fem)=" << l2Norm( u_fem )    << "\n";
                                    LOG(INFO) << "H1(fem)=" << h1Norm( u_fem )    << "\n";

                                    l2_error = l2Norm( u_error )/l2Norm( u_fem );
                                    h1_error = h1Norm( u_error )/h1Norm( u_fem );

                                    output_fem = ofem[0];
                                    time_fem = ofem[1]+time_fem_solve;

                                    if( boost::is_same<  crbmodel_type , crbmodelbilinear_type >::value )
                                    {
                                        if( solve_dual_problem )
                                        {
                                            if( boption(_name="crb.solve-fem-monolithic") )
                                            {
                                                u_dual_fem = model->solveFemDualMonolithicFormulation( mu );
                                            }
                                            else
                                            {
                                                //use affine decomposition
                                                u_dual_fem =  model->solveFemDualUsingAffineDecompositionFixedPoint( mu );
                                            }

                                            u_dual_error = model->functionSpace()->element();
                                            u_dual_error = (( u_dual_fem - u_crb_dual ).pow(2)).sqrt() ;
                                            l2_dual_error = l2Norm( u_dual_error )/l2Norm( u_dual_fem );
                                            h1_dual_error = h1Norm( u_dual_error )/h1Norm( u_dual_fem );
                                        }
                                    }

                                }//compute-fem-during-online

                                if ( crb->errorType()==2 )
                                {
                                    auto output_vector=o.template get<0>();
                                    double output_vector_size=output_vector.size();
                                    double ocrb = output_vector[output_vector_size-1];//output at last time
                                    std::vector<double> v = boost::assign::list_of( output_fem )( time_fem )( ocrb )( relative_estimated_error )( time_crb_prediction )( relative_error )( condition_number )( l2_error )( h1_error );

                                    if( proc_number == Environment::worldComm().masterRank() )
                                    {
                                        std::cout << "output=" << ocrb << " with " << o.template get<1>() << " basis functions\n";
                                        printEntry( file_summary_of_simulations, mu, v );
                                        printEntry( ostr, mu, v );
                                        //file_summary_of_simulations.close();

                                        if ( boption(_name="crb.compute-stat") && compute_fem )
                                        {
                                            relative_error_vector[curpar-1] = relative_error;
                                            l2_error_vector[curpar-1] = l2_error;
                                            h1_error_vector[curpar-1] = h1_error;
                                            time_fem_vector[curpar-1] = time_fem;
                                            time_crb_vector_prediction[curpar-1] = time_crb_prediction;
                                            time_crb_vector_error_estimation[curpar-1] = time_crb_error_estimation;
                                        }
                                        std::ofstream res(soption(_name="result-file") );
                                        res << "output="<< ocrb << "\n";
                                    }

                                }//end of crb->errorType==2
                                else
                                {
                                    //if( ! boost::is_same<  crbmodel_type , crbmodelbilinear_type >::value )
                                    //    throw std::logic_error( "ERROR TYPE must be 2 when using CRBTrilinear (no error estimation)" );

                                    auto output_vector=o.template get<0>();
                                    double output_vector_size=output_vector.size();
                                    double ocrb = output_vector[output_vector_size-1];//output at last time
                                    std::vector<double> v = boost::assign::list_of( output_fem )( time_fem )( ocrb )( relative_estimated_error )( time_crb_prediction )( relative_error )( condition_number )( l2_error )( h1_error );
                                    if( proc_number == Environment::worldComm().masterRank() )
                                    {
                                        std::cout << "output=" << ocrb << " with " << o.template get<1>() << " basis functions  (error estimation on this output : " << output_estimated_error<<") \n";
                                        //std::ofstream file_summary_of_simulations( ( boost::format( "summary_of_simulations_%d" ) % o.template get<2>() ).str().c_str() ,std::ios::out | std::ios::app );
                                        printEntry( file_summary_of_simulations, mu, v );
                                        printEntry( ostr, mu, v );
                                        //file_summary_of_simulations.close();

                                        if ( boption(_name="crb.compute-stat") && compute_fem )
                                        {
                                            relative_error_vector[curpar-1] = relative_error;
                                            l2_error_vector[curpar-1] = l2_error;
                                            h1_error_vector[curpar-1] = h1_error;
                                            time_fem_vector[curpar-1] = time_fem;
                                            time_crb_vector_prediction[curpar-1] = time_crb_prediction;
                                            time_crb_vector_error_estimation[curpar-1] = time_crb_error_estimation;
                                            relative_estimated_error_vector[curpar-1] = relative_estimated_error;
                                        }
                                        std::ofstream res(soption(_name="result-file") );
                                        res << "output="<< ocrb << "\n";
                                    }//end of proc==master
                                }//end of else (errorType==2)

                                if (boption(_name="eim.cvg-study") )
                                {
                                    bool check_name = false;
                                    std::string how_compute_unknown = soption(_name=_o( this->about().appName(),"how-compute-unkown-for-eim" ));
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

                                if (boption(_name="crb.cvg-study") && compute_fem )
                                {

                                    LOG(INFO) << "start convergence study...\n";
                                    std::map<int, boost::tuple<double,double,double,double,double,double,double> > conver;

                                    std::ofstream fileL2 ( "CrbConvergenceL2.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileH1 ( "CrbConvergenceH1.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileOutputError ( "CrbConvergenceOutputError.dat",std::ios::out | std::ios::app );
                                    std::ofstream fileOutputEstimatedError ( "CrbConvergenceOutputErrorEstimated.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileOutputErrorBoundEfficiency ( "CrbConvergenceOutputErrorBoundEfficiency.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileAbsoluteOutputErrorBoundEfficiency ( "CrbConvergenceAbsoluteOutputErrorBoundEfficiency.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileSolutionErrorBoundEfficiency ("CrbConvergencePrimalSolutionErrorBoundEfficiency.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileAbsoluteSolutionErrorBoundEfficiency ("CrbConvergencePrimalAbsoluteSolutionErrorBoundEfficiency.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileSolutionError ("CrbConvergencePrimalSolutionError.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileSolutionErrorEstimated ("CrbConvergencePrimalSolutionErrorEstimated.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileSolutionDualErrorBoundEfficiency ("CrbConvergenceDualSolutionErrorBoundEfficiency.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileAbsoluteSolutionDualErrorBoundEfficiency ("CrbConvergenceDualAbsoluteSolutionErrorBoundEfficiency.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileSolutionDualError ("CrbConvergenceDualSolutionError.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileSolutionDualErrorEstimated ("CrbConvergenceDualSolutionErrorEstimated.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream filePrimalResidualNorm ("CrbConvergencePrimalResidualNorm.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileDualResidualNorm ( "CrbConvergenceDualResidualNorm.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileCrbTimePrediction ( "CrbConvergenceCrbTimePrediction.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileCrbTimeErrorEstimation ( "CrbConvergenceCrbTimeErrorEstimation.dat" ,std::ios::out | std::ios::app );
                                    std::ofstream fileFemTime ( "CrbConvergenceFemTime.dat" ,std::ios::out | std::ios::app );

                                    int Nmax = ioption("crb.dimension-max");
                                    if( Environment::worldComm().isMasterRank() )
                                    {
                                            fileL2 << Nmax <<"\t";
                                            fileH1 << Nmax <<"\t";
                                            fileOutputError << Nmax <<"\t";
                                            fileOutputEstimatedError << Nmax << "\t";
                                            fileOutputErrorBoundEfficiency <<  Nmax << "\t";
                                            fileAbsoluteOutputErrorBoundEfficiency <<  Nmax << "\t";
                                            fileSolutionErrorBoundEfficiency << Nmax << "\t";
                                            fileAbsoluteSolutionErrorBoundEfficiency << Nmax << "\t";
                                            fileSolutionError << Nmax << "\t";
                                            fileSolutionErrorEstimated << Nmax << "\t";
                                            fileSolutionDualErrorBoundEfficiency << Nmax << "\t" ;
                                            fileAbsoluteSolutionDualErrorBoundEfficiency << Nmax << "\t" ;
                                            fileSolutionDualError << Nmax << "\t";
                                            fileSolutionDualErrorEstimated << Nmax << "\t";
                                            filePrimalResidualNorm << Nmax << "\t";
                                            fileDualResidualNorm << Nmax << "\t";
                                            fileCrbTimePrediction << Nmax << "\t";
                                            fileCrbTimeErrorEstimation << Nmax << "\t";
                                            fileFemTime << Nmax << "\t" ;
                                    }
                                    std::string str = "\t";
                                    vectorN_type crb_time;
                                    for( int N = 1; N <= Nmax ; N++ )
                                    {
                                        auto o= crb->run( mu, crb_time, online_tol , N, print_rb_matrix);

                                        auto output_vector=o.template get<0>();
                                        double output_vector_size=output_vector.size();
                                        double ocrb = output_vector[output_vector_size-1];//output at last time
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
                                        auto vector_output_estimated_error = all_upper_bounds.template get<0>();
                                        int last_time = vector_output_estimated_error.size()-1;
                                        //output estimated error for last time
                                        output_estimated_error = vector_output_estimated_error[ last_time ];
                                        solution_estimated_error = all_upper_bounds.template get<1>();
                                        dual_solution_estimated_error = all_upper_bounds.template get<2>();

                                        //auto o = crb->run( mu,  doption(_name="crb.online-tolerance") , N);
                                        double rel_err = std::abs( output_fem-ocrb ) /output_fem;
                                        double err = std::abs( output_fem-ocrb );

                                        double output_relative_estimated_error = output_estimated_error / output_fem;

                                        double primal_residual_norm = o.template get<4>();
                                        double dual_residual_norm = o.template get<5>();

                                        double solution_error=0;
                                        double dual_solution_error=0;
                                        double square_solution_error=0;
                                        double square_dual_solution_error=0;
                                        double ref_primal=0;
                                        double ref_dual=0;
                                        int qlinear=model->QLinearDecompositionA();
                                        bool symmetric = boption(_name="crb.use-symmetric-matrix");
                                        if( model->hasEim() && (qlinear > 0) )
                                        {

                                            if( model->isSteady() )
                                            {
                                                //all loops are not really necessary in the elliptic case
                                                //we could also use directly sqrt( model->scalarProduct( u_error , u_error ) ) ;
                                                //but we need to be sure that the matrix associated to scalar product
                                                //was assembled using reference parameter (or at least, to known what parameter was used).
                                                //Moreover when we deal with transient problems we need to build a( u_error^k , u_error^k ; muref )
                                                //where u_error^k means u_error at time index k
                                                //so it is not possible to do that only using model->scalarProduct()
                                                for(int q=0; q<qlinear;q++)
                                                {
                                                    for(int m=0; m<model->mMaxLinearDecompositionA(q); m++)
                                                    {
                                                        solution_error +=  ref_betaLinearDecompositionAqm[0][q][m]*model->linearDecompositionAqm(q,m,u_error,u_error) ;
                                                        ref_primal +=  ref_betaLinearDecompositionAqm[0][q][m]*model->linearDecompositionAqm(q,m,u_fem,u_fem);
                                                    }
                                                }
                                                if( solve_dual_problem )
                                                {
                                                    for(int q=0; q<model->QLinearDecompositionA();q++)
                                                    {
                                                        for(int m=0; m<model->mMaxLinearDecompositionA(q); m++)
                                                        {
                                                            dual_solution_error += ref_betaLinearDecompositionAqm[0][q][m]*model->linearDecompositionAqm(q,m,u_dual_error,u_dual_error);
                                                            ref_dual += ref_betaLinearDecompositionAqm[0][q][m]*model->linearDecompositionAqm(q,m,u_dual_fem,u_dual_fem);
                                                        }
                                                    }
                                                    square_dual_solution_error = dual_solution_error;
                                                    dual_solution_error = math::sqrt( dual_solution_error );
                                                    ref_dual = math::sqrt( ref_dual );
                                                }
                                                square_solution_error = solution_error;
                                                solution_error = math::sqrt( solution_error );
                                                ref_primal = math::sqrt( ref_primal );
                                            }//steady
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
                                                    for(int q=0; q<model->QLinearDecompositionA();q++)
                                                    {
                                                        for(int m=0; m<model->mMaxLinearDecompositionA(q); m++)
                                                        {
                                                            solution_error +=  ref_betaLinearDecompositionAqm[time_index][q][m]*model->linearDecompositionAqm(q,m,u_error,u_error) * dt;
                                                            ref_primal +=  ref_betaLinearDecompositionAqm[time_index][q][m]*model->linearDecompositionAqm(q,m,u_fem,u_fem) * dt;
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
                                                        for(int q=0; q<model->QLinearDecompositionA();q++)
                                                        {
                                                            for(int m=0; m<model->mMaxLinearDecompositionA(q); m++)
                                                            {
                                                                dual_solution_error +=  ref_betaLinearDecompositionAqm[time_index][q][m]*model->linearDecompositionAqm(q,m,u_dual_error,u_dual_error) * dt;
                                                                ref_dual +=  ref_betaLinearDecompositionAqm[time_index][q][m]*model->linearDecompositionAqm(q,m,u_dual_fem,u_dual_fem) * dt;
                                                            }
                                                        }
                                                    }
                                                    square_dual_solution_error = dual_solution_error;
                                                    dual_solution_error = math::sqrt( dual_solution_error );
                                                    ref_dual = math::sqrt( ref_dual );
                                                }//if solve-dual

                                            }//transient

                                        }//use EIM && qlinear > 0
                                        else
                                        {
                                            CHECK( symmetric ) << "Your model doesn' use a symmetric bilinear form a() so you have to implement computeLinearDecompositionA() function\n";
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
                                        }//no use EIM

                                        double l2_error = l2Norm( u_error )/l2Norm( u_fem );
                                        double h1_error = h1Norm( u_error )/h1Norm( u_fem );
                                        auto matrix_info = o.template get<3>();
                                        double condition_number = matrix_info.template get<0>();
                                        double output_error_bound_efficiency = 1;
                                        double absolute_output_error_bound_efficiency = 1;

                                        if( err > 1e-12 )
                                        {
                                            output_error_bound_efficiency = output_relative_estimated_error / rel_err;
                                            absolute_output_error_bound_efficiency = output_estimated_error / err;
                                        }
                                        double relative_primal_solution_error = solution_error / ref_primal ;
                                        double relative_primal_solution_estimated_error = solution_estimated_error / ref_primal;
                                        double relative_primal_solution_error_bound_efficiency=1;
                                        double absolute_primal_solution_error_bound_efficiency=1;
                                        if( solution_error > 1e-12 )
                                        {
                                            //compute bounds efficiency only if error is big enough
                                            relative_primal_solution_error_bound_efficiency = relative_primal_solution_estimated_error / relative_primal_solution_error;
                                            absolute_primal_solution_error_bound_efficiency = solution_estimated_error / solution_error;
                                        }

                                        // if( relative_primal_solution_error_bound_efficiency < 1 )
                                        // {
                                        //     vector_sampling_for_primal_efficiency_under_1[N-1]->push_back( mu , 1);
                                        // }

                                        double relative_dual_solution_error = 1;
                                        double relative_dual_solution_estimated_error = 1;
                                        double relative_dual_solution_error_bound_efficiency = 1;
                                        double absolute_dual_solution_error_bound_efficiency = 1;
                                        if( solve_dual_problem )
                                        {
                                            relative_dual_solution_error = dual_solution_error / ref_dual ;
                                            relative_dual_solution_estimated_error = dual_solution_estimated_error / ref_dual;
                                            if( dual_solution_error > 1e-12 )
                                            {
                                                relative_dual_solution_error_bound_efficiency = relative_dual_solution_estimated_error / relative_dual_solution_error;
                                                absolute_dual_solution_error_bound_efficiency = dual_solution_estimated_error / dual_solution_error;
                                            }
                                            // if( relative_dual_solution_error_bound_efficiency < 1 )
                                            // {
                                            //     vector_sampling_for_dual_efficiency_under_1[N-1]->push_back( mu , 0);
                                            // }
                                        }
                                        conver[N]=boost::make_tuple( rel_err, l2_error, h1_error , relative_estimated_error, condition_number , output_error_bound_efficiency , relative_primal_solution_error_bound_efficiency );

                                        //LOG(INFO) << "N=" << N << " " << rel_err << " " << l2_error << " " << h1_error << " " <<condition_number<<"\n";
                                        if ( proc_number == Environment::worldComm().masterRank() )
                                        {
                                            std::cout << "N=" << N << " Output =  "<< output_fem <<" OutputError = "<<rel_err <<" OutputErrorEstimated = "<<relative_estimated_error
                                                      <<"  L2Error = "<< l2_error << "  H1Error = " << h1_error <<std::endl;

                                            if( N == Nmax )
                                                str="\n";
                                            fileL2 << l2_error <<str;
                                            fileH1 << h1_error <<str;
                                            fileOutputError << rel_err <<str;
                                            fileOutputEstimatedError << output_relative_estimated_error << str;
                                            fileOutputErrorBoundEfficiency <<  output_error_bound_efficiency << str;
                                            fileAbsoluteOutputErrorBoundEfficiency <<  absolute_output_error_bound_efficiency << str;
                                            fileSolutionErrorBoundEfficiency << relative_primal_solution_error_bound_efficiency << str;
                                            fileAbsoluteSolutionErrorBoundEfficiency << absolute_primal_solution_error_bound_efficiency << str;
                                            fileSolutionError << relative_primal_solution_error << str;
                                            fileSolutionErrorEstimated << relative_primal_solution_estimated_error << str;
                                            fileSolutionDualErrorBoundEfficiency << relative_dual_solution_error_bound_efficiency << str ;
                                            fileAbsoluteSolutionDualErrorBoundEfficiency << absolute_dual_solution_error_bound_efficiency << str ;
                                            fileSolutionDualError << relative_dual_solution_error << str;
                                            fileSolutionDualErrorEstimated <<  relative_dual_solution_estimated_error << str;
                                            filePrimalResidualNorm << primal_residual_norm << str;
                                            fileDualResidualNorm <<  dual_residual_norm << str;
                                            fileCrbTimePrediction<< crb_time(0) << str;
                                            fileCrbTimeErrorEstimation<< crb_time(1) << str;
                                            fileFemTime<< time_fem << str;
                                        }
                                        if( boption(_name="crb.compute-matrix-information") )
                                        {
                                            auto matrix_info = o.template get<3>();// conditioning of primal reduced matrix + determinant
                                            double conditioning = matrix_info.template get<0>();
                                            double determinant = matrix_info.template get<1>();
                                            LOG( INFO ) << " primal reduced matrix information ";
                                            LOG( INFO ) << std::setprecision(15)<<"mu : \n"<<mu;
                                            LOG( INFO ) << std::setprecision(15)<<"conditioning : "<<conditioning;
                                            LOG( INFO ) << std::setprecision(15)<<"determinant : "<<determinant;
                                        }
                                        // if( relative_primal_solution_error_bound_efficiency < 1 )
                                        // {
                                        //     LOG( INFO ) << "N : "<<N;
                                        //     LOG( INFO ) << std::setprecision(15)<<"efficiency of error estimation on primal solution is "<<relative_primal_solution_error_bound_efficiency<<" ( should be >= 1 )";
                                        //     LOG( INFO ) << std::setprecision(15)<<"mu : \n"<<mu;
                                        //     LOG( INFO ) << std::setprecision(15)<<"relative_primal_solution_estimated_error : "<<relative_primal_solution_estimated_error;
                                        //     LOG( INFO ) << std::setprecision(15)<<"relative_primal_solution_error : "<<relative_primal_solution_error;
                                        //     LOG( INFO ) << std::setprecision(15)<<"primal_solution_estimated_error : "<<solution_estimated_error;
                                        //     LOG( INFO ) << std::setprecision(15)<<"primal_solution_error : "<<solution_error;
                                        //     LOG( INFO ) << std::setprecision(15)<<"square error : "<<square_solution_error;
                                        //     //LOG( INFO ) << std::setprecision(15)<<"u_crb : \n"<<u_crb[size-1];
                                        //     LOG( INFO ) << std::setprecision(15)<<"primal solution norme  : "<<uN.l2Norm();
                                        // }
                                        // if( relative_dual_solution_error_bound_efficiency < 1 )
                                        // {
                                        //     LOG( INFO ) <<std::setprecision(15)<< "efficiency of error estimation on dual solution is "<<relative_dual_solution_error_bound_efficiency<<" ( should be >= 1 )";
                                        //     LOG( INFO ) <<std::setprecision(15)<< "mu : \n"<<mu;
                                        //     LOG( INFO ) <<std::setprecision(15)<<"relative_dual_solution_estimated_error : "<<relative_dual_solution_estimated_error;
                                        //     LOG( INFO ) <<std::setprecision(15)<<"relative_dual_solution_error : "<<relative_dual_solution_error;
                                        //     LOG( INFO ) <<std::setprecision(15)<<"dual_solution_estimated_error : "<<dual_solution_estimated_error;
                                        //     LOG( INFO ) <<std::setprecision(15)<<"dual_solution_error : "<<dual_solution_error;
                                        //     LOG( INFO ) << std::setprecision(15)<<"square error : "<<square_dual_solution_error;
                                        //     //LOG( INFO ) << std::setprecision(15)<<"u_crb_du : \n"<<u_crb_du[0];
                                        //     LOG( INFO ) << std::setprecision(15)<<"dual solution norme  : "<<uNdu.l2Norm();
                                        // }
                                        // if( output_error_bound_efficiency < 1 )
                                        // {
                                        //     LOG( INFO ) <<std::setprecision(15)<<"efficiency of error estimation on output is "<<output_error_bound_efficiency<<" ( should be >= 1 )";
                                        //     LOG( INFO ) <<std::setprecision(15)<< "mu : \n"<<mu;
                                        //     LOG( INFO ) <<std::setprecision(15)<< "output_relative_estimated_error : "<<output_relative_estimated_error;
                                        //     LOG( INFO ) <<std::setprecision(15)<< "output_relative_error : "<<rel_err;
                                        //     LOG( INFO ) <<std::setprecision(15)<< "output_estimated_error : "<<output_estimated_error;
                                        //     LOG( INFO ) <<std::setprecision(15)<< "output_error : "<< std::abs( output_fem-ocrb );
                                        // }

                                        if ( boption(_name="crb.check.residual")  && solve_dual_problem  )
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
                                vectorN_type time_crb;
                                auto o = crb->run( mu, time_crb, doption(_name="crb.online-tolerance") );
                                auto output_vector=o.template get<0>();
                                double output_vector_size=output_vector.size();
                                double ocrb = output_vector[output_vector_size-1];//output at last time

                                if ( crb->errorType()==2 )
                                    {
                                        std::vector<double> v = boost::assign::list_of( ocrb )( ti.elapsed() );
                                        std::cout << "output=" << ocrb << " with " << o.template get<1>() << " basis functions\n";
                                        printEntry( ostr, mu, v );
                                    }

                                else
                                    {
                                        auto all_upper_bounds = o.template get<6>();
                                        auto vector_output_estimated_error = all_upper_bounds.template get<0>();
                                        int last_time = vector_output_estimated_error.size()-1;
                                        //output estimated error for last time
                                        double output_estimated_error = vector_output_estimated_error[ last_time ];
                                        double relative_estimated_error = output_estimated_error / output_fem;
                                        auto output_vector = o.template get<0>();
                                        double output_vector_size = output_vector.size();
                                        double output = output_vector[ output_vector_size-1 ];
                                        std::vector<double> v = boost::assign::list_of( output )( output_estimated_error )( ti.elapsed() );
                                        std::cout << "output=" << ocrb << " with " << o.template get<1>() <<
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

                                if (boption(_name="crb.scm.cvg-study")  )
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

            //one parameter change : plot
            if( vary_mu_comp0 > -1 || vary_comp_time )
            {
                parameter_type user_mu ( model->parameterSpace() );
                double mu_size = user_mu.size();
                std::vector< int > sampling_each_direction ( mu_size );
                for(int i=0; i<mu_size; i++)
                {
                    if( i == vary_mu_comp0  )
                        sampling_each_direction[i]=cutting_direction0;
                    else
                        sampling_each_direction[i]=0;
                }

                int N =  ioption(_name="crb.dimension");
                double online_tol = doption(_name="crb.online-tolerance");

                vectorN_type outputs_storage;
                vectorN_type mu0_storage;
                vectorN_type estimated_error_outputs_storage;

                //have min/max
                Sampling->equidistribute( 2 );
                vectorN_type time_crb;

                curpar=1;
                if( ! vary_comp_time )
                {
                    outputs_storage.resize( cutting_direction0 );
                    mu0_storage.resize( cutting_direction0 );
                    auto mu_=Sampling->min().template get<0>();
                    if( select_parameter_via_one_feel )
                    {
                        mu_ = user_mu_onefeel;
                    }
                    estimated_error_outputs_storage.resize( cutting_direction0 );
                    Sampling->logEquidistributeProduct( sampling_each_direction , mu_ );
                    BOOST_FOREACH( auto mu, *Sampling )
                    {
                        double x = mu(vary_mu_comp0);
                        double mu0 = mu(vary_mu_comp0);
                        auto o = crb->run( mu, time_crb, online_tol , N);
                        auto output_vector=o.template get<0>();
                        double output_vector_size=output_vector.size();
                        double ocrb = output_vector[output_vector_size-1];//output at last time
                        auto all_upper_bounds = o.template get<6>();
                        auto vector_output_estimated_error = all_upper_bounds.template get<0>();
                        int last_time = vector_output_estimated_error.size()-1;
                        //output estimated error for last time
                        double output_estimated_error = vector_output_estimated_error[ last_time ];
                        outputs_storage(curpar-1)=ocrb;
                        mu0_storage(curpar-1)=mu0;
                        estimated_error_outputs_storage(curpar-1)=output_estimated_error;
                        curpar++;
                    }
                }
                else
                {
                    vectorN_type time_crb;
                    auto mu=Sampling->min().template get<0>();
                    if( select_parameter_via_one_feel )
                    {
                        mu = user_mu_onefeel;
                    }
                    auto o = crb->run( mu, time_crb, online_tol , N);
                    auto output_vector=o.template get<0>();
                    auto all_upper_bounds = o.template get<6>();
                    auto vector_output_estimated_error = all_upper_bounds.template get<0>();
                    int size=output_vector.size();
                    outputs_storage.resize( size );
                    mu0_storage.resize( size );
                    estimated_error_outputs_storage.resize( size );
                    double dt = model->timeStep();
                    double time=0;
                    for(int i=0; i<size; i++)
                    {
                        mu0_storage(i)=time;
                        outputs_storage(i)=output_vector[i];
                        estimated_error_outputs_storage(i)=vector_output_estimated_error[i];
                        time+=dt;
                    }
                }
                model->generateGeoFileForOutputPlot( outputs_storage , mu0_storage, estimated_error_outputs_storage );
            }
            //two parameters change : response surface
            bool draw_surface=false;
            if( vary_mu_comp0 > -1 && vary_mu_comp1 > -1 )
                draw_surface=true;
            if( vary_mu_comp1 > -1 && vary_comp_time )
                draw_surface=true;
            if( draw_surface )
            {

                double Ti=model->timeInitial();
                double Tf=model->timeFinal();
                double dt=model->timeStep();

                int N =  ioption(_name="crb.dimension");
                double online_tol = doption(_name="crb.online-tolerance");
                CHECK( Environment::worldComm().globalSize() == 1 )<<"implemented only in sequential (because of dof filling)\n";
                typename crb_type::sampling_ptrtype S( new typename crb_type::sampling_type( model->parameterSpace() ) );
                bool all_procs_have_same_sampling=true;
                S->equidistribute( 2 , all_procs_have_same_sampling );
                auto min = S->min().template get<0>();
                auto max = S->max().template get<0>();
                //vector of indices of components vary
                std::vector<int> components_vary(2);
                components_vary[0]=vary_mu_comp0;
                components_vary[1]=vary_mu_comp1;
                //vector containing parameters min and max
                std::vector<parameter_type> extremums(2);
                extremums[0]=min;
                extremums[1]=max;
                //cutting in each direction
                std::vector<int> cutting(2);
                cutting[0]=cutting_direction0;
                cutting[1]=cutting_direction1;
                std::vector<double> time_cutting(3);
                time_cutting[0]=Ti;
                time_cutting[1]=Tf;
                time_cutting[2]=dt;

                auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                            _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                            _desc = model->createStructuredGrid( components_vary, extremums , cutting, time_cutting, vary_comp_time ) );

                auto Vh = Pch<1>( mesh );
                auto px = project(_space=Vh, _expr=Px() );
                auto py = project(_space=Vh, _expr=Py() );
                int ndof=Vh->nDof();
                auto u = Vh->element();
                parameter_type mu ( model->parameterSpace() );
                int time_index=0;
                vectorN_type time_crb;
                for(int i=0; i<ndof; i++)
                {
                    double mux=px(i);
                    double muy=py(i);
                    mu=min;//first, initialize all components at minimum
                    mu( vary_mu_comp0 ) = mux;

                    if( vary_comp_time )
                    {
                        //time is always represented by the X-axis
                        //look at which index mux is associated
                        time_index=(mux-Ti)/dt;
                    }
                    else
                    {
                        mu( vary_mu_comp1 ) = muy;
                    }

                    auto run = crb->run( mu, time_crb, online_tol, N );
                    auto output_vector=run.template get<0>();
                    double output=0;
                    if( vary_comp_time )
                    {
                        output = output_vector[time_index];
                    }
                    else
                    {
                        int output_vector_size=output_vector.size();
                        output = output_vector[output_vector_size-1];//output at last time
                    }
                    u(i)=output;
                }
                auto exp = exporter( _mesh=mesh );
                exp->add( "response surface", u );
                exp->save();

                /* Export geo file */
                if(soption(_name="exporter.format") == "gmsh")
                {
                    std::cout << exp->prefix() << std::endl;
                    std::ofstream of;
                    std::ostringstream oss;
                    oss.str("");
                    // To get the correct filename, we use the timeset name (The prefix is not used for exporting the data file)
                    // We might have to get the right one, if there are several ones
                    oss << exp->defaultTimeSet()->name() << "-" <<  Environment::worldComm().size() << "_" << Environment::worldComm().rank() << ".geo";
                    of.open(oss.str());

                    if(of.is_open())
                    {
                        oss.str("");
                        oss << exp->defaultTimeSet()->name() << "-" <<  Environment::worldComm().size() << "_" << Environment::worldComm().rank() << ".msh";

                        of << "Merge \"" << oss.str() << "\";" << std::endl;
                        of << "vid = PostProcessing.NbViews;" << std::endl;
                        //of << "For i In {vid-N:vid-1}" << std::endl;
                        of << "View[vid-1].Axes = 1;" << std::endl;
                        //of << "EndFor" << std::endl;
                    }

                    of.close();
                }
            }

            //model->computationalTimeEimStatistics();
            if( export_solution )
                e->save();

            if( proc_number == Environment::worldComm().masterRank() ) std::cout << ostr.str() << "\n";

            if (boption(_name="eim.cvg-study") && M_mode==CRBModelMode::CRB)
                this->doTheEimConvergenceStat( Sampling->size() );

            if (boption(_name="crb.cvg-study") && compute_fem && M_mode==CRBModelMode::CRB )
                this->doTheCrbConvergenceStat( Sampling->size() );

            if (boption(_name="crb.scm.cvg-study") && M_mode==CRBModelMode::SCM )
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
                Eigen::MatrixXf::Index index_max_time_crb_prediction;
                Eigen::MatrixXf::Index index_min_time_crb_prediction;
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
                double max_time_crb_prediction = time_crb_vector_prediction.maxCoeff(&index_max_time_crb_prediction);
                double min_time_crb_prediction = time_crb_vector_prediction.minCoeff(&index_min_time_crb_prediction);
                double mean_time_crb_prediction = time_crb_vector_prediction.mean();
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
                    file_summary_of_simulations <<"max of time CRB : "<<max_time_crb_prediction<<" at the "<<index_max_time_crb_prediction+1<<"^th simulation\n";
                    file_summary_of_simulations <<"min of time CRB : "<<min_time_crb_prediction<<" at the "<<index_min_time_crb_prediction+1<<"^th simulation\n";
                    file_summary_of_simulations <<"mean of time CRB : "<<mean_time_crb_prediction<<"\n\n";
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
        int N = ioption("crb.dimension-max");
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
        ComputeNormL2InCompositeCase compute_normL2_in_composite_case( u );
        index_vector_type index_vector;
        fusion::for_each( index_vector, compute_normL2_in_composite_case );
        return compute_normL2_in_composite_case.norm();
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
                auto mesh = u.functionSpace()->mesh();
                double norm  = normL2(_range=elements( mesh ),_expr=( idv(u) ) );
                M_vec(i)= norm ;
            }

        double norm()
            {
                return M_vec.sum();
            }

        mutable vectorN_type M_vec;
        element_type M_composite_u;
    };


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

            int max=eim->mMax();

            if( max > 1 )
            {
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

        for(int i=0; i<eim_sd_vector.size(); i++)
        {
            auto eim = eim_sd_vector[i];

            int max=eim->mMax();

            if( max > 1 )
            {
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
    }

    void doTheCrbConvergenceStat( int sampling_size )
    {
        int N = ioption("crb.dimension-max");

        std::vector< vectorN_type > L2, H1, OutputError, OutputErrorEstimated, OutputErrorBoundEfficiency, AbsoluteOutputErrorBoundEfficiency;
        std::vector< vectorN_type > PrimalSolutionError, PrimalSolutionErrorEstimated, PrimalSolutionErrorBoundEfficiency , PrimalAbsoluteSolutionErrorBoundEfficiency;
        std::vector< vectorN_type > DualSolutionError, DualSolutionErrorEstimated, DualSolutionErrorBoundEfficiency, DualAbsoluteSolutionErrorBoundEfficiency;
        std::vector< vectorN_type > CrbTimePrediction, CrbTimeErrorEstimation, FemTime;

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
        filename = "CrbConvergenceAbsoluteOutputErrorBoundEfficiency.dat";
        model->readConvergenceDataFromFile( AbsoluteOutputErrorBoundEfficiency, filename );
        filename = "CrbConvergencePrimalSolutionError.dat";
        model->readConvergenceDataFromFile( PrimalSolutionError , filename );
        filename = "CrbConvergencePrimalSolutionErrorEstimated.dat";
        model->readConvergenceDataFromFile( PrimalSolutionErrorEstimated , filename );
        filename = "CrbConvergencePrimalSolutionErrorBoundEfficiency.dat";
        model->readConvergenceDataFromFile( PrimalSolutionErrorBoundEfficiency , filename );
        filename = "CrbConvergencePrimalAbsoluteSolutionErrorBoundEfficiency.dat";
        model->readConvergenceDataFromFile( PrimalAbsoluteSolutionErrorBoundEfficiency , filename );
        filename = "CrbConvergenceDualSolutionError.dat";
        model->readConvergenceDataFromFile( DualSolutionError , filename );
        filename = "CrbConvergenceDualSolutionErrorEstimated.dat";
        model->readConvergenceDataFromFile( DualSolutionErrorEstimated , filename );
        filename = "CrbConvergenceDualSolutionErrorBoundEfficiency.dat";
        model->readConvergenceDataFromFile( DualSolutionErrorBoundEfficiency , filename );
        filename = "CrbConvergenceDualAbsoluteSolutionErrorBoundEfficiency.dat";
        model->readConvergenceDataFromFile( DualAbsoluteSolutionErrorBoundEfficiency , filename );
        filename = "CrbConvergenceCrbTimePrediction.dat";
        model->readConvergenceDataFromFile( CrbTimePrediction , filename );
        filename = "CrbConvergenceCrbTimeErrorEstimation.dat";
        model->readConvergenceDataFromFile( CrbTimeErrorEstimation , filename );
        filename = "CrbConvergenceFemTime.dat";
        model->readConvergenceDataFromFile( FemTime , filename );

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
        filename = "cvg-crb-AbsoluteOutputErrorBoundEfficiency.dat";
        model->writeConvergenceStatistics( AbsoluteOutputErrorBoundEfficiency , filename);
        filename = "cvg-crb-PrimalSolutionError.dat";
        model->writeConvergenceStatistics( PrimalSolutionError , filename);
        filename = "cvg-crb-PrimalSolutionErrorEstimated.dat";
        model->writeConvergenceStatistics( PrimalSolutionErrorEstimated , filename);
        filename = "cvg-crb-PrimalSolutionErrorBoundEfficiency.dat";
        model->writeConvergenceStatistics( PrimalSolutionErrorBoundEfficiency , filename);
        filename = "cvg-crb-PrimalAbsoluteSolutionErrorBoundEfficiency.dat";
        model->writeConvergenceStatistics( PrimalAbsoluteSolutionErrorBoundEfficiency , filename);
        filename = "cvg-crb-DualSolutionError.dat";
        model->writeConvergenceStatistics( DualSolutionError , filename);
        filename = "cvg-crb-DualSolutionErrorEstimated.dat";
        model->writeConvergenceStatistics( DualSolutionErrorEstimated , filename);
        filename = "cvg-crb-DualSolutionErrorBoundEfficiency.dat";
        model->writeConvergenceStatistics( DualSolutionErrorBoundEfficiency , filename);
        filename = "cvg-crb-DualAbsoluteSolutionErrorBoundEfficiency.dat";
        model->writeConvergenceStatistics( DualAbsoluteSolutionErrorBoundEfficiency , filename);
        filename = "cvg-crb-CrbTimePrediction.dat";
        model->writeConvergenceStatistics( CrbTimePrediction , filename, "totaltime");
        filename = "cvg-crb-CrbTimeErrorEstimation.dat";
        model->writeConvergenceStatistics( CrbTimeErrorEstimation , filename, "totaltime");
        filename = "cvg-crb-FemTime.dat";
        model->writeConvergenceStatistics( FemTime , filename , "totaltime" );

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

#if 0
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
#endif
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
        std::ifstream cfg_file( soption(_name="config-file") );
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
                std::ifstream cfg_file( soption(_name="config-file") );
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

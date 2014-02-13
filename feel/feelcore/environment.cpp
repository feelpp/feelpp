/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-04-14

  Copyright (C) 2010,2011 Université Joseph Fourier (Grenoble I)

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
   \file environment.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-04-14
 */
#include <cstdlib>
#include <pwd.h>

#include <boost/program_options.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/assign/std/vector.hpp>
#include <gflags/gflags.h>

#if defined ( FEELPP_HAS_PETSC_H )
#include <petscsys.h>
#endif

#include <feel/feelconfig.h>
#include <feel/feelcore/feel.hpp>


#include <feel/feelcore/environment.hpp>

#include <feel/feelcore/feelpetsc.hpp>
#include <feel/options.hpp>

namespace GiNaC {
extern void cleanup_ex( bool verbose );
}
namespace detail
{
class Env{
public:
    static std::string getUserName()
        {
            struct passwd *pw;
            uid_t uid;
            int c;

            uid = geteuid ();
            pw = getpwuid (uid);
            if (pw)
            {
                return std::string(pw->pw_name);
            }
            return std::string("");
        }
};
}
namespace google
{
namespace glog_internal_namespace_
{
bool IsGoogleLoggingInitialized();
}
}
namespace Feel
{

namespace detail
{
FEELPP_NO_EXPORT
std::pair<std::string, std::string>
at_option_parser_2( std::string const&s )
{
    if ( '@' == s[0] )
        return std::make_pair( std::string( "response-file" ), s.substr( 1 ) );

    else
        return std::pair<std::string, std::string>();
}

/**
   \fn freeargv -- free an argument vector

   void freeargv (char** vector)

   Free an argument vector that was built using dupargv.  Simply scans
   through the vector, freeing the memory for each argument until the
   terminating NULL is found, and then frees the vector itself.

*/

void freeargv (char** vector)
{
    char **scan;

    if (vector != NULL)
    {
        for (scan = vector; *scan != NULL; scan++)
        {
            free (*scan);
        }
        free (vector);
    }
}


/**
   \fn dupargv -- duplicate an argument vector

   char **dupargv (char** vector)

   Duplicate an argument vector.  Simply scans through the
   vector, duplicating each argument until the
   terminating NULL is found.

   \return a pointer to the argument vector if
   successful. Returns NULL if there is insufficient memory to
   complete building the argument vector.
*/
char **
dupargv (char** argv)
{
  int argc;
  char **copy;

  if (argv == NULL)
    return NULL;

  /* the vector */
  for (argc = 0; argv[argc] != NULL; argc++);
  copy = (char **) malloc ((argc + 1) * sizeof (char *));
  if (copy == NULL)
      return NULL;

  /* the strings */
  for (argc = 0; argv[argc] != NULL; argc++)
  {
      int len = strlen (argv[argc]);
      copy[argc] = (char*)malloc (sizeof (char *) * (len + 1));
      if (copy[argc] == NULL)
      {
          freeargv (copy);
          return NULL;
      }
      strcpy (copy[argc], argv[argc]);
  }
  copy[argc] = NULL;
  return copy;
}

void
Environment::generateOLFiles( int argc, char** argv, std::string const& appName)
{
  //Application path
  bool isNum = false;
  fs::path p(argv[0]);
  std::ostringstream appPath;

  appPath.str("");

  /*
  std::cout << fs::current_path() << std::endl;
  std::cout << appName << std::endl;
  std::cout << appPath.str() << std::endl;
  */

  /* get app name */
  appPath.str("");
  appPath << fs::absolute(p).string();
  //std::cout << appPath.str() << std::endl;

  std::ostringstream optionPath;

  std::ofstream ol;
  ol.open(appPath.str() + ".ol", std::ofstream::out | std::ofstream::trunc); //.ol file
  std::ofstream cfgol;
  cfgol.open(appPath.str() + ".onelab.cfg.ol", std::ofstream::out | std::ofstream::trunc); //.cfg.ol file
  const std::vector<boost::shared_ptr<po::option_description>>& listOpt = Environment::optionsDescription().options(); //Fetching options in an option_description vector

  for(boost::shared_ptr<po::option_description> option : listOpt)
  {
      //Informations about the option
      std::string optName = option->format_name().erase(0,2);//Putting the option name in a variable for easier manipulations
      std::string defVal = ""; //option->format_parameter(); //Putting the option default value in a variable for easier manipulations
      std::string desc=option->description(); // Option description

      // reset option path
      optionPath.str("");
      // reset type
      isNum = false;

      //std::cout << optName << ";" << defVal << ";" << desc << std::endl;

      std::string ens,funcName,test;
      /*
      if(defVal.size() > 3) // We want defVal.size() > 3 to exclude options without a default value
      {
          defVal.erase(0,6); // Removing "arg" and a space and the first parenthesis
          defVal.erase(defVal.size()-1,1); // Removing closing parenthesis
      }
      else
      {
          defVal = "";
      }
      */

      std::vector<std::string> strings; //Vector of the split name
      boost::split(strings,optName,boost::is_any_of(".")); //Spliting option name
      ens = "";
      if(strings.size() > 1)
      {
          ens = strings[0] + "/"; //Getting the first split element for the option set in the .cfg.ol file
          for(size_t i = 1; i < strings.size() - 1; i++) //Getting the option set
          {
            ens += strings[i] + "/";
          }
      }
      funcName = strings[strings.size() - 1]; //Raw option name

      //std::cout << ens << " " << funcName << std::endl;

      /* skip some options */
      if(funcName == "generate-ol" ||
         funcName == "config-file")
      {
          continue;
      }

      /* if an option has been set either through command line */
      /* or through the initial config file */
      /* we use its configuration */
      //if(S_vm.count(optName) && !(S_vm[optName].defaulted()))
      if(S_vm.count(optName))
      {
          std::ostringstream oss;
          oss.str("");
          //std::cout << defVal;
          
          //std::cout << "Entry for " << optName << ": ";
          if(S_vm[optName].defaulted())
          {
            //std::cout << "defaulted ";
            optionPath << "GeneralParameters/" << ens;
          }
          else
          {
            optionPath << "DefinedParameters/" << ens;
          }
          
          //optionPath << "Parameters/" << ens;

          if(optName == "licence")
          {
              const std::type_info & ti = S_vm[optName].value().type();
              std::cout << ti.name() << " " << std::endl;
          }

          if(S_vm[optName].empty())
          {
            //std::cout << "empty ";
          }
          else
          {
              const std::type_info & ti = S_vm[optName].value().type();
              //std::cout << ti.name() << " ";
              if(ti == typeid(bool))
              {
                  oss.str("");
                  oss << (S_vm[optName].as<bool>() ? "1" : "0");
                  isNum = false;
              }
              else if(ti == typeid(int))
              {
                  oss.str("");
                  oss << S_vm[optName].as<int>();
                  isNum = true;
              }
              else if(ti == typeid(size_type))
              {
                  oss.str("");
                  oss << S_vm[optName].as<size_type>();
                  isNum = true;
              }
              else if(ti == typeid(float))
              {
                  oss.str("");
                  oss << S_vm[optName].as<float>();
                  isNum = true;
              }
              else if(ti == typeid(double))
              {
                  oss.str("");
                  oss << S_vm[optName].as<double>();
                  isNum = true;
              }
              else if(ti == typeid(std::string))
              {
                  oss.str("");
                  oss <<  S_vm[optName].as<std::string>();
                  isNum = false;
              }
              else 
              {
                  std::cout << "Unknown type for parameter " << optName << "(" << typeid(void).name() << ")" << std::endl;
                  isNum = false;
              }
          }
          //std::cout << oss.str() << std::endl;
          defVal = oss.str();

      /* Force Gmsh as a the default exporter */
      /* as we are using OneLab */
      if(ens == "exporter/" && funcName == "format")
      {
          defVal = "gmsh";
      }

      /*
      if(defVal.size() != 0) //Excluding options without a default value
      {
      */
          if(isNum)
          {
              ol << funcName << ".number(" << defVal << ", " << optionPath.str() << ");" << " # "<< desc << std::endl;
              cfgol << optName << "=OL.get(" << optionPath.str() << funcName << ")" << std::endl;
          }
          else 
          {
              ol << funcName << ".string(" << defVal << ", " << optionPath.str() << ");" << " # "<< desc << std::endl;
              cfgol << optName << "=OL.get(" << optionPath.str() << funcName << ")" << std::endl;
          }
      //}
      }

  }

  ol << "" << std::endl;
  
  /* Mesher instructions */
  ol << "Mesher.register(native,/usr/bin/gmsh);" << std::endl;
  /*
  ol << "Mesher.in(feel.geo);" << std::endl;
  ol << "Mesher.out(feel.msh);" << std::endl;
  */
  ol << "OL.if(OL.get(Parameters/gmsh/filename) == untitled.geo)" << std::endl;
  ol << "OL.msg(No geo file specified. Using a default one);" << std::endl;
  ol << "OL.endif" << std::endl;
  
  /* Application instructions */
  //ol << "FeelApp.register(interfaced, mpirun -np " << worldComm().size() << " " + appPath.str() + ");" << std::endl;
  ol << "FeelApp.register(interfaced, " + appPath.str() + ");" << std::endl;
  ol << "FeelApp.in(OL.get(Arguments/FileName).onelab.cfg.ol);" << std::endl;
  ol << "FeelApp.run( --config-file OL.get(Arguments/FileName).onelab.cfg --nochdir );" << std::endl;
  ol << "FeelApp.out(" + appName + "-1_0.msh);" << std::endl;

  /* Take into account th fact that we use MPI */
  /* Generating multiple output files, so we merge them */
  //for(int i = 0; i < worldComm().size(); i++)
  //{
    ol << "FeelApp.merge(" + appName + "-1_0.msh);" << std::endl;
  //}
  
  ol.close();
  cfgol.close();
}

void
Environment::processGenericOptions()
{
//     // leave this to subclasses or users
// #if 0
//     if ( S_vm.count( "help" ) )
//         std::cout << S_desc << "\n";

// #endif


//     if ( S_vm.count( "response-file" ) )
//     {
//         using namespace std;
//         // Load the file and tokenize it
//         ifstream ifs( S_vm["response-file"].as<std::string>().c_str() );

//         if ( !ifs )
//         {
//             cout << "Could not open the response file\n";
//             return ;
//         }

//         // Read the whole file into a string
//         stringstream ss;
//         ss << ifs.rdbuf();
//         // Split the file content
//         boost::char_separator<char> sep( " \n\r" );
//         boost::tokenizer<boost::char_separator<char> > tok( ss.str(), sep );
//         vector<string> args;
//         copy( tok.begin(), tok.end(), back_inserter( args ) );

//         parseAndStoreOptions( po::command_line_parser( args ) );
//     }

    if ( worldComm().isMasterRank() )
    {

        if ( S_vm.count( "feelinfo" ) )
            std::cout << std::setw( 15 ) << std::right << "Feel Version : " << Info::versionString() << "\n"
                      << std::setw( 15 ) << std::right << "Major : " << Info::versionMajor() << "\n"
                      << std::setw( 15 ) << std::right << "Minor : " << Info::versionMinor() << "\n"
                      << std::setw( 15 ) << std::right << "Micro : " << Info::versionMicro() << "\n"
                      << std::setw( 15 ) << std::right << "Revision : " << Info::revision() << "\n"
                      << std::setw( 15 ) << std::right << "BuildId : " << Info::buildId() << "\n"
                      << std::setw( 15 ) << std::right << "Feel Prefix : " << Info::prefix() << "\n"
                      << std::setw( 15 ) << std::right << "Feel DataDir : " << Info::datadir() << "\n";

        if ( S_vm.count( "verbose" ) ||
             S_vm.count( "help" ) ||
             S_vm.count( "version" ) ||
             S_vm.count( "copyright" ) ||
             S_vm.count( "license" ) ||
             S_vm.count( "authors" ) )
        {
            std::cout << S_about.appName() << ": " << S_about.shortDescription() <<  "\n";
        }

        if ( S_vm.count( "version" ) )
        {
            std::cout << " version : " << S_about.version() << "\n";
        }

        if ( S_vm.count( "copyright" ) )
        {
            std::cout << " copyright : " << S_about.copyrightStatement() << "\n";
        }

        if ( S_vm.count( "license" ) )
        {
            std::cout << " license : " << S_about.license() << "\n";
        }

        if ( S_vm.count( "authors" ) )
        {
#if 0
            std::cout << std::setw( 30 )
                      << "Author Name"
                      << " " << std::setw( 15 )
                      << "Task"
                      << " " << std::setw( 40 )
                      << "Email Address"
                      << "\n";
            std::cout << std::setw( 85+3 ) << std::setfill( '-' ) << "\n" << std::setfill( ' ' );
            std::for_each( S_about.authors().begin(),
                           S_about.authors().end(),
                           std::cout
                           << std::setw( 30 )
                           << lambda::bind( &AboutPerson::name,
                                            lambda::_1 )
                           << " " << std::setw( 15 )
                           << lambda::bind( &AboutPerson::task,
                                            lambda::_1 )
                           << " " << std::setw( 40 )
                           << lambda::bind( &AboutPerson::emailAddress,
                                            lambda::_1 )
                           << "\n" );
#endif
        }

        if ( S_vm.count( "help" ) )
        {
            std::cout << optionsDescription() << "\n";
        }
    }
    if ( S_vm.count( "verbose" ) ||
         S_vm.count( "help" ) ||
         S_vm.count( "version" ) ||
         S_vm.count( "copyright" ) ||
         S_vm.count( "license" ) ||
         S_vm.count( "authors" ) )
    {
        if ( Environment::initialized() )
        {
            worldComm().barrier();
            MPI_Finalize();
        }
        exit(0);
    }
#if 0
    std::cout << "count = " << S_vm.count( "debug" ) << "\n"
              << "string = " << S_vm["debug"].as<std::string>() << "\n";
#endif

    VLOG(2) << "[processGenericOptions] done\n";
}

void
Environment::parseAndStoreOptions( po::command_line_parser parser, bool extra_parser )
{
    VLOG(2) << " parsing options...\n";

    boost::shared_ptr<po::parsed_options> parsed;

    if ( extra_parser )
    {
        parsed = boost::shared_ptr<po::parsed_options>( new po::parsed_options( parser
                                                                                .options( *S_desc )
                                                                                .extra_parser( at_option_parser_2 )
                                                                                .allow_unregistered()
                                                                                .run() ) );
    }

    else
    {
        parsed = boost::shared_ptr<po::parsed_options>( new po::parsed_options( parser
                                                                                .options( *S_desc )
                                                                                .allow_unregistered()
                                                                                .run() ) );
    }

    VLOG(2) << "[parseAndStoreOptions] parsing options done\n";

    S_to_pass_further = po::collect_unrecognized( parsed->options, po::include_positional );
    VLOG(2)<< " number of unrecognized options: " << ( S_to_pass_further.size() ) << "\n";

    BOOST_FOREACH( std::string const& s, S_to_pass_further )
    {
        VLOG(2)<< " option: " << s << "\n";
    }
    std::vector<po::basic_option<char> >::iterator it = parsed->options.begin();
    std::vector<po::basic_option<char> >::iterator en  = parsed->options.end();

    for ( ; it != en ; ++it )
        if ( it->unregistered )
        {
            VLOG(2)<< " remove from vector " << it->string_key << "\n";
            parsed->options.erase( it );
        }

    po::store( *parsed, S_vm );
}


void
Environment::doOptions( int argc, char** argv, po::options_description const& desc, std::string const& appName )
{
    //std::locale::global(std::locale(""));
    try
    {
        parseAndStoreOptions( po::command_line_parser( argc, argv ), true );
        processGenericOptions();

        VLOG(2) << "options parsed and stored in database";
        /**
         * parse config file if given to command line
         */
        if ( S_vm.count( "config-file" ) )
        {

            if ( fs::exists(  S_vm["config-file"].as<std::string>() ) )
            {
                LOG(INFO) << "Reading " << S_vm["config-file"].as<std::string>() << "...";

                std::ifstream ifs( S_vm["config-file"].as<std::string>().c_str() );
                po::store( parse_config_file( ifs, desc, true ), S_vm );
                po::notify( S_vm );
            }
        }
        else
        {
            using namespace boost::assign;
            std::vector<fs::path> prefixes = S_paths;
#if 0
            prefixes += boost::assign::list_of( fs::current_path() )
                ( fs::path ( Environment::localConfigRepository() ) )
                ( fs::path ( Environment::systemConfigRepository().get<0>() ) )
                ( fs::path ( "/usr/share/feel/config" ) )
                ( fs::path ( "/usr/local/share/feel/config" ) )
                ( fs::path ( "/opt/local/share/feel/config" ) );
#endif
            char* env;
            env = getenv("FEELPP_DIR");
            if (env != NULL && env[0] != '\0')
            {
                prefixes.push_back( fs::path( env ) );
            }

            VLOG(2) << "try processing cfg files...\n";
            std::string config_name;
            bool found = false;
            BOOST_FOREACH( auto prefix, prefixes )
            {
                config_name = ( boost::format( "%1%/%2%.cfg" ) % prefix.string() % appName ).str();
                VLOG(2) << " Looking for " << config_name << "\n";

                if ( fs::exists( config_name ) )
                {
                    found = true;
                    break;
                }
                else
                {
                    // try with a prefix feel_
                    config_name = ( boost::format( "%1%/feelpp_%2%.cfg" ) % prefix.string() % appName ).str();
                    VLOG(2) << " Looking for " << config_name << "\n";
                    if ( fs::exists( config_name ) )
                    {
                        found = true;
                        break;
                    }
                }
            }
            if ( found )
            {
                LOG(INFO) << "Reading  " << config_name << "...\n";
                std::ifstream ifs( config_name.c_str() );
                store( parse_config_file( ifs, desc, true ), S_vm );
                LOG(INFO) << "Reading  " << config_name << " done.\n";
                //po::store(po::parse_command_line(argc, argv, desc), S_vm);
                po::notify( S_vm );
            }
        }


        /* handle the generation of onelab files after having processed */
        /* the regular config file, so we have parsed user defined parameters */
        /* or restored a previous configuration */
        if ( worldComm().isMasterRank() )
        {
            if ( S_vm.count( "generate-ol" ) )
            {
                Environment::generateOLFiles( argc, argv, appName );
            }
        }
        if ( S_vm.count( "generate-ol" ) )
        {
            if ( Environment::initialized() )
            {
                worldComm().barrier();
                MPI_Finalize();
            }
            exit(0);
        }
        
    }
    // catches program_options exceptions
    catch ( boost::program_options::multiple_occurrences const& e )
    {
        LOG(WARNING) << "Command line or config file option parsing error: " << e.what() << "\n"
                     << "  o faulty option: " << e.get_option_name() << "\n"
                     << "Warning: the .cfg file or some options may not have been read properly\n";
    }

    catch ( boost::program_options::ambiguous_option const& e )
    {
        LOG(WARNING) << "Command line or config file option parsing error: " << e.what() << "\n"
                     << "  o faulty option: " << e.get_option_name() << "\n"
                     << "  o possible alternatives: " ;
        std::for_each( e.alternatives().begin(), e.alternatives().end(), []( std::string const& s )
                       {
                           LOG(WARNING) << s << " ";
                       } );
        LOG(WARNING) << "\n"
                     << "Warning: the .cfg file or some options may not have been read properly\n";
    }

    // catches program_options exceptions
    catch ( std::exception& e )
    {
        LOG(WARNING) << "Application option parsing: unknown option:" << e.what() << " (the .cfg file or some options may not have been read properly)\n";
    }

    catch ( ... )
    {
        LOG(WARNING) << "Application option parsing: unknown exception triggered  (the .cfg file or some options may not have been read properly)\n";
    }
}

Environment::Environment()
    :
    M_env(false)
{
    // Initialize Google's logging library.
    if ( !google::glog_internal_namespace_::IsGoogleLoggingInitialized() )
        google::InitGoogleLogging("feel++");
    google::InstallFailureSignalHandler();
#if defined( FEELPP_HAS_TBB )
    int n = tbb::task_scheduler_init::default_num_threads();
    VLOG(2) << "[Feel++] TBB running with " << n << " threads\n";
#else
    int n = 1 ;
#endif

    //tbb::task_scheduler_init init(1);

    S_worldcomm = worldcomm_type::New();
    CHECK( S_worldcomm ) << "Environment : creating worldcomm failed\n";
    S_worldcommSeq.reset( new WorldComm(S_worldcomm->subWorldCommSeq()) );

#if defined ( FEELPP_HAS_PETSC_H )
    PetscTruth is_petsc_initialized;
    PetscInitialized( &is_petsc_initialized );

    if ( !is_petsc_initialized )
    {
        i_initialized = true;
        int ierr = PetscInitializeNoArguments();

        boost::ignore_unused_variable_warning( ierr );
        CHKERRABORT( *S_worldcomm,ierr );
    }

    // make sure that petsc do not catch signals and hence do not print long
    //and often unuseful messages
    PetscPopSignalHandler();
#endif // FEELPP_HAS_PETSC_H


}

fs::path scratchdir()
{
    const char* env;
    // if scratsch dir not defined, define it
    env = getenv("FEELPP_SCRATCHDIR");
    if (env == NULL || env[0] == '\0')
    {
        env = getenv("SCRATCHDIR");
        if (env != NULL && env[0] != '\0')
        {
            std::string value = (boost::format("%1%/%2%/feelpp/") % env % ::detail::Env::getUserName()).str();
            setenv("FEELPP_SCRATCHDIR", (boost::format("%1%/%2%/feelpp/") % env % ::detail::Env::getUserName() ).str().c_str(),0);
        }
        else
        {
            env = getenv("SCRATCH");
            if (env != NULL && env[0] != '\0')
            {
                std::string value = (boost::format("%1%/%2%/feelpp/") % env % ::detail::Env::getUserName()).str();
                setenv("FEELPP_SCRATCHDIR", (boost::format("%1%/%2%/feelpp/") % env % ::detail::Env::getUserName() ).str().c_str(),0);
            }
            else
            {
                std::string value = (boost::format("/tmp/%1%/feelpp/") % ::detail::Env::getUserName()).str();
                setenv("FEELPP_SCRATCHDIR", value.c_str(),0);
            }
        }
    }
    env = getenv("FEELPP_SCRATCHDIR");
    if (env != NULL && env[0] != '\0')
    {
        return fs::path( env );
    }
    std::string value = (boost::format("/tmp/%1%/feelpp/") % ::detail::Env::getUserName()).str();
    return fs::path( value );
}

Environment::Environment( int& argc, char**& argv )
    :
    M_env( argc, argv, false )
{
    S_scratchdir = scratchdir()/fs::path(argv[0]).filename();
    if ( !fs::exists( S_scratchdir ) )
        fs::create_directories( S_scratchdir );
    FLAGS_log_dir=S_scratchdir.string();

    google::AllowCommandLineReparsing();
    google::ParseCommandLineFlags(&argc, &argv, false);

    // Initialize Google's logging library.
    if ( !google::glog_internal_namespace_::IsGoogleLoggingInitialized() )
    {
        if ( argc > 0 )
            google::InitGoogleLogging(argv[0]);
        else
            google::InitGoogleLogging("feel++");
    }
    google::InstallFailureSignalHandler();

#if defined( FEELPP_HAS_TBB )
    int n = tbb::task_scheduler_init::default_num_threads();
    //int n = 2;
    //VLOG(2) << "[Feel++] TBB running with " << n << " threads\n";
    //tbb::task_scheduler_init init(2);
#endif

    S_worldcomm = worldcomm_type::New();
    CHECK( S_worldcomm ) << "Feel++ Environment: creang worldcomm failed!";
    S_worldcommSeq.reset( new WorldComm(S_worldcomm->subWorldCommSeq()) );

#if defined ( FEELPP_HAS_PETSC_H )
    PetscTruth is_petsc_initialized;
    PetscInitialized( &is_petsc_initialized );

    if ( !is_petsc_initialized )
    {
        i_initialized = true;
#if defined( FEELPP_HAS_SLEPC )
        int ierr = SlepcInitialize( &argc,&argv, PETSC_NULL, PETSC_NULL );
#else
        int ierr = PetscInitialize( &argc, &argv, PETSC_NULL, PETSC_NULL );
#endif
        boost::ignore_unused_variable_warning( ierr );
        CHKERRABORT( *S_worldcomm,ierr );
    }

    // make sure that petsc do not catch signals and hence do not print long
    //and often unuseful messages
    PetscPopSignalHandler();
#endif // FEELPP_HAS_PETSC_H

}

void
Environment::init( int argc, char** argv, po::options_description const& desc, AboutData const& about )
{
    S_worldcomm = worldcomm_type::New();
    CHECK( S_worldcomm ) << "Feel++ Environment: creang worldcomm failed!";
    S_worldcommSeq.reset( new WorldComm(S_worldcomm->subWorldCommSeq()) );

    S_scratchdir = scratchdir();
    fs::path a0 = std::string(argv[0]);
    const int Nproc = 200;
    if ( S_worldcomm->size() > Nproc )
    {
        std::string smin = boost::lexical_cast<std::string>( Nproc*std::floor(S_worldcomm->rank()/Nproc) );
        std::string smax = boost::lexical_cast<std::string>( Nproc*std::ceil(double(S_worldcomm->rank()+1)/Nproc)-1 );
        std::string replog = smin + "-" + smax;
        S_scratchdir/= a0.filename()/replog;
    }
    else
        S_scratchdir/= a0.filename();
    // only one processor every Nproc creates the corresponding log directory
    if ( S_worldcomm->rank() % Nproc == 0 )
    {
        if ( !fs::exists( S_scratchdir ) )
            fs::create_directories( S_scratchdir );
    }
    FLAGS_log_dir=S_scratchdir.string();

    // duplicate argv before passing to gflags because gflags is going to
    // rearrange them and it screws badly the flags for PETSc/SLEPc
    char** envargv = dupargv( argv );

    google::AllowCommandLineReparsing();
    google::ParseCommandLineFlags(&argc, &argv, false);

#if 0
    std::cout << "argc=" << argc << "\n";
    for(int i = 0; i < argc; ++i )
    {
        std::cout << "argv[" << i << "]=" << argv[i] << "\n";
    }
#endif
    // Initialize Google's logging library.
    if ( !google::glog_internal_namespace_::IsGoogleLoggingInitialized() )
    {
        if ( FLAGS_no_log )
        {
            if ( S_worldcomm->rank() == 0 && FLAGS_no_log == 1 )
                FLAGS_no_log = 0;
            google::InitGoogleLogging(argv[0]);
        }
        else if ( argc > 0 )
            google::InitGoogleLogging(argv[0]);
        else
            google::InitGoogleLogging("feel++");
    }
    google::InstallFailureSignalHandler();
#if defined( FEELPP_HAS_TBB )
    int n = tbb::task_scheduler_init::default_num_threads();
    //int n = 2;
    //VLOG(2) << "[Feel++] TBB running with " << n << " threads\n";
    //tbb::task_scheduler_init init(2);
#endif

#if defined ( FEELPP_HAS_PETSC_H )
    PetscTruth is_petsc_initialized;
    PetscInitialized( &is_petsc_initialized );

    if ( !is_petsc_initialized )
    {
        i_initialized = true;
#if defined( FEELPP_HAS_SLEPC )
        int ierr = SlepcInitialize( &argc,&envargv, PETSC_NULL, PETSC_NULL );
#else
        int ierr = PetscInitialize( &argc, &envargv, PETSC_NULL, PETSC_NULL );
#endif
        boost::ignore_unused_variable_warning( ierr );
        CHKERRABORT( *S_worldcomm,ierr );
    }

    // make sure that petsc do not catch signals and hence do not print long
    //and often unuseful messages
    PetscPopSignalHandler();
#endif // FEELPP_HAS_PETSC_H


    S_about = about;
    doOptions( argc, envargv, *S_desc, about.appName() );

    // make sure that we pass the proper verbosity level to glog
    if ( S_vm.count("v") )
        FLAGS_v = S_vm["v"].as<int>();

    freeargv( envargv );

}
void
Environment::clearSomeMemory()
{
    Environment::logMemoryUsage( "Environment::clearSomeMemory before:" );

    // send signal to all deleters
    S_deleteObservers();
    google::FlushLogFiles(google::GLOG_INFO);
    VLOG(2) << "clearSomeMemory: delete signal sent" << "\n";

    Environment::logMemoryUsage( "Environment::clearSomeMemory after:" );
}
Environment::~Environment()
{
    VLOG(2) << "[~Environment] sending delete to all deleters" << "\n";

    Environment::clearSomeMemory();

    if ( i_initialized )
    {
        VLOG(2) << "clearing known paths\n";
        S_paths.clear();

        VLOG(2) << "[~Environment] cleaning up global excompiler\n";

        GiNaC::cleanup_ex( false );

        VLOG(2) << "[~Environment] finalizing slepc,petsc and mpi\n";
#if defined ( FEELPP_HAS_PETSC_H )
        PetscTruth is_petsc_initialized;
        PetscInitialized( &is_petsc_initialized );

        if ( is_petsc_initialized )
        {
#if defined( FEELPP_HAS_SLEPC )
            SlepcFinalize();
#else
            PetscFinalize();
#endif // FEELPP_HAS_SLEPC
        }

#endif // FEELPP_HAS_PETSC_H

        google::ShutdownGoogleLogging();
    }


}



bool
Environment::initialized()
{

#if defined( FEELPP_HAS_PETSC_H )
    PetscTruth is_petsc_initialized;
    PetscInitialized( &is_petsc_initialized );
    return mpi::environment::initialized() && is_petsc_initialized ;
#else
    return mpi::environment::initialized() ;
#endif
}

bool
Environment::finalized()
{
#if defined( FEELPP_HAS_PETSC_H )
    PetscTruth is_petsc_initialized;
    PetscInitialized( &is_petsc_initialized );
    return mpi::environment::finalized() && !is_petsc_initialized;
#else
    return mpi::environment::finalized();
#endif
}


std::string
Environment::rootRepository()
{
    char * senv = ::getenv( "FEELPP_REPOSITORY" );
    if ( senv != NULL && senv[0] != '\0' )
    {
        return std::string(senv);
    }
    senv = ::getenv( "FEELPP_WORKDIR" );
    if ( senv != NULL && senv[0] != '\0' )
    {
        return std::string( senv );
    }
    senv = ::getenv( "WORK" );
    if ( senv != NULL && senv[0] != '\0' )
    {
        return std::string( senv )+"/feel";
    }
    senv = ::getenv( "WORKDIR" );
    if ( senv != NULL && senv[0] != '\0' )
    {
        return std::string( senv )+"/feel";
    }
    senv = ::getenv( "HOME" );
    if ( senv != NULL && senv[0] != '\0' )
    {
        return std::string( senv ) + "/feel";
    }
    return std::string();
}
std::string
Environment::findFile( std::string const& filename )
{
    fs::path cp = fs::current_path();

    fs::path p( filename );
    if ( p.is_absolute() && fs::exists( p ) )
    {
        LOG(INFO) << "File " << filename << " found";
        return filename;
    }

#if 0
    // first try in the current path
    if ( fs::exists( cp / filename ) )
    {
        LOG(INFO) << "File " << (cp/filename) << " found";
        return ( cp/filename ).string();
    }
#endif

    // look in to paths list from end-1 to begin
    auto it = std::find_if( S_paths.rbegin(), S_paths.rend(),
                            [&filename] ( fs::path const& p ) -> bool
                            {
                                if ( fs::exists( p/filename) )
                                    return true;
                                return false;
                            } );
    if ( it != S_paths.rend() )
    {
        LOG(INFO) << "File " << (*it/filename) << " found";
        return ( *it / filename ).string();
    }

    if ( fs::exists( cp / filename ) )
    {
        LOG(INFO) << "File " << (cp/filename) << " found";
        return ( cp/filename ).string();
    }

    if ( fs::path( filename ).extension() == ".geo" || fs::path( filename ).extension() == ".msh" )
    {
        if ( fs::exists( fs::path( Environment::localGeoRepository() ) / filename ) )
        {
            LOG(INFO) << "File " << (fs::path( Environment::localGeoRepository() ) / filename ) << " found";
            return ( fs::path( Environment::localGeoRepository() ) / filename ).string();
        }

        if ( Environment::systemGeoRepository().get<1>()  &&
             fs::exists( fs::path( Environment::systemGeoRepository().get<0>() ) / filename ) )
        {
            LOG(INFO) << "File" << ( fs::path( Environment::systemGeoRepository().get<0>() ) / filename ) << " found";
            return ( fs::path( Environment::systemGeoRepository().get<0>() ) / filename ).string();
        }
    }
    LOG(INFO) << "File " << filename << " not found";
    return std::string();
}
std::vector<std::string>
Environment::geoPathList()
{
    std::vector<std::string> plist;
    plist.push_back( fs::current_path().string() );
    std::for_each( S_paths.rbegin(), S_paths.rend(),
                   [&plist] ( fs::path const& p )
                   {
                       plist.push_back( p.string() );
                   } );
    if ( fs::exists( Environment::localGeoRepository() ) )
        plist.push_back( Environment::localGeoRepository() );
    if ( Environment::systemGeoRepository().get<1>()  &&
         fs::exists( Environment::systemGeoRepository().get<0>() ) )
        plist.push_back( Environment::systemGeoRepository().get<0>() );
    return plist;
}
std::string
Environment::localGeoRepository()
{
    fs::path rep_path;

    rep_path = Environment::rootRepository();
    rep_path /= "geo";

    if ( !fs::exists( rep_path ) )
        fs::create_directory( rep_path );

    return rep_path.string();
}
boost::tuple<std::string,bool>
Environment::systemGeoRepository()
{
    fs::path rep_path;

    rep_path = Info::prefix();
    rep_path /= "share/feel/geo";
    return boost::make_tuple( rep_path.string(), fs::exists( rep_path ) );
}

std::string
Environment::localConfigRepository()
{
    fs::path rep_path;

    rep_path = Environment::rootRepository();
    rep_path /= "config";

    if ( !fs::exists( rep_path ) )
        fs::create_directory( rep_path );

    return rep_path.string();
}
boost::tuple<std::string,bool>
Environment::systemConfigRepository()
{
    fs::path rep_path;

    rep_path = Info::prefix();
    rep_path /= "share/feel/config";
    return boost::make_tuple( rep_path.string(), fs::exists( rep_path ) );
}

void
Environment::changeRepositoryImpl( boost::format fmt, std::string const& logfilename, bool add_subdir_np )
{
    if ( Environment::vm().count( "nochdir" ) )
        return;

    fs::path rep_path;
    S_paths.push_back( fs::current_path() );

    typedef std::vector< std::string > split_vector_type;

    split_vector_type dirs; // #2: Search for tokens
    std::string fmtstr = fmt.str();
    boost::split( dirs, fmtstr, boost::is_any_of( "/" ) );

    fs::path p = dirs.front();
    if ( p.relative_path() != ".")
        rep_path = Environment::rootRepository();

    if ( !fs::exists( rep_path ) )
    {
        LOG(INFO) << "Creating directory " << rep_path << "...";
        fs::create_directory( rep_path );
    }


    BOOST_FOREACH( std::string const& dir, dirs )
    {
        //VLOG(2)<< " option: " << s << "\n";
        rep_path = rep_path / dir;

        if ( !fs::exists( rep_path ) )
            fs::create_directory( rep_path );
    }
    if ( add_subdir_np )
    {
        rep_path = rep_path / (boost::format( "np_%1%" ) % Environment::numberOfProcessors() ).str();
        if ( !fs::exists( rep_path ) )
            fs::create_directory( rep_path );
        LOG(INFO) << "changing directory to " << rep_path << "\n";
    }
    ::chdir( rep_path.string().c_str() );

    setLogs( logfilename );
}

#if 0
po::variables_map
Environment::vm( po::options_description const& desc )
{
    po::variables_map vm;
    po::store( po::parse_command_line( 0, ( char** )0, desc ), vm );
    po::notify( vm );

    return vm;
}
#endif

void
Environment::setLogs( std::string const& prefix )
{


}

std::vector<WorldComm> const&
Environment::worldsComm( int n )
{
    CHECK( S_worldcomm ) << "Environment: worldcomm not allocated\n";
    return S_worldcomm->subWorlds(n);
}

std::vector<WorldComm> const&
Environment::worldsCommSeq( int n )
{
    CHECK( S_worldcommSeq ) << "Environment: worldcomm not allocated\n";
    return S_worldcommSeq->subWorlds(n);
}

std::vector<WorldComm> const&
Environment::worldsCommGroupBySubspace( int n )
{
#if 0
    std::cout << "n=" << n << "\n";
    S_worldcomm->showMe();
    S_worldcomm->masterWorld(n).showMe();
    std::cout << "size=" << S_worldcomm->subWorlds(n).size() <<  "\n";
    S_worldcomm->subWorlds(n).begin()->showMe();
#endif
    return S_worldcomm->subWorldsGroupBySubspace(n);
}


WorldComm const&
Environment::masterWorldComm( int n )
{
    return S_worldcomm->masterWorld(n);
}


MemoryUsage
Environment::logMemoryUsage( std::string const& message )
{
    MemoryUsage mem;
#if defined ( FEELPP_HAS_PETSC_H )
    PetscMemoryGetCurrentUsage( &mem.memory_usage );
    LOG(INFO) << message << " PETSC get current memory usage (resident memory): "  << mem.memory_usage/1e3 << "  KBytes "  << mem.memory_usage/1e6 << "  MBytes " << mem.memory_usage/1e9 << " GBytes" ;
    //PetscMemoryGetMaximumUsage( &mem );
    //LOG(INFO) << logMessage << " PETSC get maximum memory usag (resident memory): " << mem/1e6 << "  MBytes " << mem/1e9 << " GBytes" ;

    PetscMallocGetCurrentUsage( &mem.petsc_malloc_usage );
    LOG(INFO) << message << " PETSC get current PETSC Malloc usage: "  << mem.petsc_malloc_usage/1e3 << "  KBytes " << mem.petsc_malloc_usage/1e6 << " MBytes " << mem.petsc_malloc_usage/1e9 << " GBytes" ;
    PetscMallocGetMaximumUsage( &mem.petsc_malloc_maximum_usage );
    LOG(INFO) << message << " PETSC get maximum PETSC Malloc usage(largest memory ever used so far): "  << mem.petsc_malloc_maximum_usage/1e3 << "  KBytes " << mem.petsc_malloc_maximum_usage/1e6 << " MBytes " << mem.petsc_malloc_maximum_usage/1e9 << " GBytes" ;
#endif
    return mem;
}

AboutData Environment::S_about;
po::variables_map Environment::S_vm;
boost::shared_ptr<po::options_description> Environment::S_desc;
std::vector<std::string> Environment::S_to_pass_further;

boost::signals2::signal<void()> Environment::S_deleteObservers;

boost::shared_ptr<WorldComm> Environment::S_worldcomm;
boost::shared_ptr<WorldComm> Environment::S_worldcommSeq;

std::vector<fs::path> Environment::S_paths = { fs::current_path(),
                                               Environment::systemConfigRepository().get<0>(),
                                               Environment::systemGeoRepository().get<0>() };
fs::path Environment::S_scratchdir;

} // detail

}

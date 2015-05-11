/* -*- mode: c++; coding: utf-8 -*- */
namespace Feel {
/*! \page Environment Environment


\tableofcontents

\li \b Previous: \ref Keywords
\li \b Next: \ref Mesh

<hr>
In this section, we present some tools to initialize and manipulate \feel environment. For more information, see \ref detail::Environment "Environment".

\section Initialize Initialize Feel++
Environment class is necessary to initialize your application, as seen in \ref FirstApp. Interface is as follows:
\co
Environment env( _argc, _argv, _desc, _about );
\eco
None of those parameters are required but it is highly recommended to use the minimal declaration:
\co
  Environment env( _argc=argc, _argv=argv,
                   _desc=feel_option(),
                   _about=about(_name="name_of_your_app",
                                _author="your_name",
                                _email="your_email_adress") );
\eco

\li \c _argc and \c _argv are the arguments of your main function.
\li \c _desc is a description of your options.
\li \c _about is a brief description of your application.

\section Options Feel++ Options Description
\subsection Options_Custom Add Options
\c feel_option() returns a list of default options used in \feel.<br>

You can create a personal list of options as seen in \ref FirstApp.

You can also add a list of options, using a routine as follows:
\co
  using namespace Feel;


  inline
  po::options_description
  makeOptions()
  {
    po::options_description myappOptions( "My app options" );
    myappOptions.add_options()
      ( "option1", po::value<type1>()->default_value( value1 ), "description1" )
      ( "option2", po::value<type2>()->default_value( value2 ), "description2" )
      ( "option3", po::value<type3>()->default_value( value3 ), "description3" )
      ;
    return myappOptions.add( feel_options() ); // Add the default feel options to your list
  }
\eco
\li \c makeOptions is the usual name of this routine but you can change it
\li \c myappOptions: the name of you options list
\li \c option#: the name of parameter #
\li \c type#: the type parameter #
\li \c value#: the default value of parameter #
\li \c description#: the description of parameter #

The data returned is typically used as an argument of a Feel::Application subclass.

This routine has to be declared before your \c main function. Then you can use it to initialize \feel Environment:
\co
  Environment env( _argc=argc, _argv=argv,
                   _desc=makeOptions(),
                   _about=about(_name="myapp",
                                _author="myname",
                                _email="my@email.com") );
\eco

So you can change this parameter when you execute your app:
\verbatim
  ./myapp --option1=alpha --option2=beta --option3=gama
\endverbatim


<b>Example:</b><br>
From \c "/doc/manual/laplacian/laplacian.cpp":
\co
  using namespace Feel;

  inline
  po::options_description
  makeOptions()
  {
      po::options_description laplacianoptions( "Laplacian options" );
      laplacianoptions.add_options()
          ( "hsize", po::value<double>()->default_value( 0.2 ), "mesh size" )
          ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)" )
          ( "nu", po::value<double>()->default_value( 1 ), "grad.grad coefficient" )
          ( "weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
          ( "penaldir", Feel::po::value<double>()->default_value( 10 ),
          "penalisation parameter for the weak boundary Dirichlet formulation" )
          ( "exact1D", po::value<std::string>()->default_value( "sin(2*Pi*x)" ), "exact 1D solution" )
          ( "exact2D", po::value<std::string>()->default_value( "sin(2*Pi*x)*cos(2*Pi*y)" ), "exact 2D solution" )
          ( "exact3D", po::value<std::string>()->default_value( "sin(2*Pi*x)*cos(2*Pi*y)*cos(2*Pi*z)" ), "exact 3D solution" )
          ( "rhs1D", po::value<std::string>()->default_value( "" ), "right hand side 1D" )
          ( "rhs2D", po::value<std::string>()->default_value( "" ), "right hand side 2D" )
          ( "rhs3D", po::value<std::string>()->default_value( "" ), "right hand side 3D" )
          ;
      return laplacianoptions.add( Feel::feel_options() );
  }
\eco


\subsection Options_Accessors Accessors
<b>Options Description:</b><br>
\co Environment::optionsDescription();\eco
Returns options description data structure (\c po::options_description).<br>


<b>Variable map</b><br>
You can access to the parameters of your application environment using the following function:
\co
  Environment::vm(_name);
\eco
\c _name is the name of the parameter as seen in the previous paragraph.<br>
This function returns a \c po::variable_value .<br>
Use template methode to cast the parameter into the appropriate type.<br>

<b>Examples:</b><br>
From \c "doc/manual/solid/beam.cpp":
\co
  const double E = Environment::vm(_name="E").template as<double>();
  const double nu = Environment::vm(_name="nu").template as<double>();
\eco
From \c "doc/manual/fd/penalisation.cpp":
\co
  Tfinal =  Environment::vm( _name="test" ).template as<int>()*dt;
\eco


<a href="#" class="top">top</a>
<hr>
\section Repository Repository
\subsection Change changeRepository
You can change the default repository.

\Interface
\co
void changeRepository( _directory, _subdir, _filename );
\eco
Required Parameters:
\li <tt>_directory</tt>: new directory

Optional Parameters:
\li <tt>_subdir</tt>: Default = <tt>true</tt>
\li <tt>_filename</tt>: Default = <tt>"logfile"</tt>

You can use \c boost::format to customize the path. <br>
<b>Example:</b><br>
From \c "doc/manual/laplacian/laplacian.cpp":
\co
    Environment::changeRepository( boost::format( "doc/manual/laplacian/%1%/%2%-%3%/P%4%/h_%5%/" )
                                   % this->about().appName()
                                   % shape
                                   % Dim
                                   % Order
                                   % meshSize );
\eco
Then results will be store in: "/doc/manual/laplacian/<appName>/<shape>-<Dim>/P<Order>/h_<meshSize>/"


\subsection Find findFile
\Interface
\co
std::string findFile( std::string const& filename );\eco
Returns the string containing the filename path.

The lookup is as follows:
\li look into current path
\li look into paths that went through changeRepository(), it means that we look for example into the path from which the executable was run

If the file has an extension .geo or .msh, try also to
\li look into \c localGeoRepository() which is usually $HOME/feel/geo
\li look into \c systemGeoRepository() which is usually $FEELPP_DIR/share/feel/geo

If <tt>filename</tt> is not found, then the empty string is returned.


\subsection SetLogs setLogs
\Interface
\co
void setLogs( std::string const& prefix );
\eco
Required Parameters:
\li \c prefix: prefix for log filenames.


<hr>

\li \b Next: \ref Mesh

*/
}

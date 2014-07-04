
#include <feel/feel.hpp>
//#include "environment.hpp"

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <mpi4py/mpi4py.h>

#include <boost/parameter/keyword.hpp>
#include <boost/parameter/preprocessor.hpp>
#include <boost/parameter/python.hpp>
#include <boost/python/ptr.hpp>

using namespace boost::python;
namespace py = boost::parameter::python;
namespace mpl = boost::mpl;

class Argv {
    public :
        char** argv;
        int argc;

        Argv()
        {
            argc = 0;
            argv = NULL;    
        }

        Argv(boost::python::list arg)
        {
            
            argc=boost::python::len(arg);
            argv =new char* [argc+1];
            boost::python::stl_input_iterator<std::string> begin(arg), end;
            int i=0;
            while (begin != end)
            {
                argv[i] =strdup((*begin).c_str());
                begin++;
                i++;
            }
            argv[argc]=NULL;
            std::cout<<argc<<std::endl;
            for(int i=0;i<argc;i++)
                std::cout<<argv[i]<<std::endl;
              
        }

        std::string test()
        {
            return "test";
        }
       
       ~Argv()
       {
          if(argv != NULL)
          {
              for(int i=0;i<argc;i++)
                  delete argv[i];
              delete argv;
          }
       }
       
    

        char** getArgv () 
        {
            return argv;
        }

        static Argv* create ()
        {
            return new Argv();
        }
                   

               
};

/*
static PyObject* create(PyObject* object)
        {
            PyObject* instance= PyObject_CallObject(object,NULL);
            Argv* arg=extract<Argv*>(instance);
            return instance;
        }
        */



/*
std::auto_ptr<Feel::detail::Environment> class_arg(boost::python::list arg)
{
    int argc=boost::python::len(arg);
    char** argv =new char* [argc+1];
    boost::python::stl_input_iterator<std::string> begin(arg), end;
    int i=0;
    while (begin != end)
    { 
        argv[i] =strdup((*begin).c_str());
        begin++;
        i++;
    }
    argv[argc]=NULL;

    std::auto_ptr<Feel::detail::Environment> Env(new Feel::detail::Environment(argc,argv));

    for(int i=0;i<argc;i++)
        delete argv[i];
    delete argv;    
    
    return Env;
}
*/

// build the module, named libEnvironment, from previous methods
BOOST_PYTHON_MODULE(libEnvironment)
{
   
   if (import_mpi4py()<0) return ;
    /*
    class_<Feel::detail::Environment,boost::noncopyable> ("Environment", init<int&, char**&>())
      .def("__init__", make_constructor(class_arg)) 
        ;
    */
       
    class_<Argv>("Argv",init<boost::python::list>())
       .def("create",&Argv::create,return_value_policy<manage_new_object>())
       .def("test",&Argv::test);
        
        
        //.def("get",&Argv::getArgv,return_value_policy<manage_new_object>());
        //.def_readonly("argv",&Argv::argv);
       // .def_readonly("argc",&Argv::argc);
    /*add_property("argv",
    make_getter(&Argv::argv,return_value_policy<manage_new_object>()),
    make_setter(&Argv::argv,return_value_policy<manage_new_object>()));
*/



    class_<Feel::detail::Environment,boost::noncopyable>("Environment", init<boost::python::list>());


   /* 
       class_<Feel::Environment,boost::noncopyable, bases<Feel::detail::Environment> > ("Environment",no_init)
       .def(py::init<mpl::vector<Feel::tag::argc(int)
       ,Feel::tag::argv(char**)
       ,Feel::tag::desc*(Feel::po::options_description const&)
       ,Feel::tag::desc_lib*(Feel::po::options_description const&)
       ,Feel::tag::about*(Feel::AboutData const&)
       ,Feel::tag::directory*(std::string)
       >
       >()
       );  
    */
      
     // class_<Feel::Environment,boost::noncopyable, bases<Feel::detail::Environment> > ("Environment",init<boost::python::list>());



}

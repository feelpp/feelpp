//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 14 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
#include <vector>

#include <Eigen/Core>
#include <boost/shared_ptr.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
//#include <pybind11/eigen.h>
//#include <pybind11/numpy.h>


#include <feel/feelcore/feel.hpp>
#include <feel/feelcrb/crbenums.hpp>
#include <feel/feelcrb/crbdata.hpp>
#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelcrb/crbplugin_interface.hpp>
#include <feel/feelcrb/options.hpp>

namespace py = pybind11;
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);
using namespace Feel;

namespace Eigen {
namespace internal {
template<>
struct traits<ParameterSpaceX::Element>: traits<Eigen::VectorXd>
{
};
}
}

using VectorXr=Eigen::Matrix<Real,Eigen::Dynamic,1>;
    
template <class T>
void vector_setitem(std::vector<T>& v, int index, T value)
{
    if (index >= 0 && index < v.size()) {
        v[index] = value;
    }
    else {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        //throw_error_already_set();
    }
}

template <class T>
T vector_getitem(std::vector<T> &v, int index)
{
    if (index >= 0 && index < v.size()) {
        return v[index];
    }
    else {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        //throw_error_already_set();
    }
}
void IndexError() { PyErr_SetString(PyExc_IndexError, "Index out of range"); }
template<class T>
struct std_item
{
    typedef typename T::value_type V;
    static V const& get(T const& x, int i)
        {
            if( i<0 ) i+=x.size();
            if( i>=0 && i<x.size() ) return x[i];
            IndexError();
        }
    static void set(T & x, int i, V const& v)
        {
            if( i<0 ) i+=x.size();
            if( i>=0 && i<x.size() ) x[i]=v;
            else IndexError();
        }
    static void del(T& x, int i)
        {
            if( i<0 ) i+=x.size();
            if( i>=0 && i<x.size() ) x.erase(i);
            else IndexError();
        }
    static void add(T& x, V const& v)
        {
            x.push_back(v);
        }
    static size_t size(T const& x)
        {
            return x.size();
        }
    static std::vector<V> const& getVector(T const& v )
        {
            return v;
        }
    
};
template<class T>
struct str_item
{
    static std::string to_string( T const& t )
        {
            std::ostringstream os;
            os << t;
            return os.str();
        }
};


po::options_description
makeCRBOptions()
{
    po::options_description opts( "crb online run lib options" );
    opts.add(crbOptions())
        .add(crbSEROptions())
        .add(eimOptions())
        .add(podOptions())
        .add(backend_options("backend-primal"))
        .add(backend_options("backend-dual"))
        .add(backend_options("backend-l2"))
        ;
    return opts;
}

class CRBPluginAPIWrap : public CRBPluginAPI
{
public:
    std::string const& name() const override
        {
            PYBIND11_OVERLOAD_PURE(std::string const&, CRBPluginAPI, name, );
        }
    
    void loadDB( std::string const& db, crb::load l ) override
        {
            PYBIND11_OVERLOAD_PURE( void, CRBPluginAPI, loadDB, db, l );
        }
    void loadDBFromId( std::string const& i, crb::load l, std::string const& root  ) override
        {
            PYBIND11_OVERLOAD_PURE( void, CRBPluginAPI, loadDBFromId, i, l, root );
        }
    std::shared_ptr<ParameterSpaceX> parameterSpace() const override
        {
            PYBIND11_OVERLOAD_PURE(std::shared_ptr<ParameterSpaceX>, CRBPluginAPI, parameterSpace,  );
        }
    CRBResults run( ParameterSpaceX::Element const& mu, 
                    double eps , int N, bool print_rb_matrix ) const override
        {
            PYBIND11_OVERLOAD(CRBResults, CRBPluginAPI, run, mu, eps, N, print_rb_matrix ); 
        }
    std::vector<CRBResults> run( std::vector<ParameterSpaceX::Element> const& S, 
                                 double eps , int N, bool print_rb_matrix ) const override
        {
            PYBIND11_OVERLOAD(std::vector<CRBResults>, CRBPluginAPI, run, S, eps, N, print_rb_matrix ); 
        }
    void initExporter() override
        {
            PYBIND11_OVERLOAD_PURE(void, CRBPluginAPI, initExporter, );
        }
    void exportField( std::string const & name, CRBResults const& results ) override
        {
            PYBIND11_OVERLOAD_PURE(void, CRBPluginAPI, exportField, name, results );
            
        }
    void saveExporter() const override { PYBIND11_OVERLOAD_PURE(void, CRBPluginAPI, saveExporter,  ); }
protected:
    void setName( std::string const& n ) override { PYBIND11_OVERLOAD_PURE(void, CRBPluginAPI, setName, n ); }
};



namespace py = pybind11;
PYBIND11_MODULE( _mor, m )
{
    m.def("makeCRBOptions", &makeCRBOptions, "Create CRB Options" );
    m.def("factoryCRBPlugin", &factoryCRBPlugin, "Factory for CRB plugins",
          py::arg("name"),py::arg("libname")=std::string(""),py::arg("dirname")=Info::libdir() );

    py::enum_<crb::stage>(m, "CRBStage")
        .value("offline", crb::stage::offline)
        .value("online", crb::stage::online)
        ;

    py::enum_<crb::load>(m,"CRBLoad")
        .value("rb", crb::load::rb)
        .value("fe", crb::load::fe)
        .value("all", crb::load::all);
    
    py::class_<CRBResults>(m,"CRBResults")
        .def("setParameter", &CRBResults::setParameter)
        .def("parameter", &CRBResults::parameter)
        .def("output", &CRBResults::output)
        .def("errorBound", &CRBResults::errorbound)
        ;

    py::class_<std::vector<CRBResults>>(m,"VectorCRBResults")
        .def(py::init<>())
        .def("clear", &std::vector<CRBResults>::clear)
        .def("pop_back", &std::vector<CRBResults>::pop_back)
        .def("__len__", [](const std::vector<CRBResults> &v) { return v.size(); })
        .def("__iter__", [](std::vector<CRBResults> &v) {
                return py::make_iterator(v.begin(), v.end());
            }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */

    using ElementP = ParameterSpaceX::Element;
    py::class_<ElementP>(m,"ParameterSpaceElement")
        //.def("__repr__",[](const ParameterSpaceX::Element &a) {
        //return "<ParameterSpaceX named '" + a + "'>";
        //});
        .def("__str__", [](const ParameterSpaceX::Element &a) {
                            std::ostringstream os;
                            os << "[";
                            os.precision(2);
                            os.setf(std::ios::scientific);
                            for ( int i = 0; i < a.size(); ++i )
                            {
                                os << a[i];
                                if ( i != a.size()-1)
                                    os << ",";
                            }
                            os << "]";
                            return os.str();
                        })
        .def("view", &ElementP::view, "view the parameters names, with its values")
        .def("parameterNamed", static_cast<double& (ElementP::*)(std::string)>(&ElementP::parameterNamed), "return the parameter named", py::arg("name") )
        .def("parameterName", &ElementP::parameterName, "return the i-th name ", py::arg("i"))
        .def("size", &ElementP::size, "return the size of the parameters" )
        .def("__call__", static_cast<double& (ElementP::*)(int)>(&ElementP::coeff), "return the ith parameter", py::arg("i") )
        .def("setParameter", static_cast<void (ElementP::*)(int, double)>(&ElementP::setParameter), "set the i-th to the value", py::arg("i"), py::arg("value"))
        .def("setParameterNamed", static_cast<void (ElementP::*)(std::string, double)>(&ElementP::setParameterNamed), "set the named parameter to the value", py::arg("name"), py::arg("value"))
        .def("setParameters", static_cast<void (ElementP::*)(std::vector<double>)>(&ElementP::setParameters), "set the parameter to the given list", py::arg("values"))
        .def("setParameters", static_cast<void (ElementP::*)(std::map<std::string, double> values)>(&ElementP::setParameters), "set the parameter to the given dict with values", py::arg("values"))
        ;

    //!
    //! Sampling wrapping
    //!
    py::class_<std::vector<ParameterSpaceX::Element>>(m,"VectorParameterSpaceElement");

    py::class_<ParameterSpaceX::Sampling, std::shared_ptr<ParameterSpaceX::Sampling>>(m,"ParameterSpaceSampling")
        .def(py::init<std::shared_ptr<ParameterSpaceX>,int,std::shared_ptr<ParameterSpaceX::Sampling>>())
        .def("sampling",&ParameterSpaceX::Sampling::sampling)
        .def("__getitem__", &std_item<ParameterSpaceX::Sampling>::get,py::return_value_policy::reference)
        .def("__setitem__", &std_item<ParameterSpaceX::Sampling>::set)
        .def("getVector", &std_item<ParameterSpaceX::Sampling>::getVector)
        .def("__len__", &std_item<ParameterSpaceX::Sampling>::size)
        ;

    //!
    //! Parameter Space
    //!
    py::class_<ParameterSpaceX,std::shared_ptr<ParameterSpaceX>>(m,"ParameterSpace")
        .def( py::init<>() )
        .def("sampling", &ParameterSpaceX::sampling)
        .def("element", &ParameterSpaceX::element, "return a parameter from the space\n  - broadcast : share the parameter to all processors\n  - apply_log : log random chosen parameter",
            py::arg("broadcast")=true, py::arg("apply_log")=false)
        //.def("New", &ParameterSpaceX::New, ParameterSpaceX_New_overloads(args("dim", "WorldComm"), "New")).staticmethod("New")
        .def_static("create",&ParameterSpaceX::create )
        ;

    py::class_<CRBPluginAPI,CRBPluginAPIWrap,std::shared_ptr<CRBPluginAPI>>(m,"CRBPlugin")
        //.def(py::init<>())
        .def("name",&CRBPluginAPI::name,py::return_value_policy::reference)
        //.def("setName",pure_virtual(&CRBPluginAPI::setName))
        .def("loadDB",&CRBPluginAPI::loadDB,"load a database from filename",py::arg("filename"),py::arg("load")=crb::load::rb )
        .def("loadDBFromId",&CRBPluginAPI::loadDBFromId, "load a database from its id", py::arg(
                 "id"), py::arg("load")=crb::load::rb, py::arg("root")=Environment::rootRepository())
        .def("isReducedBasisModelDBLoaded",&CRBPluginAPI::isReducedBasisModelDBLoaded, "returns true if Reduced Basis Model DB is loaded, false otherwise")
        .def("isFiniteElementModelDBLoaded",&CRBPluginAPI::isFiniteElementModelDBLoaded, "returns true if Finite Element Model DB is loaded, false otherwise")
        .def("isAllLoaded",&CRBPluginAPI::isAllLoaded, "returns true if all DB is loaded, false otherwise")
        .def("isDBLoaded",&CRBPluginAPI::isDBLoaded, "returns true if some data from the model DB is loaded, false otherwise")
        // finite element data
        .def("meshes",&CRBPluginAPI::meshes)
        .def("doftables",&CRBPluginAPI::doftables)
        .def("feElement",&CRBPluginAPI::feElement)
        .def("feSubElements",&CRBPluginAPI::feSubElements)
        .def("primalBasis",&CRBPluginAPI::reducedBasisFunctionsPrimal)
        .def("dualBasis",&CRBPluginAPI::reducedBasisFunctionsDual)
        // exporter
        .def("initExporter",&CRBPluginAPI::initExporter)
        .def("saveExporter",&CRBPluginAPI::saveExporter)
        .def("exportField",&CRBPluginAPI::exportField)
        
        .def("parameterSpace",&CRBPluginAPI::parameterSpace)
        .def("run",py::overload_cast<ParameterSpaceX::Element const&,double,int,bool>(&CRBPluginAPI::run, py::const_),
             "run online code for a parameter mu", py::arg("mu"),py::arg("eps")=1e-6,py::arg("N")=-1,py::arg("print_reduced_matrix")=false)
        .def("run",py::overload_cast<std::vector<ParameterSpaceX::Element> const&,double,int,bool>(&CRBPluginAPI::run, py::const_),
             "run online code for a parameter mu", py::arg("mu"),py::arg("eps")=1e-6,py::arg("N")=-1,py::arg("print_reduced_matrix")=false);
        //.def("run",&CRBPluginAPI::run)

}

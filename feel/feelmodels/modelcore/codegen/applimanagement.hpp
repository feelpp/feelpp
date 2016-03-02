/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

#ifndef __APPLIMANAGEMENT_H
#define __APPLIMANAGEMENT_H 1

#include <feel/options.hpp>

#include <feel/feelmodels/modelcore/feelmodelscoreconstconfig.h>
#include <feel/feelmodels/modelcore/options.hpp>

//#include <fsi/fsidiscr/toolboxmesh.hpp>



#include "feelmodelscoreconfig.h"

//-----------------------------//
// fluid
//-----------------------------//

#if defined( FLUIDMECHANICS )
#undef APPLICATION_TYPE_NAME
#define APPLICATION_TYPE_NAME "fluid"
#if defined( FEELPP_MODELS_ENABLE_FLUIDMECHANICS )
#include "fluid/codegen_fluidmec.hpp"
#define APPLI_FLUID_UORDER FLUIDMECHANICS_ORDER_VELOCITY
#define APPLI_FLUID_PORDER FLUIDMECHANICS_ORDER_PRESSURE
#define APPLI_FLUID_GORDER FLUIDMECHANICS_ORDERGEO
#endif
#endif

//-----------------------------//
// struct
//-----------------------------//

#if defined( SOLIDMECHANICS )
#undef APPLICATION_TYPE_NAME
#define APPLICATION_TYPE_NAME "solid"
#include "solid/codegen_solidmec.hpp"
#define APPLI_SOLID_UORDER SOLIDMECHANICS_ORDER_DISPLACEMENT
#define APPLI_SOLID_GORDER SOLIDMECHANICS_ORDERGEO
#endif

//-----------------------------//
// fsi
//-----------------------------//

#if defined( FLUIDMECHANICS ) && defined( SOLIDMECHANICS )
#undef APPLICATION_TYPE_NAME
#define APPLICATION_TYPE_NAME "fsi"
#define USEFULL_FULLY_FSI 1
#include <fsi/fsidiscr/toolboxfsi.hpp>
#endif


//-----------------------------//
// thermo-dynamics model
//-----------------------------//
#undef THERMODYNAMICS
#include "feelmodelscoreconfig.h"
#if defined( THERMODYNAMICS )
#undef APPLICATION_TYPE_NAME
#define APPLICATION_TYPE_NAME "thermo-dynamics"
#include "thermodyn/codegen_thermodyn.hpp"
#endif




//-----------------------------//
// multi-fluid
//-----------------------------//
#if defined( FEELPP_MODELS_ENABLE_FLUIDMECHANICS )
// fluid0
#undef FLUIDMECHANICS
#undef FLUIDMECHANICS0
#undef FLUIDMECHANICS1
#undef FLUIDMECHANICS2
#include "feelmodelscoreconfig.h"
#if defined( FLUIDMECHANICS0 )
#undef __FLUIDMECHANICS_H
#include "fluidmec0/fluidmec.hpp"
#endif
// fluid1
#undef FLUIDMECHANICS
#undef FLUIDMECHANICS0
#undef FLUIDMECHANICS1
#undef FLUIDMECHANICS2
#include "feelmodelscoreconfig.h"
#if defined( FLUIDMECHANICS1 )
#undef __FLUIDMECHANICS_H
#undef __STOKESVARIANTS_H
#include "fluidmec1/fluidmec.hpp"
#endif
// fluid2
#undef FLUIDMECHANICS
#undef FLUIDMECHANICS0
#undef FLUIDMECHANICS1
#undef FLUIDMECHANICS2
#include "feelmodelscoreconfig.h"
#if defined( FLUIDMECHANICS2 )
#undef __FLUIDMECHANICS_H
#undef __STOKESVARIANTS_H
#include "fluidmec2/fluidmec.hpp"
#endif

#endif // if defined( FEELPP_MODELS_ENABLE_FLUIDMECHANICS )

//-----------------------------//
// multi-solid
//-----------------------------//
#undef SOLIDMECHANICS
#undef SOLIDMECHANICS0
#undef SOLIDMECHANICS1
#undef SOLIDMECHANICS2
#include "feelmodelscoreconfig.h"
#if defined( SOLIDMECHANICS0 )
#undef __SOLIDMECHANICS_H
#include "solidmec0/solidmec.hpp"
#endif
//-----------------------------//
#undef SOLIDMECHANICS
#undef SOLIDMECHANICS0
#undef SOLIDMECHANICS1
#undef SOLIDMECHANICS2
#include "feelmodelscoreconfig.h"
#if defined( SOLIDMECHANICS1 )
#undef __SOLIDMECHANICS_H
#include "solidmec1/solidmec.hpp"
#endif
//-----------------------------//
#undef SOLIDMECHANICS
#undef SOLIDMECHANICS0
#undef SOLIDMECHANICS1
#undef SOLIDMECHANICS2
#include "feelmodelscoreconfig.h"
#if defined( SOLIDMECHANICS2 )
#undef __SOLIDMECHANICS_H
#include "solidmec2/solidmec.hpp"
#endif
//-----------------------------//
// multi-thermo
//-----------------------------//
// thermo0
#undef THERMODYNAMICS
#undef THERMODYNAMICS0
#undef THERMODYNAMICS1
#undef THERMODYNAMICS2
#include "feelmodelscoreconfig.h"
#if defined( THERMODYNAMICS0 )
#undef __THERMODYNAMICS_H
#include "thermodyn0/thermodyn.hpp"
namespace Feel { namespace detail {
uint16_type getPropertyThermo0PolyOrder() { return THERMODYNAMICS_ORDERPOLY; }
uint16_type getPropertyThermo0GeoOrder() { return THERMODYNAMICS_ORDERGEO; }
} }
#endif
// thermo1
#undef THERMODYNAMICS
#undef THERMODYNAMICS0
#undef THERMODYNAMICS1
#undef THERMODYNAMICS2
#include "feelmodelscoreconfig.h"
#if defined( THERMODYNAMICS1 )
#undef __THERMODYNAMICS_H
#include "thermodyn1/thermodyn.hpp"
namespace Feel { namespace detail {
uint16_type getPropertyThermo1PolyOrder() { return THERMODYNAMICS_ORDERPOLY; }
uint16_type getPropertyThermo1GeoOrder() { return THERMODYNAMICS_ORDERGEO; }
} }
#endif
// thermo2
#undef THERMODYNAMICS
#undef THERMODYNAMICS0
#undef THERMODYNAMICS1
#undef THERMODYNAMICS2
#include "feelmodelscoreconfig.h"
#if defined( THERMODYNAMICS2 )
#undef __THERMODYNAMICS_H
#include "thermodyn2/thermodyn.hpp"
namespace Feel { namespace detail {
uint16_type getPropertyThermo2PolyOrder() { return THERMODYNAMICS_ORDERPOLY; }
uint16_type getPropertyThermo2GeoOrder() { return THERMODYNAMICS_ORDERGEO; }
} }
#endif

//-----------------------------//
// reload feelmodelscoreconfig.h
//-----------------------------//
#undef FLUIDMECHANICS
#undef FLUIDMECHANICS0
#undef FLUIDMECHANICS1
#undef FLUIDMECHANICS2
#undef SOLIDMECHANICS
#undef SOLIDMECHANICS0
#undef SOLIDMECHANICS1
#undef SOLIDMECHANICS2
#undef THERMODYNAMICS
#undef THERMODYNAMICS0
#undef THERMODYNAMICS1
#undef THERMODYNAMICS2
#include "feelmodelscoreconfig.h"

//-----------------------------//
// fbm fluid
//-----------------------------//
#if defined( FLUIDMECHANICS0 ) && defined( FLUIDMECHANICS1 ) && defined( FLUIDMECHANICS2 ) && defined( SOLIDMECHANICS0 )
#undef APPLICATION_TYPE_NAME
#define APPLICATION_TYPE_NAME "three_fluid_one_struct"
#elif defined( FLUIDMECHANICS0 ) && defined( FLUIDMECHANICS1 ) && defined( FLUIDMECHANICS2 )
#undef APPLICATION_TYPE_NAME
#define APPLICATION_TYPE_NAME "three_fluid"
#elif defined( FLUIDMECHANICS0 ) && defined( FLUIDMECHANICS1 )
#undef APPLICATION_TYPE_NAME
#define APPLICATION_TYPE_NAME "two_fluid"
#endif

#if defined( FLUIDMECHANICS0 ) && defined( FLUIDMECHANICS1 )
#if defined( FBMTOOL )
#define FBM_DONT_INCLUDE_FLUID 1
#include "FBM/fbm.hpp"
#undef FBM_DONT_INCLUDE_FLUID
#endif
#endif


//-----------------------------//
// fbm thermo
//-----------------------------//
#if defined( THERMODYNAMICS0 ) && defined( THERMODYNAMICS1 ) && defined( THERMODYNAMICS2 )
#undef APPLICATION_TYPE_NAME
#define APPLICATION_TYPE_NAME "three_thermo"
#elif defined( THERMODYNAMICS0 ) && defined( THERMODYNAMICS1 )
#undef APPLICATION_TYPE_NAME
#define APPLICATION_TYPE_NAME "two_thermo"
#endif

#if defined( THERMODYNAMICS0 ) && defined( THERMODYNAMICS1 )
#if defined( FBMTOOL )
#define FBM_DONT_INCLUDE_THERMO 1
#include "FBM/fbmthermodyn.hpp"
#undef FBM_DONT_INCLUDE_THERMO
#endif
#endif


namespace Feel
{
    namespace parameter = boost::parameter;

    BOOST_PARAMETER_NAME(buildMesh)
    BOOST_PARAMETER_NAME(subPrefix)


    namespace detail
    {

        // copy/paste from http://www.boost.org/doc/libs/1_55_0/libs/format/doc/format.html
        // specific wrapper function for building format objects with the right exceptions settings
        boost::format
        myBoostFormat(const std::string & f_string)
        {
            using namespace boost::io;
            boost::format fmter(f_string);
            fmter.exceptions( all_error_bits ^ ( too_many_args_bit | too_few_args_bit )  );
            return fmter;
        }
    } // namespace detail

} // namespace Feel

namespace Feel
{

    class AppliManagement
    {
    public :

        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//
#if defined( FEELPP_MODELS_ENABLE_FLUIDMECHANICS )
#if defined( FLUIDMECHANICS )
        typedef Feel::FeelModels::FluidMechanics FluidMechanics_type;
        typedef boost::shared_ptr<FluidMechanics_type> FluidMechanics_ptrtype;
#endif
#if defined( FLUIDMECHANICS0 )
        typedef Feel::FeelModels::FluidMechanics0 FluidMechanics0_type;
        typedef boost::shared_ptr<FluidMechanics0_type> FluidMechanics0_ptrtype;
#endif
#if defined( FLUIDMECHANICS1 )
        typedef Feel::FeelModels::FluidMechanics1 FluidMechanics1_type;
        typedef boost::shared_ptr<FluidMechanics1_type> FluidMechanics1_ptrtype;
#endif
#if defined( FLUIDMECHANICS2 )
        typedef Feel::FeelModels::FluidMechanics2 FluidMechanics2_type;
        typedef boost::shared_ptr<FluidMechanics2_type> FluidMechanics2_ptrtype;
#endif
#endif
        //-------------------------------------------------------------------------------------//
#if defined( SOLIDMECHANICS )
        typedef Feel::FeelModels::SolidMechanics SolidMechanics_type;
        typedef boost::shared_ptr<SolidMechanics_type> SolidMechanics_ptrtype;
#endif
#if defined( SOLIDMECHANICS0 )
        typedef Feel::FeelModels::SolidMechanics0 SolidMechanics0_type;
        typedef boost::shared_ptr<SolidMechanics0_type> SolidMechanics0_ptrtype;
#endif
#if defined( SOLIDMECHANICS1 )
        typedef Feel::FeelModels::SolidMechanics1 SolidMechanics1_type;
        typedef boost::shared_ptr<SolidMechanics1_type> SolidMechanics1_ptrtype;
#endif
#if defined( SOLIDMECHANICS2 )
        typedef Feel::FeelModels::SolidMechanics2 SolidMechanics2_type;
        typedef boost::shared_ptr<SolidMechanics2_type> SolidMechanics2_ptrtype;
#endif
        //-------------------------------------------------------------------------------------//
#if defined( THERMODYNAMICS )
        typedef Feel::FeelModels::ThermoDynamics ThermoDynamics_type;
        typedef boost::shared_ptr<ThermoDynamics_type> ThermoDynamics_ptrtype;
#endif
#if defined( THERMODYNAMICS0 )
        typedef Feel::FeelModels::ThermoDynamics0 ThermoDynamics0_type;
        typedef boost::shared_ptr<ThermoDynamics0_type> ThermoDynamics0_ptrtype;
#endif
#if defined( THERMODYNAMICS1 )
        typedef Feel::FeelModels::ThermoDynamics1 ThermoDynamics1_type;
        typedef boost::shared_ptr<ThermoDynamics1_type> ThermoDynamics1_ptrtype;
#endif
#if defined( THERMODYNAMICS2 )
        typedef Feel::FeelModels::ThermoDynamics2 ThermoDynamics2_type;
        typedef boost::shared_ptr<ThermoDynamics2_type> ThermoDynamics2_ptrtype;
#elif defined( THERMODYNAMICS1 )
        typedef ThermoDynamics1_type ThermoDynamics2_type;
        typedef ThermoDynamics1_ptrtype ThermoDynamics2_ptrtype;
#endif
        //-------------------------------------------------------------------------------------//
#if defined( FBMTOOL )
#if defined( FLUIDMECHANICS0 ) && defined( FLUIDMECHANICS1 ) && defined( FEELPP_MODELS_ENABLE_FLUIDMECHANICS )
        typedef Feel::FeelModels::FBM FBM_type;
        typedef boost::shared_ptr<FBM_type> FBM_ptrtype;
#endif
#if defined( THERMODYNAMICS0 ) && defined( THERMODYNAMICS1 )
        typedef Feel::FeelModels::FBMThermoDynamics FBMThermoDynamics_type;
        typedef boost::shared_ptr<FBMThermoDynamics_type> FBMThermoDynamics_ptrtype;
#endif
#endif
        //-------------------------------------------------------------------------------------//

#if defined( FLUIDMECHANICS )
        static const uint16_type nDim = FluidMechanics_type::mesh_type::nDim;
        static const uint16_type nGeoOrder = FluidMechanics_type::mesh_type::nOrder;
        static const uint16_type nRealDim = FluidMechanics_type::mesh_type::nRealDim;
#elif defined( FLUIDMECHANICS0 )
        static const uint16_type nDim = FluidMechanics0_type::mesh_type::nDim;
        static const uint16_type nGeoOrder = FluidMechanics0_type::mesh_type::nOrder;
        static const uint16_type nRealDim = FluidMechanics0_type::mesh_type::nRealDim;
#elif defined( SOLIDMECHANICS )
        static const uint16_type nDim = SolidMechanics_type::mesh_type::nRealDim;
        static const uint16_type nGeoOrder = SolidMechanics_type::mesh_type::nOrder;
        static const uint16_type nRealDim = SolidMechanics_type::mesh_type::nRealDim;
#elif defined( THERMODYNAMICS )
        static const uint16_type nDim = ThermoDynamics_type::mesh_type::nRealDim;
        static const uint16_type nGeoOrder = ThermoDynamics_type::mesh_type::nOrder;
        static const uint16_type nRealDim = ThermoDynamics_type::mesh_type::nRealDim;
#elif defined( THERMODYNAMICS0 )
        static const uint16_type nDim = ThermoDynamics0_type::mesh_type::nRealDim;
        static const uint16_type nGeoOrder = ThermoDynamics0_type::mesh_type::nOrder;
        static const uint16_type nRealDim = ThermoDynamics0_type::mesh_type::nRealDim;
#else
        static const uint16_type nDim = 3;
        static const uint16_type nGeoOrder = 1;
        static const uint16_type nRealDim = 3;
#endif

        typedef node<double>::type node_type;

        //-------------------------------------------------------------------------------------//

        AppliManagement(int argc, char** argv ,const Feel::AboutData & ad)
            :
            M_environment( _argc=argc, _argv=argv,_desc=feelmodels_options(APPLICATION_TYPE_NAME) ),
            M_verbose( boption(_prefix="master",_name="verbose") )

        {
            this->init();
        }

        //-------------------------------------------------------------------------------------//

        AppliManagement(int argc, char** argv ,const Feel::AboutData & ad, const Feel::po::options_description & od)
            :
            M_environment( _argc=argc, _argv=argv,_desc=feelmodels_options(APPLICATION_TYPE_NAME).add(od) ),
            M_verbose( boption(_prefix="master",_name="verbose") )
        {
            this->init();
        }

        //-------------------------------------------------------------------------------------//

        bool verbose() { return M_verbose; }

        double initialTime() const { return M_Ti; }
        double stepTime() const { return M_dt; }
        double finalTime() const { return M_Tf; }
        bool isStationary() const { return M_isStationary; }
        bool doRestart() const { return option(_name="bdf.restart").as<bool>(); }

        const std::list<boost::tuple<std::string,node_type> > &
        ptsToEvaluate() const { return M_listPtEval;}

        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//
#if defined( FEELPP_MODELS_ENABLE_FLUIDMECHANICS )
#if defined( FLUIDMECHANICS )
        BOOST_PARAMETER_MEMBER_FUNCTION( (FluidMechanics_ptrtype), // return type
                                         FluidMechanics, //name function
                                         tag,
                                         (optional
                                          (worldcomm, (WorldComm), Environment::worldComm() )
                                          (buildMesh     , (bool)     , true)
                                          (subPrefix, (std::string), "" ) ) //optional
                                         )
        {
            FluidMechanics_ptrtype fm(new FluidMechanics_type( /*M_isStationary,*/"fluid",
                                                               buildMesh, worldcomm,  subPrefix,
                                                               this->appliShortRepository() ) );
            return fm;
        }
#endif
        //-------------------------------------------------------------------------------------//
#if defined( FLUIDMECHANICS0 )
        BOOST_PARAMETER_MEMBER_FUNCTION( (FluidMechanics0_ptrtype), // return type
                                         FluidMechanics0, //name function
                                         tag,
                                         (optional
                                          (worldcomm, (WorldComm), Environment::worldComm() )
                                          (buildMesh     , (bool)     , true)
                                          (subPrefix, (std::string), "" ) ) //optional
                                         )
        {
            FluidMechanics0_ptrtype fm(new FluidMechanics0_type( /*M_isStationary,*/"fluid0",
                                                                 buildMesh, worldcomm, subPrefix,
                                                                 this->appliShortRepository() ) );
            return fm;
        }
#endif
        //-------------------------------------------------------------------------------------//
#if defined( FLUIDMECHANICS1 )
        BOOST_PARAMETER_MEMBER_FUNCTION( (FluidMechanics1_ptrtype), // return type
                                         FluidMechanics1, //name function
                                         tag,
                                         (optional
                                          (worldcomm, (WorldComm), Environment::worldComm() )
                                          (buildMesh     , (bool)     , true)
                                          (subPrefix, (std::string), "" ) ) //optional
                                         )
        {
            FluidMechanics1_ptrtype fm(new FluidMechanics1_type( /*M_isStationary,*/"fluid1",
                                                                 buildMesh, worldcomm, subPrefix,
                                                                 this->appliShortRepository() ) );
            return fm;
        }
#endif
        //-------------------------------------------------------------------------------------//
#if defined( FLUIDMECHANICS2 )
        BOOST_PARAMETER_MEMBER_FUNCTION( (FluidMechanics2_ptrtype), // return type
                                         FluidMechanics2, //name function
                                         tag,
                                         (optional
                                          (worldcomm, (WorldComm), Environment::worldComm() )
                                          (buildMesh     , (bool)     , true)
                                          (subPrefix, (std::string), "" ) ) //optional
                                         )
        {
            FluidMechanics2_ptrtype fm(new FluidMechanics2_type( /*M_isStationary,*/"fluid2",
                                                                 buildMesh, worldcomm, subPrefix,
                                                                 this->appliShortRepository() ) );
            return fm;
        }
#endif
#endif
        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//
#if defined( SOLIDMECHANICS )
        BOOST_PARAMETER_MEMBER_FUNCTION( (SolidMechanics_ptrtype), // return type
                                         SolidMechanics, //name function
                                         tag,
                                         (optional
                                          (worldcomm, (WorldComm), Environment::worldComm() )
                                          (buildMesh     , (bool)     , true)
                                          (subPrefix, (std::string), "" ) ) //optional
                                         )
        {
            SolidMechanics_ptrtype sm(new SolidMechanics_type( /*M_isStationary,*/"struct",
                                                               buildMesh, worldcomm, subPrefix,
                                                               this->appliShortRepository() ) );
            return sm;
        }
#endif
        //-------------------------------------------------------------------------------------//
#if defined( SOLIDMECHANICS0 )
        BOOST_PARAMETER_MEMBER_FUNCTION( (SolidMechanics0_ptrtype), // return type
                                         SolidMechanics0, //name function
                                         tag,
                                         (optional
                                          (worldcomm, (WorldComm), Environment::worldComm() )
                                          (buildMesh     , (bool)     , true)
                                          (subPrefix, (std::string), "" ) ) //optional
                                         )
        {
            SolidMechanics0_ptrtype sm(new SolidMechanics0_type( /*M_isStationary,*/"struct0",
                                                                 buildMesh, worldcomm, subPrefix,
                                                                 this->appliShortRepository() ) );
            return sm;
        }
#endif
        //-------------------------------------------------------------------------------------//
#if defined( SOLIDMECHANICS1 )
        BOOST_PARAMETER_MEMBER_FUNCTION( (SolidMechanics1_ptrtype), // return type
                                         SolidMechanics1, //name function
                                         tag,
                                         (optional
                                          (worldcomm, (WorldComm), Environment::worldComm() )
                                          (buildMesh     , (bool)     , true)
                                          (subPrefix, (std::string), "" ) ) //optional
                                         )
        {
            SolidMechanics1_ptrtype sm(new SolidMechanics1_type( /*M_isStationary,*/"struct1",
                                                                 buildMesh, worldcomm, subPrefix,
                                                                 this->appliShortRepository() ) );
            return sm;
        }
#endif
        //-------------------------------------------------------------------------------------//
#if defined( SOLIDMECHANICS2 )
        BOOST_PARAMETER_MEMBER_FUNCTION( (SolidMechanics2_ptrtype), // return type
                                         SolidMechanics2, //name function
                                         tag,
                                         (optional
                                          (worldcomm, (WorldComm), Environment::worldComm() )
                                          (buildMesh     , (bool)     , true)
                                          (subPrefix, (std::string), "" ) ) //optional
                                         )
        {
            SolidMechanics2_ptrtype sm(new SolidMechanicsé_type( /*M_isStationary,*/"struct2",
                                                                 buildMesh, worldcomm, subPrefix,
                                                                 this->appliShortRepository() ) );
            return sm;
        }
#endif
        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//

#if defined( FLUIDMECHANICS ) && defined( SOLIDMECHANICS )
      boost::shared_ptr<ToolBoxFSI<FluidMechanics_type,SolidMechanics_type> >
      toolBoxFSI(FluidMechanics_ptrtype FM,SolidMechanics_ptrtype SM)
      {
        typedef ToolBoxFSI<FluidMechanics_type,SolidMechanics_type> toolboxfsi_type;

        boost::shared_ptr<toolboxfsi_type> tbfsi(new toolboxfsi_type(FM,SM));
        return tbfsi;
      }
#endif


        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//

#if defined( THERMODYNAMICS )
        BOOST_PARAMETER_MEMBER_FUNCTION( (ThermoDynamics_ptrtype), // return type
                                         ThermoDynamics, //name function
                                         tag,
                                         (optional
                                          (worldcomm, (WorldComm), Environment::worldComm() )
                                          (buildMesh     , (bool)     , true)
                                          (prefix     , (std::string)     , "thermo" )
                                          (subPrefix, (std::string), "" ) ) //optional
                                         )
        {
            ThermoDynamics_ptrtype thermoDyn(new ThermoDynamics_type( /*M_isStationary,*/prefix,
                                                                      buildMesh, worldcomm, subPrefix,
                                                                      this->appliShortRepository() ) );
            return thermoDyn;
        }
#endif
        //-------------------------------------------------------------------------------------//
#if defined( THERMODYNAMICS0 )
        BOOST_PARAMETER_MEMBER_FUNCTION( (ThermoDynamics0_ptrtype), // return type
                                         ThermoDynamics0, //name function
                                         tag,
                                         (optional
                                          (worldcomm, (WorldComm), Environment::worldComm() )
                                          (buildMesh     , (bool)     , true)
                                          (prefix     , (std::string)     , "thermo0" )
                                          (subPrefix, (std::string), "" ) ) //optional
                                         )
        {
            ThermoDynamics0_ptrtype thermoDyn(new ThermoDynamics0_type( /*M_isStationary,*/prefix,
                                                                        buildMesh, worldcomm, subPrefix,
                                                                        this->appliShortRepository() ) );
            return thermoDyn;
        }
#endif
        //-------------------------------------------------------------------------------------//
#if defined( THERMODYNAMICS1 )
        BOOST_PARAMETER_MEMBER_FUNCTION( (ThermoDynamics1_ptrtype), // return type
                                         ThermoDynamics1, //name function
                                         tag,
                                         (optional
                                          (worldcomm, (WorldComm), Environment::worldComm() )
                                          (buildMesh     , (bool)     , true)
                                          (prefix     , (std::string)     , "thermo1" )
                                          (subPrefix, (std::string), "" ) ) //optional
                                         )
        {
            ThermoDynamics1_ptrtype thermoDyn(new ThermoDynamics1_type( /*M_isStationary,*/prefix,
                                                                        buildMesh, worldcomm, subPrefix,
                                                                        this->appliShortRepository() ) );
            return thermoDyn;
        }
#endif
        //-------------------------------------------------------------------------------------//
#if defined( THERMODYNAMICS2 )
        BOOST_PARAMETER_MEMBER_FUNCTION( (ThermoDynamics2_ptrtype), // return type
                                         ThermoDynamics2, //name function
                                         tag,
                                         (optional
                                          (worldcomm, (WorldComm), Environment::worldComm() )
                                          (buildMesh     , (bool)     , true)
                                          (prefix     , (std::string)     , "thermo2" )
                                          (subPrefix, (std::string), "" ) ) //optional
                                         )
        {
            ThermoDynamics2_ptrtype thermoDyn(new ThermoDynamics2_type( /*M_isStationary,*/prefix,
                                                                        buildMesh, worldcomm, subPrefix,
                                                                        this->appliShortRepository() ) );
            return thermoDyn;
        }
#endif
        //-------------------------------------------------------------------------------------//



        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//


        std::string rootRepository() const
        {
            return Environment::rootRepository();
        }

        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//

        std::string appliShortRepository() const
        {
            CHECK(  Environment::vm().count("exporter.directory") ) << "exporter.directory must be given\n";
          std::string exportDirectoryOption=option(_name="exporter.directory").as<std::string>();
          std::string exportDirectoryStr = exportDirectoryOption;

          //-----------------//
          // fsi
#if defined(APPLI_FLUID_UORDER) && defined(APPLI_SOLID_UORDER)
          exportDirectoryStr = ( Feel::detail::myBoostFormat(exportDirectoryOption)
                                 %APPLI_FLUID_UORDER %APPLI_FLUID_PORDER %APPLI_FLUID_GORDER
                                 %APPLI_SOLID_UORDER %APPLI_SOLID_GORDER ).str();
          //-----------------//
          // fluid
#elif defined(APPLI_FLUID_UORDER)
          //std::string exportDirectoryStr = (boost::format(exportDirectoryOption) %APPLI_FLUID_UORDER %APPLI_FLUID_PORDER %APPLI_FLUID_GORDER).str();
          exportDirectoryStr = ( Feel::detail::myBoostFormat(exportDirectoryOption)
                                 %APPLI_FLUID_UORDER %APPLI_FLUID_PORDER %APPLI_FLUID_GORDER ).str();
          //-----------------//
          // solid
#elif defined(APPLI_SOLID_UORDER)
          exportDirectoryStr = ( Feel::detail::myBoostFormat(exportDirectoryOption)
                                             %APPLI_SOLID_UORDER %APPLI_SOLID_GORDER ).str();
#endif



#if defined( THERMODYNAMICS0 ) && defined( THERMODYNAMICS1 )
          //exportDirectoryStr = ( Feel::detail::myBoostFormat(exportDirectoryOption)
          //                       %APPLI_THERMO0_TORDER %APPLI_THERMO0_GORDER
          //                       %APPLI_THERMO1_TORDER %APPLI_THERMO1_GORDER ).str();
          exportDirectoryStr = ( Feel::detail::myBoostFormat(exportDirectoryOption)
                                 %Feel::detail::getPropertyThermo0PolyOrder() %Feel::detail::getPropertyThermo0GeoOrder()
                                 %Feel::detail::getPropertyThermo1PolyOrder() %Feel::detail::getPropertyThermo1GeoOrder() ).str();
#endif

          //std::cout << "exportDirectoryOption " << exportDirectoryOption << std::endl;
          //std::cout << "exportDirectoryStr " << exportDirectoryStr << std::endl;

          return exportDirectoryStr;
        }


        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//

        WorldComm const& worldCommEnvironment()
        {
            return Environment::worldComm();
        }

        //-------------------------------------------------------------------------------------//

        bool useNonStandartCompositeSpace() const
        {
            return boption(_name="master.comm.useNonStandartCompositeSpace");
        }

        WorldComm worldCommEnvironmentNonStandartCompositeSpace()
        {
            const int VelocityWorld=0;
            const int PressureWorld=1;

            std::vector<int> MapWorld(Environment::worldComm().globalSize());
            if ( Environment::worldComm().globalSize()>1 )
                {
                    for (int proc = 0 ; proc <Environment::worldComm().globalSize(); ++proc)
                        {
                            if (proc < Environment::worldComm().globalSize()/2 ) // if (proc%2==0 )
                                MapWorld[proc] = VelocityWorld;
                            else
                                MapWorld[proc] = PressureWorld;
                        }
                    return WorldComm(MapWorld);
                }
            else return Environment::worldComm();
        }

        //-------------------------------------------------------------------------------------//


    private :

        //-------------------------------------------------------------------------------------//

        void init()
        {
#if 0
            if ( Environment::vm().count( "help" ) )
            {
                if ( Environment::worldComm().globalRank()==Environment::worldComm().masterRank() )
                    std::cout << Environment::optionsDescription() << "\n";
                exit(0);
            }
#endif
            //boost::format theAppliRep = boost::format( "%1%/" ) % option(_name="exporter.directory").as<std::string>();
            boost::format theAppliRep = boost::format( "%1%/" ) % this->appliShortRepository();
            Environment::changeRepository( _directory=theAppliRep );

            if ( option(_name="rebuildOnlyAMeshPartion",_prefix="master").as<bool>() )
            {
                this->partitioningMesh();
                exit(0);
            }


            M_Ti = option(_name="bdf.time-initial").as<double>();
            M_Tf = option(_name="bdf.time-final").as<double>();
            M_dt = option(_name="bdf.time-step").as<double>();

            if (M_Ti + M_dt == M_Tf)
                M_isStationary=true;
            else
                M_isStationary=false;

            this->loadPtsEvaluate();

        }

        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//

        void
        loadPtsEvaluate()
        {
            M_listPtEval.clear();

            if (Environment::vm().count("master.evalpt-path"))
                {
                    auto ptsfile=Environment::vm()["master.evalpt-path"].as< std::string >() ;
                    std::ifstream ifstr(ptsfile.c_str(), std::ios::in);
                    while ( ! ifstr.eof() )
                        {
                            std::string namePt;
                            ifstr>>namePt;
                            node_type pt(nRealDim);
                            for (uint16_type i=0;i<nRealDim;++i)
                                { ifstr>>pt[i]; }
                            if (! ifstr.eof()) M_listPtEval.push_back(boost::make_tuple(namePt,pt));
                        }
                    ifstr.close();
                }
            //auto it=M_listPtEval.begin();
            //auto en=M_listPtEval.end();
            //for (; it!=en;++it) std::cout << "\n "<<boost::get<0>(*it)<<" "<<boost::get<1>(*it)<<"\n";
        }

        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//
        //-------------------------------------------------------------------------------------//

        void partitioningMesh()
        {
            Feel::FeelModels::Log( "Master","partitioningMesh", "start",
                            this->worldCommEnvironment(),false);

            Gmsh gmsh( nDim, nGeoOrder, this->worldCommEnvironment() );

            int partitions = option(_name="gmsh.partitions",_prefix="master.meshpartitioning").as<int>();
            gmsh.setNumberOfPartitions( partitions );
            int partitionerType = option(_name="gmsh.partitioner",_prefix="master.meshpartitioning").as<int>();
            gmsh.setPartitioner( (GMSH_PARTITIONER)partitionerType );
            int format = option(_name="gmsh.format",_prefix="master.meshpartitioning").as<int>();
            gmsh.setFileFormat( (GMSH_FORMAT)format );

            std::string filenameInputMesh = option(_name="gmsh.filename",_prefix="master.meshpartitioning").as<std::string>();
            std::string filenameOutputMesh = option(_name="gmsh.outputfilename",_prefix="master.meshpartitioning").as<std::string>();
            gmsh.rebuildPartitionMsh(filenameInputMesh,filenameOutputMesh);

            Feel::FeelModels::Log( "Master","partitioningMesh", "finish",
                            this->worldCommEnvironment(),false);
        }

        //-------------------------------------------------------------------------------------//

    private :

        Environment M_environment;

        bool M_verbose;

        double M_Ti;
        double M_Tf;
        double M_dt;

        bool M_isStationary;

        std::list<boost::tuple<std::string,node_type> > M_listPtEval;

    };

} //end namespace Feel


#endif /* __APPLIMANAGEMENT_H */

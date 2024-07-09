#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/solid/solidmechanics.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelcore/json.hpp>
#include <fmt/core.h>
#include <Eigen/Geometry> 
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/matvec.hpp>
#include <mpi.h>
#include <feel/feelmesh/meshsupportbase.hpp>
#include <feel/feelmesh/traits.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feeldiscr/localization.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feelvf/detail/ginacmatrix.hpp>
#include <fmt/ostream.h>
#include <feel/feelmesh/bvh.hpp>
#include <feel/feelmodels/fsi/fsi.hpp>
#include <feel/feeldiscr/minmax.hpp>

#include "contactforceStatic.hpp"
#include "contactforceDynamic.hpp"

using namespace Feel;
using namespace Feel::FeelModels;
using index_type = uint32_type;

template<int nDim, std::size_t residualType, typename SolidMechanics, typename DataType>
void
contactForceModels(SolidMechanics const& t, DataType & data)
{
    bool buildCstPart = data.buildCstPart();
    if(buildCstPart)
        return;
    
    fs::path path (Environment::expand(soption(_name="solid.filename")));
    json jsonCollisionForce;

    if (fs::exists(path))
    {
        std::ifstream i(path.string().c_str());
        json j = json::parse(i);
        jsonCollisionForce = j["CollisionForce"]["body"]["setup"];
    }

    // Get parameters form json
    std::string model = jsonCollisionForce["model"].get<std::string>(); 
    std::string type = jsonCollisionForce["type"].get<std::string>();
    std::string method = jsonCollisionForce["method"].get<std::string>();
    double epsilon = jsonCollisionForce["epsilon"].get<double>();
    std::vector<double> direction = jsonCollisionForce["walls"]["directions"].get<std::vector<double>>();
    double gamma0 = jsonCollisionForce["gamma_0"].get<double>();
    double theta =  jsonCollisionForce["theta"].get<double>();
    int setA = jsonCollisionForce["updateA"].get<int>();
    std::vector<double> pressurePoint = jsonCollisionForce["pressurePoint"].get<std::vector<double>>();
    int withMarker = jsonCollisionForce["withMarker"].get<int>();
    double tolerance = jsonCollisionForce["tolerance"].get<double>();
    int fixedPoint = jsonCollisionForce["fixedPoint"].get<int>();            
    int dispCondition = jsonCollisionForce["dispCondition"].get<int>();
    std::vector<double> initCondition = jsonCollisionForce["initCondition"].get<std::vector<double>>();
    int raytracing = jsonCollisionForce["raytracing"].get<int>();
    double distance = 0;
    distance = jsonCollisionForce["walls"]["distances"].get<double>(); 
    int nbr_walls = jsonCollisionForce["walls"]["nbr"].get<int>();

    // Init distance
    if (t.timeStepBase()->iteration() == 1)
        initG(t, direction, initCondition);

    // Model : steady
    if (model.compare("steady") == 0)
    {
        // Type case : unilateral
        if (type.compare("unilateral") == 0)
        {
            if (nbr_walls == 1)
            {                
                std::cout << " Steady collision algorithm " << std::endl;
                SingleUnilateralSteady<nDim, residualType>(t, data, direction, gamma0, theta, setA, pressurePoint, dispCondition, withMarker, raytracing, distance);                 
            }   
            else if (nbr_walls > 1)
            {
                // TODO
            }        
        }
        else if (type.compare("bilateral") == 0)
        {
            // TODO
        }
    }    

    // Model : dynamic
    if (model.compare("dynamic") == 0)
    {
        if (type.compare("unilateral") == 0)
        {
            if (nbr_walls == 1)
            {                
                if (fixedPoint == 1)
                {
                    std::cout << " Dynamic Fixed point collision algorithm " << std::endl;
                    SingleUnilateralDynamicFixedPoint<nDim, residualType>(t, data, direction, gamma0, theta, tolerance, setA,  dispCondition, withMarker, raytracing, distance);
                }
                else 
                {
                    if (method.compare("penalty") == 0)
                    {
                        std::cout << " Dynamic collision algorithm " << std::endl;
                        SingleUnilateralDynamicPenalty<nDim, residualType>(t, data, direction, epsilon, distance);
                    }
                    else if (method.compare("nitsche") == 0)
                    {
                        std::cout << " Dynamic collision algorithm " << std::endl;
                        SingleUnilateralDynamic<nDim, residualType>(t, data, direction, gamma0, theta,setA, dispCondition, withMarker, raytracing, distance);
                    }

                }
                
            }   
            else if (nbr_walls > 1)
            {
                // TODO
            }

        }
        else if (type.compare("bilateral") == 0)
        {
            // TODO
        }
    }
}


template<int nDim, std::size_t residualType, typename FSIModel, typename SolidMechanics, typename DataType>
void
contactForceModelsFSI(SolidMechanics const& t, typename FSIModel::element_solid_normalstressfromfluid_ptrtype sigmafn, DataType & data)
{
    bool buildCstPart = data.buildCstPart();
    if(buildCstPart)
        return;
    
    std::cout << "apply fsi contact model" << std::endl;

    // Get function sapces
    auto Xh = t.functionSpaceDisplacement();
   
    auto const& u = t.fieldDisplacement();
    
    auto nw = vec(cst(0.),cst(-1.0));
    double g = 0.00025;
    //double g = 0.0;
    
    auto dist = Pch<1>(t.mesh())->element();
    dist = project( _space=Pch<1>(t.mesh()), _range=markedfaces( t.mesh(),"contact" ), _expr = trans(nw)*(idv(u)));


    //auto [maxU,arg_maxU] = maxelt(_range=markedfaces(t.mesh(),"contact"), _element=dist);
    auto [maxU,arg_maxU] = maxelt(_range= boundaryfaces(t.mesh()), _element=dist);
   
    std::cout << "maxU : " << maxU << std::endl;
    
    // maxU - g
    if (maxU - g >= 0)
    {
        std::cout << "Contact forces have to be applied" << std::endl;
    
        /* 
        On doit trouver l'ensemble des faces en contact
        */ 

        typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );
        
        int nbrFaces = 0;
        auto const& trialDofIdToContainerId =  form2( _test=Xh,_trial=Xh).dofIdToContainerIdTest();
        for (auto const& theface : boundaryfaces(t.mesh()) )
        {                
            auto & face = boost::unwrap_ref( theface );
            int contactDof = 0;
            for( auto const& ldof : Xh->dof()->faceLocalDof( face.id() ) )
            {
                index_type thedof = ldof.index();
                thedof = trialDofIdToContainerId[ thedof ];

                if (dist[thedof] - g >= 0 )
                    contactDof++;
                    
                if (contactDof == 2)
                {
                    nbrFaces++;
                    myelts->push_back( boost::cref( face ) );
                }      
            }
        }
    
        myelts->shrink_to_fit();

        std::cout << "Number of faces in contact : " << nbrFaces << std::endl;
        auto myfaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), myelts->begin(), myelts->end(), myelts );


        sparse_matrix_ptrtype& A = data.matrix();
        vector_ptrtype& F = data.rhs();

        auto bilinearFormDD = form2( _test=Xh,_trial=Xh,_matrix=A );
        auto linearFormDisp = form1( _test=Xh, _vector=F);

        double E = 5.6e6;
        double nu = 0.4;
        double gamma = 5e5;

        //double E = 1.5e8;
        //double nu = 0.49;
        //double gamma = 1.5e8;

        double lambda = E*nu/( (1+nu)*(1-2*nu) );
        double mu =  E/(2*(1+nu));
     
        auto const Id = eye<nDim,nDim>();
        auto epst = sym(gradt(u));
        auto eps = sym(grad(u));
        auto epsv = sym(gradv(u));

        auto sigmat = (lambda*trace(epst)*Id + 2*mu*epst)*N();
        auto sigma = (lambda*trace(eps)*Id + 2*mu*eps)*N();
        auto sigmav = (lambda*trace(epsv)*Id + 2*mu*epsv)*N();

        std::cout << "add contact terms using Nitsche theta = 0" << std::endl;

        //bilinearFormDD += integrate (_range=markedfaces(t.mesh(), "contact"),_expr= inner(cst(gamma)*trans(nw)*idt(u) - trans(nw)*sigmat, trans(nw)*id(u)) ,_geomap=t.geomap() );
        //linearFormDisp += integrate (_range=markedfaces(t.mesh(), "contact"),_expr= inner(cst(gamma)*cst(g) + trans(nw)*idv(sigmafn), trans(nw)*id(u)) ,_geomap=t.geomap() ); 
        
        //linearFormDisp += integrate (_range=markedfaces(t.mesh(), "contact"),_expr= inner(trans(nw)*idv(sigmafn), trans(nw)*id(u)) ,_geomap=t.geomap() );
    
        bilinearFormDD += integrate (_range=myfaces,_expr= inner(cst(gamma)*trans(nw)*idt(u) - trans(nw)*sigmat, trans(nw)*id(u)) ,_geomap=t.geomap() );
        linearFormDisp += integrate (_range=myfaces,_expr= inner(cst(gamma)*cst(g) + trans(nw)*idv(sigmafn), trans(nw)*id(u)) ,_geomap=t.geomap() ); 
       

        auto exp =  exporter(_mesh = t.mesh(), _name = "Collision" );
        exp->addRegions();
    
        auto Vh = Pch<1>(t.mesh());
        //auto a = project(_space= Vh, _range=markedfaces(t.mesh(), "contact"), _expr = cst(gamma)*trans(nw)*idv(u) - trans(nw)*sigmav);
        auto a = project(_space= Vh, _range=myfaces, _expr = cst(gamma)*trans(nw)*idv(u) - trans(nw)*sigmav);

        //std::cout << "a : " << a << std::endl;

        //auto f = project(_space= Vh, _range=markedfaces(t.mesh(), "contact"), _expr = cst(gamma)*cst(g) + trans(nw)*idv(sigmafn));
        auto f = project(_space= Vh, _range=myfaces, _expr = cst(gamma)*cst(g) + trans(nw)*idv(sigmafn));

        //std::cout << "f : " << f << std::endl;

        //auto pressure = project(_space= Vh, _range=markedfaces(t.mesh(), "contact"), _expr = trans(nw)*sigmav);
        auto pressure = project(_space= Vh, _range=myfaces, _expr = trans(nw)*sigmav);

        auto faces = project(_space= Vh, _range=myfaces, _expr = cst(1.0));

        //std::cout << "pressure : " << pressure << std::endl;

        //t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "faces", cst(1), myfaces, "Pch1" );


        exp->add( "a", a);
        exp->add( "f", f);
        exp->add( "pressure", pressure);
        exp->add( "sigmafn", idv(sigmafn));
        exp->add( "faces", faces);
        exp->add( "dist", dist);
    
        exp->save();
    
        A->close();
        F->close();
    }
    
}

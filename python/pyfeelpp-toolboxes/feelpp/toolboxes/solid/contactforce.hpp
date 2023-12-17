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

using namespace Feel;
using namespace Feel::FeelModels;
using index_type = uint32_type;
using size_type_ = index_type;



template<int nDim, typename SolidMechanics>
void
setWalls(SolidMechanics const& t)
{  
    std::string filename = "$cfgdir/wall.geo";
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<nDim>>("Wall"), _filename = filename);
    auto exp = exporter(_mesh = mesh, _name = fmt::format("Wall"));
    exp->addRegions();
    exp->save();
}

template<int nDim, std::size_t residualType,typename SolidMechanics,typename DataType>
void
SingleUnilateralSteady( SolidMechanics t, DataType & data, const double g, std::vector<double> const& direction, const double gamma0, const double theta, int & it, double & volume ) 
{
    auto Xh = t.functionSpaceDisplacement();
    auto const& u = t.fieldDisplacement();

    auto nwall = [&direction]() 
    { 
        if constexpr(nDim == 2) 
            return vec(cst(direction[0]),cst(direction[1]));
        else if constexpr(nDim == 3)
            return vec(cst(direction[0]),cst(direction[1]),cst(direction[2]));  
    };
    
    auto nw = nwall();
    auto epst = sym(gradt(u));
    auto eps = sym(grad(u));
    auto epsv = sym(gradv(u));

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    size_type rowStartInMatrix = t.rowStartInMatrix();
    size_type colStartInMatrix = t.colStartInMatrix();

    auto bilinearFormDD = form2( _test=Xh,_trial=Xh,_matrix=A );
    auto linearFormDisp = form1( _test=Xh, _vector=F);
    
    auto const Id = eye<nDim,nDim>();
    auto Vh = Pch<1>(t.mesh());
    
    for ( auto const& [physicName,physicData] : t.physicsFromCurrentType() )
    {
        for ( std::string const& matName : t.materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& matProperties = t.materialsProperties()->materialProperties( matName );
            
            auto lameFirstExpr = matProperties.property( "Lame-first-parameter" ).exprScalar().evaluate()(0,0);
            auto lameSecondExpr = matProperties.property( "Lame-second-parameter" ).exprScalar().evaluate()(0,0);

            auto sigmat = (lameFirstExpr*trace(epst)*Id + 2*lameSecondExpr*epst)*N();
            auto sigma = (lameFirstExpr*trace(eps)*Id + 2*lameSecondExpr*eps)*N();
            auto sigmav = (lameFirstExpr*trace(epsv)*Id + 2*lameSecondExpr*epsv)*N();

            auto contactregion = project(_space=Vh, _range=boundaryfaces(t.mesh()), _expr =  (cst(gamma0)/hFace())*(trans(nw)*idv(u) - abs(cst(g) - Py())) - trans(nw)*sigmav );

            typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );
            
            it = 0;
            auto const& trialDofIdToContainerId = bilinearFormDD.dofIdToContainerIdTest();
            for (auto const& theface : boundaryfaces(t.mesh()) )
            {                
                auto & face = boost::unwrap_ref( theface );
                for( auto const& ldof : Vh->dof()->faceLocalDof( face.id() ) )
                {
                    index_type thedof = ldof.index();
                    thedof = trialDofIdToContainerId[ thedof ];

                    if (contactregion[thedof] >= -1e-8 )
                    {
                        it++;
                        myelts->push_back( boost::cref( face ) );
                    }   
                    break;
                }
            }
            
            myelts->shrink_to_fit();
            auto myfaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), myelts->begin(), myelts->end(), myelts );
                     
            if (it > 0)
            {
                bilinearFormDD += integrate (_range=myfaces,_expr= - cst(theta)/(cst(gamma0)/hFace())*inner(trans(nw)*sigmat, trans(nw)*sigma) ,_geomap=t.geomap() ); 
                bilinearFormDD += integrate (_range=myfaces,_expr= cst(1.)/(cst(gamma0)/hFace())*inner((cst(gamma0)/hFace())*trans(nw)*idt(u) - trans(nw)*sigmat, (cst(gamma0)/hFace())*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() );
                linearFormDisp += integrate (_range=myfaces,_expr= cst(1.)/(cst(gamma0)/hFace())*inner((cst(gamma0)/hFace())*abs(cst(g) - Py()), cst(gamma0)/hFace()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() ); 
            }

            // Exports 
            t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactFaces", cst(0.), elements(t.mesh()), "Pch1" );
            t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactFaces", cst(1.), myfaces, "Pch1" );
            t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactPressure", trans(nw)*sigmav, myfaces, "Pch1" );
        }
    }

    A->close();
    F->close();

    // Volume
    volume = integrate (_range = elements(t.mesh()), _expr = cst(1.0)).evaluate()( 0, 0 );
}

template<int nDim, std::size_t residualType,typename SolidMechanics,typename DataType>
void
SingleUnilateralDynamic( SolidMechanics t, DataType & data, const double g, std::vector<double> const& direction, const double gamma0, const double theta, double & Eh,double & Eh_1, double & Eh_2, double & R, double & E_tot, double & Ep, double & E, int & it, double & volume ) 
{
    auto Xh = t.functionSpaceDisplacement();
    auto const& u = t.fieldDisplacement();

    auto nwall = [&direction]() 
    { 
        if constexpr(nDim == 2) 
            return vec(cst(direction[0]),cst(direction[1]));
        else if constexpr(nDim == 3)
            return vec(cst(direction[0]),cst(direction[1]),cst(direction[2]));  
    };
    
    auto nw = nwall();
    auto epst = sym(gradt(u));
    auto eps = sym(grad(u));
    auto epsv = sym(gradv(u));

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    size_type rowStartInMatrix = t.rowStartInMatrix();
    size_type colStartInMatrix = t.colStartInMatrix();

    auto bilinearFormDD = form2( _test=Xh,_trial=Xh,_matrix=A );
    auto linearFormDisp = form1( _test=Xh, _vector=F);
    
    auto const Id = eye<nDim,nDim>();
    auto Vh = Pch<1>(t.mesh());

    
    for ( auto const& [physicName,physicData] : t.physicsFromCurrentType() )
    {
        for ( std::string const& matName : t.materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& matProperties = t.materialsProperties()->materialProperties( matName );
            
            auto lameFirstExpr = matProperties.property( "Lame-first-parameter" ).exprScalar().evaluate()(0,0);
            auto lameSecondExpr = matProperties.property( "Lame-second-parameter" ).exprScalar().evaluate()(0,0);

            auto sigmat = (lameFirstExpr*trace(epst)*Id + 2*lameSecondExpr*epst)*N();
            auto sigma = (lameFirstExpr*trace(eps)*Id + 2*lameSecondExpr*eps)*N();
            auto sigmav = (lameFirstExpr*trace(epsv)*Id + 2*lameSecondExpr*epsv)*N();

            auto contactregion = project(_space=Vh, _range=boundaryfaces(t.mesh()), _expr =  (cst(gamma0)/hFace())*(trans(nw)*idv(u) - abs(cst(g) - Py())) - trans(nw)*sigmav );

            typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );
            
            it = 0;
            auto const& trialDofIdToContainerId = bilinearFormDD.dofIdToContainerIdTest();
            for (auto const& theface : boundaryfaces(t.mesh()) )
            {                
                auto & face = boost::unwrap_ref( theface );
                for( auto const& ldof : Vh->dof()->faceLocalDof( face.id() ) )
                {
                    index_type thedof = ldof.index();
                    thedof = trialDofIdToContainerId[ thedof ];

                    if (contactregion[thedof] >= -1e-8 )
                    {
                        it++;
                        myelts->push_back( boost::cref( face ) );
                    }   
                    break;
                }
            }
            
            myelts->shrink_to_fit();
            auto myfaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), myelts->begin(), myelts->end(), myelts );
                    
            if (it > 0)
            {
                bilinearFormDD += integrate (_range=myfaces,_expr= - cst(theta)/(cst(gamma0)/hFace())*inner(trans(nw)*sigmat, trans(nw)*sigma) ,_geomap=t.geomap() ); 
                bilinearFormDD += integrate (_range=myfaces,_expr= cst(1.)/(cst(gamma0)/hFace())*inner(cst(gamma0)/hFace()*trans(nw)*idt(u) - trans(nw)*sigmat, cst(gamma0)/hFace()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() );
                linearFormDisp += integrate (_range=myfaces,_expr= cst(1.)/(cst(gamma0)/hFace())*inner(cst(gamma0)/hFace()*abs(cst(g) - Py()), cst(gamma0)/hFace()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() ); 
            }

            // Exports
            t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactFaces", cst(0.), elements(t.mesh()), "Pch1" );
            t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactFaces", cst(1.), myfaces, "Pch1" );
            t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactPressure", trans(nw)*sigmav, myfaces, "Pch1" );
            
            // Energy
            auto const& v = t.fieldVelocity();
            Eh_1 = 0.5*normL2Squared(_range = elements(t.mesh()), _expr = idv(v));
            Eh_2 = 0.5*integrate( _range= elements( t.mesh() ), _expr= inner(lameFirstExpr*trace(sym(gradv(u)))*Id + 2*lameSecondExpr*sym(gradv(u)),gradv(u)) ).evaluate()( 0,0 );
            Eh = Eh_1 + Eh_2;

            R = theta/2 * (normL2Squared( _range=myfaces, _expr= cst(1.)/sqrt(cst(gamma0)/hFace()) * trans(nw)*(lameFirstExpr*trace(sym(gradv(u)))*Id + 2*lameSecondExpr*sym(gradv(u)))*N() ) - normL2Squared( _range=myfaces, _expr= sqrt(cst(gamma0)/hFace()) * (cst(gamma0)/hFace() *trans(nw)*idv(u) - trans(nw)*(lameFirstExpr*trace(sym(gradv(u)))*Id + 2*lameSecondExpr*sym(gradv(u)))*N() ))); 
            E_tot = Eh + R;

            auto ctx = Xh->context();
            node_type t1(nDim);
            t1(0)=0.0; t1(1)=0.0;
            ctx.add( t1 );
            auto evaluateX = evaluateFromContext( _context=ctx, _expr= idv(u) );
            
            auto physicSolidData = std::static_pointer_cast<ModelPhysicSolid<nDim>>(physicData);
            for ( auto const& bodyForce : physicSolidData->bodyForces() )
            {
                auto const& theExpr = bodyForce.expr();
                //auto resToNum = toNumericValues( theExpr );
                // std::cout << resToNum[0].second << std::endl;
                Ep = integrate( _range= elements( t.mesh() ), _expr = cst(-0.1)).evaluate()( 0,0 ) ;
                Ep*=abs(evaluateX[1]);
                std::cout << "force : " << Ep * abs(evaluateX[1]) << std::endl;
                E = Eh + Ep;
            }
        }
    }

    A->close();
    F->close();

    // Volume
    volume = integrate (_range = elements(t.mesh()), _expr = cst(1.0)).evaluate()( 0, 0 );
}

template<typename SolidMechanics>
void
storeData( SolidMechanics const& t ) 
{
    auto const& u = t.fieldDisplacement();
    u.saveHDF5("solution_ref.h5");

    auto const& mesh = t.mesh();
    saveGMSHMesh(_mesh=mesh,_filename= "mesh_ref.msh" );
}


template<int nDim, typename SolidMechanics>
void
error( SolidMechanics const& t) 
{
    auto meshref = loadMesh(_mesh=new typename SolidMechanics::mesh_type, _filename = "/data/scratch/vanlandeghem/feel/convergence/np_1/mesh_ref.msh");
    auto uref = SolidMechanics::space_displacement_type::New(meshref)->elementPtr() ;
    uref->loadHDF5("/data/scratch/vanlandeghem/feel/convergence/np_1/solution_ref.h5");

    auto op_inter = opInterpolation(_domainSpace =  SolidMechanics::space_displacement_type::New(t.mesh()), _imageSpace = SolidMechanics::space_displacement_type::New(meshref) );
    auto uinter =  SolidMechanics::space_displacement_type::New(meshref)->element(); 
    op_inter->apply(t.fieldDisplacement(), uinter);

    // Test
    auto exp = exporter(_mesh = meshref, _name = fmt::format("Test_error"));
    exp->addRegions();
    exp->add("uref",*uref);
    exp->add("uinter",uinter);
    exp->save();

    auto h1norm = normH1( _range=elements(meshref), _expr=idv(*uref), _grad_expr=gradv(*uref) );
    auto h1err = normH1( _range=elements(meshref), _expr = idv(uinter) - idv(*uref), _grad_expr = gradv(uinter) - gradv(*uref));
    std::cout << "H1 relative error : " << h1err/h1norm << std::endl;

    // Get gamma0
    fs::path path (Environment::expand(soption(_name="solid.filename")));
    json jsonCollisionForce;

    if (fs::exists(path))
    {
        std::ifstream i(path.string().c_str());
        json j = json::parse(i);
        jsonCollisionForce = j["CollisionForce"]["body"]["setup"];
    }
    double gamma0 = jsonCollisionForce["gamma_0"].get<double>();

    // Get nw
    std::vector<double> direction = jsonCollisionForce["walls"]["directions"].get<std::vector<double>>();
    auto nwall = [&direction]() 
    { 
        if constexpr(nDim == 2) 
            return vec(cst(direction[0]),cst(direction[1]));
        else if constexpr(nDim == 3)
            return vec(cst(direction[0]),cst(direction[1]),cst(direction[2]));  
    };
    auto nw = nwall();
    auto const Id = eye<nDim,nDim>();

    auto epsvref = sym(gradv(*uref));
    auto epsvinter = sym(gradv(uinter));

    for ( auto const& [physicName,physicData] : t.physicsFromCurrentType() )
    {
        for ( std::string const& matName : t.materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& matProperties = t.materialsProperties()->materialProperties( matName );
            
            auto lameFirstExpr = matProperties.property( "Lame-first-parameter" ).exprScalar().evaluate()(0,0);
            auto lameSecondExpr = matProperties.property( "Lame-second-parameter" ).exprScalar().evaluate()(0,0);

            auto sigmavref = (lameFirstExpr*trace(epsvref)*Id + 2*lameSecondExpr*epsvref)*N();
            auto sigmavinter = (lameFirstExpr*trace(epsvinter)*Id + 2*lameSecondExpr*epsvinter)*N();

            auto l2norm = normL2( _range=markedfaces(meshref,"contact_wall"), _expr= cst(1.)/(cst(gamma0)/hFace()) * (cst(gamma0)/hFace()*trans(nw)*idv(*uref) - trans(nw)*sigmavref) );
            auto l2err = normL2( _range=markedfaces(meshref,"contact_wall"), _expr= sqrt(cst(gamma0)/hFace())*(cst(1.)/(cst(gamma0)/hFace()) * (cst(gamma0)/hFace()*trans(nw)*idv(*uref) - trans(nw)*sigmavref) - cst(1.)/(cst(gamma0)/hFace()) * (cst(gamma0)/hFace()*trans(nw)*idv(uinter) - trans(nw)*sigmavinter)) );
            std::cout << "Contact stress relative error : " << l2err/l2norm << std::endl;
        }
    }
}


class Measures
{
    public:
        static int nbr;
        static std::vector<int> it;
        static std::vector<double> Eh;
        static std::vector<double> Eh_1;
        static std::vector<double> Eh_2;
        static std::vector<double> R;
        static std::vector<double> E_tot;
        static std::vector<double> Ep;
        static std::vector<double> E;
        static std::vector<double> volume;
};

int Measures::nbr = 0;
std::vector<int> Measures::it = {0};
std::vector<double> Measures::Eh = {0.};
std::vector<double> Measures::Eh_1 = {0.};
std::vector<double> Measures::Eh_2 = {0.};
std::vector<double> Measures::R = {0.};
std::vector<double> Measures::E_tot = {0.};
std::vector<double> Measures::Ep = {0.};
std::vector<double> Measures::E = {0.};
std::vector<double> Measures::volume = {0.};

template<typename SolidMechanics>
void
init_Measures(SolidMechanics const& t)
{
    Measures::nbr = 0;
    Measures::it.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::Eh.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::Eh_1.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::Eh_2.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::R.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::E_tot.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::Ep.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::E.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::volume.resize(t.timeStepBase()->iterationNumber()-1);
}

template<typename SolidMechanics>
void
export_Measures(SolidMechanics const& t)
{
    std::ofstream myFile("measures.csv");
    
    myFile << "iter,Eh,Eh_1,Eh_2,R,EhR,Ep,E,it,volume\n";
    
    for(int i = 0; i < Measures::nbr; ++i)
    {
        myFile << i << "," << Measures::Eh[i] << "," << Measures::Eh_1[i] << "," << Measures::Eh_2[i] << "," << Measures::R[i] << "," << Measures::E_tot[i] << "," << Measures::Ep[i] << "," << Measures::E[i] << "," << Measures::it[i] << "," << Measures::volume[i] << "\n";
    }
    
    myFile.close();
    
    // Reset
    Measures::nbr = 0;
    Measures::it.resize(1);
    Measures::Eh.resize(1);
    Measures::Eh_1.resize(1);
    Measures::Eh_2.resize(1);
    Measures::R.resize(1);
    Measures::E_tot.resize(1);
    Measures::Ep.resize(1);
    Measures::E.resize(1);
    Measures::volume.resize(1);
}


/*
Read JSON and execute function
*/

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

    std::string model = jsonCollisionForce["model"].get<std::string>(); // steady or dynamic
    std::string type = jsonCollisionForce["type"].get<std::string>(); // unilateral or bilateral
    
    // Model : steady
    if (model.compare("steady") == 0)
    {
        double volume = 0.;
        int it = 0;

        // Type case : unilateral
        if (type.compare("unilateral") == 0)
        {
            int nbr_walls = jsonCollisionForce["walls"]["nbr"].get<int>(); // 1 or more

            if (nbr_walls == 1)
            {
                double g = jsonCollisionForce["walls"]["distances"].get<double>();
                std::vector<double> direction = jsonCollisionForce["walls"]["directions"].get<std::vector<double>>();
                double gamma0 = jsonCollisionForce["gamma_0"].get<double>();
                double theta =  jsonCollisionForce["theta"].get<double>();
                
                SingleUnilateralSteady<nDim, residualType>(t, data, g, direction, gamma0, theta, it, volume);

                // Store measures
                Measures::it[Measures::nbr] = it;
                Measures::volume[Measures::nbr] = volume;
                Measures::nbr += 1;

            }   
            else if (nbr_walls > 1)
            {
                std::vector<double> g = jsonCollisionForce["walls"]["distances"].get<std::vector<double>>();
                std::vector<std::vector<double>> directions = jsonCollisionForce["walls"]["directions"].get<std::vector<std::vector<double>>>();
                double gamma0 = jsonCollisionForce["gamma_0"].get<double>();
                double theta =  jsonCollisionForce["theta"].get<double>();
            }        
        }
        else if (type.compare("bilateral") == 0)
        {

        }
    }    

    // Model : dynamic
    if (model.compare("dynamic") == 0)
    {
        double Eh = 0.;
        double Eh_1 = 0.;
        double Eh_2 = 0.;
        double R = 0;
        double E_tot = 0;
        double Ep = 0;
        double E = 0;
        int it = 0; 
        double volume = 0;

        if (type.compare("unilateral") == 0)
        {
            int nbr_walls = jsonCollisionForce["walls"]["nbr"].get<int>();

            if (nbr_walls == 1)
            {
                double g = jsonCollisionForce["walls"]["distances"].get<double>();
                std::vector<double> direction = jsonCollisionForce["walls"]["directions"].get<std::vector<double>>();
                double gamma0 = jsonCollisionForce["gamma_0"].get<double>();
                double theta =  jsonCollisionForce["theta"].get<double>();
                
                SingleUnilateralDynamic<nDim, residualType>(t, data, g, direction, gamma0, theta, Eh, Eh_1, Eh_2, R, E_tot, Ep ,E,it, volume);

                // Store measures
                Measures::Eh[Measures::nbr] = Eh;
                Measures::Eh_1[Measures::nbr] = Eh_1;
                Measures::Eh_2[Measures::nbr] = Eh_2;
                Measures::R[Measures::nbr] = R;
                Measures::Ep[Measures::nbr] = Ep;
                Measures::E_tot[Measures::nbr] = E_tot;
                Measures::E[Measures::nbr] = E;
                Measures::it[Measures::nbr] = it;
                Measures::volume[Measures::nbr] = volume;
                Measures::nbr += 1;

            }   
            else if (nbr_walls > 1)
            {
                std::vector<double> g = jsonCollisionForce["walls"]["distances"].get<std::vector<double>>();
                std::vector<std::vector<double>> directions = jsonCollisionForce["walls"]["directions"].get<std::vector<std::vector<double>>>();
                double gamma0 = jsonCollisionForce["gamma_0"].get<double>();
                double theta =  jsonCollisionForce["theta"].get<double>();
            }

        }
        else if (type.compare("bilateral") == 0)
        {

        }
    }
    
}


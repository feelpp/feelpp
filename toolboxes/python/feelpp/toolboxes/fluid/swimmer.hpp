#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelvf/vf.hpp>


using namespace Feel;
using namespace Feel::FeelModels;
using json = nl::json;

typedef Eigen::Matrix<double, 1, Eigen::Dynamic> RowVectord;

namespace ns {
    struct MagnetoParam
    {
        std::string trajectory;
        double freq;
        double amp;
        double mx;
        double my;
        double mz;
        double bx;
        double by;
        double bz;
    };
    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(MagnetoParam,trajectory,freq,amp,mx,my,mz,bx,by,bz);
}

template<int nDim, std::size_t residualType, typename FSIModel, typename FluidMechanics, typename DataType>
void
magnetoTorqueModelFSI(FluidMechanics const& t, DataType & data)
{
    //std::cout << "apply fsi magneto torque" << std::endl;

    // Get parameters from JSON
    fs::path path (Environment::expand(soption(_name="fsi.filename")));
    json jsonMagneto;

    if (fs::exists(path))
    {
        std::ifstream i(path.string().c_str());
        json j = json::parse(i,nullptr,true,true);
        jsonMagneto = j["MagnetoTorque"]["body"]["setup"];
    }

    ns::MagnetoParam torqueParam = jsonMagneto["torqueParam"].get<ns::MagnetoParam>();
    std::string traj = torqueParam.trajectory;
    double freq = torqueParam.freq;
    double amp = torqueParam.amp;
    double bx = torqueParam.bx;
    double by = torqueParam.by;
    double bz = torqueParam.bz;
    double mx = torqueParam.mx;
    double my = torqueParam.my;
    double mz = torqueParam.mz;


    
    bool buildCstPart = data.buildCstPart();
    if(buildCstPart)
        return;
    
    // Get current orientation
    double orientation = 0;
    double x_curr = 0;
    
    for ( auto const& [bpname,bpbc] : t.bodySetBC() )
    {
        auto angle = bpbc.body().rigidRotationAngles();

        if constexpr(nDim == 2)  
            orientation = angle(0,0);
        else if constexpr(nDim == 3)
            orientation = angle(2,0);

        auto massCenter = bpbc.body().massCenter();
        x_curr = massCenter(0,0);
        //std::cout << "Current orientation : " << orientation << std::endl;
    }
    
    // Compute true theta
    
    //double theta_true = amp * sin( 2 * M_PI * freq * t.currentTime() ); 
    //double T_head = integrate( _range = markedelements( t.mesh(), "Head" ), _expr = cst(mx) * std::cos(orientation) * cst(by) * std::sin(theta_true) - cst(my) * std::sin(orientation) * cst(bx) * std::cos(theta_true)).evaluate()(0,0);
    double T_head = 0; 
    //std::cout << "Applied torue to head : " << T_head << std::endl;
    
    if (traj.compare("droite") == 0)
    {
        T_head = integrate( _range = markedelements( t.mesh(), "Head" ), _expr = cst(mx) * std::cos(orientation) * cst(by) * std::sin(2 * M_PI * freq * t.currentTime()) - cst(my) * std::sin(orientation) * cst(bx)).evaluate()(0,0);
    }
    else if (traj.compare("pipe") == 0)
    {
        T_head = 0;
    }
    else if (traj.compare("cos") == 0)
    {
        // n = (sin(x), 1), t = (1, - sin(x))

        T_head = integrate( _range = markedelements( t.mesh(), "Head" ), _expr = cst(mx) * std::cos(orientation) * (cst(bx) - cst(by) * std::sin(2 * M_PI * freq * t.currentTime()) * std::sin(x_curr)) - cst(my) * std::sin(orientation) * (cst(bx)*std::sin(x_curr) + cst(by) * std::sin(2 * M_PI * freq * t.currentTime()))).evaluate()(0,0);
    }
    else if (traj.compare("circular") == 0)
    {
        double thetaTraj = (2 * M_PI * t.currentTime())/40.;
        T_head = integrate( _range = markedelements( t.mesh(), "Head" ), _expr = cst(mx) * std::cos(orientation) * (cst(bx)*std::sin(thetaTraj) + cst(by) * std::sin(2 * M_PI * freq * t.currentTime()) * std::cos(thetaTraj)) - cst(my) * std::sin(orientation) * (cst(bx)*std::cos(thetaTraj) - cst(by) * std::sin(2 * M_PI * freq * t.currentTime())*std::sin(thetaTraj))).evaluate()(0,0);
    }

    // Add torque to newton eq
    auto r = [&data]() 
    { 
        if constexpr(residualType == 1) 
            return data.residual(); 
        else if constexpr(residualType == 0)
            return data.rhs();
    };
    
    auto rowStartInVector = t.rowStartInVector();
    
    // Add torque
    for (auto const& [bpname,bpbc] : t.bodySetBC())
    {
        size_type startBlockIndexAngularVelocity = t.startSubBlockSpaceIndex("body-bc."+bpbc.name()+".angular-velocity");
        int nLocalDofAngularVelocity = bpbc.spaceAngularVelocity()->nLocalDofWithoutGhost();
    
        if (bpbc.spaceAngularVelocity()->nLocalDofWithoutGhost() > 0)
        {
            auto const& basisToContainerGpAngularVelocityVector = r()->map().dofIdToContainerId(rowStartInVector+startBlockIndexAngularVelocity);
    
            if (residualType == 1) 
            {
                if constexpr(nDim == 2)  
                    r()->add(basisToContainerGpAngularVelocityVector[0],-T_head);
                else if constexpr(nDim == 3)  
                    r()->add(basisToContainerGpAngularVelocityVector[2],-T_head);
            }
                  
            else if (residualType == 0)
            {
                if constexpr(nDim == 2)  
                    r()->add(basisToContainerGpAngularVelocityVector[0],T_head);
                else if constexpr(nDim == 3)  
                    r()->add(basisToContainerGpAngularVelocityVector[2],T_head);
            }
        }       
                     
    }    
}



/*
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
    
    auto dist = Pch<1>(t.mesh())->element();
    dist = project( _space=Pch<1>(t.mesh()), _range=markedfaces( t.mesh(),"contact" ), _expr = trans(nw)*(idv(u)));


    //auto [maxU,arg_maxU] = maxelt(_range=markedfaces(t.mesh(),"contact"), _element=dist);
    auto [maxU,arg_maxU] = maxelt(_range= boundaryfaces(t.mesh()), _element=dist);
   
    std::cout << "maxU : " << maxU << std::endl;
    
    // maxU - g
    if (maxU - g >= 0)
    {
        std::cout << "Contact forces have to be applied" << std::endl;
    
         
        //On doit trouver l'ensemble des faces en contact
         

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

        bilinearFormDD += integrate (_range=myfaces,_expr= inner(cst(gamma)*trans(nw)*idt(u) - trans(nw)*sigmat, trans(nw)*id(u)) ,_geomap=t.geomap() );
        linearFormDisp += integrate (_range=myfaces,_expr= inner(cst(gamma)*cst(g) + trans(nw)*idv(sigmafn), trans(nw)*id(u)) ,_geomap=t.geomap() ); 
    
        A->close();
        F->close();
    }
    
}
*/
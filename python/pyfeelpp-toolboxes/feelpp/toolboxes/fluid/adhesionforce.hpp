#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feells/distancetorange.hpp>
#include <feel/feeldiscr/minmax.hpp>
#include <feel/feelvf/vf.hpp>

typedef Eigen::Matrix<double, 1, Eigen::Dynamic> RowVectord;
typedef Feel::node<double>::type node_type;

struct Collision_param
{
    RowVectord coord_C;
    RowVectord coord_S;
    double radius_C;
    double dist;        
};

struct Adhesion_param
{
    node_type coord_R;
    node_type coord_S;
    int id_R;
    //int id_S;
    double dist;        
};

class Bonds
{
    public:
        static int nbr;
        static std::vector<node_type> ligandNodes;
        static std::vector<int> receptorIDs;
        static std::vector<std::string> cellIDs;
};

int Bonds::nbr = 0;
std::vector<node_type> Bonds::ligandNodes(1000,node_type());
std::vector<int> Bonds::receptorIDs(1000,0);
std::vector<std::string> Bonds::cellIDs(1000,"");

template<typename FluidMechanics>
void
reset_bonds(FluidMechanics const& t)
{
    Bonds::nbr = 0;
    
    for (int i=0; i<100; i++)
    {
        Bonds::ligandNodes[i] = node_type();
        Bonds::receptorIDs[i] = 0;
        Bonds::cellIDs[i] = "";
    }
}

/*
double
distanceNodes( const node_type & p1, const node_type & p2 )
{
    double res=0.0;

    for ( uint16_type i=0; i<p1.size(); ++i )
    {
        res+=( p2( i )-p1( i ) )*( p2( i )-p1( i ) );
    }

    return res;
}
*/

RowVectord
unitVector( const node_type & p1, const node_type & p2, int dim, double dist )
{
    RowVectord res = RowVectord::Zero(dim);
    for ( uint16_type i=0; i<p1.size(); ++i )
    {
        res[i] = (p1( i ) - p2( i ))/dist;
    }

    return res;
}

template<typename FluidMechanics>
void
index_boundary_points(FluidMechanics const& t)
{
    std::set<int> pointIDs;

    for (auto & sface : markedfaces(t.mesh(),"Circle"))
    {
        auto & face = boost::unwrap_ref( sface );
        
        pointIDs.insert(face.point(0).id());
        pointIDs.insert(face.point(1).id());
    }

    LOG(INFO) << "IDs : " << pointIDs << std::endl;

}

template<std::size_t residualType,typename FluidMechanics,typename DataType>
void
adhesionModel(FluidMechanics const& t, DataType & data)
{
    // Get parameters
    int const dim = t.nDim;
    //auto HC = 0.04;
    auto HC = 0.04*math::pow(10.,-6);
    //auto A = 0.00005;
    auto A = 5*math::pow(10.,-20);
    auto lambda = 0.02*math::pow( 10., -6 );
    auto sigma = 0.002;
    auto T = 310.;
    auto kb = 1.38*math::pow( 10., -23 );
    auto kf0 = 0.01;
    auto kr0 = 0.01*math::pow( 10., -12 );
    auto sigmaTS = 0.001;
    auto deltat = 0.00005;
    auto N = 47*math::pow( 10., 12 );
    auto tmp = math::pow(10.,5);
    //auto tmp = 1;
    std::vector<std::string> marker_surface = {"walls"};
    
    // Get for each cell all the receptors IDs
    std::map<std::string, std::set<int>> CellReceptors;
    std::map<std::string, RowVectord> CellG;

    for ( auto const& [bpname,bpbc] : t.bodySetBC() )
    {
        std::string marker = bpbc.name();
        
        for (auto & bface : markedfaces(t.mesh(),marker))
        {
            auto & face = boost::unwrap_ref( bface ); 
            
            CellReceptors[marker].insert(face.point(0).id());
            CellReceptors[marker].insert(face.point(1).id());

            if (dim == 3)
                CellReceptors[marker].insert(face.point(2).id());
        }
        CellG[marker] = bpbc.body().massCenter();         
    }

    LOG(INFO) << "CellReceptors : " << CellReceptors << std::endl;
    
    // Compute distance to surface
    int nbr_surface = marker_surface.size();
    auto Vh = Pch<1>(t.mesh()); 
    std::vector<decltype( Vh->element() )> distToSurface;

    for (int s_i = 0;s_i < nbr_surface; s_i++)
    {
        distToSurface.push_back(distanceToRange(_space=Vh, _range=markedfaces(t.mesh(),marker_surface[s_i]),_max_distance=HC));
    }

    // Collision and adhesion detection
    std::map<std::string,Collision_param> CollisionParam;
    std::map<std::string,std::vector<Adhesion_param>> AdhesionParam;

    for (auto c_i : CellReceptors)
    {
        for (int s_i = 0;s_i < nbr_surface; s_i++)
        {
            auto [minToCi,arg_minToCi] = minelt(_range=markedfaces(t.mesh(),c_i.first), _element=distToSurface[s_i]);

            LOG(INFO) << "Distance : " << minToCi << std::endl;

            if (minToCi <= HC)
            {   
                // Collision parameter
                std::string id = c_i.first;

                Collision_param col;
                col.coord_C = arg_minToCi;
                col.dist = minToCi;

                auto distToCi = distanceToRange(_space=Vh, _range=markedfaces(t.mesh(),c_i.first),_max_distance=HC);                
                auto [minToSi,arg_minToSi] = minelt(_range=markedfaces(t.mesh(),marker_surface[s_i]), _element=distToCi);
                col.coord_S = arg_minToSi;

                double meas = integrate( _range=markedfaces(t.mesh(),c_i.first), _expr=cst(1.0)).evaluate()(0,0);
                if (dim == 2)
                    col.radius_C = meas/(2*Pi);
                else if (dim == 3)
                    col.radius_C = sqrt(meas/(4*Pi));

                CollisionParam[id] = col;


                // Adhesion parameter

                // Get possible ligands points 
                std::vector<node_type> ligands;
                //std::vector<int> ligandsIDs;
                
                for (auto & sface : markedfaces(t.mesh(),marker_surface[s_i]))
                {
                    auto & face = boost::unwrap_ref( sface );
                    
                    std::vector<unsigned int> pointIDs;
                    if ( dim == 2 ) 
                        pointIDs = {face.point(0).id(), face.point(1).id()};
                        
                    if ( dim == 3 )
                        pointIDs = {face.point(0).id(), face.point(1).id(), face.point(2).id()};
                        
                    int loc = 0;
                    while (std::find(pointIDs.begin(), pointIDs.end(), face.element(0).point(loc).id()) == pointIDs.end())
                        loc++;
                    
                    auto DOF = Vh->dof()->localToGlobal( face.element(0).id(), loc, 0 ).index();
                    auto dist = distToCi[DOF];
                    
                    //if ((dist < HC) && (std::find(Bonds::ligandIDs.begin(), Bonds::ligandIDs.end(), face.element(0).point(loc).id()) == Bonds::ligandIDs.end()))
                    if (dist < HC)
                    {
                        ligands.push_back(face.element(0).point(loc).node());
                        //ligandsIDs.push_back(face.element(0).point(loc).id());
                    }
                }
                
                for (auto & cface : markedfaces(t.mesh(),c_i.first))
                {
                    auto & face = boost::unwrap_ref( cface );

                    std::vector<unsigned int> pointIDs;
                    if ( dim == 2 ) 
                        pointIDs = {face.point(0).id(), face.point(1).id()};
                        
                    if ( dim == 3 )
                        pointIDs = {face.point(0).id(), face.point(1).id(), face.point(2).id()};
                        
                    int loc = 0;
                    while (std::find(pointIDs.begin(), pointIDs.end(), face.element(0).point(loc).id()) == pointIDs.end())
                        loc++;

                    auto DOF = Vh->dof()->localToGlobal( face.element(0).id(), loc, 0 ).index();
                    auto dist = distToSurface[s_i][DOF];

                    if ((dist < HC) && (std::find(CellReceptors[c_i.first].begin(), CellReceptors[c_i.first].end(), face.element(0).point(loc).id()) != CellReceptors[c_i.first].end()) && (std::find(Bonds::receptorIDs.begin(), Bonds::receptorIDs.end(), face.element(0).point(loc).id()) == Bonds::receptorIDs.end()))
                    {
                        Adhesion_param adh;
                        adh.coord_R = face.element(0).point(loc).node();
                        adh.id_R = face.element(0).point(loc).id();
                        
                        // Search closest ligand
                        auto closest = [&ligands]( node_type p ) 
                        {
                            auto const it = std::min_element(ligands.begin(), ligands.end(),[p] (node_type a, node_type b) {
                                return distanceNodes(p,a) < distanceNodes(p,b);
                            });
                            return *it; 
                        };

                        adh.coord_S =  closest(adh.coord_R);

                        /*
                        for (int j=0;j<ligands.size();j++)
                        {
                            if ((adh.coord_S(0) == ligands[j](0)) && (adh.coord_S(1) == ligands[j](1)))
                            {
                                adh.id_S = ligandsIDs[j];
                                break;   
                            }     
                        }
                        */

                        adh.dist = math::sqrt(distanceNodes(adh.coord_R,adh.coord_S));
                        AdhesionParam[c_i.first].push_back(adh);
                    }
                }   
            }
        } 
    }

    // Compute forces and torques
    std::map<std::string,RowVectord> forces;
    std::map<std::string,Eigen::Vector3d> torques;

    // Collision forces
    for (const auto& [id, data] : CollisionParam)
    {       
        RowVectord F = A/(8*math::sqrt(2.)) * math::sqrt(data.radius_C/data.dist) * (data.coord_C - data.coord_S)/data.dist;
        
        if (forces.find(id) == forces.end()) 
            forces[id] = F*10;
        else
            forces[id] = forces[id] + F*10;
    }


    // Adhesion forces
    
    // Old bonds
    std::vector<int> itErase;
    for (int i = 0;i < Bonds::nbr; i++)
    {
        // Compute proba of bond rupture
        auto dist = math::sqrt(distanceNodes(t.mesh()->point(Bonds::receptorIDs[i]).node(),Bonds::ligandNodes[i]));
        auto rupture = kr0*math::exp((sigma-sigmaTS)*(dist-lambda)*(dist-lambda)/(2*kb*T));
        auto probaRupture = 1-math::exp(-rupture*deltat);
        
        LOG(INFO) << "Proba rupture : " << probaRupture << std::endl;
                
        // If rupture
        if (probaRupture > ((double) rand() / (RAND_MAX)))
            itErase.push_back(i);
        else 
        {
            RowVectord F = sigma*(dist - lambda)*unitVector(Bonds::ligandNodes[i], t.mesh()->point(Bonds::receptorIDs[i]).node(), dim, dist);
            if (forces.find(Bonds::cellIDs[i]) == forces.end()) 
                forces[Bonds::cellIDs[i]] = F*tmp;
            else
                forces[Bonds::cellIDs[i]] = forces[Bonds::cellIDs[i]] + F*tmp;
            
            if (dim == 2)
            {   
                Eigen::Vector3d Gx(t.mesh()->point(Bonds::receptorIDs[i]).node()(0) - CellG[Bonds::cellIDs[i]][0],t.mesh()->point(Bonds::receptorIDs[i]).node()(1) - CellG[Bonds::cellIDs[i]][1],0.);
                Eigen::Vector3d F_(F[0],F[1],0.);

                if (torques.find(Bonds::cellIDs[i]) == torques.end()) 
                    torques[Bonds::cellIDs[i]] = -F_.cross(Gx)*tmp;
                else
                    torques[Bonds::cellIDs[i]] = torques[Bonds::cellIDs[i]] - F_.cross(Gx)*tmp;

            }    
            else if (dim == 3)
            {
                Eigen::Vector3d Gx(t.mesh()->point(Bonds::receptorIDs[i]).node()(0) - CellG[Bonds::cellIDs[i]][0],t.mesh()->point(Bonds::receptorIDs[i]).node()(1) - CellG[Bonds::cellIDs[i]][1],t.mesh()->point(Bonds::receptorIDs[i]).node()(2)- CellG[Bonds::cellIDs[i]][2]); 
                Eigen::Vector3d F_ = F;

                if (torques.find(Bonds::cellIDs[i]) == torques.end()) 
                    torques[Bonds::cellIDs[i]] = -F_.cross(Gx)*tmp;
                else
                    torques[Bonds::cellIDs[i]] = torques[Bonds::cellIDs[i]] - F_.cross(Gx)*tmp;                  
            } 
        }
    }

    // Erase 
    int addit=0;
    for (auto it : itErase)
    {
        LOG(INFO) << "Apply erase" << std::endl;
        Bonds::ligandNodes.erase(Bonds::ligandNodes.begin()+it-addit);
        Bonds::receptorIDs.erase(Bonds::receptorIDs.begin()+it-addit);
        Bonds::cellIDs.erase(Bonds::cellIDs.begin()+it-addit);
        Bonds::nbr--;
        addit++;
    }  

    std::cout << "Nbr bonds : " << Bonds::nbr << std::endl;

    // New bonds     
    for (const auto& [id, bonds] : AdhesionParam)
    {      
        for (auto data: bonds)
        {
            // Compute proba of bond construction
            auto construction = kf0*math::exp(-(sigmaTS)*(data.dist-lambda)*(data.dist-lambda)/(2*kb*T));
            auto probaConstruction = 1-math::exp(-N*construction*deltat);

            LOG(INFO) << "Proba construction : " << probaConstruction << std::endl;

            if (probaConstruction > ((double) rand() / (RAND_MAX)))
            {
                RowVectord F = sigma*(data.dist - lambda)*unitVector(data.coord_S, data.coord_R, dim, data.dist);
            
                Bonds::cellIDs[Bonds::nbr] = id;
                Bonds::ligandNodes[Bonds::nbr] = data.coord_S;
                Bonds::receptorIDs[Bonds::nbr] = data.id_R;
                Bonds::nbr +=1;
            
                if (forces.find(id) == forces.end()) 
                    forces[id] = F*tmp;
                else
                    forces[id] = forces[id] + F*tmp;
            
            
                if (dim == 2)
                {   
                    Eigen::Vector3d Gx(data.coord_R(0) - CellG[id][0],data.coord_R(1) - CellG[id][1],0.);
                    Eigen::Vector3d F_(F[0],F[1],0.);

                    if (torques.find(id) == torques.end()) 
                        torques[id] = -F_.cross(Gx)*tmp;
                    else
                        torques[id] = torques[id] - F_.cross(Gx)*tmp;

                }    

                else if (dim == 3)
                {
                    Eigen::Vector3d Gx(data.coord_R(0) - CellG[id][0],data.coord_R(1) - CellG[id][1],data.coord_R(2)- CellG[id][2]); 
                    Eigen::Vector3d F_ = F;

                    if (torques.find(id) == torques.end()) 
                        torques[id] = - F_.cross(Gx)*tmp;
                    else
                        torques[id] = torques[id] - F_.cross(Gx)*tmp;                  
                } 
            }
        }
    }

    // Insert forces
    auto r = [&data]() 
    { 
        if constexpr(residualType == 1) 
            return data.residual(); 
        else if constexpr(residualType == 0)
            return data.rhs();
    };

    auto rowStartInVector = t.rowStartInVector();

    for ( auto const& [bpname,bpbc] : t.bodySetBC() )
    {
        if (forces.find(bpbc.name()) != forces.end()) 
        {
            size_type startBlockIndexTranslationalVelocity = t.startSubBlockSpaceIndex("body-bc."+bpbc.name()+".translational-velocity");
            r()->setIsClosed(false);
            if ( bpbc.spaceTranslationalVelocity()->nLocalDofWithoutGhost() > 0 )
            {
                auto const& basisToContainerGpTranslationalVelocityVector = r()->map().dofIdToContainerId(rowStartInVector+startBlockIndexTranslationalVelocity);    

                std::cout << "ADDED FORCE : " << forces[bpbc.name()] << std::endl;
                for (int d=0;d<dim;++d)
                {   
                    if (residualType == 1) 
                        r()->add(basisToContainerGpTranslationalVelocityVector[d],-forces[bpbc.name()][d]);  
                    else if (residualType == 0)
                        r()->add(basisToContainerGpTranslationalVelocityVector[d],forces[bpbc.name()][d]);        
                }
            } 
        }
        if (torques.find(bpbc.name()) != torques.end())
        {
            size_type startBlockIndexAngularVelocity = t.startSubBlockSpaceIndex("body-bc."+bpbc.name()+".angular-velocity");
            int nLocalDofAngularVelocity = bpbc.spaceAngularVelocity()->nLocalDofWithoutGhost();

            if (bpbc.spaceAngularVelocity()->nLocalDofWithoutGhost() > 0)
            {
                auto const& basisToContainerGpAngularVelocityVector = r()->map().dofIdToContainerId(rowStartInVector+startBlockIndexAngularVelocity);

                std::cout << "ADDED TORQUE : " << torques[bpbc.name()] << std::endl;
                if (dim == 2)
                {
                    if (residualType == 1) 
                        r()->add(basisToContainerGpAngularVelocityVector[0],-torques[bpbc.name()][2]);  
                    else if (residualType == 0)
                        r()->add(basisToContainerGpAngularVelocityVector[0],torques[bpbc.name()][2]); 
                }        
                else if (dim == 3)
                {
                    for (int d=0;d<nLocalDofAngularVelocity;++d)
                    {       
                        if (residualType == 1) 
                            r()->add(basisToContainerGpAngularVelocityVector[d],-torques[bpbc.name()][d]);  
                        else if (residualType == 0)
                            r()->add(basisToContainerGpAngularVelocityVector[d],torques[bpbc.name()][d]);     
                    }         
                }     
            }
        }              
    }    
}
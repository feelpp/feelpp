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
    int id_S;
    double dist;        
};

class Bonds
{
    public:
        static int nbr;
        static std::vector<int> ligandIDs;
        static std::vector<int> receptorIDs;
        static std::vector<std::string> cellIDs;
};

int Bonds::nbr = 0;
std::vector<int> Bonds::ligandIDs(100,0);
std::vector<int> Bonds::receptorIDs(100,0);
std::vector<std::string> Bonds::cellIDs(100,"");

template<typename FluidMechanics>
void
reset_bonds(FluidMechanics const& t)
{
    Bonds::nbr = 0;
    
    for (int i=0; i<100; i++)
    {
        Bonds::ligandIDs[i] = 0;
        Bonds::receptorIDs[i] = 0;
        Bonds::cellIDs[i] = "";
    }
}

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

template<std::size_t residualType,typename FluidMechanics,typename DataType>
void
adhesionModel(FluidMechanics const& t, DataType & data)
{
    // Get parameters
    int const dim = t.nDim;
    auto HC = 0.04;
    auto A = 0.00005;
    auto lambda = 0.02;
    auto sigma = 2.;
    auto T = 310.;
    auto kb = 0.0000000138;
    auto kf0 = 0.01;
    auto kr0 = 0.01;
    auto sigmaTS = 1.;
    auto deltat = 0.001;
    std::vector<std::string> marker_surface = {"walls"};
    
    // Get for each cell all the receptors IDs
    std::map<std::string, std::set<int>> CellReceptors;
    std::map<std::string, RowVectord> CellG;
    int nbr_cells;

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
        nbr_cells ++;
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
                std::vector<int> ligandsIDs;
                
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
                    
                    if ((dist < HC) && (std::find(Bonds::ligandIDs.begin(), Bonds::ligandIDs.end(), face.element(0).point(loc).id()) == Bonds::ligandIDs.end()))
                    {
                        ligands.push_back(face.element(0).point(loc).node());
                        ligandsIDs.push_back(face.element(0).point(loc).id());
                    }
                }

                LOG(INFO) << "Ligands : " << ligands << std::endl;
                LOG(INFO) << "LigandsIDs : " << ligandsIDs << std::endl;

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

                        for (int j=0;j<ligands.size();j++)
                        {
                            if ((adh.coord_S(0) == ligands[j](0)) && (adh.coord_S(1) == ligands[j](1)))
                            {
                                adh.id_S = ligandsIDs[j];
                                break;   
                            }     
                        }

                        adh.dist = math::sqrt(distanceNodes(adh.coord_R,adh.coord_S));
                        AdhesionParam[c_i.first].push_back(adh);
                    }
                }   
            }
        } 
    }


    for (auto m : CollisionParam)
    {
        LOG(INFO) << "ID : " << m.first << std::endl;
        LOG(INFO) << "coord c : " << m.second.coord_C << std::endl;
        LOG(INFO) << "coord s : " << m.second.coord_S << std::endl;
        LOG(INFO) << "rayon c : " << m.second.radius_C << std::endl;
        LOG(INFO) << "dist : " << m.second.dist << std::endl;
    }    

    for (auto m : AdhesionParam)
    {
        std::cout << "ID : " << m.first << std::endl;
        for (auto n : m.second)
        {
            std::cout << "coord r : " << n.coord_R << std::endl;
            std::cout << "coord s : " << n.coord_S << std::endl;
            std::cout << "id r : " << n.id_R << std::endl;
            std::cout << "id s : " << n.id_S << std::endl;
            std::cout << "dist : " << n.dist << std::endl;
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
            forces[id] = F;
        else
            forces[id] = forces[id] + F;
    }


    // Adhesion forces
    
    // Old bonds
    std::vector<int> itErase;
    for (int i = 0;i < Bonds::nbr; i++)
    {
        // Compute proba of bond rupture
        auto dist = math::sqrt(distanceNodes(t.mesh()->point(Bonds::receptorIDs[i]).node(),t.mesh()->point(Bonds::ligandIDs[i]).node()));
        auto rupture = kr0*math::exp((sigma-sigmaTS)*(dist-lambda)*(dist-lambda)/(2*kb*T));
        auto probaRupture = 1-math::exp(-rupture*deltat);
        
        std::cout << "Proba rupture : " << probaRupture << std::endl;
                
        // If rupture -> use Marco Carlo
        if (probaRupture > 0.5)
            itErase.push_back(i);
        else 
        {
            RowVectord F = sigma*(dist - lambda)*unitVector(t.mesh()->point(Bonds::ligandIDs[i]).node(), t.mesh()->point(Bonds::receptorIDs[i]).node(), dim, dist);
            if (forces.find(Bonds::cellIDs[i]) == forces.end()) 
                forces[Bonds::cellIDs[i]] = F;
            else
                forces[Bonds::cellIDs[i]] = forces[Bonds::cellIDs[i]] + F;
            
            if (dim == 2)
            {   
                Eigen::Vector3d Gx(t.mesh()->point(Bonds::receptorIDs[i]).node()(0) - CellG[Bonds::cellIDs[i]][0],t.mesh()->point(Bonds::receptorIDs[i]).node()(1) - CellG[Bonds::cellIDs[i]][1],0.);
                Eigen::Vector3d F_(F[0],F[1],0.);

                if (torques.find(Bonds::cellIDs[i]) == torques.end()) 
                    torques[Bonds::cellIDs[i]] = F_.cross(Gx);
                else
                    torques[Bonds::cellIDs[i]] = torques[Bonds::cellIDs[i]] + F_.cross(Gx);

            }    
            else if (dim == 3)
            {
                Eigen::Vector3d Gx(t.mesh()->point(Bonds::receptorIDs[i]).node()(0) - CellG[Bonds::cellIDs[i]][0],t.mesh()->point(Bonds::receptorIDs[i]).node()(1) - CellG[Bonds::cellIDs[i]][1],t.mesh()->point(Bonds::receptorIDs[i]).node()(2)- CellG[Bonds::cellIDs[i]][2]); 
                Eigen::Vector3d F_ = F;

                if (torques.find(Bonds::cellIDs[i]) == torques.end()) 
                    torques[Bonds::cellIDs[i]] = F_.cross(Gx);
                else
                    torques[Bonds::cellIDs[i]] = torques[Bonds::cellIDs[i]] + F_.cross(Gx);                  
            } 
        }
    }

    // Erase 
    int addit=0;
    for (auto it : itErase)
    {
        Bonds::ligandIDs.erase(Bonds::ligandIDs.begin()+it-addit);
        Bonds::receptorIDs.erase(Bonds::receptorIDs.begin()+it-addit);
        Bonds::cellIDs.erase(Bonds::cellIDs.begin()+it-addit);
        Bonds::nbr--;
        addit++;
    }  

    // New bonds     
    for (const auto& [id, bonds] : AdhesionParam)
    {      
        for (auto data: bonds)
        {
            // Compute proba of bond construction
            auto construction = kf0*math::exp(-(sigmaTS)*(data.dist-lambda)*(data.dist-lambda)/(2*kb*T));
            auto probaConstruction = 1-math::exp(-construction*deltat);
            std::cout << "Proba construction : " << probaConstruction << std::endl;
            
            if (probaConstruction > 0.0000000001)
            {
                RowVectord F = sigma*(data.dist - lambda)*unitVector(data.coord_S, data.coord_R, dim, data.dist);
            
                Bonds::cellIDs[Bonds::nbr] = id;
                Bonds::ligandIDs[Bonds::nbr] = data.id_S;
                Bonds::receptorIDs[Bonds::nbr] = data.id_R;
                Bonds::nbr +=1;
            
                if (forces.find(id) == forces.end()) 
                    forces[id] = F;
                else
                    forces[id] = forces[id] + F;
            
            
                if (dim == 2)
                {   
                    Eigen::Vector3d Gx(data.coord_R(0) - CellG[id][0],data.coord_R(1) - CellG[id][1],0.);
                    Eigen::Vector3d F_(F[0],F[1],0.);

                    if (torques.find(id) == torques.end()) 
                        torques[id] = F_.cross(Gx);
                    else
                        torques[id] = torques[id] + F_.cross(Gx);

                }    

                else if (dim == 3)
                {
                    Eigen::Vector3d Gx(data.coord_R(0) - CellG[id][0],data.coord_R(1) - CellG[id][1],data.coord_R(2)- CellG[id][2]); 
                    Eigen::Vector3d F_ = F;

                    if (torques.find(id) == torques.end()) 
                        torques[id] = F_.cross(Gx);
                    else
                        torques[id] = torques[id] + F_.cross(Gx);                  
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
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
    RowVectord GC;
    double dist;        
};

class Bonds
{
    public:
        static int nbr;
        static std::vector<int> ligandIDs;
        static std::vector<int> receptorIDs;
};

int Bonds::nbr = 0;
std::vector<int> Bonds::ligandIDs(100,0.);
std::vector<int> Bonds::receptorIDs(100,0.);

template<typename FluidMechanics>
void
reset_bonds(FluidMechanics const& t)
{
    Bonds::nbr = 0;
    
    for (int i=0; i<100; i++)
    {
        Bonds::ligandIDs[i] = 0;
        Bonds::receptorIDs[i] = 0;
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
    double HC = 0.5;
    double A = 0.00005;
    double lambda = 0.02;
    double sigma = 2.;
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
                        adh.GC = CellG[c_i.first];
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

    for (const auto& [id, data] : CollisionParam)
    {       
        RowVectord F = A/(8*math::sqrt(2.)) * math::sqrt(data.radius_C/data.dist) * (data.coord_C - data.coord_S)/data.dist;
        
        if (forces.find(id) == forces.end()) 
            forces[id] = F;
        else
            forces[id] = forces[id] + F;
    }


    // New bonds     
    for (const auto& [id, bonds] : AdhesionParam)
    {      
        for (auto data: bonds)
        {
            // Compute proba of bond construction
            
            std::cout << "Found Bonds" << std::endl;

            RowVectord F = sigma*(data.dist - lambda)*unitVector(data.coord_S, data.coord_R, dim, data.dist);
            
            Bonds::ligandIDs[Bonds::nbr] = data.id_S;
            Bonds::receptorIDs[Bonds::nbr] = data.id_R;
            Bonds::nbr +=1;
            
            if (forces.find(id) == forces.end()) 
                forces[id] = F;
            else
                forces[id] = forces[id] + F;
            
            
            if (dim == 2)
            {   
                Eigen::Vector3d Gx(data.coord_R(0) - data.GC[0],data.coord_R(1) - data.GC[1],0.);
                Eigen::Vector3d F_(F[0],F[1],0.);

                if (torques.find(id) == torques.end()) 
                    torques[id] = F_.cross(Gx);
                else
                    torques[id] = torques[id] + F_.cross(Gx);

            }    

            else if (dim == 3)
            {
                Eigen::Vector3d Gx(data.coord_R(0) - data.GC[0],data.coord_R(1) - data.GC[1],data.coord_R(2)- data.GC[2]); 
                Eigen::Vector3d F_ = F;

                if (torques.find(id) == torques.end()) 
                    torques[id] = F_.cross(Gx);
                else
                    torques[id] = torques[id] + F_.cross(Gx);                  
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
        
                for (int d=0;d<dim;++d)
                {   
                    std::cout << "ADDED FORCE : " << forces[bpbc.name()][d] << std::endl;
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


                if (dim == 2)
                {
                    std::cout << "ADDED TORQUE : " << torques[bpbc.name()][2] << std::endl;
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
#include <feel/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/inner.hpp>
#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feelvf/functionexpr.hpp>
#include <feel/feells/distancetorange.hpp>
#include <fmt/color.h>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelfilters/filters.hpp>
#include <feel/feelmesh/meshmover.hpp>
#include <feel/feelmesh/remesher.hpp>
#include <feel/feeldiscr/minmax.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelmesh/concatenate.hpp>
#include <feel/feelmesh/submeshdata.hpp>
#if defined( FEELPP_HAS_MMG ) && defined( FEELPP_HAS_PARMMG )
#include <mmg/libmmg.h>
#include <parmmg/libparmmg.h>
#endif
#include <variant>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace Feel;

template<uint16_type Dim>
class InsertingObjects
{
    public:

        using mesh_type = Mesh<Simplex<Dim>>;
        using mesh_ptrtype = std::shared_ptr<mesh_type>;
        using space_type = FunctionSpace<mesh_type, bases<Lagrange<1>>> ;
        using space_ptrtype =  std::shared_ptr<space_type> ;
        using index_type = uint32_type;
        using size_type_ = index_type;

        /*
            - getData : Read json file and store data.
            - setAmbiantMesh : Set the ambiant mesh.
            - setObjectMesh : Set the object mesh.
            - moveObjectMesh : Adapt (translate, rotate, scale) object mesh.
            - remeshAmbiantMesh : Execute local remesh.
            - setMarkers : Set map (MarkerName, MarkerID) for all volume and boundary markers.
            - getMinElts : Get min coordinates of each volume marker.
            - setDistanceFields : Set distance fields to all object boundary markers.
            - setLevelsets : Set levelset functions for all object boundary markers.
            - setPhysicalGroups : Set markers identifiers and names. 
            - testBoundaryMarkers : Test if the boundary markers are well recovered.
            - run : Execute entire algorithm.
        */

        InsertingObjects() = default;
        void getData(std::string path);
        void setAmbiantMesh();
        void setObjectMesh();
        void moveObjectMesh();
        void setMarkers();
        void setDistanceFields();
        void setLevelsets();
        void remeshAmbiantMesh();
        void getMinElts();
        void setPhysicalGroups();
        void testBoundaryMarkers();
        int run();

        
    
    private:

        // Parameters of json file
        std::string M_fileNameAmbiant;
        std::string M_fileNameObject;
        double M_hObject;
        int M_multiPhysics;
        std::string M_markerAmbiant;
        std::string M_boundarymarkerAmbiant;
        int M_position;
        std::vector<double> M_translationVector;
        double M_scalingFactor;
        double M_rotationAngle;
        int M_remesh;
        double M_hRemesh;
        int M_export;
        int M_checkers;
        bool M_articulated = false;

        // Meshes and spaces
        mesh_ptrtype M_ObjectMesh;
        mesh_ptrtype M_AmbiantMesh;
        space_ptrtype M_VhObjectMesh;
        space_ptrtype M_VhAmbiantMesh;
        mesh_ptrtype M_ResultMesh;

        // Markers
        std::map<int,std::string> M_VolumeMarkersObject;
        std::map<int,std::string> M_VolumeMarkersAmbiant;
        std::map<int,std::string> M_BoundaryMarkersObject;
        std::map<int,std::string> M_BoundaryMarkersAmbiant;
        std::map<std::string, std::vector<size_type_> > M_AllMarkerNames;
        std::vector<std::string> M_BoundaryObjectNames;
        std::vector<std::string> M_VolumeObjectNames;
        int M_nbrObjects;

        // Markers Multi-physics
        std::vector<int> M_MultiPhysicsVolumeMarkersSolid;
        std::vector<int> M_MultiPhysicsVolumeMarkersFluid;
        std::vector<std::string> M_MultiPhysicsSolidNames;
        std::map<int,std::string> M_MultiPhysicsBoundaryMarkers;

        // Distances and levelsets
        std::vector<decltype( M_VhObjectMesh->element() )> M_distToObjectBoundaries;
        std::vector<decltype( M_VhAmbiantMesh->element() )> M_levelsetsObjectBoundaries;
        decltype( M_VhObjectMesh->element() ) M_distToObject;
        decltype( M_VhAmbiantMesh->element() ) M_levelsetObject;

        // Remesh metric
        decltype( M_VhAmbiantMesh->element() ) M_metric;

        // Physical groups recovery
        std::map<int,std::vector<double>> M_minEltsMarker;

};


template<uint16_type Dim>
void 
InsertingObjects<Dim>::getData(std::string filename)
{
    fs::path path = fs::path(Environment::findFile(filename));
    std::ifstream i(path.string().c_str());
    json jsonData = json::parse(i);

    M_fileNameAmbiant = jsonData["meshes"]["ambiant_mesh"].get<std::string>();
    M_fileNameObject = jsonData["meshes"]["object_mesh"].get<std::string>();
    M_hObject = jsonData["meshes"]["h_object"].get<double>();
    M_multiPhysics = jsonData["multiPhysics"]["has_multiPhysics"].get<double>();
    M_markerAmbiant = jsonData["multiPhysics"]["markerAmiant"].get<std::string>();
    M_boundarymarkerAmbiant = jsonData["multiPhysics"]["boundarymarkerAmbiant"].get<std::string>();
    M_position = jsonData["position"]["update_position"].get<int>();
    M_translationVector = jsonData["position"]["translation_vector"].get<std::vector<double>>();
    M_scalingFactor = jsonData["position"]["scaling_factor"].get<double>();
    M_rotationAngle = jsonData["position"]["rotation_angle"].get<double>();
    M_remesh = jsonData["remesh"]["local_remesh"].get<int>();
    M_hRemesh = jsonData["remesh"]["h_remesh"].get<double>();
    M_export = jsonData["postProcess"]["add_exports"].get<int>();
    M_checkers = jsonData["postProcess"]["add_checkers"].get<int>();
}

template<uint16_type Dim>
void 
InsertingObjects<Dim>::setAmbiantMesh()
{
    auto create_mesh = [&]() {
        if constexpr ( Dim == 2 )
        {
            return loadMesh( _mesh = new Mesh<Simplex<2>>( "mesh" ), _filename = M_fileNameAmbiant );
        }
        else if constexpr ( Dim == 3 )
        {
            return loadMesh( _mesh = new Mesh<Simplex<3, 1>>( "mesh" ), _filename = M_fileNameAmbiant );
        }
    };
    
    M_AmbiantMesh = create_mesh();
    M_VhAmbiantMesh = space_type::New(M_AmbiantMesh); 

    /*    
    auto r = remesher(M_AmbiantMesh);
    auto met = M_VhAmbiantMesh->element();

    auto [ havg, hmin, hmax ] = hMeasures( M_AmbiantMesh );
    met.on(_range=elements(M_AmbiantMesh), _expr=cst(havg));
    
    r.setMetric(met);

    M_AmbiantMesh = r.execute();
    saveGMSHMesh( _mesh=M_AmbiantMesh, _filename="remeshedsubmesh.msh" );  
    M_VhAmbiantMesh = space_type::New(M_AmbiantMesh); 

    auto exp1 = exporter( _mesh = M_AmbiantMesh , _name = fmt::format( "Remesh"));
    exp1->addRegions();
    exp1->save();
    */

}

template<uint16_type Dim>
void
InsertingObjects<Dim>::setObjectMesh()
{
    auto create_mesh = [&]() {
        if constexpr ( Dim == 2 )
        {
            return loadMesh( _mesh = new Mesh<Simplex<2>>( "mesh" ), _filename = M_fileNameObject, _h = M_hObject );
        }
        else if constexpr ( Dim == 3 )
        {
            return loadMesh( _mesh = new Mesh<Simplex<3, 1>>( "mesh" ), _filename = M_fileNameObject, _h = M_hObject );
        }
    };
    
    M_ObjectMesh = create_mesh();
    M_VhObjectMesh = space_type::New(M_ObjectMesh); 
}

template<uint16_type Dim>
void
InsertingObjects<Dim>::moveObjectMesh()
{
    auto Vhv = Pchv<1>( M_ObjectMesh );
    auto disp = Vhv->element();

    if constexpr (Dim == 2 )
        disp.on(_range=elements(M_ObjectMesh),_expr=vec(cst(M_translationVector[0])+cst(M_scalingFactor)*(Px()*cst(std::cos(M_rotationAngle)) - Py()*cst(std::sin(M_rotationAngle))) - Px(),cst(M_translationVector[1])+cst(M_scalingFactor)*(Px()*cst(std::sin(M_rotationAngle)) + Py()*cst(std::cos(M_rotationAngle))) - Py()));
    else if constexpr (Dim == 3 )
        disp.on(_range=elements(M_ObjectMesh),_expr=vec(cst(M_translationVector[0])+cst(M_scalingFactor)*(Px()*cst(std::cos(M_rotationAngle)) - Py()*cst(std::sin(M_rotationAngle))),cst(M_translationVector[1])+cst(M_scalingFactor)*(Px()*cst(std::sin(M_rotationAngle)) + Py()*cst(std::cos(M_rotationAngle))),cst(M_translationVector[2])+cst(M_scalingFactor)*Pz()) );

    meshMove( M_ObjectMesh , disp );
}


template<uint16_type Dim>
void 
InsertingObjects<Dim>::setMarkers()
{
    // Get object markers (name and id)
    for (auto m : M_ObjectMesh->markerNames())
    {
        if constexpr (Dim == 2)
        {
            if (m.second[1] == 1) 
                M_BoundaryMarkersObject.insert(std::pair<int,std::string>(m.second[0],m.first));
            else if (m.second[1] == 2) 
                M_VolumeMarkersObject.insert(std::pair<int,std::string>(m.second[0],m.first));
            
            M_AllMarkerNames[m.first] = m.second;   
        }
        else if constexpr ( Dim == 3 )
        {
            if (m.second[1] == 2)
                M_BoundaryMarkersObject.insert(std::pair<int,std::string>(m.second[0],m.first));
            else if (m.second[1] == 3) 
                M_VolumeMarkersObject.insert(std::pair<int,std::string>(m.second[0],m.first));
        
            M_AllMarkerNames[m.first] = m.second; 
        }
    }

    // Get names of object markers
    for (auto m : M_BoundaryMarkersObject)
        M_BoundaryObjectNames.push_back(m.second);

    M_nbrObjects = M_BoundaryObjectNames.size();

    for (auto m : M_VolumeMarkersObject)
        M_VolumeObjectNames.push_back(m.second);

    // Print 
    LOG(INFO) << "Object boundary markers : " << M_BoundaryObjectNames << std::endl;
    LOG(INFO) << "Object volume markers : " << M_VolumeObjectNames << std::endl;


    // Get ambiant markers (name + id)
    for (auto m : M_AmbiantMesh->markerNames())
    {
        if constexpr (Dim == 2)
        {
            if (m.second[1] == 1) 
            {
                if (M_multiPhysics == 1)
                {
                    if (M_boundarymarkerAmbiant.compare(m.first) != 0)
                    {
                        M_MultiPhysicsBoundaryMarkers.insert(std::pair<int,std::string>(m.second[0],m.first));
                    }   
                }

                M_BoundaryMarkersAmbiant.insert(std::pair<int,std::string>(m.second[0],m.first));
            }
                
            else if (m.second[1] == 2) 
            {
                if (M_multiPhysics == 1)
                {
                    if (M_markerAmbiant.compare(m.first) == 0)
                    {
                        M_VolumeMarkersAmbiant.insert(std::pair<int,std::string>(m.second[0],m.first));
                        M_MultiPhysicsVolumeMarkersFluid.push_back(m.second[0]);
                    }
                    else
                    {
                        M_MultiPhysicsVolumeMarkersSolid.push_back(m.second[0]);
                        M_MultiPhysicsSolidNames.push_back(m.first);
                    }
                        
                }
                else
                    M_VolumeMarkersAmbiant.insert(std::pair<int,std::string>(m.second[0],m.first));
            }
               
            M_AllMarkerNames[m.first] = m.second; 
        }
        else if constexpr ( Dim == 3 )
        {
            if (m.second[1] == 2)   
            {
                if (M_multiPhysics == 1)
                {
                    if (M_boundarymarkerAmbiant.compare(m.first) != 0)
                    {
                        M_MultiPhysicsBoundaryMarkers.insert(std::pair<int,std::string>(m.second[0],m.first));
                    }   
                }
                M_BoundaryMarkersAmbiant.insert(std::pair<int,std::string>(m.second[0],m.first));
            }
                
            else if (m.second[1] == 3)
            {
                if (M_multiPhysics == 1)
                {
                    if (M_markerAmbiant.compare(m.first) == 0)
                    {
                        M_VolumeMarkersAmbiant.insert(std::pair<int,std::string>(m.second[0],m.first));
                        M_MultiPhysicsVolumeMarkersFluid.push_back(m.second[0]);
                    }
                    else 
                    {
                        M_MultiPhysicsVolumeMarkersSolid.push_back(m.second[0]);
                        M_MultiPhysicsSolidNames.push_back(m.first);
                    }
                }
                else
                    M_VolumeMarkersAmbiant.insert(std::pair<int,std::string>(m.second[0],m.first));
            }
                
            M_AllMarkerNames[m.first] = m.second; 
        }
    }

    // Print 
    LOG(INFO) << "Ambiant boundary markers : " << M_BoundaryMarkersAmbiant << std::endl;
    LOG(INFO) << "Ambiant volume markers : " << M_VolumeMarkersAmbiant << std::endl;
    LOG(INFO) << "Multi-physics boundary markers : " << M_MultiPhysicsBoundaryMarkers << std::endl;
    LOG(INFO) << "Multi-physics volume markers : " << M_MultiPhysicsSolidNames << std::endl;
    LOG(INFO) << "All markers : " << M_AllMarkerNames << std::endl;

}

template<uint16_type Dim>
void 
InsertingObjects<Dim>::setDistanceFields()
{
    auto exp = exporter( _mesh = M_ObjectMesh , _name = fmt::format( "distanceFields"));
    exp->addRegions();

    for (int i = 0;i<M_nbrObjects;i++)
    {
        auto dist = distanceToRange(_space=M_VhObjectMesh, _range=markedfaces(M_ObjectMesh,M_BoundaryObjectNames[i]));

        for (int j=0;j<M_nbrObjects;j++)
        {
            if (j != i)
                dist.on(_range=markedelements(M_ObjectMesh,M_VolumeObjectNames[j]), _expr = cst(-1.));    
        }

        exp->add("dist_" + std::to_string(i), dist);
        M_distToObjectBoundaries.push_back(dist);
    }
        
    if (M_nbrObjects > 1)
    {
        M_articulated = true;
        M_distToObject = distanceToRange(_space=M_VhObjectMesh, _range=markedfaces(M_ObjectMesh,M_BoundaryObjectNames));
        exp->add("M_distToObject", M_distToObject);
    }   
    exp->save();
}

template<uint16_type Dim>
void 
InsertingObjects<Dim>::setLevelsets()
{
    auto exp = exporter( _mesh = M_AmbiantMesh , _name = fmt::format( "levelSetFunctions"));
    exp->addRegions();

    auto op_inter = opInterpolation(_domainSpace = M_VhObjectMesh, _imageSpace = M_VhAmbiantMesh );

    int it = 0;
    for (auto dist : M_distToObjectBoundaries)
    {
        auto lambda = M_VhAmbiantMesh->element(); 
        op_inter->apply(-dist, lambda);
        M_levelsetsObjectBoundaries.push_back(lambda);
        exp->add("Levelset" + std::to_string(it), lambda);
        it++;
    }

    if (M_articulated == true)
    {
        auto lambda = M_VhAmbiantMesh->element(); 
        op_inter->apply(-M_distToObject, lambda);
        M_levelsetObject = lambda;
    }
    else 
        M_levelsetObject = M_levelsetsObjectBoundaries[0];

    exp->add("M_levelsetObject", M_levelsetObject);
    exp->save();
}

template<uint16_type Dim>
void
InsertingObjects<Dim>::remeshAmbiantMesh()
{
    auto met = M_VhAmbiantMesh->element();

    auto [ havg, hmin, hmax ] = hMeasures( M_AmbiantMesh );
    met.on(_range=elements(M_AmbiantMesh), _expr=cst(hmax));

    for (int it = 0; it < M_nbrObjects; it++)
        met.on(_range=elements(M_AmbiantMesh, idv(M_levelsetsObjectBoundaries[it]), _selector=select_elements_from_expression::with_negative_values), _expr = cst(M_hRemesh));
    
    if (M_multiPhysics == 1)
    {
        for (auto name : M_MultiPhysicsSolidNames)
            met.on(_range=markedelements(M_AmbiantMesh, name), _expr = cst(M_hRemesh));
    }

    if (M_export == 1)
    {
        auto exp = exporter( _mesh = M_AmbiantMesh , _name = fmt::format( "metric"));
        exp->addRegions();
        exp->add("met", met);
        exp->save();
    }

    M_metric = met;
}

template<uint16_type Dim>
void 
InsertingObjects<Dim>::getMinElts()
{
    for (auto item :  M_VolumeMarkersObject)
    {
        auto ux = project(_space=M_VhObjectMesh, _range=elements(M_ObjectMesh), _expr= Px());
        auto [xmin,arg_xmin] = minelt(_range=markedelements(M_ObjectMesh,item.second), _element=ux);
        M_minEltsMarker[item.first].push_back(xmin);
        
        auto uy = project(_space=M_VhObjectMesh, _range=elements(M_ObjectMesh), _expr= Py());
        auto [ymin,arg_ymin] = minelt(_range=markedelements(M_ObjectMesh,item.second), _element=uy);
        M_minEltsMarker[item.first].push_back(ymin);
        
        if constexpr ( Dim == 3 )
        {
            auto uz = project(_space=M_VhObjectMesh, _range=elements(M_ObjectMesh), _expr= Pz());
            auto [zmin,arg_zmin] = minelt(_range=markedelements(M_ObjectMesh,item.second), _element=uz);
            M_minEltsMarker[item.first].push_back(zmin);
        }
    }

    if (M_multiPhysics == 1)
    {
        for (int i = 0; i < M_MultiPhysicsSolidNames.size();i++)
        {
            auto ux = project(_space=M_VhAmbiantMesh, _range=elements(M_AmbiantMesh), _expr= Px());
            auto [xmin,arg_xmin] = minelt(_range=markedelements(M_AmbiantMesh,M_MultiPhysicsSolidNames[i]), _element=ux);
            M_minEltsMarker[M_MultiPhysicsVolumeMarkersSolid[i]].push_back(xmin);
        
            auto uy = project(_space=M_VhAmbiantMesh, _range=elements(M_AmbiantMesh), _expr= Py());
            auto [ymin,arg_ymin] = minelt(_range=markedelements(M_AmbiantMesh,M_MultiPhysicsSolidNames[i]), _element=uy);
            M_minEltsMarker[M_MultiPhysicsVolumeMarkersSolid[i]].push_back(ymin);
        
            if constexpr ( Dim == 3 )
            {
                auto uz = project(_space=M_VhAmbiantMesh, _range=elements(M_AmbiantMesh), _expr= Pz());
                auto [zmin,arg_zmin] = minelt(_range=markedelements(M_AmbiantMesh,M_MultiPhysicsSolidNames[i]), _element=uz);
                M_minEltsMarker[M_MultiPhysicsVolumeMarkersSolid[i]].push_back(zmin);
            }
        }
    }
    // Print
    LOG(INFO) << "M_minEltsMarker : " << M_minEltsMarker << std::endl;
}


template<uint16_type Dim>
void 
InsertingObjects<Dim>::setPhysicalGroups()
{
    
    std::map<int,int> Mtest;
    std::vector<int> BodiesElts;

    std::cout << "Get fluid volume" << std::endl;
    // Distinguish solid and fluid elements
    for (auto & elem : M_ResultMesh->elements())
    {
        auto & [key, elt] = boost::unwrap_ref( elem );
        
        if (M_multiPhysics == 1)
        {
            if (elt.marker().value() == 2000) //Fluid
            {
                int idF = 0;
                for (auto m : M_VolumeMarkersAmbiant)
                    idF = m.first; // Fluid only has one marker
                elt.setMarker(idF);
                Mtest[elt.id()] =idF;
            } 
            else  //Solid
                BodiesElts.push_back(elt.id());
        }
        else
        {
            if (elt.marker().value() == 2) //Fluid
            {
                int idF = 0;
                for (auto m : M_VolumeMarkersAmbiant)
                    idF = m.first; // Fluid only has one marker
                elt.setMarker(idF);
                Mtest[elt.id()] =idF;
            } 
            else  //Solid
                BodiesElts.push_back(elt.id());
        }
    }

    LOG(INFO) << "Nbr BodiesElts : " << BodiesElts.size();

    std::cout << "Get solid groups" << std::endl;
    // Distinguish the different objects and get marker
    int nb_foundGroups = 0;
    while (nb_foundGroups < M_VolumeObjectNames.size() + M_MultiPhysicsSolidNames.size())
    {
        int GroupFull = 0;
        std::vector<int> GroupElts;
        int it = 0;

        double minx;
        double miny;
        double minz = 0;

        while (GroupFull == 0)
        {
            int tmp;

            if (it == 0)
            {
                tmp = BodiesElts[it];
                minx = M_ResultMesh->element(tmp).point(0).node()[0];
                miny = M_ResultMesh->element(tmp).point(0).node()[1];
                if constexpr ( Dim == 3 )
                    minz =  M_ResultMesh->element(tmp).point(0).node()[2];

                for (int i=1;i<Dim+1;i++)
                {
                    if (M_ResultMesh->element(tmp).point(i).node()[0] < minx)
                        minx = M_ResultMesh->element(tmp).point(i).node()[0];

                    if (M_ResultMesh->element(tmp).point(i).node()[1] < miny)
                        miny = M_ResultMesh->element(tmp).point(i).node()[1];

                    if (Dim == 3 && M_ResultMesh->element(tmp).point(i).node()[2] < minz)
                        minz = M_ResultMesh->element(tmp).point(i).node()[2];
                }
                
            }
            else 
                tmp = GroupElts[it];
            

            for (int i = 0; i < M_ResultMesh->element(tmp).nNeighbors();i++)
            {
                int id = M_ResultMesh->element(tmp).neighbor(i);
                if (std::find(BodiesElts.begin(), BodiesElts.end(), id) != BodiesElts.end())
                {
                    GroupElts.push_back(id);

                    for (int i=0;i<Dim+1;i++)
                    {
                        if (M_ResultMesh->element(id).point(i).node()[0] < minx)
                            minx = M_ResultMesh->element(id).point(i).node()[0];

                        if (M_ResultMesh->element(id).point(i).node()[1] < miny)
                            miny = M_ResultMesh->element(id).point(i).node()[1];

                        if (Dim == 3 && M_ResultMesh->element(id).point(i).node()[2] < minz)
                            minz = M_ResultMesh->element(id).point(i).node()[2];
                    }
                    
                    for (int j=0;j<BodiesElts.size();j++)
                    {
                        if (BodiesElts[j] == id)
                            BodiesElts.erase(BodiesElts.begin()+j);
                    }
                }  
            }

            it++;
            if (it >= GroupElts.size())
            {
                GroupFull = 1;
                LOG(INFO) << "GroupElts : " << GroupElts.size();

                int marker;

                if constexpr(Dim == 3)
                {
                    for (auto item : M_minEltsMarker)
                    {
                        if (math::abs(item.second[0]-minx) < 1.0 && math::abs(item.second[1]-miny) < 1.0 && math::abs(item.second[2]-minz) < 1.0)
                            marker = item.first;
                    }
                }
                else
                {
                    for (auto item : M_minEltsMarker)
                    {
                        if (math::abs(item.second[0]-minx) < 0.5 && math::abs(item.second[1]-miny) < 0.5)
                            marker = item.first;
                    }
                }
                
                std::cout << "Marker : " << marker << std::endl;
                for ( auto & elem : M_ResultMesh->elements() )
                {
                    auto & [key, elt] = boost::unwrap_ref( elem );
                    if (std::find(GroupElts.begin(), GroupElts.end(), key) != GroupElts.end()) 
                    {
                        elt.setMarker(marker); 
                        Mtest[elt.id()] = marker;
                    }                   
                }    
            }      
        }
        nb_foundGroups ++;
    }

    
    std::cout << "Get solid boundary markers" << std::endl;
    // Recover boundary markers
    std::map<int,int> boundaryfacesIdMarker;
    std::vector<int> eltIDs;
    
    for (auto & elem : M_ResultMesh->elements())
    {
        auto & [key, elt] = boost::unwrap_ref( elem );

        auto markerID = Mtest[elt.id()];

        for (int i=0;i<elt.nNeighbors();i++)
        {
            auto markerIDneighbor = Mtest[elt.neighbor(i)];

            //if ( M_VolumeMarkersObject.find(markerIDneighbor) == M_VolumeMarkersObject.end() && std::find(M_MultiPhysicsVolumeMarkersSolid.begin(), M_MultiPhysicsVolumeMarkersSolid.end(), markerIDneighbor) == M_MultiPhysicsVolumeMarkersSolid.end()) 
                //LOG(INFO) << "Not found Neighbor : " << markerIDneighbor << std::endl;
            //else 
            if ( M_VolumeMarkersObject.find(markerIDneighbor) != M_VolumeMarkersObject.end() || std::find(M_MultiPhysicsVolumeMarkersSolid.begin(), M_MultiPhysicsVolumeMarkersSolid.end(), markerIDneighbor) != M_MultiPhysicsVolumeMarkersSolid.end())
            {
                if (markerIDneighbor != markerID) // Boundary 
                {
                    if (std::find(eltIDs.begin(), eltIDs.end(), elt.id())==eltIDs.end() && std::find(eltIDs.begin(), eltIDs.end(), elt.neighbor(i))==eltIDs.end()) 
                    {
                        eltIDs.push_back(elt.id());
                        eltIDs.push_back(elt.neighbor(i));

                        // Get corresponding boundary marker id
                        std::string name;
                        int markerIDO;

                        if (std::find(M_MultiPhysicsVolumeMarkersSolid.begin(), M_MultiPhysicsVolumeMarkersSolid.end(), markerIDneighbor) != M_MultiPhysicsVolumeMarkersSolid.end())
                        {
                            for (int i=0;i<M_MultiPhysicsVolumeMarkersSolid.size();i++)
                            {
                                if (M_MultiPhysicsVolumeMarkersSolid[i] == markerIDneighbor)
                                {
                                    int j = 0;
                                    for (auto item : M_MultiPhysicsBoundaryMarkers)
                                    {
                                        if (j==i)
                                            markerIDO = item.first;
                                        j++;
                                    }
                                    
                                    break;
                                } 
                            }
                        }
                        else 
                        {
                            if (M_VolumeMarkersObject.find(markerIDneighbor) != M_VolumeMarkersObject.end())
                                name = M_VolumeMarkersObject[markerIDneighbor];
                            else    
                                name = M_VolumeMarkersObject[markerID];
                        
                            int pos=0;
                            while (M_VolumeObjectNames[pos].compare(name) != 0)
                                pos++;
                        
                            for (auto item : M_BoundaryMarkersObject)
                            {
                                if (M_BoundaryObjectNames[pos].compare(item.second) == 0 )
                                    markerIDO = item.first;
                            }
                        }
                        

                        // Find common face
                        for (auto & face : elt.faces())
                        {
                            for (auto & faceN : M_ResultMesh->element(elt.neighbor(i)).faces())
                            {
                                if (face->id() == faceN->id())
                                {
                                    boundaryfacesIdMarker[face->id()] = markerIDO;
                                    break;
                                }
                            }    
                        }
                    }
                }
            }       
        } 
    }

    
    std::cout << "Get fluid boundary markers" << std::endl;
    for (auto & bfaceF : boundaryfaces(M_ResultMesh))
    {
        auto & faceF = boost::unwrap_ref( bfaceF ); 

        if constexpr ( Dim == 2 ) 
        {
            std::vector<double> centerF = {(faceF.point(0).node()[0] + faceF.point(1).node()[0])/2.,(faceF.point(0).node()[1] + faceF.point(1).node()[1])/2.};

            double minDist = 0;
            int minMarker = 0;
            int first = 1;

            for (auto & bfaceI : boundaryfaces(M_AmbiantMesh))
            {
                auto & faceI = boost::unwrap_ref( bfaceI ); 
                std::vector<double> centerI = {(faceI.point(0).node()[0] + faceI.point(1).node()[0])/2.,(faceI.point(0).node()[1] + faceI.point(1).node()[1])/2.};
                auto dist = sqrt((centerI[0]-centerF[0])*(centerI[0]-centerF[0])+(centerI[1]-centerF[1])*(centerI[1]-centerF[1]));

                if (first == 1)
                {
                    minDist = dist;
                    minMarker = faceI.marker().value();
                    first = 0;
                }
                else if ((first == 0) && (dist < minDist))
                {
                    minDist = dist;
                    minMarker = faceI.marker().value();
                }
                
            }
            boundaryfacesIdMarker[faceF.id()] = minMarker;
        }
        
        
        else if constexpr ( Dim == 3 ) 
        {
            std::vector<double> centerF = {(faceF.point(0).node()[0] + faceF.point(1).node()[0] + faceF.point(2).node()[0])/3.,(faceF.point(0).node()[1] + faceF.point(1).node()[1]+ faceF.point(2).node()[1])/3.,(faceF.point(0).node()[2] + faceF.point(1).node()[2]+ faceF.point(2).node()[2])/3.};

            double minDist = 0;
            int minMarker = 0;
            int first = 1;

            for (auto & bfaceI : boundaryfaces(M_AmbiantMesh))
            {
                auto & faceI = boost::unwrap_ref( bfaceI ); 
                std::vector<double> centerI = {(faceI.point(0).node()[0] + faceI.point(1).node()[0] + faceI.point(2).node()[0])/3.,(faceI.point(0).node()[1] + faceI.point(1).node()[1]+ faceI.point(2).node()[1])/3.,(faceI.point(0).node()[2] + faceI.point(1).node()[2]+ faceI.point(2).node()[2])/3.};
                auto dist = sqrt((centerI[0]-centerF[0])*(centerI[0]-centerF[0])+(centerI[1]-centerF[1])*(centerI[1]-centerF[1])+(centerI[2]-centerF[2])*(centerI[2]-centerF[2]));

                if (first == 1)
                {
                    minDist = dist;
                    minMarker = faceI.marker().value();
                    first = 0;
                }
                else if ((first == 0) && (dist < minDist))
                {
                    minDist = dist;
                    minMarker = faceI.marker().value();
                }
                
            }
            boundaryfacesIdMarker[faceF.id()] = minMarker;
        }

    }
   

    std::cout << "Set boundary markers" << std::endl;
    for ( auto & bface : M_ResultMesh->faces() )
    {
        auto & [key, face] = boost::unwrap_ref( bface );
        if (boundaryfacesIdMarker.find(key) != boundaryfacesIdMarker.end()) 
        {
            face.setMarker(boundaryfacesIdMarker[key]); 
        }                   
    }
    

    // Set all markers
    M_ResultMesh->setMarkerNames( M_AllMarkerNames );
    M_ResultMesh->updateMeshFragmentation();
    M_ResultMesh->updateForUse();

    // Save mesh
    saveGMSHMesh( _mesh=M_ResultMesh, _filename="finalmesh.msh" );        
}

template<uint16_type Dim>
void 
InsertingObjects<Dim>::testBoundaryMarkers()
{
    auto expfinal = exporter( _mesh = M_ResultMesh , _name = fmt::format( "TestMarkers"));
    expfinal->addRegions();

    for (auto m : M_BoundaryMarkersObject)
    {
        auto test = project(_space=Pch<1>(M_ResultMesh), _range=markedfaces(M_ResultMesh,m.second),_expr=cst(m.first));
        expfinal->add(m.second,test);
    }

    for (auto m : M_BoundaryMarkersAmbiant)
    {
        auto test = project(_space=Pch<1>(M_ResultMesh), _range=markedfaces(M_ResultMesh,m.second),_expr=cst(m.first));
        expfinal->add(m.second,test);
    }

    if (M_multiPhysics == 1)
    {
        for (auto m : M_MultiPhysicsBoundaryMarkers)
        {
            auto test = project(_space=Pch<1>(M_ResultMesh), _range=markedfaces(M_ResultMesh,m.second),_expr=cst(m.first));
            expfinal->add(m.second,test);
        }
    }

    expfinal->save();
}

template<uint16_type Dim>
int
InsertingObjects<Dim>::run()
{
    // Get data
    std::cout << "Read Json file and store data" << std::endl;
    this->getData("cases/mesh_from_ls.json");

    // Set meshes
    std::cout << "Load and store meshes" << std::endl;
    this->setObjectMesh();
    this->setAmbiantMesh();

    // Move object mesh
    std::cout << "Update position of object mesh" << std::endl;
    if (M_position == 1)
        this->moveObjectMesh();
    
    // Set markers
    std::cout << "Set boundary and volume markers" << std::endl;
    this->setMarkers();

    // Set distance functions 
    std::cout << "Set the distance fields" << std::endl;
    this->setDistanceFields();

    // Set levelset functions
    std::cout << "Set the levelset functions" << std::endl;
    this->setLevelsets();

    // Set metric field
    std::cout << "Set remesh metric field" << std::endl;
    if (M_remesh == 1)
        this->remeshAmbiantMesh();
    
    // Get minimal coordinates of each solid for physical group recovery
    std::cout << "Get minimal coordinates" << std::endl;
    this->getMinElts();

    // Run remesh
    std::cout << "Run remesh" << std::endl;

    auto expfinal = exporter( _mesh = M_AmbiantMesh , _name = fmt::format( "M_AmbiantMesh"));
    expfinal->addRegions();
    expfinal->save();

    auto r = remesher(M_AmbiantMesh);
    if (M_remesh == 1)
        M_ResultMesh = r.mesh_from_ls_met(M_levelsetObject,M_metric,M_MultiPhysicsVolumeMarkersSolid, M_MultiPhysicsVolumeMarkersFluid);
    else 
        M_ResultMesh = r.mesh_from_ls(M_levelsetObject,M_MultiPhysicsVolumeMarkersSolid, M_MultiPhysicsVolumeMarkersFluid);

    // Set physical groups
    std::cout << "Physical names" << std::endl;
    this->setPhysicalGroups();

    // Checkers
    if (M_checkers == 1)
        this->testBoundaryMarkers();
        
    return 0;
}


int main( int argc, char** argv)
{
    // Initialize Feel++ Environment
    Environment env( _argc=argc, _argv=argv,
                     _about=about( _name="Inserting_objects"));
    
    InsertingObjects<2> inserting_object;
    inserting_object.run();  


    return 0;
}
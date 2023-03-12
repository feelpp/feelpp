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
            - run : Execute entire algorithm.
            - testBoundaryMarkers : Test if the boundary markers are well recovered.
            - checkers : Compute area and perimeter of each object.
        */

        InsertingObjects() = default;
        void getData(std::string path);
        void setAmbiantMesh();
        void setObjectMesh();
        void moveObjectMesh();
        void remeshAmbiantMesh();
        void setMarkers();
        void getMinElts();
        void setDistanceFields();
        void setLevelsets();
        void setPhysicalGroups();
        int run();
        void testBoundaryMarkers();
        void checkers();
        
    
    private:

        // Parameters of json file
        std::string M_fileNameAmbiant;
        std::string M_fileNameObject;
        double M_hObject;
        int M_multiPhysics;
        std::string M_markerAmbiant;
        std::string M_boundarymarkerAmbiant;
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
        std::map<int,std::string> M_BoundaryMarkersAmbiant_;
        std::map<std::string, std::vector<size_type_> > M_AllMarkerNames;
        std::vector<std::string> M_BoundaryMarkersONames;
        std::vector<std::string> M_VolumeMarkersONames;
        int M_nbrObjects;
        std::map<int,std::vector<double>> M_minEltsMarker;
        std::vector<int> M_VolumeMarkersSolidAmbiant;
        std::vector<int> M_VolumeMarkersFluidAmbiant;
        std::vector<std::string> M_VolumeMarkersSolidAmbiantNames;
    
        

        // Distances and levelsets
        std::vector<decltype( M_VhObjectMesh->element() )> M_distToObjectBoundaries;
        std::vector<decltype( M_VhAmbiantMesh->element() )> M_levelsetsObjectBoundaries;
        decltype( M_VhObjectMesh->element() ) M_distToObject;
        decltype( M_VhAmbiantMesh->element() ) M_levelsetObject;
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
InsertingObjects<Dim>::remeshAmbiantMesh()
{
    auto r = remesher(M_AmbiantMesh);
    auto met = M_VhAmbiantMesh->element();

    auto [ havg, hmin, hmax ] = hMeasures( M_AmbiantMesh );
    met.on(_range=elements(M_AmbiantMesh), _expr=cst(hmax));

    for (int it = 0; it < M_nbrObjects; it++)
        met.on(_range=elements(M_AmbiantMesh, idv(M_distToObjectBoundaries[it]), _selector=select_elements_from_expression::with_positive_values), _expr = cst(M_hRemesh));
    
    if (M_export == 1)
    {
        auto exp = exporter( _mesh = M_AmbiantMesh , _name = fmt::format( "metric"));
        exp->addRegions();
        exp->add("met", met);
        exp->save();
    }

    r.setMetric(met);
    M_AmbiantMesh = r.execute();
    M_VhAmbiantMesh = space_type::New(M_AmbiantMesh); 
}

template<uint16_type Dim>
void 
InsertingObjects<Dim>::setMarkers()
{
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
                        M_BoundaryMarkersAmbiant_.insert(std::pair<int,std::string>(m.second[0],m.first));
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
                        M_VolumeMarkersFluidAmbiant.push_back(m.second[0]);
                    }
                    else
                    {
                        M_VolumeMarkersSolidAmbiant.push_back(m.second[0]);
                        M_VolumeMarkersSolidAmbiantNames.push_back(m.first);
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
                        M_BoundaryMarkersAmbiant_.insert(std::pair<int,std::string>(m.second[0],m.first));
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
                        M_VolumeMarkersFluidAmbiant.push_back(m.second[0]);
                    }
                    else 
                    {
                        M_VolumeMarkersSolidAmbiant.push_back(m.second[0]);
                        M_VolumeMarkersSolidAmbiantNames.push_back(m.first);
                    }
                }
                else
                    M_VolumeMarkersAmbiant.insert(std::pair<int,std::string>(m.second[0],m.first));
            }
                
            M_AllMarkerNames[m.first] = m.second; 
        }
    }

    // Get names of object markers
    for (auto m : M_BoundaryMarkersObject)
        M_BoundaryMarkersONames.push_back(m.second);

    M_nbrObjects = M_BoundaryMarkersONames.size();

    for (auto m : M_VolumeMarkersObject)
        M_VolumeMarkersONames.push_back(m.second);
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

    for (int i = 0; i < M_VolumeMarkersSolidAmbiantNames.size();i++)
    {
        auto ux = project(_space=M_VhAmbiantMesh, _range=elements(M_AmbiantMesh), _expr= Px());
        auto [xmin,arg_xmin] = minelt(_range=markedelements(M_AmbiantMesh,M_VolumeMarkersSolidAmbiantNames[i]), _element=ux);
        M_minEltsMarker[M_VolumeMarkersSolidAmbiant[i]].push_back(xmin);
        
        auto uy = project(_space=M_VhAmbiantMesh, _range=elements(M_AmbiantMesh), _expr= Py());
        auto [ymin,arg_ymin] = minelt(_range=markedelements(M_AmbiantMesh,M_VolumeMarkersSolidAmbiantNames[i]), _element=uy);
        M_minEltsMarker[M_VolumeMarkersSolidAmbiant[i]].push_back(ymin);
        
        if constexpr ( Dim == 3 )
        {
            auto uz = project(_space=M_VhAmbiantMesh, _range=elements(M_AmbiantMesh), _expr= Pz());
            auto [zmin,arg_zmin] = minelt(_range=markedelements(M_AmbiantMesh,M_VolumeMarkersSolidAmbiantNames[i]), _element=uz);
            M_minEltsMarker[M_VolumeMarkersSolidAmbiant[i]].push_back(zmin);
        }
    }

}

template<uint16_type Dim>
void 
InsertingObjects<Dim>::setDistanceFields()
{
    
    for (int i = 0;i<M_nbrObjects;i++)
    {
        auto dist = distanceToRange(_space=M_VhObjectMesh, _range=markedfaces(M_ObjectMesh,M_BoundaryMarkersONames[i]));

        for (int j=0;j<M_nbrObjects;j++)
        {
            if (j != i)
                dist.on(_range=markedelements(M_ObjectMesh,M_VolumeMarkersONames[j]), _expr = cst(-1.));    
        }

        M_distToObjectBoundaries.push_back(dist);
    }
        
    if (M_nbrObjects > 1)
    {
        M_articulated = true;
        M_distToObject = distanceToRange(_space=M_VhObjectMesh, _range=markedfaces(M_ObjectMesh,M_BoundaryMarkersONames));
    }

}

template<uint16_type Dim>
void 
InsertingObjects<Dim>::setLevelsets()
{
    auto op_inter = opInterpolation(_domainSpace = M_VhObjectMesh, _imageSpace = M_VhAmbiantMesh );

    int it = 0;
    for (auto dist : M_distToObjectBoundaries)
    {
        auto lambda = M_VhAmbiantMesh->element(); 
        op_inter->apply(-dist, lambda);
        M_levelsetsObjectBoundaries.push_back(lambda);
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
}



template<uint16_type Dim>
void 
InsertingObjects<Dim>::setPhysicalGroups()
{
    
    std::map<int,int> Mtest;
    std::vector<int> BodiesElts;
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
            else //(elt.marker().value() == 3000) //Objects
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
            else //(elt.marker().value() == 3000) //Objects
                BodiesElts.push_back(elt.id());
        }
        
    }

    LOG(INFO) << "Nbr BodiesElts : " << BodiesElts.size();

    int nb_foundGroups = 0;
    while (nb_foundGroups < M_VolumeMarkersONames.size() + M_VolumeMarkersSolidAmbiant.size())
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
                        if (math::abs(item.second[0]-minx) < 0.1 && math::abs(item.second[1]-miny) < 0.1 && math::abs(item.second[2]-minz) < 0.1)
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

    // Recover faces markers
    std::map<int,int> boundaryfacesIdMarker;
    std::vector<int> eltIDs;

    for (auto & elem : M_ResultMesh->elements())
    {
        auto & [key, elt] = boost::unwrap_ref( elem );

        auto markerID = Mtest[elt.id()];

        for (int i=0;i<elt.nNeighbors();i++)
        {
            
            auto markerIDneighbor = Mtest[elt.neighbor(i)];

            if ( M_VolumeMarkersObject.find(markerIDneighbor) == M_VolumeMarkersObject.end() && std::find(M_VolumeMarkersSolidAmbiant.begin(), M_VolumeMarkersSolidAmbiant.end(), markerIDneighbor) == M_VolumeMarkersSolidAmbiant.end()) 
                LOG(INFO) << "Not found Neighbor : " << markerIDneighbor << std::endl;
                    
            else 
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

                        if (std::find(M_VolumeMarkersSolidAmbiant.begin(), M_VolumeMarkersSolidAmbiant.end(), markerIDneighbor) != M_VolumeMarkersSolidAmbiant.end())
                        {
                            for (int i=0;i<M_VolumeMarkersSolidAmbiant.size();i++)
                            {
                                if (M_VolumeMarkersSolidAmbiant[i] == markerIDneighbor)
                                {
                                    int j = 0;
                                    for (auto item : M_BoundaryMarkersAmbiant_)
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
                            while (M_VolumeMarkersONames[pos].compare(name) != 0)
                                pos++;
                        
                            for (auto item : M_BoundaryMarkersObject)
                            {
                                if (M_BoundaryMarkersONames[pos].compare(item.second) == 0 )
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

    expfinal->save();
}

template<uint16_type Dim>
void 
InsertingObjects<Dim>::checkers()
{
    int size = M_BoundaryMarkersONames.size();
    for (int i = 0;i<size;i++)
    {
        auto area = integrate( _range=markedelements(M_ResultMesh,M_VolumeMarkersONames[i]), _expr=cst(1.) ).evaluate()( 0,0 );
        auto meas = integrate( _range=markedfaces(M_ResultMesh,M_BoundaryMarkersONames[i]), _expr=cst(1.)).evaluate()(0,0);

        std::cout << "Object : " << M_VolumeMarkersONames[i] << std::endl;
        std::cout << " Area : " << area << std::endl;
        std::cout << " Permieter : " << meas << std::endl;
    }
}

template<uint16_type Dim>
int
InsertingObjects<Dim>::run()
{
    // Get data
    this->getData("cases/mesh_from_ls.json");

    // Set meshes
    this->setObjectMesh();
    this->setAmbiantMesh();

    // Move object mesh
    this->moveObjectMesh();
    
    // Set markers
    this->setMarkers();

    this->getMinElts();
    
    // Set distance functions 
    this->setDistanceFields();

    // Local remesh
    if (M_remesh == 1)
        this->remeshAmbiantMesh();

    // Set levelset
    this->setLevelsets();
    
    // Run remesh
    auto r = remesher(M_AmbiantMesh);
    M_ResultMesh = r.mesh_from_ls(M_levelsetObject,M_VolumeMarkersSolidAmbiant, M_VolumeMarkersFluidAmbiant);

    // Set markers name
    this->setPhysicalGroups();

    // Checkers
    if (M_checkers == 1)
    {
        this->checkers();
        this->testBoundaryMarkers();
    }
        
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
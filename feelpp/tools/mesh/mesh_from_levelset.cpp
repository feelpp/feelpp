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
            - getAABB : Get bounding box of the object.
            - remeshAmbiantMesh : Execute local remesh.
            - setMarkers : Set map (MarkerName, MarkerID) for all volume and boundary markers.
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
        void getAABB();
        void remeshAmbiantMesh();
        void setMarkers();
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
        std::vector<double> M_translationVector;
        double M_scalingFactor;
        double M_rotationAngle;
        int M_remesh;
        double M_tolRemesh;
        double M_hRemesh;
        double M_tolLs;
        int M_filter;
        int M_export;
        int M_checkers;
        bool M_articulated = false;

        // Meshes and spaces
        mesh_ptrtype M_ObjectMesh;
        mesh_ptrtype M_AmbiantMesh;
        space_ptrtype M_VhObjectMesh;
        space_ptrtype M_VhAmbiantMesh;
        mesh_ptrtype M_ResultMesh;

        // Coordinates of AABB
        std::vector<double> M_AABBx;
        std::vector<double> M_AABBy;
        std::vector<double> M_AABBz;

        // Markers
        std::map<int,std::string> M_VolumeMarkersObject;
        std::map<int,std::string> M_VolumeMarkersAmbiant;
        std::map<int,std::string> M_BoundaryMarkersObject;
        std::map<int,std::string> M_BoundaryMarkersAmbiant;
        std::map<std::string, std::vector<size_type_> > M_AllMarkerNames;
        std::vector<std::string> M_BoundaryMarkersONames;
        std::vector<std::string> M_VolumeMarkersONames;

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
    M_translationVector = jsonData["position"]["translation_vector"].get<std::vector<double>>();
    M_scalingFactor = jsonData["position"]["scaling_factor"].get<double>();
    M_rotationAngle = jsonData["position"]["rotation_angle"].get<double>();
    M_remesh = jsonData["remesh"]["local_remesh"].get<int>();
    M_tolRemesh = jsonData["remesh"]["tol_remesh"].get<double>();
    M_hRemesh = jsonData["remesh"]["h_remesh"].get<double>();
    M_tolLs = jsonData["ls"]["tol_ls"].get<double>();
    M_filter = jsonData["filter"]["closest_filter"].get<int>();
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
InsertingObjects<Dim>::getAABB()
{
    auto x = M_VhObjectMesh->element();
    auto y = M_VhObjectMesh->element();

    x = project(_space=M_VhObjectMesh, _range=elements(M_ObjectMesh),_expr=Px());
    y = project(_space=M_VhObjectMesh, _range=elements(M_ObjectMesh),_expr=Py());

    M_AABBx.push_back(x.min());
    M_AABBx.push_back(x.max());

    M_AABBy.push_back(y.min());
    M_AABBy.push_back(y.max());

    if constexpr(Dim == 3)
    {
        auto z = M_VhObjectMesh->element();

        z = project(_space=M_VhObjectMesh, _range=elements(M_ObjectMesh),_expr=Pz());

        M_AABBz.push_back(z.min());
        M_AABBz.push_back(z.max());
    }
}

template<uint16_type Dim>
void 
InsertingObjects<Dim>::remeshAmbiantMesh()
{
    auto r = remesher(M_AmbiantMesh);
    auto met = M_VhAmbiantMesh->element();

    met.on(_range=elements(M_AmbiantMesh), _expr=cst(M_hRemesh));
    auto [ havg, hmin, hmax ] = hMeasures( M_AmbiantMesh );
    
    met.on(_range=elements(M_AmbiantMesh, Px() - M_AABBx[0] + M_tolRemesh, _selector=select_elements_from_expression::with_negative_values), _expr=cst(hmax));
    met.on(_range=elements(M_AmbiantMesh, Py() - M_AABBy[0] + M_tolRemesh, _selector=select_elements_from_expression::with_negative_values), _expr=cst(hmax));
    met.on(_range=elements(M_AmbiantMesh, Px() - M_AABBx[1] - M_tolRemesh, _selector=select_elements_from_expression::with_positive_values), _expr=cst(hmax));
    met.on(_range=elements(M_AmbiantMesh, Py() - M_AABBy[1] - M_tolRemesh, _selector=select_elements_from_expression::with_positive_values), _expr=cst(hmax));

    if constexpr ( Dim == 3 )
    {
        met.on(_range=elements(M_AmbiantMesh, Pz() - M_AABBz[0] + M_tolRemesh,_selector=select_elements_from_expression::with_negative_values), _expr=cst(hmax));
        met.on(_range=elements(M_AmbiantMesh, Pz() - M_AABBz[1] - M_tolRemesh,_selector=select_elements_from_expression::with_positive_values), _expr=cst(hmax));
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
                M_BoundaryMarkersAmbiant.insert(std::pair<int,std::string>(m.second[0],m.first));

            else if (m.second[1] == 2) 
                M_VolumeMarkersAmbiant.insert(std::pair<int,std::string>(m.second[0],m.first));
            
            M_AllMarkerNames[m.first] = m.second; 
        }
        else if constexpr ( Dim == 3 )
        {
            if (m.second[1] == 2)
                M_BoundaryMarkersAmbiant.insert(std::pair<int,std::string>(m.second[0],m.first));
            else if (m.second[1] == 3) 
                M_VolumeMarkersAmbiant.insert(std::pair<int,std::string>(m.second[0],m.first));
        
            M_AllMarkerNames[m.first] = m.second; 
        }
    }

    // Get names of object markers
    for (auto m : M_BoundaryMarkersObject)
        M_BoundaryMarkersONames.push_back(m.second);
    for (auto m : M_VolumeMarkersObject)
        M_VolumeMarkersONames.push_back(m.second);
}

template<uint16_type Dim>
void 
InsertingObjects<Dim>::setDistanceFields()
{
    int size = M_BoundaryMarkersONames.size();

    for (int i = 0;i<size;i++)
    {
        auto dist = distanceToRange(_space=M_VhObjectMesh, _range=markedfaces(M_ObjectMesh,M_BoundaryMarkersONames[i]));

        for (int j=0;j<size;j++)
        {
            if (j != i)
                dist.on(_range=markedelements(M_ObjectMesh,M_VolumeMarkersONames[j]), _expr = cst(-1.));    
        }

        M_distToObjectBoundaries.push_back(dist);
    }
        
    if (M_BoundaryMarkersONames.size() > 1)
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
    auto VhMarkers = Pdh<0>(M_AmbiantMesh);
    //auto VhBoundaryMarkers = Pdh<1>(M_AmbiantMesh);
    
    auto VolumeMarkers = VhMarkers->element();    
    //auto BoundaryMarkers = VhBoundaryMarkers->element();   

    // Volume markers
    for (auto m : M_VolumeMarkersAmbiant)
        VolumeMarkers.on(_range=markedelements(M_AmbiantMesh,m.second), _expr = cst(m.first));
    
    int it = 0;
    for (auto m : M_VolumeMarkersObject)
    {
        VolumeMarkers.on(_range=elements(M_AmbiantMesh, idv(M_levelsetsObjectBoundaries[it]), _selector=select_elements_from_expression::with_negative_values), _expr = cst(m.first));
        it++;
    }
    
    // Fluid boundary markers
    /*
    for (auto m : M_BoundaryMarkersAmbiant)
        BoundaryMarkers.on(_range=markedfaces(M_AmbiantMesh,m.second), _expr = cst(m.first));
    */

    // Export markers on init mesh
    if (M_export == 1)
    {
        auto exp = exporter( _mesh = M_AmbiantMesh , _name = fmt::format( "InitMesh"));
        exp->addRegions();
        exp->add("VolumeMarkers", VolumeMarkers);
        //exp->add("BoundaryMarkers_P0", BoundaryMarkers);
        //exp->add("BoundaryMarkers_P1", BoundaryMarkers,"nodal");
        exp->save();
    }

    // Interpolation
    auto VhFinalMarkers = Pdh<0>(M_ResultMesh);
    //auto VhBoundaryFinalMarkers = Pdh<1>(M_ResultMesh);

    auto op_ = opInterpolation(_domainSpace = VhMarkers, _imageSpace = VhFinalMarkers);
    //auto op_boundary = opInterpolation(_domainSpace = VhBoundaryMarkers, _imageSpace = VhBoundaryFinalMarkers,_range= boundaryfaces(M_ResultMesh));

    auto Finalmarkers = VhFinalMarkers->element();
    //auto finalboundarymarkerstmp = VhBoundaryFinalMarkers->element();

    op_->apply(VolumeMarkers, Finalmarkers);
    //op_boundary->apply(BoundaryMarkers, finalboundarymarkerstmp);

    // Add closest filter
    /*
    auto finalboundarymarkers = VhBoundaryFinalMarkers->element();

    if (M_filter == 1)
    {
        std::vector<int> keys;
        for (auto it = M_BoundaryMarkersAmbiant.begin(); it != M_BoundaryMarkersAmbiant.end(); it++) 
            keys.push_back(it->first);
        std::sort(keys.begin(),keys.end());
    
        auto filter = [&keys]( double p ) 
        {
            auto const it = std::min_element(keys.begin(), keys.end(),[p] (double a, double b) {
                return std::abs(p - a) < std::abs(p - b);
            });
            return *it; 
        };

        auto filterexpr = functionExpr( filter, idv(finalboundarymarkerstmp) );
        finalboundarymarkers = project( _space=VhBoundaryFinalMarkers, _range=elements(M_ResultMesh), _expr=filterexpr );
    }
    else 
        finalboundarymarkers = finalboundarymarkerstmp;
    */

    if (M_export == 1)
    {
        auto exp_ = exporter( _mesh = M_ResultMesh , _name = fmt::format( "FinalMesh"));
        exp_->addRegions();
        exp_->add("VolumeMarkers",Finalmarkers);
        //exp_->add("Boundarymarkers_P0",finalboundarymarkers);
        //exp_->add("Boundarymarkers_P1",finalboundarymarkers,"nodal");
        exp_->save();
    }

    // Recover volume markers
    for (auto & elem : M_ResultMesh->elements())
    {
        auto & [key, elt] = boost::unwrap_ref( elem );
        auto DOF = VhFinalMarkers->dof()->localToGlobal( elt.id(), 0, 0 ).index();
        auto markerID = (int)Finalmarkers[DOF];
        elt.setMarker(markerID);
    }

    // Recover faces markers
    std::map<int,int> boundaryfacesIdMarker;
    std::vector<int> eltIDs;

    for (auto & elem : M_ResultMesh->elements())
    {
        auto & [key, elt] = boost::unwrap_ref( elem );
        auto markerID = elt.marker().value();
        
        for (int i=0;i<elt.nNeighbors();i++)
        {
            auto DOFneighbor = VhFinalMarkers->dof()->localToGlobal( elt.neighbor(i), 0, 0 ).index();
            auto markerIDneighbor = (int)Finalmarkers[DOFneighbor];

            if (markerIDneighbor != markerID) // Boundary 
            {
                if (std::find(eltIDs.begin(), eltIDs.end(), elt.id())==eltIDs.end() && std::find(eltIDs.begin(), eltIDs.end(), elt.neighbor(i))==eltIDs.end()) 
                {
                    eltIDs.push_back(elt.id());
                    eltIDs.push_back(elt.neighbor(i));

                    // Get corresponding boundary marker id
                    std::string name;
                    if (M_VolumeMarkersObject.find(markerIDneighbor) != M_VolumeMarkersObject.end())
                        name = M_VolumeMarkersObject[markerIDneighbor];
                    else    
                        name = M_VolumeMarkersObject[markerID];
                    int pos=0;
                    while (M_VolumeMarkersONames[pos].compare(name) != 0)
                        pos++;
                    int markerIDO;
                    for (auto item : M_BoundaryMarkersObject)
                    {
                        if (M_BoundaryMarkersONames[pos].compare(item.second) == 0 )
                            markerIDO = item.first;
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

    /* Using interpolation and filter
    int loc;
    std::vector<unsigned int> pointIDs;

    for (auto & bface : boundaryfaces(M_ResultMesh))
    {
        auto & face = boost::unwrap_ref( bface );

        if constexpr ( Dim == 2 ) 
            pointIDs = {face.point(0).id(), face.point(1).id()};
        if constexpr ( Dim == 3 )
            pointIDs = {face.point(0).id(), face.point(1).id(), face.point(2).id()};

        loc = 0;
        while (std::find(pointIDs.begin(), pointIDs.end(), face.element(0).point(loc).id()) == pointIDs.end())
            loc++;

        auto DOF = VhBoundaryFinalMarkers->dof()->localToGlobal( face.element(0).id(), loc, 0 ).index();
        int markerID = finalboundarymarkers[DOF];
        
        boundaryfacesIdMarker[face.id()] = markerID;
    }
    */

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

    auto testwall = project(_space=Pch<1>(M_ResultMesh), _range=markedfaces(M_ResultMesh,"wall"),_expr=cst(2));
    auto testinlet = project(_space=Pch<1>(M_ResultMesh), _range=markedfaces(M_ResultMesh,"inlet"),_expr=cst(3));
    auto testoutlet = project(_space=Pch<1>(M_ResultMesh), _range=markedfaces(M_ResultMesh,"outlet"),_expr=cst(4));
    
    expfinal->add("wall",testwall);
    expfinal->add("inlet",testinlet);
    expfinal->add("outlet",testoutlet);
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
    
    // Local remesh
    if (M_remesh == 1)
    {
        this->getAABB();
        this->remeshAmbiantMesh();
    }

    // Set markers
    this->setMarkers();
    
    // Set distance functions and levelsets
    this->setDistanceFields();
    this->setLevelsets();
    
    // Run remesh
    auto r = remesher(M_AmbiantMesh);
    M_ResultMesh = r.mesh_from_ls(M_levelsetObject);

    // Set markers name
    this->setPhysicalGroups();

    // Checkers
    if (M_checkers == 1)
        this->checkers();
    
    // Test boundary markers
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
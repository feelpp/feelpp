#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/inner.hpp>
#include <feel/feelvf/cross.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feelfilters/geotool.hpp>

namespace Feel
{

po::options_description
biotsavart_options()
{
    po::options_description bo( "BiotSavart options" );
    bo.add_options()
	( "biotsavart.unit", po::value<std::string>()->default_value("m"), "unit of the mesh (m or mm)" )
	( "biotsavart.hsize", po::value<double>()->default_value(0.1), "characteristic mesh size")
	;
    return bo;
}

template<int DimB>
class BiotSavart
{
public:
    using mesh_conductor_type = Mesh<Simplex<3> >;
    using mesh_conductor_ptrtype = std::shared_ptr<mesh_conductor_type>;
    
    using convex_box_type = Simplex<DimB,1,3>;
    using mesh_box_type = Mesh<convex_box_type>;
    using mesh_box_ptrtype = std::shared_ptr<mesh_box_type>;

    using magneticfield_space_type = FunctionSpace<mesh_box_type,
                                                   bases<Lagrange<1, Vectorial,
                                                                  Continuous, PointSetFekete> > >;
    using magneticfield_space_ptrtype = std::shared_ptr<magneticfield_space_type>;
    using magneticfield_element_type = typename magneticfield_space_type::element_type;

    using dof_point_type = boost::tuple<node_type, size_type, uint16_type >;
    using dof_points_type = typename std::vector<dof_point_type>;

private:
    mesh_conductor_ptrtype M_meshCond;
    mesh_box_ptrtype M_meshMgn;

    std::set<std::string> M_markersCond;
    std::set<std::string> M_markersMgn;

    magneticfield_space_ptrtype M_Xh;
    magneticfield_element_type M_A;
    magneticfield_element_type M_B;

    std::map<int, dof_points_type> M_dofMgn;

    double M_h;
    std::string M_unit;
public:
    BiotSavart(mesh_conductor_ptrtype const& mesh, std::set<std::string> markers = std::set<std::string>());
    template<int D = DimB>
    BiotSavart(mesh_conductor_ptrtype const& mesh, GeoTool::Node p1, GeoTool::Node p2,
	       typename std::enable_if<D==1>::type* = nullptr);
    template<int D = DimB>
    BiotSavart(mesh_conductor_ptrtype const& mesh, GeoTool::Node center, double r,
	       typename std::enable_if<D==2>::type* = nullptr);
    template<int D = DimB>
    BiotSavart(mesh_conductor_ptrtype const& mesh, GeoTool::Node center, double r,
	       typename std::enable_if<D==3>::type* = nullptr);
    template<int D = DimB>
    BiotSavart(mesh_conductor_ptrtype const& mesh,
	       GeoTool::Node center, GeoTool::Node direction, double r, double length,
	       typename std::enable_if<D==3>::type* = nullptr);
    void init();
    void createXh(std::set<std::string> markers, hana::int_<3>);
    void createXh(std::set<std::string> markers, hana::int_<2>);
    void createXh(std::set<std::string> markers, hana::int_<1>);
    template<typename Expr>
    void compute(Expr j, bool computeB = true, bool computeA = true, std::set<std::string> markers = std::set<std::string>());
    Eigen::MatrixXd magneticFieldValue(double x, double y, double z);
    Eigen::MatrixXd magneticPotentialValue(double x, double y, double z);

    mesh_box_ptrtype mesh() const { return M_meshMgn; }
    magneticfield_space_ptrtype const mgnSpace() { return M_Xh; }
    magneticfield_element_type const& magneticPotential() const { return M_A; }
    magneticfield_element_type const& magneticField() const { return M_B; }
    double hsize() const { return M_h; }
    void setHsize(double d) { M_h = d; }
    std::string unit() const { return M_unit; }
    void setUnit(std::string u) { M_unit = u == "mm" ? u : "m"; }
};

template<int DimB>
BiotSavart<DimB>::BiotSavart(mesh_conductor_ptrtype const& mesh, std::set<std::string> markers)
    : M_meshCond(mesh),
      M_markersMgn(markers),
      M_h(doption("biotsavart.hsize")),
      M_unit(soption("biotsavart.unit"))
{
    tic();
    this->createXh(markers, hana::int_<DimB>());

    M_A = M_Xh->element();
    M_B = M_Xh->element();
    toc("BiotSavart constructor");        
}

template<int DimB>
template<int D>
BiotSavart<DimB>::BiotSavart(mesh_conductor_ptrtype const& mesh, GeoTool::Node p1, GeoTool::Node p2,
			     typename std::enable_if<D==1>::type*)
    : M_meshCond(mesh),
      M_h(doption("biotsavart.hsize")),
      M_unit(soption("biotsavart.unit"))
{
    tic();
    GeoTool::Line l(M_h,"BOX", p1, p2);
    l.setMarker(_type="line", _name="box", _markerAll=true);
    M_meshMgn = l.createMesh(_mesh=new mesh_box_type, _name="box");
    this->M_Xh = magneticfield_space_type::New(M_meshMgn);

    M_A = M_Xh->element();
    M_B = M_Xh->element();
    toc("BiotSavart constructor");        
}

template<int DimB>
template<int D>
BiotSavart<DimB>::BiotSavart(mesh_conductor_ptrtype const& mesh, GeoTool::Node center, double r,
			     typename std::enable_if<D==2>::type*)
    : M_meshCond(mesh),
      M_h(doption("biotsavart.hsize")),
      M_unit(soption("biotsavart.unit"))
{
    tic();
    GeoTool::Sphere s(M_h,"BOX", center, r);
    s.setMarker(_type="surface", _name="box", _markerAll=true);
    M_meshMgn = s.createMesh(_mesh=new mesh_box_type, _name="box");
    this->M_Xh = magneticfield_space_type::New(M_meshMgn);

    M_A = M_Xh->element();
    M_B = M_Xh->element();
    toc("BiotSavart constructor");            
}

template<int DimB>
template<int D>
BiotSavart<DimB>::BiotSavart(mesh_conductor_ptrtype const& mesh, GeoTool::Node center, double r,
			     typename std::enable_if<D==3>::type*)
    : M_meshCond(mesh),
      M_h(doption("biotsavart.hsize")),
      M_unit(soption("biotsavart.unit"))
{
    tic();
    GeoTool::Sphere s(M_h,"BOX", center, r);
    s.setMarker(_type="volume", _name="box", _markerAll=true);
    M_meshMgn = s.createMesh(_mesh=new mesh_box_type, _name="box");
    this->M_Xh = magneticfield_space_type::New(M_meshMgn);

    M_A = M_Xh->element();
    M_B = M_Xh->element();
    toc("BiotSavart constructor");            
}

template<int DimB>
template<int D>
BiotSavart<DimB>::BiotSavart(mesh_conductor_ptrtype const& mesh,
			     GeoTool::Node center, GeoTool::Node direction, double r, double length,
			     typename std::enable_if<D==3>::type*)
    : M_meshCond(mesh),
      M_h(doption("biotsavart.hsize")),
      M_unit(soption("biotsavart.unit"))
{
    tic();
    GeoTool::Cylindre c(M_h,"BOX", center, direction, r, length);
    c.setMarker(_type="volume", _name="box", _markerAll=true);
    M_meshMgn = c.createMesh(_mesh=new mesh_box_type, _name="box");
    this->M_Xh = magneticfield_space_type::New(M_meshMgn);

    M_A = M_Xh->element();
    M_B = M_Xh->element();
    toc("BiotSavart constructor");            
}

template<int DimB>
void BiotSavart<DimB>::init()
{
    tic();
    auto nDofs = M_Xh->nDof();
    Feel::cout << "Computing on box with " << nDofs << " dofs" << std::endl;
    M_dofMgn.clear();
    int isIn = M_Xh->nLocalDof() > 0 ? 1 : 0;
    std::vector<int> isInGlob( Environment::worldComm().size() );
    mpi::all_gather( Environment::worldComm(), isIn, isInGlob );
    for(int i=0; i<Environment::worldComm().size(); i++)
    {
	if( isInGlob[i] )
	{
	    int dofSize;
	    dof_points_type dofM; //dofs
	    if( Environment::rank() == i )
	    {
		for ( size_type dof_id = 0; dof_id < this->M_Xh->nLocalDofWithGhost() ; ++dof_id )
		{
		    auto dofpoint = this->M_Xh->dof()->dofPoint(dof_id);
		    dofM.push_back( dofpoint );
		}
		dofSize = dofM.size();
	    }

	    mpi::broadcast( Environment::worldComm(), dofSize, i);
	    dofM.resize( dofSize );
	    mpi::broadcast( Environment::worldComm(), dofM.data(), dofSize, i);
	    M_dofMgn.insert(std::make_pair(i,dofM));
	}
    }
    toc("BiotSavart init");
}

template<int DimB>
template<typename Expr>
void BiotSavart<DimB>::compute(Expr j, bool computeB, bool computeA, std::set<std::string> markers)
{
    tic();
    // auto coeff = 1/(4*M_PI);
    // auto mu0 = 4*M_PI*1e-4; //SI unit : H.m-1 = m.kg.s-2.A-2
    auto unit = M_unit == "mm" ? 1e-4 : 1e-7;

    decltype(elements(M_meshCond)) range;
    if( markers.empty() )
        range = elements(M_meshCond);
    else
        range = markedelements(M_meshCond, markers);
    auto nElts = nelements(range,true);
    Feel::cout << "Integral on " << nElts << " elements" << std::endl;

    for( auto const& [rank,dofPoints] : M_dofMgn )
    {
	int dofSize = dofPoints.size();
	int pointSize = dofSize/3;
	std::vector<Eigen::Matrix<double,3,1>> coords( pointSize );
	for ( size_type d = 0; d < pointSize ; d++ )
	{
	    auto dofCoord = dofPoints[d*3].template get<0>();
	    Eigen::Matrix<double,3,1> coord;
	    coord << dofCoord[0], dofCoord[1], dofCoord[2];
	    coords[d] = coord;
	}

	std::vector<Eigen::Matrix<double,3,1> > mgnFields, mgnPot;
	auto dist = inner( _e1v-P(), _e1v-P(),
			   mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
    
	if( computeB )
	{
	    mgnFields = integrate(_range=range,
				  _expr=unit*cross(j, _e1v-P())/(dist*dist*dist),
				  _quad=_Q<1>()
		).template evaluate(coords);
	}
	if( computeA )
	{
	    mgnPot = integrate(_range=range,
			       _expr=unit*j/dist,
			       _quad=_Q<1>()
		).template evaluate(coords);
	}
	if( Environment::rank() == rank )
	{
	    for( int d = 0; d < dofSize; ++d )
	    {
		auto dofComp = dofPoints[d].template get<2>();
		if( computeB )
		    M_B.set(dofPoints[d].template get<1>(), mgnFields[d/3](dofComp,0));
		if( computeA )
		    M_A.set(dofPoints[d].template get<1>(), mgnPot[d/3](dofComp,0));
	    }
	}
    }
    M_B.close();
    M_A.close();
    toc("BiotSavart compute");
}

template<int DimB>
Eigen::MatrixXd BiotSavart<DimB>::magneticFieldValue(double x, double y, double z)
{
    tic();
    auto ctx = this->M_Xh->context();
    node_type t(3);
    t(0) = x;
    t(1) = y;
    t(2) = z;
    ctx.add(t);
    auto bEval = evaluateFromContext(_context=ctx,_expr=idv(this->M_B));
    toc("evaluate");
    return bEval;
}

template<int DimB>
Eigen::MatrixXd BiotSavart<DimB>::magneticPotentialValue(double x, double y, double z)
{
    tic();
    auto ctx = this->M_Xh->context();
    node_type t(3);
    t(0) = x;
    t(1) = y;
    t(2) = z;
    ctx.add(t);
    auto aEval = evaluateFromContext(_context=ctx,_expr=idv(this->M_A));
    toc("evaluate");
    return aEval;    
}

template<int DimB>
void BiotSavart<DimB>::createXh(std::set<std::string> markers, hana::int_<3>)
{
    M_meshMgn = M_meshCond;
    this->M_Xh = magneticfield_space_type::New(_mesh=M_meshCond,
                                               _range=markedelements(M_meshCond, markers));
}

template<int DimB>
void BiotSavart<DimB>::createXh(std::set<std::string> markers, hana::int_<2>)
{
    M_meshMgn = createSubmesh(_mesh=M_meshCond,
                              _range=markedfaces(M_meshCond, markers));
    this->M_Xh = magneticfield_space_type::New(M_meshMgn);
}

template<int DimB>
void BiotSavart<DimB>::createXh(std::set<std::string> markers, hana::int_<1>)
{
    M_meshMgn = createSubmesh(_mesh=M_meshCond,
                              _range=markededges(M_meshCond, markers));
    this->M_Xh = magneticfield_space_type::New(M_meshMgn);
}

template<typename Expr, typename Range>
Eigen::MatrixXd biotsavartComputeA(Expr j, Range r, double x, double y, double z, std::string unit = "m")
{
    auto coeff = unit == "m" ? 1e-7 : 1e-4;
    auto coord = vec(cst(x),cst(y),cst(z));
    auto dist = inner( coord-P(), coord-P(),
		       mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
    auto mgnPot = integrate(_range=r,
			    _expr=coeff*j/dist,
			    _quad=_Q<1>() ).evaluate();
    return mgnPot;
}

template<typename Expr, typename Range>
Eigen::MatrixXd biotsavartComputeB(Expr j, Range r, double x, double y, double z, std::string unit = "m")
{
    auto coeff = unit == "m" ? 1e-7 : 1e-4;
    auto coord = vec(cst(x),cst(y),cst(z));
    auto dist = inner( coord-P(), coord-P(),
		       mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
    auto mgnField = integrate(_range=r,
			       _expr=coeff*cross(j,coord-P())/(dist*dist*dist),
			       _quad=_Q<1>() ).evaluate();
    return mgnField;
}

} // namespace Feel

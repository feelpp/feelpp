#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/filters.hpp>
#include <feel/feeldiscr/projector.hpp>

using namespace Feel;
using namespace Feel::vf;

namespace Feel
{
  inline
    Feel::po::options_description
    makeOptions()
    {
      return Feel::feel_options() ;
    }

  template<int Dim, int Order>
    class cv : public Simget
  {
    typedef Simget super;
    public:
    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    //** Mesh **
    typedef Mesh<Simplex<Dim>> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //** Basis **
    //Scalar
    typedef bases<Lagrange<Order,Scalar>> basis_type;
    typedef boost::shared_ptr<basis_type> basis_ptrtype;

    //Scalar
    typedef FunctionSpace<mesh_type,basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;

    /**
     * Constructor
     */
    cv() : super()
    {}

    void run();
  };

  template<int Dim, int Order>
    void
    cv<Dim,Order>::run()
    {
      auto mesh   = loadMesh(_mesh = new mesh_type ,_update = MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );
      auto   Xh  = space_type::New( mesh );
      auto phi = Xh->element();
      phi = vf::project(Xh,elements(mesh),Px()+Py());
      auto   e = exporter(_mesh=     mesh,_name= "myExporter");
      e->add( "phi", phi );
      e->save();
    }
} // Feel

int main(int argc, char**argv )
{
  Environment env( _argc=argc, _argv=argv,
      _desc=makeOptions(),
      _about=about(_name="test_exporter",
        _author="Cemosis",
        _email="vincent.huber@cemosis.fr"));

  Application app1;app1.add(new cv<2,2>);app1.run();
}


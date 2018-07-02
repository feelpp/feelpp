
#ifndef _FEELPP_DICT_GEIM
#define _FEELPP_DICT_GEIM 1



#include <feel/feeldiscr/fsfunctionallinear.hpp>
#include <vector>
#include <feel/feelvf/vf.hpp>

namespace Feel
{

// Linear form for sensor modelising
template<typename Space>
class GeimFunctionalModel: public FsFunctionalLinear<Space>
{
public:

    // -- TYPEDEFS --
    typedef GeimFunctionalModel<Space> this_type;
    typedef FsFunctionalLinear<Space> super_type;

    typedef Space space_type;

    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename super_type::element_type element_type;

    typedef typename super_type::value_type value_type;

    typedef typename super_type::backend_type backend_type;
    typedef typename super_type::backend_ptrtype backend_ptrtype;

    typedef typename super_type::vector_type vector_type;
    typedef typename super_type::vector_ptrtype vector_ptrtype;

    GeimFunctionalModel( space_ptrtype space, std::vector<double> const& center = std::vector<double>(), double radius =0. ):
        super_type(space),
        M_center(center),
        M_radius(radius),
        M_space(space)
    {
    }

    virtual ~GeimFunctionalModel(){}

    void setCenter( std::vector<double> const center )
    {
        M_center=center;
    }

    void setRadius( double const radius )
    {
        M_radius=radius;
    }

    std::vector<double> center() const
    {
        return M_center;
    }

    double radius() const
    {
        return M_radius;
    }

    void defineForm()
    {   auto v=M_space->element();
        //auto phi=exp((Px()-M_center)*(Px()-M_center[0])/(2*M_radius*M_radius));
        auto phi= this->phiExpr( mpl::int_< space_type::nDim >() );
        auto expr=integrate(_range=elements(M_space->mesh()), _expr=id(v)*phi);
        //this->operator=( expr );
        this->close();
    }

private:
    auto phiExpr( mpl::int_<1> /**/ )
    {
         return exp( pow(Px()-M_center[0],2)/(2*std::pow(M_radius,2)));
    }
    auto phiExpr( mpl::int_<2> /**/ )
    {
         return exp(( pow(Px()-M_center[0],2)+pow(Py()-M_center[1],2))/(2*std::pow(M_radius,2)));
    }
    auto phiExpr( mpl::int_<3> /**/ )
    {
         return exp(( pow(Px()-M_center[0],2)+pow(Py()-M_center[1],2)+pow(Pz()-M_center[2],2))/(2*std::pow(M_radius,2)));
    }
    std::vector<double> M_center;
    double M_radius;
    space_ptrtype M_space;

};//class GeimFunctionalModel

template<typename Space>
class DictionnaryGeim : public std::vector<GeimFunctionalModel<Space>>
{
public:

    // -- TYPEDEFS --
    typedef DictionnaryGeim<Space> this_type;
    typedef GeimFunctionalModel<Space> super_type;

    typedef Space space_type;

    typedef boost::shared_ptr<space_type> space_ptrtype;

    DictionnaryGeim( space_type space, std::vector<std::vector<double>> Centers, std::vector<double> Radiuss, int size):
        M_Centers(Centers),
        M_Radiuss(Radiuss),
        M_size(size)
    {
    }

    virtual ~DictionnaryGeim() {}

    void run()
    {
        this->resize( M_size, super_type( M_space ) );
        for ( int i = 0; i < this->size(); ++i )
        {
            this->operator[](i).setCenter( M_Centers[i] );
            this->operator[](i).setRadius( M_Radiuss[i] );
            this->operator[](i).defineForm();
        }
        this->close();
    }

    void addElement( space_type space, std::vector<double> center, double radius )
    {
        M_size++;
        super_type newElement( space, center, radius );
        this->push_back( newElement );
        M_Centers.push_back( center );
        M_Radiuss.push_back( radius );
        this->operator[](M_size-1).defineForm();
        this->close();
    }

private:

  int M_size;
  std::vector< std::vector<double> > M_Centers;
  std::vector<double> M_Radiuss;
  boost::shared_ptr<Space> M_space;

};//class DictionnaryGeim

}//Feel

#endif

#ifndef FEEL_OMMODEL_HPP
#define FEEL_OMMODEL_HPP

#include <boost/property_tree/xml_parser.hpp>

namespace Feel
{


class OMModel
{
public :
    OMModel() :
        M_input_edited( true )
    {}

    virtual int run()=0;

protected :
    bool M_input_edited;
    boost::property_tree::ptree M_input_ptree;
};


}//namespace Feel
#endif

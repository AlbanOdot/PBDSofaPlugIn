#ifndef ORIENTEDCONSTRAINT_HPP
#define ORIENTEDCONSTRAINT_HPP
#include "PBDFEMConstraint.hpp"
#include "../../InternalData/Data/PBDOrientation.hpp"
class OrientedConstraint : public PBDFEMConstraint<sofa::defaulttype::Rigid3Types>
{
    typedef sofa::defaulttype::Vec3 vec3;
public:
    OrientedConstraint() : PBDFEMConstraint(){}
    PBDOrientation& orientation() { return  m_orientation;}
protected:
    PBDOrientation m_orientation;
};
#endif

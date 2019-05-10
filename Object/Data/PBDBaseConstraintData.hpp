#ifndef PBDBASECONSTRAINTDATA_HPP
#define PBDBASECONSTRAINTDATA_HPP

#include <SofaBaseMechanics/MechanicalObject.h>

class PBDBaseConstraintData
{
public:
    typedef sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > Mech;
    typedef sofa::core::topology::BaseMeshTopology  Topo;

    PBDBaseConstraintData(Mech * m = nullptr, Topo* t = nullptr) : m_mechanicalObject(m),m_sofa_topology(t){}
    virtual void init() = 0;
    virtual void update() = 0;

protected:
    Mech * m_mechanicalObject;
    Topo * m_sofa_topology; //<<Basic sofa topology usefull for a lot of thing
};

#endif // PBDBASECONSTRAINTDATA_HPP

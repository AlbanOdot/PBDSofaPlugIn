#ifndef PBDBASECONSTRAINTDATA_HPP
#define PBDBASECONSTRAINTDATA_HPP

#include <SofaBaseMechanics/MechanicalObject.h>

template< class T>
class PBDBaseConstraintData
{
public:
    typedef sofa::component::container::MechanicalObject< T > Mech;
    typedef sofa::core::topology::BaseMeshTopology  Topo;
    typedef sofa::helper::Quater<SReal> Quaternion;

    PBDBaseConstraintData(Mech * m = nullptr, Topo* t = nullptr) : m_mechanicalObject(m),m_sofa_topology(t){}

    /*
     * Create and init all of the data needed to solve a defined constraint.
     */
    virtual void init() = 0;

    /*
     * Reinit all of the data according to the current context
     */
    virtual void update() = 0;

protected:
    Mech * m_mechanicalObject;
    Topo * m_sofa_topology; //<<Basic sofa topology usefull for a lot of thing
};

#endif // PBDBASECONSTRAINTDATA_HPP

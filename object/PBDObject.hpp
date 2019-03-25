#ifndef PBDOBJECT_HPP
#define PBDOBJECT_HPP

#include <SofaBaseMechanics/MechanicalObject.h>

class PBDObject
{
    typedef sofa::defaulttype::Vec3Types::Coord       Coord;
    typedef sofa::helper::vector<Coord>               VecCoord;
    typedef sofa::core::objectmodel::Data<VecCoord>   Coordinates;
    typedef sofa::helper::ReadAccessor  <Coordinates> ReadCoord;
    typedef sofa::helper::WriteAccessor <Coordinates> WriteCoord;

    typedef sofa::defaulttype::Vec3Types::Deriv       Deriv;
    typedef sofa::helper::vector<Deriv>               VecDeriv;
    typedef sofa::core::objectmodel::Data<VecDeriv>   Derivatives;
    typedef sofa::helper::ReadAccessor  <Derivatives> ReadDeriv;
    typedef sofa::helper::WriteAccessor <Derivatives> WriteDeriv;

    typedef  std::vector<std::vector<std::pair<uint,SReal>>> Topology;

public:
    PBDObject(sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * mobj,
              sofa::core::topology::BaseMeshTopology * topo = nullptr);

    inline const ReadCoord & rest() const   {return m_rest[0];}
    inline WriteCoord position()            {return m_mechanicalObject->writePositions ();}
    inline WriteDeriv velocity()            {return m_mechanicalObject->writeVelocities ();}
    inline const Topology& topology() {return m_topology;}
    inline sofa::core::topology::BaseMeshTopology * sofaTopology() { return m_sofa_topology;}
    inline sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * object() {return m_mechanicalObject;}

    inline void setTopology(sofa::core::topology::BaseMeshTopology * topology);
    void optimizeTopology();


protected:
    sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * m_mechanicalObject;
    sofa::core::topology::BaseMeshTopology * m_sofa_topology; //<<Basic sofa topology usefull for a lot of thing
    Topology m_topology;//<<Topology optimized to apply constraints on vertex
    std::vector<ReadCoord> m_rest;

};
#endif

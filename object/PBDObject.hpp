#ifndef PBDOBJECT_HPP
#define PBDOBJECT_HPP

#include <SofaBaseMechanics/MechanicalObject.h>
#include "BeamElement.hpp"

//            indices of
//                 x0,x1         Q
typedef std::pair<uint[2],Eigen::Matrix4d> bendingStruct;


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

    typedef std::vector<std::vector<std::pair<uint,SReal>>> VertexTopology;
    typedef std::vector<std::vector<bendingStruct>>         BendingTopology;
    typedef std::vector<std::pair<float,Eigen::Matrix3d>> TetrahedronBasis; //<< Based described by the 4 points of a tetrahedron
    typedef std::vector<BeamElement> Beam;
private:

public:
    PBDObject(sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * mobj,
              sofa::core::topology::BaseMeshTopology * topo = nullptr);

    inline const SReal& invMass() const                                                             { return m_invVertexMass;}
    inline const SReal& mass() const                                                                { return m_mass;}
    inline const ReadCoord & rest() const                                                           {return m_rest[0];}
    inline WriteCoord position()                                                                    {return m_mechanicalObject->writePositions ();}
    inline WriteDeriv velocity()                                                                    {return m_mechanicalObject->writeVelocities ();}
    inline const VertexTopology& topology()                                                         {return m_stretch_topology;}
    inline  BendingTopology& bend_topology()                                                        {return m_bending_topology;}
    inline const std::vector<std::vector<std::pair<float,float>>>& areas()                          { return m_triangle_rest_area;}
    inline sofa::core::topology::BaseMeshTopology * sofaTopology()                                  { return m_sofa_topology;}
    inline sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * object()  {return m_mechanicalObject;}
    inline TetrahedronBasis& tetrahedraBases()                                                      {return m_tetra_bases;}
    inline Beam& beam()                                                                             {return m_beam;}
    inline std::vector<Eigen::Quaterniond>& freeOrientation()                                         {return m_freeOrientation;}
    inline std::vector<Eigen::Vector3d>& angularSpeed()                                             {return m_angularSpeed;}
    inline std::vector<Eigen::Vector3d>& torque()                                                   {return m_torque;}
    inline std::vector<Eigen::Matrix3d>& inertia()                                                  {return m_intertia;}

    inline void setTopology(sofa::core::topology::BaseMeshTopology * topology);
    void optimizeTopology();
    void computeStretchTopology();
    void computeBendingTopology();
    void computeTetrahedraBasis();
    void computeBeam();
    void computeQ(const sofa::defaulttype::Vec3 *x[4],Eigen::Matrix4d& Q, std::pair<float,float>& area);
    void computeStiffRod();
    void setupAngularVelocity();
    void applyFixedPoint(const std::vector<uint>& idx);




protected:
    sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * m_mechanicalObject;
    sofa::core::topology::BaseMeshTopology * m_sofa_topology; //<<Basic sofa topology usefull for a lot of thing
    VertexTopology m_stretch_topology;//<<Topology optimized to apply constraints on vertex
    BendingTopology m_bending_topology;
    TetrahedronBasis m_tetra_bases;
    std::vector<std::vector<std::pair<float,float>>> m_triangle_rest_area;
    std::vector<ReadCoord> m_rest;

    std::vector<Eigen::Quaterniond> m_freeOrientation;
    std::vector<Eigen::Vector3d> m_angularSpeed;
    std::vector<Eigen::Vector3d> m_torque;
    std::vector<Eigen::Matrix3d> m_intertia;
    Beam m_beam;

    SReal m_invVertexMass;
    SReal m_mass;
};
#endif

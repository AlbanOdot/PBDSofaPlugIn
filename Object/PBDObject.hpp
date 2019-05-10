#ifndef PBDOBJECT_HPP
#define PBDOBJECT_HPP

#include <SofaBaseMechanics/MechanicalObject.h>

#include "./Data/PBDVertexMass.hpp"

#include "./Data/Elastic/PBDVertexTopology.hpp"
#include "./Data/Elastic/PBDBendingTopology.hpp"

#include "./Data/FEM/BeamElement.hpp"
#include "./Data/FEM/PBDTetrahedronBasis.hpp"

namespace PBDDataType
{
#define STRETCH     1
#define BENDING     2
#define TETRAHEDRON 4
#define MASS        8
#define ORIENTED    16
#define BEAM        32
#define PBDMASS     64
};

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

    typedef std::vector<BeamElement> Beam;

public:
    PBDObject(sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * mobj,
              sofa::core::topology::BaseMeshTopology * topo = nullptr);

    inline       VertexMassData& mass()                                                                             {return m_mass.m();}
    inline       VertexMassData& invMass()                                                                          {return m_mass.w();}
    inline       SReal& mass(uint i)                                                                                {return m_mass.m(i);}
    inline       SReal& invMass(uint i)                                                                             {return m_mass.w(i);}
    inline const ReadCoord & rest() const                                                                           {return m_rest[0];}
    inline       WriteCoord position()                                                                              {return m_mechanicalObject->writePositions ();}
    inline       WriteDeriv velocity()                                                                              {return m_mechanicalObject->writeVelocities ();}
    inline       VertexTopologyData& topology()                                                                     {return m_stretch_topology.data();}
    inline       BendingTopologyData& bendTopology()                                                                {return m_bending_topology.bendingData ();}
    inline const TriangleAreaData& areas()                                                                          {return m_bending_topology.triangleAreaData();}
    inline       sofa::core::topology::BaseMeshTopology * sofaTopology()                                            {return m_sofa_topology;}
    inline       sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * object()            {return m_mechanicalObject;}
    inline       TetrahedronBasisData& tetrahedraBases()                                                            {return m_tetra_bases.data();}
    inline       Beam& beam()                                                                                       {return m_beam;}
    inline       std::vector<Eigen::Quaterniond>& freeOrientation()                                                 {return m_freeOrientation;}
    inline       std::vector<Eigen::Vector3d>& angularSpeed()                                                       {return m_angularSpeed;}
    inline       std::vector<Eigen::Vector3d>& torque()                                                             {return m_torque;}
    inline       std::vector<Eigen::Matrix3d>& inertia()                                                            {return m_intertia;}
    inline       unsigned long int dataType()                                                               const   {return m_dataType;}

    inline       void setTopology(sofa::core::topology::BaseMeshTopology * topology);
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
    sofa::core::topology::BaseMeshTopology *                                       m_sofa_topology; //<<Basic sofa topology usefull for a lot of thing

    VertexTopology m_stretch_topology;//<<Topology optimized to apply constraints on vertex
    BendingTopology m_bending_topology;
    TetrahedronBasis m_tetra_bases;
    std::vector<ReadCoord> m_rest;
    VertexMass m_mass;
    unsigned long int m_dataType;

    std::vector<Eigen::Quaterniond> m_freeOrientation;
    std::vector<Eigen::Vector3d> m_angularSpeed;
    std::vector<Eigen::Vector3d> m_torque;
    std::vector<Eigen::Matrix3d> m_intertia;
    Beam m_beam;

};
#endif

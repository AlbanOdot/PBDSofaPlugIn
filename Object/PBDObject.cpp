#include "PBDObject.hpp"
#include <Eigen/MatrixFunctions>
#include <sofa/core/objectmodel/BaseContext.h>
#include <SofaBaseMechanics/UniformMass.h>

PBDObject::PBDObject(sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * mobj,
                     sofa::core::topology::BaseMeshTopology * topo ): m_mechanicalObject(mobj), m_sofa_topology(topo)
{
    if(m_mechanicalObject && m_sofa_topology)
    {
        m_rest.emplace_back(m_mechanicalObject->readRestPositions ());
        m_mass = VertexMass(m_mechanicalObject,m_sofa_topology);
    }
    m_dataType = 0;
}

inline void PBDObject::setTopology(sofa::core::topology::BaseMeshTopology *topology)
{
    m_sofa_topology = topology;
    if(m_mechanicalObject && m_sofa_topology)
    {
        m_rest.emplace_back(m_mechanicalObject->readRestPositions ());
        m_mass = VertexMass(m_mechanicalObject,m_sofa_topology);
    }
    m_dataType = 0;
}


void PBDObject::computeStretchTopology()
{
    m_stretch_topology  = VertexTopology(m_mechanicalObject,m_sofa_topology);
    m_dataType |= STRETCH;
    m_integration_type |= NORMAL;
}

void PBDObject::computeBendingTopology()
{
    m_bending_topology = PBDBendingTopology(m_mechanicalObject,m_sofa_topology);
    m_dataType |= BENDING;
    m_integration_type |= NORMAL;
}

void PBDObject::computeTetrahedraBasis()
{
    m_tetra_bases = PBDTetrahedronBasis(m_mechanicalObject,m_sofa_topology);
    m_dataType |= TETRAHEDRON;
    m_integration_type |= NORMAL;
}

void PBDObject::computeOrientation()
{
    m_orientation = PBDOrientation(m_mechanicalObject,m_sofa_topology);
    m_dataType |= ORIENTED;
    m_integration_type |= NORMAL;
}
void PBDObject::computeElasticRod()
{
    m_elasticRod = PBDElasticRodData(m_mechanicalObject,m_sofa_topology);
    m_dataType |= ELASTICROD;
    m_integration_type |= ANGULAR;

}

void PBDObject::computeStiffRod()
{
    m_stiffRod = PBDStiffRodData(m_mechanicalObject,m_sofa_topology);
    m_dataType |= STIFFROD;
    m_integration_type |= NORMAL;
}

void PBDObject::setupAngularVelocity(const std::vector<Vector3r>& as)
{
    if(m_dataType & ORIENTED)
        m_orientation.setAngularVelocity(as);
}


void PBDObject::applyFixedPoint(const std::vector<uint>& idx)
{
    if(m_dataType & ELASTICROD)
        m_elasticRod.setToZero (idx);
}

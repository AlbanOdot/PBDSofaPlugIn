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

void PBDObject::optimizeTopology()
{
    computeStretchTopology ();
    computeBendingTopology ();
    computeTetrahedraBasis ();
}


void PBDObject::computeStretchTopology()
{
    m_stretch_topology  = VertexTopology(m_mechanicalObject,m_sofa_topology);
    m_dataType |= STRETCH;
}

void PBDObject::computeBendingTopology()
{
    m_bending_topology = PBDBendingTopology(m_mechanicalObject,m_sofa_topology);
    m_dataType |= BENDING;
}

void PBDObject::computeTetrahedraBasis()
{
    m_tetra_bases = PBDTetrahedronBasis(m_mechanicalObject,m_sofa_topology);
    m_dataType |= TETRAHEDRON;
}

void PBDObject::computeBeam()
{
    //Compute the vertice oriented topology
    const auto& edges = m_sofa_topology->getEdges ();
    const auto& rest = m_mechanicalObject->readRestPositions ();
    for(uint e = 0; e < edges.size(); ++e)
    {
        uint edge[2] = {edges[e][0],edges[e][1]};
        uint a = e+1 == edges.size () ? e : e+1;
        uint edgeN[2] = {edges[a][0],edges[a][1]};
        SReal l = 0.5* ( (rest[edge[0]] - rest[edge[1]]).norm () + (rest[edgeN[0]] - rest[edgeN[1]]).norm ());
        m_beam.emplace_back(BeamElement(m_mechanicalObject,edge,e,l));
    }
    //Post init of darboux vectors
    for(uint e = 0; e < edges.size() - 1; ++e)
    {
        m_beam[e].initRestDarboux(m_beam[e+1].m_q);
    }
    m_beam[m_beam.size() - 1].initRestDarboux (m_beam[m_beam.size() - 1].q ());
    m_beam.emplace_back(m_beam[m_beam.size() - 1]);
    auto& last = m_beam[m_beam.size() - 1];
    last.m_segmentIndicies.first = last.m_segmentIndicies.second;
    setupAngularVelocity ();
    m_dataType |= BEAM;
    m_dataType |= ORIENTED;
}

void PBDObject::computeStiffRod()
{
}

void PBDObject::setupAngularVelocity()
{
    const Eigen::Matrix3d id; id.Identity ();
    for(uint i = 0; i< m_beam.size (); ++i)
    {
        m_angularSpeed.emplace_back(Eigen::Vector3d(0,0,0));
        m_torque.emplace_back(Eigen::Vector3d(0,0,0));
        m_freeOrientation.emplace_back(m_beam[i].q());
        m_intertia.emplace_back(id);
    }
    for(uint i = 0; i< m_beam.size (); ++i)
    {
        m_intertia[i].setIdentity ();
    }
}


void PBDObject::applyFixedPoint(const std::vector<uint>& idx)
{
    for(auto& a : idx)
    {
        m_beam[a].setwq(0.0);
    }
}

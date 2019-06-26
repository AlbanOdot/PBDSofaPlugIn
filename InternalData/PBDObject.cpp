#include "PBDObject.hpp"
#include <Eigen/MatrixFunctions>
#include <sofa/core/objectmodel/BaseContext.h>
#include <SofaBaseMechanics/UniformMass.h>
template < class T >
PBDObject<T>::PBDObject(sofa::component::container::MechanicalObject< T > * mobj,
                     sofa::core::topology::BaseMeshTopology * topo ): m_mechanicalObject(mobj), m_sofa_topology(topo)
{
    if(m_mechanicalObject && m_sofa_topology)
    {
        m_rest.emplace_back(m_mechanicalObject->readRestPositions ());
        m_mass = PBDVertexMass<T>(m_mechanicalObject,m_sofa_topology);
    }
    m_dataType = 0;
    m_integration_type = 0;
}

template < class T >
inline void PBDObject<T>::setTopology(sofa::core::topology::BaseMeshTopology *topology)
{
    m_sofa_topology = topology;
    if(m_mechanicalObject && m_sofa_topology)
    {
        m_rest.emplace_back(m_mechanicalObject->readRestPositions ());
        m_mass = PBDVertexMass<T>(m_mechanicalObject,m_sofa_topology);
    }
    m_dataType = 0;
}

template < class T >
void PBDObject<T>::computeStretchTopology()
{
    m_stretch_topology  = PBDVertexTopology<T>(m_mechanicalObject,m_sofa_topology);
    m_dataType |= STRETCH;
    if(m_integration_type == 0)
        m_integration_type |= NORMAL;
}

template < >
void PBDObject<sofa::defaulttype::Vec3Types>::computeBendingTopology()
{
    m_bending_topology = PBDBendingTopology(m_mechanicalObject,m_sofa_topology);
    m_dataType |= BENDING;
    if(m_integration_type == 0)
        m_integration_type |= NORMAL;
}

template <>
void PBDObject<sofa::defaulttype::Vec3Types>::computeTetrahedraBasis()
{
    m_tetra_bases = PBDTetrahedronBasis(m_mechanicalObject,m_sofa_topology);
    m_dataType |= TETRAHEDRON;
    if(m_integration_type == 0)
        m_integration_type |= NORMAL;
}

template < >
void PBDObject<sofa::defaulttype::RigidTypes>::computeOrientation()
{
    m_orientation = PBDOrientation(m_mechanicalObject,m_sofa_topology);
    m_dataType |= ORIENTED;
    if(m_integration_type == 0)
        m_integration_type |= NORMAL;
}

template <>
void PBDObject<sofa::defaulttype::RigidTypes>::computeElasticRod()
{
    m_elasticRod = PBDElasticRodData(m_mechanicalObject,m_sofa_topology);
    m_dataType |= ELASTICROD;
    if(m_integration_type == 0 || m_integration_type == NORMAL)
        m_integration_type |= ANGULAR;

}

template <  >
void PBDObject<sofa::defaulttype::RigidTypes>::computeStiffRod()
{
    m_stiffRod = PBDStiffRodData(m_mechanicalObject,m_sofa_topology);
    m_dataType |= STIFFROD;
    if(m_integration_type == 0)
        m_integration_type |= NORMAL;
}


template < >
void PBDObject<sofa::defaulttype::RigidTypes>::computeCosseratRod()
{
    m_PD_CosseratRod = PDCosseratRodData(m_mechanicalObject,m_sofa_topology);
    m_dataType  |= COSSERATROD;
    if(m_integration_type == 0 || m_integration_type == NORMAL)
        m_integration_type |= ANGULAR;
}

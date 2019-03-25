#include "PBDObject.hpp"

PBDObject::PBDObject(sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * mobj,
          sofa::core::topology::BaseMeshTopology * topo ): m_mechanicalObject(mobj), m_sofa_topology(topo)
{
    if(m_mechanicalObject && m_sofa_topology)
    {
        m_rest.emplace_back(m_mechanicalObject->readRestPositions ());
        optimizeTopology ();
    }
}

inline void PBDObject::setTopology(sofa::core::topology::BaseMeshTopology *topology)
{
    m_sofa_topology = topology;
    if(m_mechanicalObject && m_sofa_topology)
    {
        m_rest.emplace_back(m_mechanicalObject->readRestPositions ());
        optimizeTopology ();
    }
}

void PBDObject::optimizeTopology()
{
    //We assume 1 topology per object and vice versa
    for(uint i = 0; i < m_rest[0].size(); ++i)
    {
        //Get the neighbors of point I
        const auto& neighbors = m_sofa_topology->getVerticesAroundVertex (i);
        std::vector<std::pair<uint,SReal>> neighborhood;
        for(uint j = 0; j < neighbors.size(); ++j)
        {
            if( neighbors[j] > i )//Unidirectionnal neighborhood
            {
                SReal d = (m_rest[0][i] - m_rest[0][neighbors[j]]).norm();
                neighborhood.emplace_back(std::pair<uint,SReal>(neighbors[j],d));
            }
        }
        m_topology.emplace_back(neighborhood);
    }
}

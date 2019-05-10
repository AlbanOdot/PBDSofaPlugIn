#include "PBDVertexTopology.hpp"

PBDVertexTopology::PBDVertexTopology(Mech * m, Topo * t) : PBDBaseConstraintData (m,t)
{
    if( m && t )
        init ();
}

void PBDVertexTopology::init()
{
    const auto& rest = m_mechanicalObject->readRestPositions ();
    //Compute the vertice oriented topology
    for(uint i = 0; i < rest.size(); ++i)
    {
        //Get the neighbors of point I
        const auto& neighbors = m_sofa_topology->getVerticesAroundVertex (i);
        std::vector<std::pair<uint,SReal>> neighborhood;
        for(uint j = 0; j < neighbors.size(); ++j)
        {
            if( neighbors[j] < i )//Unidirectionnal neighborhood
            {
                SReal d = (rest[i] - rest[neighbors[j]]).norm();
                neighborhood.emplace_back(std::pair<uint,SReal>(neighbors[j],d));
            }
        }
        m_data.emplace_back(neighborhood);
    }
}

void PBDVertexTopology::update()
{
    m_data.clear ();
    if( m_mechanicalObject && m_sofa_topology )
        init ();
}

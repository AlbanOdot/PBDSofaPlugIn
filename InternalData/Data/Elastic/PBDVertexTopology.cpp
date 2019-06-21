#include "PBDVertexTopology.hpp"

template <class T>
PBDVertexTopology<T>::PBDVertexTopology(Mech * m, Topo * t) : PBDBaseConstraintData<T> (m,t)
{
    if( m && t )
        init ();
}

template <class T>
void PBDVertexTopology<T>::init()
{
    const auto& rest = PBDBaseConstraintData<T>::m_mechanicalObject->readRestPositions ();
    //Compute the vertice oriented topology
    for(uint i = 0; i < rest.size(); ++i)
    {
        //Get the neighbors of point I
        const auto& neighbors = PBDBaseConstraintData<T>::m_sofa_topology->getVerticesAroundVertex (i);
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


template <class T>
void PBDVertexTopology<T>::update()
{
    m_data.clear ();
    if( PBDBaseConstraintData<T>::m_mechanicalObject && PBDBaseConstraintData<T>::m_sofa_topology )
        init ();
}

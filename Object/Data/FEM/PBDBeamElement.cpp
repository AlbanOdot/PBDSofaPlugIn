#include "./PBDBeamElement.hpp"

PBDBeamElement::PBDBeamElement(Mech * m, Topo * t) : PBDBaseConstraintData (m,t)
{
    if( m && t )
        init ();
}

void PBDBeamElement::init()
{
    const auto& rest = m_mechanicalObject->readRestPositions ();
    const auto& edges = m_sofa_topology->getEdges ();
    for(uint i = 0; i < edges.size (); ++i)
    {
        uint a = i+1 == edges.size () ? i : i+1;
        infos inf;
        inf.info[0] = edges[i][0];
        inf.info[1] = edges[i][1];
        inf.info[2] = i;
        m_indices.emplace_back(inf);
        m_averageLength.emplace_back(0.5 * ( (rest[edges[i][0]]-rest[edges[i][1]]).norm () + (rest[edges[a][0]]-rest[edges[a][1]]).norm ()));
    }
    m_indices.emplace_back(m_indices[m_indices.size () - 1]);
    m_averageLength.emplace_back(m_averageLength[m_averageLength.size () - 1]);

}

void PBDBeamElement::update()
{

    m_indices.clear ();
    m_averageLength.clear ();
    if( m_mechanicalObject && m_sofa_topology )
        init ();
}

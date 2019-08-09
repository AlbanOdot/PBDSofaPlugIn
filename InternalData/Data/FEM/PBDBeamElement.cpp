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
    bool color = false;
    for(uint i = 0; i < edges.size () - 1; ++i)
    {
        uint a = i+1 == edges.size () ? i : i+1;
        m_averageLength.emplace_back(0.5 * ( (rest[edges[i][0]].getCenter ()-rest[edges[i][1]].getCenter ()).norm () + (rest[edges[a][0]].getCenter ()-rest[edges[a][1]].getCenter ()).norm ()));
        m_color.emplace_back(color);
        color = !color;
    }
    m_color.emplace_back(color);
    m_averageLength.emplace_back(m_averageLength[m_averageLength.size () - 1]);
}

void PBDBeamElement::update()
{

    m_averageLength.clear ();
    m_color.clear ();
    if( m_mechanicalObject && m_sofa_topology )
        init ();
}

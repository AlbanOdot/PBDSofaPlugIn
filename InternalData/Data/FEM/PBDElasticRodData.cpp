#include "./PBDElasticRodData.hpp"

PBDElasticRodData::PBDElasticRodData(Mech * m, Topo * t) : PBDBeamElement(m,t)
{
    if( m && t )
        init ();
}

void PBDElasticRodData::init()
{

    const auto& edges = m_sofa_topology->getEdges ();
    for(const auto& r : edges)
    {
        m_wq.emplace_back(1.0);
    }
}

void PBDElasticRodData::update()
{

    m_wq.clear ();
    if( m_mechanicalObject && m_sofa_topology )
        init ();
}

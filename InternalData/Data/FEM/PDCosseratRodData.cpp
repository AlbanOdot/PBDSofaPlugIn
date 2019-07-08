#include "./PDCosseratRodData.hpp"

void PDCosseratRodData::setupW(const SReal E, const SReal nu, const SReal r)
{
    m_wbt.clear ();
    m_ws.clear ();
    SReal G = E / (2.0 * (1 + nu));
    SReal ws = E * M_PI * r * r;
    SReal wbt = 2.0 * M_PI * std::pow(r,4) * G;
    for(auto & l : m_averageLength)
    {
        m_wbt.emplace_back(0.5 * wbt / l);
        m_ws.emplace_back(ws*l);
    }
}

void PDCosseratRodData::update()
{

    m_wbt.clear ();
    m_ws.clear ();
    if( m_mechanicalObject && m_sofa_topology )
        init ();
}

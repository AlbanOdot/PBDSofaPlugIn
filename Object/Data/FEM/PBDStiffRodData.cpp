#include "./PBDStiffRodData.hpp"

PBDStiffRodData::PBDStiffRodData(Mech * m, Topo * t) : PBDBeamElement(m,t)
{
    if( m && t )
        init ();
}

void PBDStiffRodData::init()
{
    const auto& rest = m_mechanicalObject->readRestPositions ();
    for(const auto& r : rest)
    {
        m_invMass.emplace_back(Mat6());
        m_invMass[m_invMass.size() - 1].setZero();
        m_lambda.emplace_back(Vec6());
        m_lambda[m_lambda.size() - 1].setZero();
    }

}

void PBDStiffRodData::update()
{

    m_invMass.clear();
    m_lambda.clear();
    if( m_mechanicalObject && m_sofa_topology )
        init ();
}

void PBDStiffRodData::setInvMass(std::vector<Vec6> invMassDiag)
{
    if(invMassDiag.size () ==  m_invMass.size ())
    {
        for(uint i = 0; i < invMassDiag.size (); ++i){
            m_invMass[i] = invMassDiag[i].asDiagonal ();
        }
    }
    else
    {
        for(auto& m : m_invMass)
            m += invMassDiag[0].asDiagonal();
    }

}

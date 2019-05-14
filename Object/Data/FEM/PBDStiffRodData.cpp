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
        m_invMass.emplace_back(Vec6());
        m_invMass[m_invMass.size() - 1].setZero();
        m_invAlpha.emplace_back(Vec6());
        m_invAlpha[m_invAlpha.size() - 1].setZero();
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

void PBDStiffRodData::setInvMass(const std::vector<Vec6> &invMassDiag)
{
    if(invMassDiag.size () ==  m_invMass.size ())
    {
        for(uint i = 0; i < invMassDiag.size (); ++i){
            m_invMass[i] = invMassDiag[i];
        }
    }
    else if( invMassDiag.size () == 1)
    {
        for(auto& m : m_invMass)
            m = invMassDiag[0];
    }
}

void PBDStiffRodData::setInvAlpha(const std::vector<Vec6> &invalpha)
{
    if(invalpha.size () ==  m_invMass.size ())
    {
        for(uint i = 0; i < invalpha.size (); ++i){
            m_invAlpha[i] = invalpha[i];
        }
    }
    else if( invalpha.size () == 1)
    {
        for(auto& m : m_invAlpha)
            m = invalpha[0];
    }
}

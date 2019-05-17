#include "./PBDStiffRodData.hpp"

PBDStiffRodData::PBDStiffRodData(Mech * m, Topo * t) : PBDBeamElement(m,t)
{
    if( m && t )
        init ();
}

void PBDStiffRodData::init()
{
    const auto& rest = m_mechanicalObject->readRestPositions ();
    Vector3r zero;
    zero.setZero ();
    for(const auto& r : rest)
    {
        m_invInertia.emplace_back(zero);
        m_invAlpha.emplace_back(zero);
        m_lambdar.emplace_back(zero);
        m_lambdabt.emplace_back(zero);
        m_stiffness.emplace_back(0.0);
    }

}

void PBDStiffRodData::update()
{

    m_invInertia.clear();
    m_lambdar.clear();
    m_invAlpha.clear();
    m_lambdabt.clear();
    m_stiffness.clear();
    if( m_mechanicalObject && m_sofa_topology )
        init ();
}

void PBDStiffRodData::setInvInertia(const std::vector<Vector3r> &invInertiaDiag)
{
    if(invInertiaDiag.size () ==  m_invInertia.size ())
    {
        for(uint i = 0; i < invInertiaDiag.size (); ++i){
            m_invInertia[i] = invInertiaDiag[i];
        }
    }
    else if( invInertiaDiag.size () == 1)
    {
        for(auto& m : m_invInertia)
            m = invInertiaDiag[0];
    }
}

void PBDStiffRodData::setInvAlpha(const std::vector<Vector3r> &invalpha)
{
    if(invalpha.size () ==  m_invAlpha.size ())
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

void PBDStiffRodData::setStiffness(const std::vector<SReal> &stiffness)
{
    if(stiffness.size () ==  m_stiffness.size ())
    {
        for(uint i = 0; i < stiffness.size (); ++i){
            m_stiffness[i] = stiffness[i];
        }
    }
    else if( stiffness.size () == 1)
    {
        for(auto& m : m_stiffness)
            m = stiffness[0];
    }
}

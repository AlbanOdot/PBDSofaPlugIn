#include "./PBDStiffRodData.hpp"

PBDStiffRodData::PBDStiffRodData(Mech * m, Topo * t) : PBDBeamElement(m,t)
{
    if( m && t )
        init ();
}

void PBDStiffRodData::init()
{
    const auto& rest = m_mechanicalObject->readRestPositions ();
    Vector6r zero;
    zero.setZero ();
    for(const auto& r : rest)
    {
        m_massMatrix.emplace_back(zero);
        m_alpha.emplace_back(zero);
        m_lambda.emplace_back(zero);
    }

}

void PBDStiffRodData::update()
{

    m_massMatrix.clear();
    m_alpha.clear();
    m_lambda.clear();
    if( m_mechanicalObject && m_sofa_topology )
        init ();
}

void PBDStiffRodData::setMassMatrix(const std::vector<Vector6r> &massMatrixDiag)
{
    if(massMatrixDiag.size () ==  m_massMatrix.size ())
    {
        for(uint i = 0; i < massMatrixDiag.size (); ++i){
            m_massMatrix[i] = massMatrixDiag[i];
        }
    }
    else if( massMatrixDiag.size () == 1)
    {
        for(auto& m : m_massMatrix)
            m = massMatrixDiag[0];
    }
}

void PBDStiffRodData::setAlpha(const std::vector<Vector6r> &alpha)
{
    if(alpha.size () ==  m_alpha.size ())
    {
        for(uint i = 0; i < alpha.size (); ++i){
            m_alpha[i] = alpha[i];
        }
    }
    else if( alpha.size () == 1)
    {
        for(auto& m : m_alpha)
            m = alpha[0];
    }
}

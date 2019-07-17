#include "PBDTetrahedronBasis.hpp"
#include <Eigen/MatrixFunctions>

PBDTetrahedronBasis::PBDTetrahedronBasis(Mech * m, Topo * t) : PBDBaseConstraintData (m,t)
{
    if( m && t )
        init ();
}

void PBDTetrahedronBasis::init()
{
    static SReal one_over_6 = 1.0/6.0;

    const auto& rest = m_mechanicalObject->readRestPositions ();
    const auto& tetrahedra = m_sofa_topology->getTetrahedra ();
    for(uint i = 0; i < tetrahedra.size(); ++i)
    {
        const auto& t = tetrahedra[i];
        const auto& p3 = rest[t[3]];
        Matrix3 Dm;
        Dm.x() = rest[t[0]] - p3;
        Dm.y() = rest[t[1]] - p3;
        Dm.z() = rest[t[2]] - p3;
        SReal volume = one_over_6 * dot(Dm.x(),Dm.y().cross(Dm.z()));
        m_data.emplace_back(std::pair<float,Matrix3>(volume,Dm.inverted ()));

    }

}

void PBDTetrahedronBasis::update()
{
    m_data.clear ();
    if( m_mechanicalObject && m_sofa_topology )
        init ();
}

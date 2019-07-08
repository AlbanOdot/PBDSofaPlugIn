#include "PBDTetrahedronBasis.hpp"
#include <Eigen/MatrixFunctions>

PBDTetrahedronBasis::PBDTetrahedronBasis(Mech * m, Topo * t) : PBDBaseConstraintData (m,t)
{
    if( m && t )
        init ();
}

void PBDTetrahedronBasis::init()
{

    const auto& rest = m_mechanicalObject->readRestPositions ();
    const auto& tetrahedra = m_sofa_topology->getTetrahedra ();
    for(uint i = 0; i < tetrahedra.size(); ++i)
    {
        const auto& x = tetrahedra[i];
        const auto& p3 = rest[x[3]];
        Matrix3 Dm;
        Dm.x() = rest[x[0]] - p3;
        Dm.y() = rest[x[1]] - p3;
        Dm.z() = rest[x[2]] - p3;
        SReal volume = 0.166666666666666666 * dot(Dm.x(),Dm.y().cross(Dm.z()));
        m_data.emplace_back(std::pair<float,Matrix3>(volume,Dm.inverted ()));

    }

}

void PBDTetrahedronBasis::update()
{
    m_data.clear ();
    if( m_mechanicalObject && m_sofa_topology )
        init ();
}

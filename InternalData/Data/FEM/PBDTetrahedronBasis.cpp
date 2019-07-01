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
    m_data.resize (m_sofa_topology->getNbTetrahedra ());
    for(uint i = 0; i < tetrahedra.size(); ++i)
    {
        const auto& x = tetrahedra[i];
        const auto& p0 = rest[x[0]];
        const auto& p1 = rest[x[1]];
        const auto& p2 = rest[x[2]];
        const auto& p3 = rest[x[3]];
        Matrix3 Dm;
        Dm.x() = p0 - p3;
        Dm.y() = p1 - p3;
        Dm.z() = p2 - p3;

        m_data[i] = std::pair<float,Matrix3>(0.5*determinant(Dm),Dm.inverted ());
    }
}

void PBDTetrahedronBasis::update()
{
    m_data.clear ();
    if( m_mechanicalObject && m_sofa_topology )
        init ();
}

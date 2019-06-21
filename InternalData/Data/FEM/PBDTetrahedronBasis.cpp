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
        const auto& r1 = rest[x[0]] - rest[x[3]];
        const auto& r2 = rest[x[1]] - rest[x[3]];
        const auto& r3 = rest[x[2]] - rest[x[3]];
        Eigen::Matrix3d Dm; Dm << r1[0],r1[1],r1[2],
                                  r2[0],r2[1],r2[2],
                                  r3[0],r3[1],r3[2];
        //We take 0.5 since the tetrahedron volume is half the det
        m_data[i] = std::pair<float,Eigen::Matrix3d>(0.5*Dm.determinant(),Dm.inverse());
    }
}

void PBDTetrahedronBasis::update()
{
    m_data.clear ();
    if( m_mechanicalObject && m_sofa_topology )
        init ();
}

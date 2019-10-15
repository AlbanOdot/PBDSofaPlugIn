#include "PBDTriDynamicBasis.hpp"
#include <unsupported/Eigen/MatrixFunctions>

PBDTriDynamicBasis::PBDTriDynamicBasis(Mech * m, Topo * t) : PBDBaseConstraintData (m,t)
{
    if( m && t )
        init ();
}

void PBDTriDynamicBasis::init()
{
    const auto& rest = m_mechanicalObject->readRestPositions ();
    const auto& triangles = m_sofa_topology->getTriangles ();
    for(uint i = 0; i < triangles.size(); ++i)
    {
        const auto& t = triangles[i];

        SReal a = rest[t[1]][0] - rest[t[0]][0]; SReal b = rest[t[2]][0] - rest[t[0]][0];
        SReal c = rest[t[1]][1] - rest[t[0]][1]; SReal d = rest[t[2]][1] - rest[t[0]][1];

        // inverse
        SReal det = a*d - b*c;

        SReal s = static_cast<SReal>(1.0) / det;
        Matrix2 Pinv;
        Pinv(0,0) =  d*s; Pinv(0,1) = -b*s;
        Pinv(1,0) = -c*s;  Pinv(1,1) =  a*s;
        m_data.emplace_back(std::pair<float,Matrix2>(det/2.0,Pinv));
    }
}

void PBDTriDynamicBasis::update()
{
    m_data.clear ();
    if( m_mechanicalObject && m_sofa_topology )
        init ();
}



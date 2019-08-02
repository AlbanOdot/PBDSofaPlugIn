#include "PBDTriBasis.hpp"
#include <Eigen/MatrixFunctions>

PBDTriBasis::PBDTriBasis(Mech * m, Topo * t) : PBDBaseConstraintData (m,t)
{
    if( m && t )
        init ();
}

void PBDTriBasis::init()
{
    const auto& rest = m_mechanicalObject->readRestPositions ();
    const auto& triangles = m_sofa_topology->getTriangles ();
    for(uint i = 0; i < triangles.size(); ++i)
    {
        const auto& t = triangles[i];

        Vec3 normal0 = (rest[t[1]] - rest[t[0]]).cross(rest[t[2]] - rest[t[0]]);
        SReal area = normal0.norm() * static_cast<SReal>(0.5);

        Vec3 axis0_1 = rest[t[1]] - rest[t[0]];
        axis0_1.normalize();
        Vec3 axis0_2 = normal0.cross(axis0_1);
        axis0_2.normalize();

        Vec2 p[3];
        p[0] = Vec2(dot(rest[t[0]],axis0_2), dot(rest[t[0]],axis0_1));
        p[1] = Vec2(dot(rest[t[1]],axis0_2), dot(rest[t[1]],axis0_1));
        p[2] = Vec2(dot(rest[t[2]],axis0_2), dot(rest[t[2]],axis0_1));

        Matrix2 P;
        P(0, 0) = p[0][0] - p[2][0];
        P(1, 0) = p[0][1] - p[2][1];
        P(0, 1) = p[1][0] - p[2][0];
        P(1, 1) = p[1][1] - p[2][1];
        m_data.emplace_back(std::pair<float,Matrix2>(area,P.inverted()));
    }
}

void PBDTriBasis::update()
{
    m_data.clear ();
    if( m_mechanicalObject && m_sofa_topology )
        init ();
}

//        SReal a = rest[t[1]][0] - rest[t[0]][0]; SReal b = rest[t[2]][0] - rest[t[0]][0];
//        SReal c = rest[t[1]][1] - rest[t[0]][1]; SReal d = rest[t[2]][1] - rest[t[0]][1];
//        SReal det = a*d - b*c;
//        Matrix2 invRestMat;
//        SReal s = 1.0 / det;
//        invRestMat(0,0) =  d*s;  invRestMat(0,1) = -b*s;
//        invRestMat(1,0) = -c*s;  invRestMat(1,1) =  a*s;
//        m_data.emplace_back(std::pair<float,Matrix2>(det,invRestMat));

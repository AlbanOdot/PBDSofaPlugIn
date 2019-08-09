#include "PBDBendingTopology.hpp"
#include "../../../Common/MathFunctions.hpp"

PBDBendingTopology::PBDBendingTopology(Mech * m, Topo * t) : PBDBaseConstraintData (m,t), m_triangle_rest_area(m,t)
{
    if( m && t )
        init ();
}

void PBDBendingTopology::init()
{
    const auto& rest = m_mechanicalObject->readRestPositions ();
    //Compute the vertice oriented topology
    const auto& triangles = m_sofa_topology->getTriangles ();
    m_bending_topology.resize(m_sofa_topology->getNbEdges ());
    std::vector<std::vector<uint>> topology;
    m_triangle_rest_area.initTopology (topology);

    for(uint i = 0; i < topology.size(); ++i)
    {
        const auto& voisins = topology[i];
        const sofa::defaulttype::Vec3 * x[4];
        x[0] = &(rest[i]);
        for( const auto& voisin : voisins)
        {
            if( voisin < i )//Unidirectionnal neighborhood
            {
                //BENDING TOPOLOGY
                x[1] = &(rest[voisin]);
                uint edge_ID = m_sofa_topology->getEdgeIndex(i,voisin);
                const auto& triangles_ID = m_sofa_topology->getTrianglesAroundEdge(edge_ID);
                for(uint t = 0; t < triangles_ID.size() - 1; ++t)
                {
                    std::pair<uint[2],std::valarray<SReal>> bs;

                    //Third point of the first triangle
                    const auto& t0 = triangles[triangles_ID[t]];
                    bs.first[0] = t0[0] == i || t0[0] == voisin ?
                                      (t0[1] == i || t0[1] == voisin ? t0[2]: t0[1] ): t0[0];
                    x[2] = &(rest[bs.first[0]]);

                    //Third point of the second triangle
                    const auto& t1 = triangles[triangles_ID[t+1]];
                    bs.first[1] = t1[0] == i || t1[0] == voisin ?
                                      (t1[1] == i || t1[1] == voisin ? t1[2]: t1[1] ): t1[0];
                    x[3] = &(rest[bs.first[1]]);
                    computeQ(x,bs.second);
                    m_bending_topology[edge_ID].emplace_back(bs);
                }
            }
        }
    }
}

void PBDBendingTopology::update()
{
    m_bending_topology.clear ();
    m_triangle_rest_area.update ();
    if( m_mechanicalObject && m_sofa_topology )
        init ();
}

void PBDBendingTopology::computeQ(const sofa::defaulttype::Vec3 *x[], std::valarray<SReal> &Q)
{

    const auto& cot_angle = [](const sofa::defaulttype::Vec3& v, const sofa::defaulttype::Vec3& w) -> SReal
    {
        const Real cosTheta = dot(v,w);
        const Real sinTheta = (v.cross(w)).norm();
        return (cosTheta / sinTheta);
    };

    const auto& x0 = *x[0];
    const auto& x1 = *x[1];
    const auto& x2 = *x[2];
    const auto& x3 = *x[3];

    //This is e_s in the paper
    sofa::defaulttype::Vec3 e[5] = {(x1-x0).normalized (),
                                    (x2-x1).normalized (),
                                    (x0-x2).normalized (),
                                    (x3-x0).normalized (),
                                    (x1-x3).normalized ()};

    //zero padding to keep up with the notation of the paper
    SReal c01 = cot_angle(e[0],e[1]);
    SReal c02 = cot_angle(e[0],e[2]);
    SReal c03 = cot_angle(e[0],e[3]);
    SReal c04 = cot_angle(e[0],e[4]);

    sofa::defaulttype::Vec4 K( c01+c04,c02+c03,-c01-c02,-c03-c04);

    //Q is a symetric matrix (Looks a lot like the Loop decimation algorithm matrix ((Just sayin')) )
    Q = {K[0] * K[0], K[0] * K[1], K[0] * K[2], K[0] * K[3],
         K[1] * K[1], K[1] * K[2], K[1] * K[3],
         K[2] * K[2], K[2] * K[3],
         K[3] * K[3]};

    const Real A0 = static_cast<Real>(0.5) * (e[0].cross(e[1])).norm();
    const Real A1 = static_cast<Real>(0.5) * (e[0].cross(e[2])).norm();
    Q *= (3.0/(A1 + A0));

}


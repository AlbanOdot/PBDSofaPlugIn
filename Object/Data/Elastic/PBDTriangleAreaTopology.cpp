#include "PBDTriangleAreaTopology.hpp"

PBDTriangleAreaTopology::PBDTriangleAreaTopology(Mech * m, Topo * t) : PBDBaseConstraintData (m,t)
{
    if( m && t )
        init ();
}

void PBDTriangleAreaTopology::initTopology(std::vector<std::vector<uint>>& t)
{
    const auto& rest = m_mechanicalObject->readRestPositions ();
    //Compute the vertice oriented topology
    for(uint i = 0; i < rest.size(); ++i)
    {
        //Get the neighbors of point I
        const auto& neighbors = m_sofa_topology->getVerticesAroundVertex (i);
        std::vector<uint> neighborhood;
        for(const auto& n : neighbors)
        {
            if( n < i )//Unidirectionnal neighborhood
            {
                neighborhood.emplace_back(n);
            }
        }
        t.emplace_back(neighborhood);
    }
}

void PBDTriangleAreaTopology::init()
{
    const auto& area_triangle = [] (const sofa::defaulttype::Vec3& A,const sofa::defaulttype::Vec3& B, const sofa::defaulttype::Vec3& C ) -> SReal
    {
        return 0.5 * ((B-A).cross(C-A)).norm();
    };
    const auto& rest = m_mechanicalObject->readRestPositions ();
    std::vector<std::vector<uint>> topology;
    initTopology(topology);
    //Compute the vertice oriented topology
    const auto& triangles = m_sofa_topology->getTriangles ();
    m_data.resize(m_sofa_topology->getNbEdges ());
    for(uint i = 0; i < topology.size (); ++i)
    {
        const auto& voisins = topology[i];
        for(const auto& voisin : voisins)
        {
            if( voisin < i )//Unidirectionnal neighborhood
            {
                //TriangleArea TOPOLOGY
                uint edge_ID = m_sofa_topology->getEdgeIndex(i,voisin);
                const auto& triangles_ID = m_sofa_topology->getTrianglesAroundEdge(edge_ID);
                for(uint t = 0; t < triangles_ID.size() - 1; ++t)
                {
                    std::pair<uint[2],Eigen::Matrix4d> bs;

                    //Third point of the first triangle
                    const auto& t0 = triangles[triangles_ID[t]];
                    uint a = t0[0] == i || t0[0] == voisin ?
                                      (t0[1] == i || t0[1] == voisin ? t0[2]: t0[1] ): t0[0];

                    //Third point of the second triangle
                    const auto& t1 = triangles[triangles_ID[t+1]];
                    uint b = t1[0] == i || t1[0] == voisin ?
                                      (t1[1] == i || t1[1] == voisin ? t1[2]: t1[1] ): t1[0];

                    std::pair<float,float> area;
                    area.first  = area_triangle(rest[i],rest[voisin],rest[a]);
                    area.second = area_triangle(rest[i],rest[voisin],rest[b]);
                    m_data[edge_ID].emplace_back(area);
                }
            }
        }
    }
}

void PBDTriangleAreaTopology::update()
{
    m_data.clear ();
    if( m_mechanicalObject && m_sofa_topology )
        init ();
}


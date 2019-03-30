#include "PBDObject.hpp"
PBDObject::PBDObject(sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * mobj,
                     sofa::core::topology::BaseMeshTopology * topo ): m_mechanicalObject(mobj), m_sofa_topology(topo)
{
    if(m_mechanicalObject && m_sofa_topology)
    {
        m_rest.emplace_back(m_mechanicalObject->readRestPositions ());
        optimizeTopology ();
    }
}

inline void PBDObject::setTopology(sofa::core::topology::BaseMeshTopology *topology)
{
    m_sofa_topology = topology;
    if(m_mechanicalObject && m_sofa_topology)
    {
        m_rest.emplace_back(m_mechanicalObject->readRestPositions ());
        optimizeTopology ();
    }
}

void PBDObject::optimizeTopology()
{

    //Compute the vertice oriented topology
    for(uint i = 0; i < m_rest[0].size(); ++i)
    {
        //Get the neighbors of point I
        const auto& neighbors = m_sofa_topology->getVerticesAroundVertex (i);
        std::vector<std::pair<uint,SReal>> neighborhood;
        std::cout << "["<<i<<"] :-> {";
        for(uint j = 0; j < neighbors.size(); ++j)
        {
            if( neighbors[j] < i )//Unidirectionnal neighborhood
            {
                std::cout << neighbors[j] <<" ";
                SReal d = (m_rest[0][i] - m_rest[0][neighbors[j]]).norm();
                neighborhood.emplace_back(std::pair<uint,SReal>(neighbors[j],d));
            }
        }
        std::cout << "}"<<std::endl;
        m_topology.emplace_back(neighborhood);
    }


    //Compute the vertice oriented topology
    const auto& triangles = m_sofa_topology->getTriangles ();
    m_bending_topology.resize(m_sofa_topology->getNbEdges ());
    m_triangle_rest_area.resize (m_sofa_topology->getNbEdges ());
    for(uint i = 0; i < m_rest[0].size(); ++i)
    {
        const auto& voisins = m_topology[i];

        const sofa::defaulttype::Vec3 * x[4];
        x[0] = &(m_rest[0][i]);
        for( const auto& voisin : voisins)
        {
            //BENDING TOPOLOGY
            uint edge_ID = m_sofa_topology->getEdgeIndex(i,voisin.first);
            const auto& triangles_ID = m_sofa_topology->getTrianglesAroundEdge(edge_ID);
            x[1] = &(m_rest[0][voisin.first]);
            for(uint t = 0; t < triangles_ID.size() - 1; ++t)
            {
                bendingStruct bs;

                //Third point of the first triangle
                const auto& t0 = triangles[triangles_ID[t]];
                bs.first[0] = t0[0] == i || t0[0] == voisin.first ?
                             (t0[1] == i || t0[1] == voisin.first ? t0[2]: t0[1] ): t0[0];
                x[2] = &(m_rest[0][bs.first[0]]);

                //Third point of the second triangle
                const auto& t1 = triangles[triangles_ID[t+1]];
                bs.first[1] = t1[0] == i || t1[0] == voisin.first ?
                             (t1[1] == i || t1[1] == voisin.first ? t1[2]: t1[1] ): t1[0];
                x[3] = &(m_rest[0][bs.first[1]]);

                std::pair<float,float> area;
                computeQ(x,bs.second,area);
                m_triangle_rest_area[edge_ID].emplace_back(area);
                m_bending_topology[edge_ID].emplace_back(bs);
            }
        }
    }
}


void PBDObject::computeQ(const sofa::defaulttype::Vec3 *x[], Eigen::Matrix4d &Q, std::pair<float,float>& area)
{

    const auto& area_triangle = [] (const sofa::defaulttype::Vec3& A,const sofa::defaulttype::Vec3& B, const sofa::defaulttype::Vec3& C ) -> SReal
    {
        return 0.5 * ((B-A).cross(C-A)).norm();
    };

    const auto& cot_angle = [](const sofa::defaulttype::Vec3& e0, const sofa::defaulttype::Vec3& e1) -> SReal
    {
        return 1.0/std::tan( std::acos(e0*e1) );
    };

    const auto& x0 = *x[0];
    const auto& x1 = *x[1];
    const auto& x2 = *x[2];
    const auto& x3 = *x[3];

    //Let's hope there is no degenerated triangles
    area.first  = area_triangle(x0,x1,x2);
    area.second = area_triangle(x0,x1,x3);

    //This is e_s in the paper
    sofa::defaulttype::Vec3 e[5] = {(x1-x0).normalized (),
                                    (x2-x1).normalized (),
                                    (x0-x2).normalized (),
                                    (x3-x0).normalized (),
                                    (x1-x3).normalized ()};

    //zero padding to keep up with the notation of the paper
    SReal c0[5] = { 0.0,
                    cot_angle(e[0],e[1]),
                    cot_angle(e[0],e[2]),
                    cot_angle(e[0],e[3]),
                    cot_angle(e[0],e[4])};

    sofa::defaulttype::Vec4 K( c0[1]+c0[4],
                               c0[2]+c0[3],
                              -c0[1]-c0[2],
                              -c0[3]-c0[4]);

    //Q is a symetric matrix (Looks a lot like the Loop decimation algorithm matrix ((Just sayin')) )
    Q(0,0) = K[0] * K[0]; Q(0,1) = K[0] * K[1]; Q(0,2) = K[0] * K[2]; Q(0,3) = K[0] * K[3];
    Q(1,0) = Q(0,1)     ; Q(1,1) = K[1] * K[1]; Q(1,2) = K[1] * K[2]; Q(1,3) = K[1] * K[3];
    Q(2,0) = Q(0,2)     ; Q(2,1) = Q(1,2)     ; Q(2,2) = K[2] * K[2]; Q(2,3) = K[2] * K[3];
    Q(3,0) = Q(0,3)     ; Q(3,1) = Q(1,3)     ; Q(3,2) = Q(2,3)     ; Q(3,3) = K[3] * K[3];

    Q *= (3.0/(area.first+area.second));
}

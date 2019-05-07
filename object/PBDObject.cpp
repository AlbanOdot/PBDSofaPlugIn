#include "PBDObject.hpp"
#include <Eigen/MatrixFunctions>
#include <sofa/core/objectmodel/BaseContext.h>
#include <SofaBaseMechanics/UniformMass.h>

PBDObject::PBDObject(sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * mobj,
                     sofa::core::topology::BaseMeshTopology * topo ): m_mechanicalObject(mobj), m_sofa_topology(topo)
{
    if(m_mechanicalObject && m_sofa_topology)
    {
        m_rest.emplace_back(m_mechanicalObject->readRestPositions ());
        m_mass = static_cast<sofa::component::mass::UniformMass<sofa::defaulttype::Vec3Types,SReal> *>(m_mechanicalObject->getContext ()->getMass())->getVertexMass ();
        m_invVertexMass = 1.0/m_mass;
    }
}

inline void PBDObject::setTopology(sofa::core::topology::BaseMeshTopology *topology)
{
    m_sofa_topology = topology;
    if(m_mechanicalObject && m_sofa_topology)
    {
        m_rest.emplace_back(m_mechanicalObject->readRestPositions ());
        m_mass = static_cast<sofa::component::mass::UniformMass<sofa::defaulttype::Vec3Types,SReal> *>(m_mechanicalObject->getContext ()->getMass())->getVertexMass ();
        m_invVertexMass = 1.0/m_mass;

    }

}

void PBDObject::optimizeTopology()
{
    computeStretchTopology ();
    computeBendingTopology ();
    computeTetrahedraBasis ();
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


void PBDObject::computeStretchTopology()
{
    //Compute the vertice oriented topology
    for(uint i = 0; i < m_rest[0].size(); ++i)
    {
        //Get the neighbors of point I
        const auto& neighbors = m_sofa_topology->getVerticesAroundVertex (i);
        std::vector<std::pair<uint,SReal>> neighborhood;
        for(uint j = 0; j < neighbors.size(); ++j)
        {
            if( neighbors[j] < i )//Unidirectionnal neighborhood
            {
                SReal d = (m_rest[0][i] - m_rest[0][neighbors[j]]).norm();
                neighborhood.emplace_back(std::pair<uint,SReal>(neighbors[j],d));
            }
        }
        m_stretch_topology.emplace_back(neighborhood);
    }

}
void PBDObject::computeBendingTopology()
{
    //Compute the vertice oriented topology
    const auto& triangles = m_sofa_topology->getTriangles ();
    m_bending_topology.resize(m_sofa_topology->getNbEdges ());
    m_triangle_rest_area.resize (m_sofa_topology->getNbEdges ());
    for(uint i = 0; i < m_rest[0].size(); ++i)
    {
        const auto& voisins = m_stretch_topology[i];

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

void PBDObject::computeTetrahedraBasis()
{
    const auto& tetrahedra = m_sofa_topology->getTetrahedra ();
    m_tetra_bases.resize (m_sofa_topology->getNbTetrahedra ());
    for(uint i = 0; i < tetrahedra.size(); ++i)
    {
        const auto& x = tetrahedra[i];
        const auto& r1 = m_rest[0][x[0]] - m_rest[0][x[3]];
        const auto& r2 = m_rest[0][x[1]] - m_rest[0][x[3]];
        const auto& r3 = m_rest[0][x[2]] - m_rest[0][x[3]];
        Eigen::Matrix3d Dm; Dm << r1[0],r1[1],r1[2],
                                  r2[0],r2[1],r2[2],
                                  r3[0],r3[1],r3[2];
        //We take 0.5 since the tetrahedron volume is half the det
        m_tetra_bases[i] = std::pair<float,Eigen::Matrix3d>(0.5*Dm.determinant(),Dm.inverse());
    }
}

void PBDObject::computeBeam()
{
    //Compute the vertice oriented topology
    const auto& edges = m_sofa_topology->getEdges ();
    std::vector<SReal> averageLength;
    const auto& rest = m_mechanicalObject->readRestPositions ();
    for(uint e = 0; e < edges.size(); ++e)
    {
        uint edge[2] = {edges[e][0],edges[e][1]};
        uint a = e+1 == edges.size () ? e : e+1;
        uint edgeN[2] = {edges[a][0],edges[a][1]};
        SReal l = 0.5* ( (rest[edge[0]] - rest[edge[1]]).norm () + (rest[edgeN[0]] - rest[edgeN[1]]).norm ());
        m_beam.emplace_back(BeamElement(m_mechanicalObject,edge,e,l));
    }
    //Post init of darboux vectors
    for(uint e = 0; e < edges.size() - 1; ++e)
    {
        m_beam[e].initRestDarboux(m_beam[e+1].m_q);
    }
    m_beam[m_beam.size() - 1].initRestDarboux (m_beam[m_beam.size() - 1].q ());
    m_beam.emplace_back(m_beam[m_beam.size() - 1]);
    auto& last = m_beam[m_beam.size() - 1];
    last.m_segmentIndicies.first = last.m_segmentIndicies.second;
    setupAngularVelocity ();
}

void PBDObject::computeStiffRod()
{
}

void PBDObject::setupAngularVelocity()
{
    const Eigen::Matrix3d id; id.Identity ();
    for(uint i = 0; i< m_beam.size (); ++i)
    {
        m_angularSpeed.emplace_back(Eigen::Vector3d(0,0,0));
        m_torque.emplace_back(Eigen::Vector3d(0,0,0));
        m_freeOrientation.emplace_back(m_beam[i].q());
        m_intertia.emplace_back(id);
    }
    for(uint i = 0; i< m_beam.size (); ++i)
    {
        m_intertia[i].setIdentity ();
    }
}


void PBDObject::applyFixedPoint(const std::vector<uint>& idx)
{
    for(auto& a : idx)
    {
        m_beam[a].setwq(0.0);
    }
}

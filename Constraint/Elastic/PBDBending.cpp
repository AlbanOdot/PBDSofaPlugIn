#include "PBDBending.hpp"
#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>

int PBDBendingClass = sofa::core::RegisterObject("Constraint that correct the bending.")
                      .add< PBDBending >();

void PBDBending::bwdInit ()
{
    auto node = dynamic_cast<sofa::simulation::Node*>(this->getContext());
    m_K = m_k.getValue ();
    m_mass = PBDVertexMass<sofa::defaulttype::Vec3Types>(m_mechanicalObject.getValue (), m_topology.getValue ());
    m_stretch_topology = PBDVertexTopology<sofa::defaulttype::Vec3Types>(m_mechanicalObject.getValue (),m_topology.getValue ());
    m_bending_topology = PBDBendingTopology(m_mechanicalObject.getValue (),m_topology.getValue ());
}

void PBDBending::solve(sofa::simulation::Node * node)
{

    WriteCoord p = m_pbdObject->getFreePosition ();
    uint pointCount = p.size();
    if(m_indices.getValue().empty())
    {
        for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
        {
            for( uint a = 0; a < pointCount; ++a)
            {
                //Get the edge of the corresponding neighbors
                for( const auto& voisin : m_stretch_topology.data ()[a])
                {
                    correction(a,voisin.first,p);
                }
            }
        }
    }
    else
    {
        const auto& idx = m_indices.getValue ();
        for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
        {
            for( const auto& a : idx)
            {
                //Get the edge of the corresponding neighbors
                for( const auto& voisin : m_stretch_topology.data ()[a])
                {
                    correction(a,voisin.first,p);
                }
            }
        }
    }
}


void PBDBending::correction ( uint a, uint b,WriteCoord&p)
{
    static SReal eps = 1e-6;
    uint edge_ID = m_topology.getValue()->getEdgeIndex(a,b);
    const auto& hessian_and_idx = m_bending_topology.bendingData ()[edge_ID];
    //Apply correction on the points
    for(const auto& data : hessian_and_idx)
    {
        const Vec3 *x[4] = { &p[a], &p[b], &p[data.first[0]], &p[data.first[1]] };
        Real invMass[4] = { m_mass.w(a),m_mass.w(b), m_mass.w(data.first[0]), m_mass.w(data.first[1])};
        const auto& Q = data.second;
        //           Real energy = 0.0;
        //           for (unsigned char k = 0; k < 4; k++)
        //               for (unsigned char j = k; j < 4; j++)
        //                   energy += data.second(j, k)*(dot(*x[k],*x[j]));
        //We only take half of the diagonal for obvious reasons
        //E(bending) = 0.5 * x Q x^T
        Real energy =  0.5 * Q(0, 0)*(dot(p[a],p[a])) + Q(1, 0)*(dot(p[b],p[a])) + Q(2, 0)*(dot(p[data.first[0]],p[a])) + Q(3, 0)*(dot(p[data.first[1]],p[a]))
                + 0.5 * Q(1, 1)*(dot(p[b],p[b])) + Q(2, 1)*(dot(p[data.first[0]],p[b])) + Q(3, 1)*(dot(p[data.first[1]],p[b]))
                + 0.5 * Q(2, 2)*(dot(p[data.first[0]],p[data.first[0]])) + Q(3, 2)*(dot(p[data.first[1]],p[data.first[0]]))
                + 0.5 * Q(3, 3)*(dot(p[data.first[1]],p[data.first[1]]));
        Vec3 gradC[4];
        gradC[0]  =  Q(0,0) * *x[0];
        gradC[0] +=  Q(0,1) * *x[1];
        gradC[0] +=  Q(0,2) * *x[2];
        gradC[0] +=  Q(0,3) * *x[3];

        gradC[1]  =  Q(0,1) * *x[0];
        gradC[1] +=  Q(1,1) * *x[1];
        gradC[1] +=  Q(1,2) * *x[2];
        gradC[1] +=  Q(1,3) * *x[3];

        gradC[2]  =  Q(0,2) * *x[0];
        gradC[2] +=  Q(1,2) * *x[1];
        gradC[2] +=  Q(2,2) * *x[2];
        gradC[2] +=  Q(2,3) * *x[3];

        gradC[3]  =  Q(0,3) * *x[0];
        gradC[3] +=  Q(1,3) * *x[1];
        gradC[3] +=  Q(2,3) * *x[2];
        gradC[3] +=  Q(3,3) * *x[3];

        Real sum_normGradC =  invMass[0] * gradC[0].norm2()
                              +invMass[1] * gradC[1].norm2()
                              +invMass[2] * gradC[2].norm2()
                              +invMass[3] * gradC[3].norm2();

        // exit early if required
        if (sum_normGradC > eps)
        {
            // compute impulse-based scaling factor
            const Real s = energy / sum_normGradC;

            p[a]             += (m_K * s * invMass[0]) * gradC[0];
            p[b]             += (m_K * s * invMass[1]) * gradC[1];
            p[data.first[0]] += (m_K * s * invMass[2]) * gradC[2];
            p[data.first[1]] += (m_K * s * invMass[3]) * gradC[3];


        }
    }
}

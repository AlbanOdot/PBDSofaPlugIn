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
    m_bending_topology.dampHighFrequencies (m_alpha_too.getValue ());
}

void PBDBending::solve(sofa::simulation::Node * node)
{

    WriteCoord p = m_pbdObject->getFreePosition ();
    WriteDeriv v = m_pbdObject->getFreeVelocity ();
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
                    correction(a,voisin.first,p,v);
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
                    correction(a,voisin.first,p,v);
                }
            }
        }
    }
}


void PBDBending::correction ( uint a, uint b,WriteCoord&p,WriteDeriv& v)
{
    static SReal eps = 1e-6;
    uint edge_ID = m_topology.getValue()->getEdgeIndex(a,b);
    const auto& hessian_and_idx = m_bending_topology.bendingData ()[edge_ID];
    //Apply correction on the points
    for(const auto& data : hessian_and_idx)
    {
        Real w[4] = { m_mass.w(a),m_mass.w(b), m_mass.w(data.first[0]), m_mass.w(data.first[1])};
        Real m[4] = { m_mass.m(a),m_mass.m(b), m_mass.m(data.first[0]), m_mass.m(data.first[1])};
        const auto& Q = data.second;
        //We only take half of the diagonal for obvious reasons
        //E(bending) = 0.5 * x Q x^T
        Real energy = 0.5 * Q[0]*(dot(p[a],p[a]))                         + Q[1]*(dot(p[b],p[a]))                        + Q[2]*(dot(p[data.first[0]],p[a])) + Q[3]*(dot(p[data.first[1]],p[a]))
                + 0.5 * Q[4]*(dot(p[b],p[b]))                         + Q[5]*(dot(p[data.first[0]],p[b]))            + Q[6]*(dot(p[data.first[1]],p[b]))
                + 0.5 * Q[7]*(dot(p[data.first[0]],p[data.first[0]])) + Q[8]*(dot(p[data.first[1]],p[data.first[0]]))
                + 0.5 * Q[9]*(dot(p[data.first[1]],p[data.first[1]]));

        Vec3 gradC[4]  =    {Q[0] * p[a] + Q[1] * p[b] + Q[2] * p[data.first[0]] + Q[3] * p[data.first[1]] - m_alpha_wann.getValue () * m[0]  * v[a],
                             Q[1] * p[a] + Q[4] * p[b] + Q[5] * p[data.first[0]] + Q[6] * p[data.first[1]] - m_alpha_wann.getValue () * m[1]  * v[b],
                             Q[2] * p[a] + Q[5] * p[b] + Q[7] * p[data.first[0]] + Q[8] * p[data.first[1]] - m_alpha_wann.getValue () * m[2]  * v[data.first[0]],
                             Q[3] * p[a] + Q[6] * p[b] + Q[8] * p[data.first[0]] + Q[9] * p[data.first[1]] - m_alpha_wann.getValue () * m[3]  * v[data.first[1]]};

        Real sum_normGradC =   w[0] * gradC[0].norm2()
                               +w[1] * gradC[1].norm2()
                               +w[2] * gradC[2].norm2()
                               +w[3] * gradC[3].norm2();

        // exit early if required
        if (sum_normGradC > eps)
        {
            // compute impulse-based scaling factor
            const Real s = energy / sum_normGradC;
            p[a]             += (m_K * s * w[0]) * gradC[0];
            p[b]             += (m_K * s * w[1]) * gradC[1];
            p[data.first[0]] += (m_K * s * w[2]) * gradC[2];
            p[data.first[1]] += (m_K * s * w[3]) * gradC[3];
        }
    }
}

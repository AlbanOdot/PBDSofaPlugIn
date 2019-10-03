#include "PBDIsometricBending.hpp"

#include <sofa/core/ObjectFactory.h>

int PBDIsometricBendingClass = sofa::core::RegisterObject("Constraint that correct the bending.")
                               .add< PBDIsometricBending >();

void PBDIsometricBending::bwdInit ()
{
    m_K = m_k.getValue();
    auto node = dynamic_cast<sofa::simulation::Node*>(this->getContext());
    m_K = 1.0-std::pow(1.0-m_k.getValue (),1.0 / ((double)m_nbIter.getValue ()));
    m_mass = PBDVertexMass<sofa::defaulttype::Vec3Types>(m_mechanicalObject.getValue (), m_topology.getValue ());
    m_stretch_topology = PBDVertexTopology<sofa::defaulttype::Vec3Types>(m_mechanicalObject.getValue (),m_topology.getValue ());
    m_bending_topology = PBDBendingTopology(m_mechanicalObject.getValue (),m_topology.getValue ());
    m_bending_topology.dampHighFrequencies (m_alpha_too.getValue ());
}

bool PBDIsometricBending::solve(sofa::simulation::Node* node)
{
    WriteCoord p = m_pbdObject->getFreePosition ();
    WriteDeriv v = m_pbdObject->getFreeVelocity ();
    uint pointCount = p.size();
    bool modification = false;

    for( uint a = 0; a < pointCount; ++a)
    {
        //Get the edge of the corresponding neighbors
        SReal w0 = m_mass.w(a);
        for( const auto& voisin : m_stretch_topology.data ()[a])
        {
            SReal w1 = m_mass.w(voisin.first);
            SReal wSum = w0 + w1;
            const sofa::defaulttype::Vec3& p_ij = p[a] - p[voisin.first];
            SReal l = p_ij.norm();
            const auto& displacement = ((l-voisin.second)/(l * wSum)) * p_ij;
            p[a]            -= w0 * displacement;
            p[voisin.first] += w1 * displacement;
            modification &= correction(a,voisin.first,p,v);
        }
    }
    return modification;
}


//We only take half of the diagonal for obvious reasons
//E(bending) = 0.5 * x Q x^T
//                    Real energy =  0.5 * Q(0, 0)*(dot(p[a],p[a])) + Q(1, 0)*(dot(p[voisin],p[a])) + Q(2, 0)*(dot(p[data.first[0]],p[a])) + Q(3, 0)*(dot(p[data.first[1]],p[a]))
//                            + 0.5 * Q(1, 1)*(dot(p[voisin],p[voisin])) + Q(2, 1)*(dot(p[data.first[0]],p[voisin])) + Q(3, 1)*(dot(p[data.first[1]],p[voisin]))
//                            + 0.5 * Q(2, 2)*(dot(p[data.first[0]],p[data.first[0]])) + Q(3, 2)*(dot(p[data.first[1]],p[data.first[0]]))
//                            + 0.5 * Q(3, 3)*(dot(p[data.first[1]],p[data.first[1]]));

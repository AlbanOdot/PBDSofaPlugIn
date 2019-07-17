#include "PBDIsometricBending.hpp"

#include <sofa/core/ObjectFactory.h>

int PBDIsometricBendingClass = sofa::core::RegisterObject("Constraint that correct the bending.")
                               .add< PBDIsometricBending >();

void PBDIsometricBending::bwdInit ()
{
    auto node = dynamic_cast<sofa::simulation::Node*>(this->getContext());
    m_K = 1.0-std::pow(1.0-m_k.getValue (),1.0 / ((double)m_nbIter.getValue ()));
    m_mass = PBDVertexMass<sofa::defaulttype::Vec3Types>(m_mechanicalObject.getValue (), m_topology.getValue ());
    m_stretch_topology = PBDVertexTopology<sofa::defaulttype::Vec3Types>(m_mechanicalObject.getValue (),m_topology.getValue ());
    m_bending_topology = PBDBendingTopology(m_mechanicalObject.getValue (),m_topology.getValue ());
}

void PBDIsometricBending::solve(sofa::simulation::Node* node)
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
            for( uint a : idx)
            {
                //Get the edge of the corresponding neighbors
                for( const auto& voisin : m_stretch_topology.data ()[a])
                {
                    //Enforce the isometric property
                    const sofa::defaulttype::Vec3& p_ij = p[a] - p[voisin.first];
                    SReal l = p_ij.norm();
                    const auto& dp = (0.5*(l-voisin.second)/l) * p_ij;
                    p[a]            -= dp;
                    p[voisin.first] += dp;
                    correction (a,voisin.first,p);//PBDBending
                }
            }
        }
    }
}


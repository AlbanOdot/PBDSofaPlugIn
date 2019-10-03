#include "PBDStretch.hpp"
#include <sofa/core/ObjectFactory.h>

int PBDStretchClass = sofa::core::RegisterObject("Constraint that correct the streching.")
                      .add< PBDStretch >();

void PBDStretch::bwdInit ()
{
    m_K = 1.0-std::pow(1.0-m_k.getValue (),1.0 / ((double)m_nbIter.getValue ()));
    m_stretch_topology = PBDVertexTopology<sofa::defaulttype::Vec3Types>(m_mechanicalObject.getValue (),m_topology.getValue ());
    m_mass = PBDVertexMass<sofa::defaulttype::Vec3Types>(m_mechanicalObject.getValue (), m_topology.getValue ());
}

bool PBDStretch::solve(sofa::simulation::Node * node)
{
    WriteCoord p = m_pbdObject->getFreePosition ();
    uint pointCount = p.size();
    if(m_indices.getValue().empty())
    {
        for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
        {
            for( uint i = 0; i < pointCount; ++i)
            {
                SReal w0 = m_mass.w(i);
                for( std::pair<uint,double>& voisin : m_stretch_topology.data()[i])
                {
                    correction(i,voisin,p,w0);
                }
            }
        }
    }
    else
    {
        const auto& idx = m_indices.getValue ();
        //#pragma omp parallel for
        for( uint i = 0; i < idx.size(); ++i)
        {
            SReal w0 = m_mass.w(i);
            for( std::pair<uint,double>& voisin : m_stretch_topology.data()[i])
            {
                correction(i,voisin,p,w0);
            }
        }
    }
    return true;
}


void PBDStretch::correction (uint i, const std::pair<uint,double>& voisin, WriteCoord& p, SReal w0)
{
    SReal w1 = m_mass.w(voisin.first);
    SReal wSum = w0 + w1;
    const sofa::defaulttype::Vec3& p_ij = p[i] - p[voisin.first];
    SReal l = p_ij.norm();
    const auto& displacement = (m_K*(l-voisin.second)/(l * wSum)) * p_ij;
    p[i]            -= w0 * displacement;
    p[voisin.first] += w1 * displacement;
}

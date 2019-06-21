#include "PBDStretch.hpp"
#include <sofa/core/ObjectFactory.h>

int PBDStretchClass = sofa::core::RegisterObject("Constraint that correct the streching.")
                      .add< PBDStretch >();

void PBDStretch::bwdInit ()
{
    m_K = 1.0-std::pow(1.0-m_k.getValue (),1.0 / ((double)m_nbIter.getValue ()));
}

void PBDStretch::solve(PBDObject<sofa::defaulttype::Vec3Types> &object, WriteCoord &p)
{

    if(!object.hasDataType(PBDObject<sofa::defaulttype::Vec3Types>::STRETCH))
        object.computeStretchTopology ();

    uint pointCount = p.ref().size();
    if(m_indices.getValue().empty())
    {
        for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
        {
            for( uint i = 0; i < pointCount; ++i)
            {
                for( const auto& voisin : object.topology().data ()[i])
                {
                    correction(i,voisin,p);
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
            for( const auto& voisin : object.topology().data ()[i])
            {
                correction(i,voisin,p);
            }
        }
    }
}


void PBDStretch::correction (uint i, const std::pair<uint,double>& voisin, WriteCoord &p)
{
    const sofa::defaulttype::Vec3& p_ij = p[i] - p[voisin.first];
    SReal l = p_ij.norm();
    const auto& displacement = (0.5*m_K*(l-voisin.second)/l) * p_ij;
    p[i]            -= displacement;
    p[voisin.first] += displacement;
}

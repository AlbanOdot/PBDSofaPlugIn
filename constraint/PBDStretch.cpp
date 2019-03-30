#include "PBDStretch.hpp"

#include <sofa/core/ObjectFactory.h>

int PBDStretchClass = sofa::core::RegisterObject("Constraint that correct the streching.")
                            .add< PBDStretch >();
PBDStretch::PBDStretch(sofa::simulation::Node *gnode)
    : PBDBaseConstraint(true),
      m_k(initData(&m_k,(SReal)1.0,"stretch","Stretching factor"))
{

}

PBDStretch::PBDStretch(unsigned int objectSize)
    : PBDBaseConstraint(true),
      m_k(initData(&m_k,(SReal)1.0,"stretch","Stretching factor"))
{

}

sofa::defaulttype::BaseMatrix * PBDStretch::getConstraintMatrix ()
{
    return nullptr;
}

void PBDStretch::solve(PBDObject &object, WriteCoord &p)
{
    uint pointCount = p.ref().size();
    float k = m_k.getValue ();
    auto correction = [&object,&p,&k](uint i)
    {
        const auto& voisins = object.topology()[i];

        for( const auto& voisin : voisins)
        {
            const sofa::defaulttype::Vec3& p_ij = p[i] - p[voisin.first];
            SReal l = p_ij.norm();
            const auto& displacement = (0.5*k*(l-voisin.second)/l) * p_ij;
            p[i]            -= displacement;
            p[voisin.first] += displacement;
        }
    };

    if(m_indices.getValue().empty())
    {
        for( uint i = 0; i < pointCount; ++i)
        {
            correction(i);
        }
    }
    else
    {
        const auto& idx = m_indices.getValue ();
        //#pragma omp parallel for
        for( uint i = 0; i < idx.size(); ++i)
        {
            correction(idx[i]);
        }
    }
}



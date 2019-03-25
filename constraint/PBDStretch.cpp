#include "PBDStretch.hpp"

#include <sofa/core/ObjectFactory.h>

int PBDStretchClass = sofa::core::RegisterObject("Constraint that fixe a point.")
                            .add< PBDStretch >();

PBDStretch::PBDStretch(unsigned int objectSize)
    : PBDBaseConstraint(true),
      m_k(initData(&m_k,1.0,"stretch","Stretching factor"))
{

}

sofa::defaulttype::BaseMatrix * PBDStretch::getConstraintMatrix ()
{
    return nullptr;
}

void PBDStretch::solve(PBDObject &object, WriteCoord &p)
{
    uint pointCount = p.ref().size();
    float k = 0.5;//m_k.getValue ();
    if(m_indices.getValue().empty())
    {
        for( uint i = 0; i < pointCount; ++i)
        {
            const auto& voisins = object.topology()[i];
            for( const auto& voisin : voisins)
            {
                const sofa::defaulttype::Vec3& p_ij = p[i] - p[voisin.first];
                SReal l = p_ij.norm();
                const auto& displacement = (0.5*k*(l-voisin.second)/l) * p_ij;
                p[i] -= displacement;
                p[voisin.first] += displacement;
            }
        }
    }
    else
    {
        for( const auto& i : m_indices.getValue())
        {
            const auto& voisins = object.topology()[i];
            for( const auto& voisin : voisins)
            {
                const sofa::defaulttype::Vec3& p_ij = p[i] - p[voisin.first];
                SReal l = p_ij.norm();
                SReal d = voisin.second;
                const auto& displacement = (0.5*k*(l-d)/l) * p_ij;
                p[i] -= displacement;
                p[voisin.first] += displacement;
            }
        }
    }
}



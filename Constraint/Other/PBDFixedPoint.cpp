#include "PBDFixedPoint.hpp"
#include <sofa/core/ObjectFactory.h>

int PBDFixedPointClass = sofa::core::RegisterObject("Constraint that fixe a point.")
                         .add< PBDFixedPoint<sofa::defaulttype::Vec3Types> >(true)
                         .add < PBDFixedPoint<sofa::defaulttype::Rigid3Types>>();

template < class T>
void PBDFixedPoint<T>::solve(PBDObject<T> &object, WriteCoord &p)
{
    for(const auto& idx : PBDBaseConstraint<T>::m_indices.getValue ())
    {
        p[idx] = object.rest()[idx];
    }
}

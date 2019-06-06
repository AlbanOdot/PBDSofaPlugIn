#include "PBDFixedPoint.hpp"
#include <sofa/core/ObjectFactory.h>

int PBDFixedPointClass = sofa::core::RegisterObject("Constraint that fixe a point.")
                            .add< PBDFixedPoint >();

PBDFixedPoint::PBDFixedPoint(unsigned int objectSize) : PBDBaseConstraint()
{
}


void PBDFixedPoint::solve(PBDObject &object, WriteCoord &p)
{
    for(const auto& idx : m_indices.getValue ())
    {
        p[idx] = object.rest()[idx];
    }
}

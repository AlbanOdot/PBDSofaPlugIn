#include "PBDFixedPoint.hpp"
#include <sofa/core/ObjectFactory.h>


int PBDFixedPointClass = sofa::core::RegisterObject("Constraint that fixes a point")
                         .add< PBDFixedPoint >();
void PBDFixedPoint::solve(sofa::simulation::Node* node)
{
    const ReadCoord& r = m_mechanicalObject.getValue ()->readRestPositions ();
    WriteCoord p = m_pbdObject->getFreePosition ();
    for(const auto& idx : m_indices.getValue ())
    {
        p[idx] = r[idx];
    }
}


int PBDFixedRigidPointClass = sofa::core::RegisterObject("Constraint that fixes a rigid point")
                              .add< PBDFixedRigidPoint >();
void PBDFixedRigidPoint::solve(sofa::simulation::Node* node)
{
    const ReadCoordR& r = m_mechanicalObject.getValue ()->readRestPositions ();
    WriteCoordR p = m_pbdObject->getFreePosition ();;
    for(const auto& idx : m_indices.getValue ())
    {
        p[idx] = r[idx];
    }
}

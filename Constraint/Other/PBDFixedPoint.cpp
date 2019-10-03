#include "PBDFixedPoint.hpp"
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/defaulttype/VecTypes.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/core/behavior/ForceField.inl>
#include <SofaBaseTopology/TopologyData.inl>

using namespace sofa::defaulttype;
using namespace sofa::defaulttype;
using namespace	sofa::component::topology;
using namespace sofa;
using namespace core::topology;

int PBDFixedPointClass = sofa::core::RegisterObject("Constraint that fixes a point")
                         .add< PBDFixedPoint >();
bool PBDFixedPoint::solve(sofa::simulation::Node* node)
{
    const ReadCoord& r = m_mechanicalObject.getValue ()->readRestPositions ();
    WriteCoord p = m_pbdObject->getFreePosition ();
    WriteDeriv v = m_mechanicalObject.getValue ()->writeVelocities ();
    static Vec3 v0(0,0,0);
    for(const auto& idx : m_indices.getValue ())
    {
        p[idx] = r[idx];
        v[idx] = v0;
    }
    return false;
}

void PBDFixedPoint::draw(const sofa::core::visual::VisualParams *vparams)
{
    vparams->drawTool()->saveLastState();
    const VecCoord& x = m_mechanicalObject.getValue ()->read(core::ConstVecCoordId::position())->getValue();
    std::vector< sofa::defaulttype::Vector3 > points;
    sofa::defaulttype::Vector3 point;
    for(const auto& idx : m_indices.getValue ())
    {
        point = Vec3Types::getCPos(x[idx]);
        points.push_back(point);
    }
    vparams->drawTool()->drawPoints(points, 5, sofa::defaulttype::Vec<4,float>(0.5,0.5,0.5,1));
}

int PBDFixedRigidPointClass = sofa::core::RegisterObject("Constraint that fixes a rigid point")
                              .add< PBDFixedRigidPoint >();
bool PBDFixedRigidPoint::solve(sofa::simulation::Node* node)
{
    const ReadCoordR& r = m_mechanicalObject.getValue ()->readRestPositions ();
    WriteCoordR p = m_pbdObject->getFreePosition ();
    WriteDerivR v = m_mechanicalObject.getValue ()->writeVelocities ();
    for(const auto& idx : m_indices.getValue ())
    {
        p[idx] = r[idx];
        v[idx] *= 0.0;
    }
    return false;
}

void PBDFixedRigidPoint::draw(const sofa::core::visual::VisualParams *vparams)
{
    vparams->drawTool()->saveLastState();
    const VecCoordR& x = m_mechanicalObject.getValue ()->read(core::ConstVecCoordId::position())->getValue();
    std::vector< sofa::defaulttype::Vector3 > points;
    sofa::defaulttype::Vector3 point;
    for(const auto& idx : m_indices.getValue ())
    {
        point = Rigid3Types::getCPos(x[idx]);
        points.push_back(point);
    }
    vparams->drawTool()->drawPoints(points, 5, sofa::defaulttype::Vec<4,float>(0.5,0.5,0.5,1));
}

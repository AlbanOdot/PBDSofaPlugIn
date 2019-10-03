#include "PBDContactManager.hpp"

#include <sofa/core/collision/NarrowPhaseDetection.h>
#include <sofa/core/collision/ContactManager.h>

#include <sofa/core/CollisionModel.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/simulation/CollisionBeginEvent.h>
#include <sofa/simulation/CollisionEndEvent.h>

using namespace sofa;
using namespace core;
using namespace collision;

PBDCollisionHandler::PBDCollisionHandler(  sofa::simulation::Node* n)
    :  m_narrowPhase(nullptr), m_node(n)
{
}

PBDCollisionHandler::~PBDCollisionHandler()
{
}

void PBDCollisionHandler::init(void)
{
    m_narrowPhase = m_node->getContext()->get<core::collision::NarrowPhaseDetection>();

}

void PBDCollisionHandler::handleCollision( )
{

    m_solver->clearCollisions();

    const NarrowPhaseDetection::DetectionOutputMap& detectionOutputsMap = m_narrowPhase->getDetectionOutputs();

    if ( detectionOutputsMap.size() == 0 )
        return;

    //// Solver create all of the collision constraint then solve it
    for (core::collision::NarrowPhaseDetection::DetectionOutputMap::const_iterator it = detectionOutputsMap.begin(); it!=detectionOutputsMap.end(); ++it )
    {
        m_solver->generateCollisionConstraint(it->first,it->second);
    }

}

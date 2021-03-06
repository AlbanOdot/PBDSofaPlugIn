/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2018 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include "./PBDInitVisitor.hpp"
#include "../Solver/PBDSolver.hpp"

#include <sofa/simulation/MechanicalVisitor.h>
#include <sofa/simulation/CollisionVisitor.h>

#include <sofa/simulation/PropagateEventVisitor.h>
#include <sofa/simulation/CollisionBeginEvent.h>
#include <sofa/simulation/CollisionEndEvent.h>
#include <sofa/simulation/IntegrateBeginEvent.h>
#include <sofa/simulation/IntegrateEndEvent.h>
#include <sofa/simulation/PropagateEventVisitor.h>

#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/VectorOperations.h>
#include <sofa/core/behavior/MultiVec.h>

#include <sofa/helper/AdvancedTimer.h>
#include "../AnimationLoop/PBDAnimationLoop.hpp"

//#include "MechanicalIntegration.h"

using namespace sofa::core;
using namespace sofa::simulation;
using namespace sofa;


PBDInitVisitor::PBDInitVisitor(const core::ExecParams* params, SReal dt)
    : Visitor(params)
    , dt(dt)
    , firstNodeVisited(false)
{
}

PBDInitVisitor::PBDInitVisitor(const core::ExecParams* params)
    : Visitor(params)
    , dt(0)
    , firstNodeVisited(false)
{
}

void PBDInitVisitor::processBehaviorModel(simulation::Node*, core::BehaviorModel* obj)
{
    obj->updatePosition(getDt());
}

void PBDInitVisitor::fwdInteractionForceField(simulation::Node*, core::behavior::BaseInteractionForceField* obj)
{

    MultiVecDerivId   ffId      = VecDerivId::externalForce();
    MechanicalParams mparams; // = MechanicalParams::defaultInstance();
    mparams.setDt(this->dt);
    obj->addForce(&mparams, ffId);
}

void PBDInitVisitor::initCollisionPipeline(simulation::Node* node, core::collision::Pipeline* obj)
{
    CollisionBeginEvent evBegin;
    PropagateEventVisitor eventPropagation( params, &evBegin);
    eventPropagation.execute(node->getContext());

    {
        CollisionResetVisitor act(this->params);
        node->execute(&act);
    }
    {
        CollisionDetectionVisitor act(this->params);
        node->execute (&act);
    }
}

void PBDInitVisitor::processOdeSolver(simulation::Node* node, core::behavior::OdeSolver* solver)
{
    solver->solve(params, getDt());
}

Visitor::Result PBDInitVisitor::processNodeTopDown(simulation::Node* node)
{
    if (!node->isActive()) return Visitor::RESULT_PRUNE;
    if (node->isSleeping()) return Visitor::RESULT_PRUNE;
    if (!firstNodeVisited)
    {
        firstNodeVisited=true;

        sofa::core::ConstraintParams cparams(*this->params);
        MechanicalResetConstraintVisitor resetConstraint(&cparams);
        node->execute(&resetConstraint);
    }

    if (dt == 0) setDt(node->getDt());
    else node->setDt(dt);

    if (node->collisionPipeline != NULL)
    {
        initCollisionPipeline(node, node->collisionPipeline);
    }
    if (!node->solver.empty() )
    {
        sofa::helper::AdvancedTimer::StepVar timer("Mechanical",node);
        SReal nextTime = node->getTime() + dt;
        {
            IntegrateBeginEvent evBegin;
            PropagateEventVisitor eventPropagation( this->params, &evBegin);
            eventPropagation.execute(node);
        }

        MechanicalBeginIntegrationVisitor beginVisitor(this->params, dt);
        node->execute(&beginVisitor);

        sofa::core::MechanicalParams m_mparams(*this->params);
        m_mparams.setDt(dt);

        for( unsigned i=0; i<node->solver.size(); i++ )
        {
            ctime_t t0 = begin(node, node->solver[i]);
            node->solver[i]->solve(params, getDt());
            end(node, node->solver[i], t0);
        }

        MechanicalPropagateOnlyPositionAndVelocityVisitor(&m_mparams, nextTime,
                                                          VecCoordId::position(),
                                                          VecDerivId::velocity(), true).execute( node );

        MechanicalEndIntegrationVisitor endVisitor(this->params, dt);
        node->execute(&endVisitor);

        {
            IntegrateEndEvent evBegin;
            PropagateEventVisitor eventPropagation(this->params, &evBegin);
            eventPropagation.execute(node);
        }

        return RESULT_PRUNE;
    }
    {
        // process InteractionForceFields
        for_each(this, node, node->interactionForceField, &PBDInitVisitor::fwdInteractionForceField);
        return RESULT_CONTINUE;
    }

}

#include "PBDAnimationLoop.hpp"
#include <sofa/core/ObjectFactory.h>


#include <SofaBaseMechanics/MechanicalObject.h>

//Visitors
#include <sofa/simulation/AnimateVisitor.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/PropagateEventVisitor.h>
#include <sofa/simulation/BehaviorUpdatePositionVisitor.h>
#include <sofa/simulation/SolveVisitor.h>
#include <sofa/simulation/AnimateEndEvent.h>
#include <sofa/simulation/UpdateMappingEndEvent.h>
#include <sofa/simulation/UpdateBoundingBoxVisitor.h>
#include <sofa/simulation/UpdateMappingVisitor.h>




using namespace sofa::component::container;
using namespace sofa::core::objectmodel;
using namespace sofa::defaulttype;
using namespace sofa::core::behavior;
using namespace sofa::simulation;

int PBDAnimationLoopClass = sofa::core::RegisterObject("Simulation loop to use in scene without constraints nor contact.")
                            .add< PBDAnimationLoop >()
                            .addDescription(R"(
                                            This loop do the following steps:
                                            - build and solve all linear systems in the scene : collision and time integration to compute the new values of the dofs
                                            - update the context (dt++)
                                            - update the mappings
                                            - update the bounding box (volume covering all objects of the scene))");

PBDAnimationLoop::PBDAnimationLoop(sofa::simulation::Node* _gnode)
    : Inherit()
    , gnode(_gnode)
{
}

PBDAnimationLoop::~PBDAnimationLoop()
{
    delete m_restPositions;
    delete m_freePosition;
}


void PBDAnimationLoop::init()
{
    if (!gnode)
        gnode = dynamic_cast<sofa::simulation::Node*>(this->getContext());
    m_context = gnode->getContext ();
    m_solver = gnode->solver[0];
    m_mechanicalObject = m_context->getObjects< MechanicalObject< sofa::defaulttype::Vec3Types > >(BaseContext::SearchDown)[0];
    // m_restPositions = new ReadCoord(m_mechanicalObject->readRestPositions());
    //m_freePosition = new WriteCoord(m_mechanicalObject->writePositions());
}

void PBDAnimationLoop::setNode( sofa::simulation::Node* n )
{
    gnode=n;
}

void PBDAnimationLoop::step(const sofa::core::ExecParams* params, SReal dt)
{
    if (dt == 0)
    {
        dt = gnode->getDt();
    }


    sofa::simulation::common::VectorOperations vop(params, getContext());
    sofa::simulation::common::MechanicalOperations mop(params, getContext());

    MultiVecCoord pos(&vop, sofa::core::VecCoordId::position() );
    MultiVecDeriv vel(&vop, sofa::core::VecDerivId::velocity() );
    MultiVecCoord freePos(&vop, sofa::core::VecCoordId::freePosition() );

    //HACKY_HACKS_DO_NOT_REPRODUCE
    MultiVecCoord zero(&vop, sofa::core::VecCoordId::position() );
    //zero = pos + pos * (-1)
    sofa::simulation::MechanicalVOpVisitor zeroProducer(params, zero, pos, pos, -1.0);
    zeroProducer.setMapped(true);
    getContext()->executeVisitor(&zeroProducer);





    sofa::core::ConstraintParams cparams(*params);
    cparams.setX(freePos);

    MultiVecDeriv dx(&vop, sofa::core::VecDerivId::dx()); dx.realloc(&vop, false, true);
    MultiVecDeriv df(&vop, sofa::core::VecDerivId::dforce()); df.realloc(&vop, false, true);

    sofa::simulation::MechanicalVInitVisitor< sofa::core::V_COORD >(params, sofa::core::VecCoordId::freePosition(), sofa::core::ConstVecCoordId::position(), true).execute(gnode);

    //AnimateBeginEvent
    sofa::simulation::AnimateBeginEvent ev ( dt );
    sofa::simulation::PropagateEventVisitor act ( params, &ev );
    gnode->execute ( act );

    //Update the position
    sofa::simulation::BehaviorUpdatePositionVisitor beh(params , dt);
    gnode->execute(&beh);


    //Prepare the system for an integration step
    sofa::simulation::MechanicalBeginIntegrationVisitor beginVisitor(params, dt);
    gnode->execute(&beginVisitor);


    //Compute the new velocity and update free position
    //vel = vel + m * sum(F_ext)
    sofa::simulation::SolveVisitor freeMotion(params, dt, true);
    gnode->execute(&freeMotion);
    mop.projectResponse(vel);
    mop.propagateDx(vel, true);

    // freePos = pos + vel * dt
    sofa::simulation::MechanicalVOpVisitor freePosVis(params, freePos, pos, vel, dt);
    freePosVis.setMapped(true);
    getContext()->executeVisitor(&freePosVis);

    /*
     *
     *
     * There should be a lot of code here to solve constraints
     *
     *
     */

    //vel = (freePos-pos)/dt
    //freePos_minus_pos = freePos + pos*(-1)
    MultiVecCoord freePos_minus_pos(&vop, sofa::core::VecCoordId::freePosition() );
    sofa::simulation::MechanicalVOpVisitor freePos_m_pos(params, freePos_minus_pos, freePos, pos, -1.0);
    freePos_m_pos.setMapped(true);
    getContext()->executeVisitor(&freePos_m_pos);
    //vel = zero + freePos_minus_pos * (1/dt)
    sofa::simulation::MechanicalVOpVisitor newVel(params, vel, zero, freePos_minus_pos, 1./dt);
    newVel.setMapped(true);
    getContext()->executeVisitor(&newVel);
    //pos = freePos
    //pos = zero + freePos * 1
    sofa::simulation::MechanicalVOpVisitor newPos(params, pos, zero, freePos, 1.0);
    newPos.setMapped(true);
    getContext()->executeVisitor(&newPos);

    mop.propagateXAndV(pos, vel);


    AnimateEndEvent evEnd ( dt );
    PropagateEventVisitor propEvEnd ( params, &evEnd );
    gnode->execute ( propEvEnd );

    gnode->execute<UpdateMappingVisitor>(params);

    UpdateMappingEndEvent evMapEnd ( dt );
    PropagateEventVisitor propEvMapEnd ( params , &evMapEnd );
    gnode->execute ( propEvMapEnd );

    gnode->execute<UpdateBoundingBoxVisitor>(params);





    /*
               //Apply external forces on velocity only

               sofa::helper::ReadAccessor<Data<VecDeriv>> externalForces = mechanicalObjects[currentMechObj]->readForces ();
               m_integrator.integrateExternalForces(velocities,externalForces,dt);

               //Temporary integration
               Data<VecCoord> newPositions;
               m_integrator.integrateTmp(newPositions,positions,velocities,dt);


               //Generate collision constraint
               //TODO ^^^^



               //Project Constraints
               //At the moment we create constraint on the fly because it's an alpha version of the alpha version

               //Save the fixed points coordinates
               //We set point 3 ,39 and 64 to be fixed this had to be done last
               std::vector<uint> fixedPointIndices = {3,39,64};

               //We set point 3 ,39 and 64 to be fixed this had to be done first
               //m_integrator.solveFixedPointConstraint(tmpPos,positions,fixedPointIndices);

               //We define rigid constraint and solve constraint in O(n) -> F(a->b) = - F(b->a)
               m_integrator.solveDistanceConstraint(tmpPos,positions);

               //Final integration
               m_integrator.PBDUpdate(tmpPos,velocities,positions,dt);
               //Find a way to put the new values into the mechanical object*/


}

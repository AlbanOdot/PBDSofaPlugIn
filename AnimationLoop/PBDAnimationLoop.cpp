#include "PBDAnimationLoop.hpp"
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/visual/VisualParams.h>
#include <SofaConstraint/LCPConstraintSolver.h>

//Mecha stuffs
#include <SofaBaseMechanics/MechanicalObject.h>

#include <sofa/helper/AdvancedTimer.h>

//Visitors
#include <sofa/simulation/PropagateEventVisitor.h>
#include <sofa/simulation/UpdateMappingEndEvent.h>
#include <sofa/simulation/UpdateMappingVisitor.h>
#include <sofa/simulation/BehaviorUpdatePositionVisitor.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/SolveVisitor.h>
#include <sofa/simulation/VectorOperations.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/AnimateVisitor.h>
#include <sofa/simulation/AnimateEndEvent.h>
#include <sofa/simulation/UpdateContextVisitor.h>
#include <sofa/simulation/UpdateBoundingBoxVisitor.h>
#include <sofa/simulation/IntegrateBeginEvent.h>
#include <sofa/simulation/IntegrateEndEvent.h>
#include <SofaConstraint/LCPConstraintSolver.h>
#include "../Visitor/PBDVisitor.hpp"

#include <omp.h>




template class SOFA_CORE_API sofa::helper::WriteAccessor<sofa::helper::vector<sofa::defaulttype::RigidCoord<3,SReal>>>;


using namespace sofa::component::container;
using namespace sofa::core::objectmodel;
using namespace sofa::defaulttype;
using namespace sofa::core::behavior;
using namespace sofa::simulation;
using namespace sofa::component;
using namespace sofa;
using helper::system::thread::CTime;
using sofa::helper::ScopedAdvancedTimer;

int PBDAnimationLoopClass = sofa::core::RegisterObject("Simulation loop to use in scene without constraints nor contact.")
                            .add< PBDAnimationLoop >()
                            .addDescription(R"(
                                            This loop do the following steps:
                                            - build and solve all linear systems in the scene : collision and time integration to compute the new values of the dofs
                                            - update the context (dt++)
                                            - update the mappings))");

PBDAnimationLoop::PBDAnimationLoop(sofa::simulation::Node* _gnode)
    :  Inherit1(_gnode),
      m_nbIter(initData(&m_nbIter,1,"iter","Number of iteration for the solver"))
    , m_solveVelocityConstraintFirst(initData(&m_solveVelocityConstraintFirst , false, "solveVelocityConstraintFirst", "solve separately velocity constraint violations before position constraint violations"))
    , d_threadSafeVisitor(initData(&d_threadSafeVisitor, false, "threadSafeVisitor", "If true, do not use realloc and free visitors in fwdInteractionForceField.")),
      f_rayleighMass( initData(&f_rayleighMass,(SReal)0.0,"rayleighMass","Rayleigh damping coefficient related to mass"))

{
}

PBDAnimationLoop::~PBDAnimationLoop()
{

}

void PBDAnimationLoop::init()
{
    m_context = gnode->getContext();
}

void PBDAnimationLoop::bwdInit ()
{
    //On récupère les topologies
    auto topologies = m_context->getObjects<sofa::core::topology::BaseMeshTopology>(BaseContext::SearchDown);
    m_solver.setupSolver(gnode,m_nbIter.getValue (),f_rayleighMass.getValue ());
    m_firststep = true;

}

void PBDAnimationLoop::setNode( sofa::simulation::Node* n )
{
    gnode=n;
}


static inline void beginEventAndComputeSofaPhysics(const sofa::core::ExecParams* params,
                                            sofa::simulation::Node* gnode,
                                            const SReal dt)
{

}

static inline void EndEventAndUpdateVisitors(const sofa::core::ExecParams* params,
                                      sofa::simulation::Node* gnode,
                                      const SReal dt){

    gnode->execute< UpdateSimulationContextVisitor >(params);

    {
        AnimateEndEvent ev ( dt );
        PropagateEventVisitor act ( params, &ev );
        gnode->execute ( act );
    }

    //Visual Information update: Ray Pick add a MechanicalMapping used as VisualMapping
    gnode->execute< UpdateMappingVisitor >(params);
    {
        UpdateMappingEndEvent ev ( dt );
        PropagateEventVisitor act ( params , &ev );
        gnode->execute ( act );
    }

    if (!SOFA_NO_UPDATE_BBOX)
    {
        sofa::helper::ScopedAdvancedTimer timer("UpdateBBox");
        gnode->execute< UpdateBoundingBoxVisitor >(params);
    }
}
void PBDAnimationLoop::step(const sofa::core::ExecParams* params,
                            SReal dt)
{
    if (dt == 0)
        dt = gnode->getDt();

    double startTime = gnode->getTime();
    {
        AnimateBeginEvent ev ( dt );
        PropagateEventVisitor act ( params, &ev );
        gnode->execute ( act );
    }

    BehaviorUpdatePositionVisitor beh(params , dt);
    gnode->execute(&beh);

    //update the orientations
    m_solver.setupOrientations (dt);
    //Solve each constraints one by one
    m_solver.solvePBDConstraints(params);

    //Integrate as said in the PBD method
    m_solver.integrate (dt);

    AnimateVisitor act(params, dt);
    gnode->execute ( act );

    gnode->setTime ( startTime + dt );

    gnode->execute< UpdateSimulationContextVisitor >(params);

    {
        AnimateEndEvent ev ( dt );
        PropagateEventVisitor act ( params, &ev );
        gnode->execute ( act );
    }

    //Visual Information update: Ray Pick add a MechanicalMapping used as VisualMapping
    gnode->execute< UpdateMappingVisitor >(params);
    {
        UpdateMappingEndEvent ev ( dt );
        PropagateEventVisitor act ( params , &ev );
        gnode->execute ( act );
    }

    if (!SOFA_NO_UPDATE_BBOX)
    {
        sofa::helper::ScopedAdvancedTimer timer("UpdateBBox");
        gnode->execute< UpdateBoundingBoxVisitor >(params);
    }
}

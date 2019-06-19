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
#include <SofaConstraint/LCPConstraintSolver.h>

#include <omp.h>





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
                                            - update the mappings))");

PBDAnimationLoop::PBDAnimationLoop(sofa::simulation::Node* _gnode)
    :  Inherit(),
      gnode(_gnode),
      m_nbIter(initData(&m_nbIter,(int)1,"iter","Number of iteration for the solver"))
{
}

PBDAnimationLoop::~PBDAnimationLoop()
{
}


void PBDAnimationLoop::init()
{
    if (!gnode)
        gnode = dynamic_cast<sofa::simulation::Node*>(this->getContext());
    m_context = gnode->getContext();
}

void PBDAnimationLoop::bwdInit ()
{
    //On récupère les topologies
    auto topologies = m_context->getObjects<sofa::core::topology::BaseMeshTopology>(BaseContext::SearchDown);
    auto mechanicalObjects = m_context->getObjects< MechanicalObject< sofa::defaulttype::Vec3Types > >(BaseContext::SearchDown);
    m_integrator.setUpIntegrator(gnode,m_nbIter.getValue ());
    for(uint i = 0; i < mechanicalObjects.size (); ++i)
    {
        m_objects.emplace_back(PBDObject(mechanicalObjects[i],topologies[i]));
    }

}

void PBDAnimationLoop::setNode( sofa::simulation::Node* n )
{
    gnode=n;
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

    AnimateVisitor act(params, dt);
    gnode->execute ( act );

    //Solve PBDConstraints
    sofa::core::MechanicalParams mparams(*params);
    static const Vec3 zero(0,0,0);
    const float inv_dt = 1.0/dt;

    for(auto& object : m_objects)
    {
        //Object parameters
        WriteCoord x = object.position();
        WriteDeriv v = object.velocity();

        //We will compute constrainst on p
        Coordinates freeCoord(x.ref());
        WriteCoord p = freeCoord;

        //Solve all of the constraints
        m_integrator.solveConstraint(object,p);

        //Apply torque and angular velocity
        m_integrator.integrateAngularVelocity(object,dt);

        //Integrate using PBD method
        m_integrator.updatePosAndVel(object,p,x,v,inv_dt);

    }

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

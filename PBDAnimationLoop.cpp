#include "PBDAnimationLoop.hpp"
#include <sofa/core/ObjectFactory.h>

//Mecha stuffs
#include <SofaBaseMechanics/MechanicalObject.h>
#include<sofa/core/behavior/MultiVec.h>

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
}


void PBDAnimationLoop::init()
{
    if (!gnode)
        gnode = dynamic_cast<sofa::simulation::Node*>(this->getContext());
    m_context = gnode->getContext();
    m_mechanicalObjects = m_context->getObjects< MechanicalObject< sofa::defaulttype::Vec3Types > >(BaseContext::SearchDown);
}

void PBDAnimationLoop::setNode( sofa::simulation::Node* n )
{
    gnode=n;
}

void PBDAnimationLoop::step(const sofa::core::ExecParams* params, SReal dt)
{
    std::cout << "COMPUTING FRAME [ "<< frame++<<" ]"<<std::endl;
    if (dt == 0)
    {
        dt = gnode->getDt();
    }
    sofa::core::MechanicalParams mparams(*params);
    for(auto& Obj : m_mechanicalObjects)
    {
        std::cout << "<<INIT BEGIN | ...";
        WriteCoord positions = Obj->writePositions ();
        WriteDeriv velocities = Obj->writeVelocities ();
        Vec3 zero;zero.set(0,0,0);
        VecDeriv zeros(positions.ref().size(),zero);
        Derivatives dforces(zeros);

        std::cout << " | INIT END >>" << std::endl;

        std::cout << "<<COMPUTE FORCES BEGIN | ...";
        //Accumulates all external forces (Gravities,interactions etc)
        for(uint i = 0;  i < gnode->forceField.size(); ++i)
        {
            auto ff = dynamic_cast<sofa::core::behavior::ForceField< sofa::defaulttype::Vec3Types > *>(gnode->forceField.get(i));
            ff->addForce(&mparams,dforces,positions.ref(),velocities.ref ());
        }
        std::cout << " | COMPUTE FORCES END>>" << std::endl;

        WriteDeriv forces = dforces;
        //std::cout << "Positions size : "<< positions.size() << std::endl;
        //std::cout << "Velocities size : "<< velocities.size() << std::endl;
        //std::cout << "Force size : "<< forces.size() << std::endl;

        std::cout << "<<INTEGRATION BEGIN | ...";
        for(uint i = 0; i < velocities.ref().size(); ++i)
        {
            velocities[i] = velocities[i] + dt * forces[i];
            positions[i] = positions[i] + dt * velocities[i];
        }
        std::cout << " | INTEGRATION END>>" << std::endl;
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;

    gnode->execute<UpdateMappingVisitor>(params);
    {
        UpdateMappingEndEvent ev ( dt );
        PropagateEventVisitor act ( params , &ev );
        gnode->execute ( act );
    }


}

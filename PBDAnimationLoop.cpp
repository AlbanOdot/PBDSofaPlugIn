#include "PBDAnimationLoop.hpp"
#include <sofa/core/ObjectFactory.h>
#include <SofaBaseMechanics/MechanicalObject.h>

using namespace sofa;
using namespace simulation;
using namespace defaulttype;
using namespace component;
using namespace container;

int PBDAnimationLoopClass = core::RegisterObject("Simulation loop to use in scene without constraints nor contact.")
                            .add< PBDAnimationLoop >()
                            .addDescription(R"(
                                            This loop do the following steps:
                                            - build and solve all linear systems in the scene : collision and time integration to compute the new values of the dofs
                                            - update the context (dt++)
                                            - update the mappings
                                            - update the bounding box (volume covering all objects of the scene))");

PBDAnimationLoop::PBDAnimationLoop(simulation::Node* _gnode)
    : Inherit()
    , gnode(_gnode)
{
    //assert(gnode);
}

PBDAnimationLoop::~PBDAnimationLoop()
{
}

void PBDAnimationLoop::init()
{
    if (!gnode)
        gnode = dynamic_cast<simulation::Node*>(this->getContext());
}

void PBDAnimationLoop::setNode( simulation::Node* n )
{
    gnode=n;
}

void PBDAnimationLoop::step(const core::ExecParams* params, SReal dt)
{
    if ( dt == 0 )
    {
        dt = this->gnode->getDt();
    }


    BaseContext* context = gnode->getContext ();
    //Get all of the mechanicals object in the scene
    std::vector<MechanicalObject<Vec3Types>*> mechanicalObjects = context->getObjects<MechanicalObject<Vec3Types>>(BaseContext::SearchDown);

    for(uint currentMechObj = 0; currentMechObj < mechanicalObjects.size (); ++currentMechObj)
    {
        //Get all of the positions and velocities
        std::vector<SReal> positions;
        std::vector<SReal> velocities;
        getVelocitiesAndPosition(mechanicalObjects[currentMechObj],velocities,positions);

        //Apply external forces on velocity only
        std::vector<SReal> extForces = {0.0,-9.8,0.0};
        //ATM there is only gravity = (0,-9.8,0)
        //computeUniformExternalForces(extForces);

        m_integrator.integrateUniformExternalForces(velocities,extForces,dt);
        //computeExternalForces(extForces);
        //m_integrator.integrate(velocities,forces,dt);

        //Temporary integration
        std::vector<SReal> tmpPos;
        m_integrator.integrateTmp(tmpPos,positions,velocities,dt);


        //Generate collision constraint
            //TODO ^^^^



        //Project Constraints
            //At the moment we create constraint on the fly because it's an alpha version of the alpha version

        //Save the fixed points coordinates
        //We set point 3 ,39 and 64 to be fixed this had to be done last
        std::vector<SReal> fixedPointPositions;
        std::vector<uint> fixedPointIndices = {3,39,64};
        for( auto index : fixedPointIndices)
        {
            fixedPointPositions.push_back (positions[3*index]);
            fixedPointPositions.push_back (positions[3*index+1]);
            fixedPointPositions.push_back (positions[3*index+2]);
        }

        //We set point 3 ,39 and 64 to be fixed this had to be done first
        m_integrator.solveFixedPointConstraint(tmpPos,positions,fixedPointIndices);

        //We define rigid constraint and solve constraint in O(n) -> F(a->b) = - F(b->a)
        m_integrator.solveDistanceConstraint(tmpPos,positions);

        //Final integration
        m_integrator.PBDUpdate(tmpPos,velocities,positions,dt);

        //Find a way to put the new values into the mechanical object


    }


}

void PBDAnimationLoop::getVelocitiesAndPosition(const MechanicalObject<Vec3Types> * mechanicalObject, std::vector<SReal>& velocities, std::vector<SReal>& positions)
{
    //Resize for parallele writting
    uint nbCoord = mechanicalObject->getSize ();
    positions.reserve (3*nbCoord);
    velocities.reserve (3*nbCoord);

    #pragma omp for
    for(uint currCoord = 0; currCoord < nbCoord; currCoord += 3)
    {
        //TODO Changer ca pour diviser par 3 le coût
        //Aller voir le code de getPX
        positions[currCoord] = mechanicalObject->getPX (currCoord);
        positions[currCoord+1] = mechanicalObject->getPY (currCoord+1);
        positions[currCoord+2] = mechanicalObject->getPZ (currCoord+2);

        velocities[currCoord] = mechanicalObject->getVX (currCoord);
        velocities[currCoord+1] = mechanicalObject->getVY (currCoord+1);
        velocities[currCoord+2] = mechanicalObject->getVZ (currCoord+2);

    }
}

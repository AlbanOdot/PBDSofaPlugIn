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

void PBDAnimationLoop::step(const sofa::core::ExecParams* params,
                            SReal dt)
{
    if (dt == 0)
    {
        dt = gnode->getDt();
    }
    sofa::core::MechanicalParams mparams(*params);

    for(auto& Obj : m_mechanicalObjects)
    {
        //Object parameters
        WriteCoord x = Obj->writePositions ();
        WriteDeriv v = Obj->writeVelocities ();
        ReadDeriv rest = Obj->readRestPositions ();
        uint pointCount = x.ref ().size();

        //External forces
        Vec3 zero;zero.set(0,0,0);
        VecDeriv zerosDeriv(x.ref().size(),zero);
        Derivatives dFext(zerosDeriv);

        //We will compute constrainst on p
        VecCoord freeVecCoord(x.ref());
        Coordinates freeCoord(freeVecCoord);
        WriteCoord p = freeCoord;

        //Apply external forces on p
        extForces (&mparams,dFext,p,x,v,dt);

        //Solve all of the constraints
        solveConstraints(rest,p);

        //Integrate using PBD method
        for(uint i = 0; i < pointCount; ++i)
        {
            v[i] = (p[i] - x[i])/dt;
            x[i] = p[i];
        }

    }


    gnode->execute<UpdateMappingVisitor>(params);
    {
        UpdateMappingEndEvent ev ( dt );
        PropagateEventVisitor act ( params , &ev );
        gnode->execute ( act );
    }


}


void PBDAnimationLoop::extForces (const sofa::core::MechanicalParams * mparams,
                                  Derivatives& f,
                                  WriteCoord& p,
                                  const WriteCoord& x,
                                  WriteDeriv& v,
                                  SReal dt)
{

    //Accumulates all external forces (Gravities,interactions etc)
    for(uint i = 0;  i < gnode->forceField.size(); ++i)
    {
        auto ff = dynamic_cast<sofa::core::behavior::ForceField< sofa::defaulttype::Vec3Types > *>(gnode->forceField.get(i));
        ff->addForce(mparams,f,x.ref(),v.ref ());
    }

    WriteDeriv Fext = f;
    //#pragma omp parallel for
    for(uint i = 0; i < v.ref().size(); ++i)
    {
        v[i] = v[i] + dt * Fext[i];
        p[i] = x[i] + dt * v[i];
    }

}

void PBDAnimationLoop::solveConstraints(ReadCoord& rest, WriteCoord& p)
{
    //From here we solve all of the constraints -> solve on p
    uint max_iter = 1;
    for(uint nbIter = 0; nbIter < max_iter; ++nbIter)
    {
        solveStretch(rest,p);

        solveFixedPoint(rest,p);
    }
}

void PBDAnimationLoop::solveStretch(ReadCoord& rest, WriteCoord& p)
{

    //Stretching constraint
    //O(0.5(n+1)n) : F(a->b) = -F(b->a)
    //Solving Gauss-Seidel's style
    //Use an accumulator to put it as Jacobi
    //TODO utiliser la topology
    uint pointCount =  p.ref().size();
    for( uint i = 0; i < pointCount; ++i)
    {
        for( uint j = i+1; j < pointCount; ++j)
        {
            const auto& p_ij = p[i] - p[j];
            SReal l = p_ij.norm();
            SReal d = (rest[i]-rest[j]).norm();
            const auto& displacement = 0.5*(l-d)/l;
            const auto& vdisplacement = displacement * p_ij;
            p[i] -= vdisplacement;
            p[j] += vdisplacement;
        }
    }
}

void PBDAnimationLoop::solveFixedPoint(ReadCoord& rest, WriteCoord& p)
{
    //Constrainst related stuff
    std::vector<uint> fixedIdx = {3};

    //Fixed point constraint
    for(const auto& idx : fixedIdx )
    {
        p[idx] = rest[idx];
    }
}

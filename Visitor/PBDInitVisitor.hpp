#ifndef PBDINITVISITOR_HPP
#define PBDINITVISITOR_HPP

#include <sofa/simulation/simulationcore.h>
#include <sofa/simulation/Visitor.h>
#include <sofa/simulation/Node.h>
#include <sofa/core/VecId.h>
#include <sofa/core/MultiVecId.h>
#include <sofa/core/ExecParams.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/core/BehaviorModel.h>
#include <sofa/core/behavior/BaseInteractionForceField.h>
#include <sofa/core/behavior/OdeSolver.h>
#include <sofa/core/behavior/BaseAnimationLoop.h>
#include <sofa/core/collision/Pipeline.h>

namespace sofa
{

namespace simulation
{

class PBDInitVisitor : public Visitor
{

protected :
    SReal dt;
    bool firstNodeVisited;
public:
    PBDInitVisitor(const core::ExecParams* params = core::ExecParams::defaultInstance());
    PBDInitVisitor(const core::ExecParams* params, SReal dt);

    void setDt(SReal v) { dt = v; }
    SReal getDt() const { return dt; }

    virtual void initCollisionPipeline(simulation::Node* node, core::collision::Pipeline* obj);
    virtual void processBehaviorModel(simulation::Node* node, core::BehaviorModel* obj);
    virtual void fwdInteractionForceField(simulation::Node* node, core::behavior::BaseInteractionForceField* obj);
    virtual void processOdeSolver(simulation::Node* node, core::behavior::OdeSolver* obj);

    virtual Result processNodeTopDown(simulation::Node* node);
    //virtual void processNodeBottomUp(simulation::Node* node);

    /// Specify whether this action can be parallelized.
    virtual bool isThreadSafe() const { return true; }

    /// Return a category name for this action.
    /// Only used for debugging / profiling purposes
    virtual const char* getCategoryName() const { return "animate"; }
    virtual const char* getClassName() const { return "AnimateVisitor"; }
};

} // namespace simulation

} // namespace sofa

#endif

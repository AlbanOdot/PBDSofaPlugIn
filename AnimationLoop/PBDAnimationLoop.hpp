#ifndef PBDANIMATIONLOOP_HPP
#define PBDANIMATIONLOOP_HPP

#include <sofa/simulation/CollisionAnimationLoop.h>
#include <sofa/simulation/Node.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/BaseConstraint.h>
#include "../Solver/PBDSolver.hpp"
#include "../InternalData/PBDObject.hpp"
#include "../Constraint/PBDBaseConstraint.hpp"

#include <sofa/simulation/PropagateEventVisitor.h>
#include <sofa/simulation/CollisionVisitor.h>
#include <sofa/simulation/CollisionBeginEvent.h>
#include <sofa/simulation/CollisionEndEvent.h>

#include <sofa/simulation/CollisionAnimationLoop.h>
#include <SofaConstraint/LCPConstraintSolver.h>

class PBDAnimationLoop : public sofa::simulation::CollisionAnimationLoop
{
    typedef sofa::core::objectmodel::Data<sofa::helper::vector<sofa::defaulttype::Vec3Types::Coord>> Coordinates;
    typedef sofa::core::objectmodel::Data<sofa::helper::vector<sofa::defaulttype::Rigid3Types::Coord>> RCoordinates;
    typedef sofa::helper::WriteAccessor <Coordinates> WriteCoord;
    typedef sofa::helper::WriteAccessor <RCoordinates> RWriteCoord;

    typedef sofa::helper::WriteAccessor <sofa::core::objectmodel::Data<sofa::helper::vector<sofa::defaulttype::Rigid3Types::Deriv>>> RWriteDeriv;
    typedef sofa::helper::WriteAccessor <sofa::core::objectmodel::Data<sofa::helper::vector<sofa::defaulttype::Vec3Types::Deriv>>> WriteDeriv;

public:
    PBDAnimationLoop(sofa::simulation::Node* gnode = NULL);
    virtual ~PBDAnimationLoop();
    typedef sofa::core::behavior::BaseAnimationLoop Inherit;
    typedef sofa::core::objectmodel::BaseContext BaseContext;
    typedef sofa::core::objectmodel::BaseObjectDescription BaseObjectDescription;
    SOFA_CLASS(PBDAnimationLoop,sofa::simulation::CollisionAnimationLoop);

    /// Set the simulation node this animation loop is controlling
    virtual void setNode( sofa::simulation::Node* );
    /*
     * Create and init all of the data needed to solve a defined constraint.
     */
    virtual void init() override;


    /*
     * Init function of sofa. It's called after the first init of the tree.
     */
    virtual void bwdInit () override;

    /// perform one animation step
    /*
     * Inputs : ExecParams *    -> Execution context
     *          SReal           -> Time Step
     *
     * Output : Compute a single step of the simulation
     */
    virtual void step(const sofa::core::ExecParams* params, SReal dt) override;

    int getNbIter() {return m_nbIter.getValue ();}
    /// Construction method called by ObjectFactory.
    template<class T>
    static typename T::SPtr create(T*, BaseContext* context, BaseObjectDescription* arg)
    {
        sofa::simulation::Node* node = dynamic_cast<sofa::simulation::Node*>(context);
        typename T::SPtr obj = sofa::core::objectmodel::New<T>(node);
        if (context) context->addObject(obj);
        if (arg) obj->parse(arg);
        return obj;
    }

protected :

    //Context and scene hierachy
    BaseContext* m_context;
    sofa::core::objectmodel::Data<bool> d_threadSafeVisitor;
    sofa::core::objectmodel::Data<bool> m_solveVelocityConstraintFirst;
    sofa::core::objectmodel::Data<SReal> f_rayleighMass;
    //Solvers
    PBDSolver m_solver;
    sofa::core::objectmodel::Data<int> m_nbIter;

    bool m_firststep;

    SReal m_minTime = 1e308;

    //Remove this later
    std::ofstream myfile;
    std::vector<SReal> e34,e629,e347,moyenne;
    bool fileok = false;

};


#endif //PBDANIMATIONLOOP_HPP

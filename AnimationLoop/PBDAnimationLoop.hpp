#ifndef PBDANIMATIONLOOP_HPP
#define PBDANIMATIONLOOP_HPP

#include <sofa/core/behavior/BaseAnimationLoop.h>
#include <sofa/simulation/Node.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/core/behavior/ForceField.h>
#include "../Solver/PBDExplicitIntegrator.hpp"
#include "../InternalData/PBDObject.hpp"
#include "../Constraint/PBDBaseConstraint.hpp"

#include <sofa/simulation/PropagateEventVisitor.h>
#include <sofa/simulation/CollisionVisitor.h>
#include <sofa/simulation/CollisionBeginEvent.h>
#include <sofa/simulation/CollisionEndEvent.h>

#include <sofa/simulation/CollisionAnimationLoop.h>
#include <SofaConstraint/LCPConstraintSolver.h>

class PBDAnimationLoop : public sofa::core::behavior::BaseAnimationLoop
{
    typedef sofa::defaulttype::Vec3Types::Coord       Coord;
    typedef sofa::helper::vector<Coord>               VecCoord;
    typedef sofa::core::objectmodel::Data<VecCoord>   Coordinates;
    typedef sofa::helper::WriteAccessor <Coordinates> WriteCoord;

    typedef sofa::defaulttype::Vec3Types::Deriv       Deriv;
    typedef sofa::helper::vector<Deriv>               VecDeriv;
    typedef sofa::core::objectmodel::Data<VecDeriv>   Derivatives;
    typedef sofa::helper::WriteAccessor <Derivatives> WriteDeriv;

    //Eigen
    typedef Eigen::Vector3f Vec;
    typedef Eigen::Vector4f Vec4;
    typedef std::vector<Vec> VecList;
    typedef std::vector<Vec4> Vec4List;
    typedef Eigen::Matrix3f Mat;
    typedef Eigen::Matrix4f Mat4;

public:
    PBDAnimationLoop(sofa::simulation::Node* gnode = NULL);
    virtual ~PBDAnimationLoop();
    typedef sofa::core::behavior::BaseAnimationLoop Inherit;
    typedef sofa::core::objectmodel::BaseContext BaseContext;
    typedef sofa::core::objectmodel::BaseObjectDescription BaseObjectDescription;
    SOFA_CLASS(PBDAnimationLoop,sofa::core::behavior::BaseAnimationLoop);

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

    /// Construction method called by ObjectFactory.
    template<class T>
    static typename T::SPtr create(T*, BaseContext* context, BaseObjectDescription* arg)
    {
        sofa::simulation::Node* gnode = dynamic_cast<sofa::simulation::Node*>(context);
        typename T::SPtr obj = sofa::core::objectmodel::New<T>(gnode);
        if (context) context->addObject(obj);
        if (arg) obj->parse(arg);
        return obj;
    }

protected :

    //Context and scene hierachy
    BaseContext* m_context;
    sofa::simulation::Node* gnode; ///< the node controlled by the loop

    //Solvers
    PBDExplicitIntegrator<sofa::defaulttype::Vec3Types> m_integrator_v;
    PBDExplicitIntegrator<sofa::defaulttype::RigidTypes> m_integrator_r;

    //Datas and transformations
    std::vector<PBDObject<sofa::defaulttype::Vec3Types>> m_objects_v;
    std::vector<PBDObject<sofa::defaulttype::RigidTypes>> m_objects_r;

    sofa::core::objectmodel::Data<int> m_nbIter;

    sofa::core::behavior::ConstraintSolver *constraintSolver;
    sofa::component::constraintset::LCPConstraintSolver::SPtr defaultSolver;

};
#endif //PBDANIMATIONLOOP_HPP

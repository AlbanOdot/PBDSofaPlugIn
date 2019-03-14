#ifndef PBDANIMATIONLOOP_HPP
#define PBDANIMATIONLOOP_HPP

#include <sofa/core/behavior/BaseAnimationLoop.h>
#include <sofa/simulation/Node.h>
#include <SofaBaseMechanics/MechanicalObject.h>

#include "solver/GaussSeidelSolver.hpp"
#include "solver/PBDExplicitIntegrator.hpp"

class PBDAnimationLoop : public sofa::core::behavior::BaseAnimationLoop
{
    typedef sofa::defaulttype::Vec3Types::Coord       Coord;
    typedef sofa::helper::vector<Coord>               VecCoord;
    typedef sofa::core::objectmodel::Data<VecCoord>   Coordinates;
    typedef sofa::helper::ReadAccessor  <Coordinates> ReadCoord;
    typedef sofa::helper::WriteAccessor <Coordinates> WriteCoord;
    typedef sofa::defaulttype::Vec3Types::Deriv       Deriv;
    typedef sofa::helper::vector<Deriv>               VecDeriv;
    typedef sofa::core::objectmodel::Data<VecDeriv>   Derivatives;
    typedef sofa::helper::ReadAccessor  <Derivatives> ReadDeriv;
    typedef sofa::helper::WriteAccessor <Derivatives> WriteDeriv;

protected:
    PBDAnimationLoop(sofa::simulation::Node* gnode = NULL);
    virtual ~PBDAnimationLoop();
public:
    typedef sofa::core::behavior::BaseAnimationLoop Inherit;
    typedef sofa::core::objectmodel::BaseContext BaseContext;
    typedef sofa::core::objectmodel::BaseObjectDescription BaseObjectDescription;
    SOFA_CLASS(PBDAnimationLoop,sofa::core::behavior::BaseAnimationLoop);

    /// Set the simulation node this animation loop is controlling
    virtual void setNode( sofa::simulation::Node* );

    /// Set the simulation node to the local context if not specified previously
    virtual void init() override;

//    virtual void bwdInit () override;

    /// perform one animation step
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
    //PBDExplicitIntegrator m_integrator;
    sofa::core::behavior::OdeSolver * m_solver;
    BaseContext* m_context;
    sofa::component::container::MechanicalObject< sofa::defaulttype::Vec3Types > * m_mechanicalObject;
    ReadCoord* m_restPositions;
    WriteCoord* m_freePosition;
    sofa::simulation::Node* gnode; ///< the node controlled by the loop
};
#endif //PBDANIMATIONLOOP_HPP

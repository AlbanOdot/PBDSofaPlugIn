#ifndef PBDANIMATIONLOOP_HPP
#define PBDANIMATIONLOOP_HPP

#include <sofa/core/behavior/BaseAnimationLoop.h>

class PBDAnimationLoop : public sofa::core::behavior::BaseAnimationLoop
{
public:
    typedef sofa::core::behavior::BaseAnimationLoop Inherit;
    typedef sofa::core::objectmodel::BaseContext BaseContext;
    typedef sofa::core::objectmodel::BaseObjectDescription BaseObjectDescription;
    SOFA_CLASS(PBDAnimationLoop,sofa::core::behavior::BaseAnimationLoop);
protected:
    PBDAnimationLoop(simulation::Node* gnode = NULL);

    virtual ~PBDAnimationLoop();
public:
    /// Set the simulation node this animation loop is controlling
    virtual void setNode( simulation::Node* );

    /// Set the simulation node to the local context if not specified previously
    virtual void init() override;

    /// perform one animation step
    virtual void step(const sofa::core::ExecParams* params, SReal dt) override;


    /// Construction method called by ObjectFactory.
    template<class T>
    static typename T::SPtr create(T*, BaseContext* context, BaseObjectDescription* arg)
    {
        simulation::Node* gnode = dynamic_cast<simulation::Node*>(context);
        typename T::SPtr obj = sofa::core::objectmodel::New<T>(gnode);
        if (context) context->addObject(obj);
        if (arg) obj->parse(arg);
        return obj;
    }

protected :

    simulation::Node* gnode;  ///< the node controlled by the loop

}

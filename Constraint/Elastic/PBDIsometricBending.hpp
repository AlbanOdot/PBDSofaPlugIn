#ifndef PBDISOMETRICBENDING_HPP
#define PBDISOMETRICBENDING_HPP

#include "PBDBending.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>

class PBDIsometricBending : public PBDBending
{
public:
    PBDIsometricBending(sofa::simulation::Node* gnode = NULL): PBDBending(){}
    /*
     * Inputs : PBDObject   -> Object on wich we will solve the constraint
     *          WriteCoord  -> Free positions on wich we apply the dispalcement
     *
     * Output : Solve the constraint adding in WriteCoord the computed displacement
     */
    virtual bool solve(sofa::simulation::Node * node);

    virtual void bwdInit() override;

    /// Construction method called by ObjectFactory.
    template<class T>
    static typename T::SPtr create(T*, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg)
    {
        typename T::SPtr obj = sofa::core::objectmodel::New<T>();
        if (context) context->addObject(obj);
        if (arg) obj->parse(arg);
        return obj;
    }

};

#endif // PBDISOOMETRICBENDING_HPP

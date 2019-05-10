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
    virtual void solve(PBDObject& object, WriteCoord& p);

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

#ifndef PBDBENDING_HPP
#define PBDBENDING_HPP

#include "PBDElasticConstraint.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>

class PBDBending : public PBDElasticConstraint
{
public:
    PBDBending(sofa::simulation::Node* gnode = NULL) : PBDElasticConstraint(){}
    virtual void solve(PBDObject& object, WriteCoord& p);
    virtual void bwdInit () override;
    /// Construction method called by ObjectFactory.
    template<class T>
    static typename T::SPtr create(T*, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg)
    {
        typename T::SPtr obj = sofa::core::objectmodel::New<T>();
        if (context) context->addObject(obj);
        if (arg) obj->parse(arg);
        return obj;
    }

protected:
    void correction(PBDObject &object, uint a, uint b, WriteCoord& x);
protected:
    SReal coeff;
};

#endif // PBDBENDING_HPP

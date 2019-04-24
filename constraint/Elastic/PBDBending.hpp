#ifndef PBDBENDING_HPP
#define PBDBENDING_HPP

#include "PBDElasticConstraint.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>

class PBDBending : public PBDElasticConstraint
{
public:
    PBDBending(sofa::simulation::Node* gnode = NULL) : PBDElasticConstraint(),
    alpha_wann(initData(&alpha_wann,(SReal)1e-2,"a1","Low frequency ondulation")),
    alpha_too(initData(&alpha_too,(SReal)4e-4,"a2","High frequency ondulation")){}
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
    void correction(PBDObject &object, uint a, uint b, WriteCoord& x, const ReadDeriv& vel);
protected:
    sofa::core::objectmodel::Data<SReal> alpha_wann;
    sofa::core::objectmodel::Data<SReal> alpha_too;
    SReal coeff;
    SReal m_K;
};

#endif // PBDBENDING_HPP

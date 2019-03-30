#ifndef PBDISOMETRICBENDING_HPP
#define PBDISOMETRICBENDING_HPP

#include "PBDBaseConstraint.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>

class PBDIsometricBending : public PBDBaseConstraint
{
public:
    PBDIsometricBending(sofa::simulation::Node* gnode = NULL);
    PBDIsometricBending(uint objectSize);
    virtual Matrix* getConstraintMatrix();
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
protected:
    sofa::component::linearsolver::FullMatrix<float> m_constraint;
    sofa::core::objectmodel::Data<SReal> m_k;
    SReal dt2;
};

#endif // PBDSTRECH_HPP

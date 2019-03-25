#ifndef PBDFIXEDPOINT_HPP
#define PBDFIXEDPOINT_HPP
#include "PBDBaseConstraint.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>

class PBDFixedPoint : public PBDBaseConstraint
{
public:

    PBDFixedPoint(sofa::simulation::Node* gnode = NULL){}
    PBDFixedPoint(uint objectSize);
    virtual Matrix* getConstraintMatrix();
    virtual void solve(const PBDObject& object, WriteCoord& p);

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
};

#endif // PBDFIXEDPOINT_HPP

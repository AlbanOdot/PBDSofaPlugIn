#ifndef PBDSTRECH_HPP
#define PBDSTRECH_HPP

#include "PBDElasticConstraint.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>

class PBDStretch : public PBDElasticConstraint
{
public:
    PBDStretch(sofa::simulation::Node* gnode = NULL) : PBDElasticConstraint(){}
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
    void correction (uint i, const std::pair<uint,double>& voisin, WriteCoord &p);
protected:
    double m_K;
};

#endif // PBDSTRECH_HPP

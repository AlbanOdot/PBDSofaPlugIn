#ifndef PBDSTRAINSHAPE_HPP
#define PBDSTRAINSHAPE_HPP


#include "PBDFEMConstraint.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>

/**
 * @brief The PBDStrainShape class
 *  This class implement https://sci-hub.tw/10.1016/j.cag.2014.07.004
*/

class PBDStrainShape : public PBDFEMConstraint
{
public:
    PBDStrainShape(sofa::simulation::Node* gnode = NULL):PBDFEMConstraint(){}
    virtual void solve(PBDObject& object, WriteCoord& p) override;
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
};

#endif // PBDSTRAINSHAPE_HPP

#ifndef PBDSTRAINSHAPE_HPP
#define PBDSTRAINSHAPE_HPP


#include "PBDBaseConstraint.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>

/**
 * @brief The PBDStrainShape class
 *  This class implement https://sci-hub.tw/10.1016/j.cag.2014.07.004
*/





class PBDStrainShape : public PBDBaseConstraint
{
public:
    PBDStrainShape(sofa::simulation::Node* gnode = NULL);
    PBDStrainShape(uint objectSize);
    virtual Matrix* getConstraintMatrix();
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
    sofa::component::linearsolver::FullMatrix<float> m_constraint;
    sofa::core::objectmodel::Data<SReal> m_young_modulus;
    sofa::core::objectmodel::Data<SReal> m_poisson_ratio;
    Eigen::Matrix3d m_C;//<<Elasticity Tensor
    SReal dt2;
};

#endif // PBDSTRAINSHAPE_HPP

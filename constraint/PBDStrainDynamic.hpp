#ifndef PBDSTRAINDYNAMIC_HPP
#define PBDSTRAINDYNAMIC_HPP


#include "PBDBaseConstraint.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>

/**
 * @brief The PBDStrainDynamic class
 *  This class implement http://matthias-mueller-fischer.ch/publications/strainBasedDynamics.pdf
*/

class PBDStrainDynamic : public PBDBaseConstraint
{
public:
    PBDStrainDynamic(sofa::simulation::Node* gnode = NULL);
    PBDStrainDynamic(uint objectSize);
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

#endif // PBDSTRAINDynamic_HPP

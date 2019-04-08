#ifndef PBDFEMCONSTRAINT_HPP
#define PBDFEMCONSTRAINT_HPP

#include "../PBDBaseConstraint.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>

class PBDFEMConstraint : public PBDBaseConstraint
{
public:
    PBDFEMConstraint(): PBDBaseConstraint(true),
        m_young_modulus(initData(&m_young_modulus,1.7,"young","Young modulus")), //Polypropylene (1.5-2)
        m_poisson_ratio(initData(&m_poisson_ratio,0.45,"poisson","Poisson ratio")) //Polymère (0.3-0.5) -> 0.45 c'ets du plexiglas)
    {}
protected:
    sofa::core::objectmodel::Data<SReal> m_young_modulus;
    sofa::core::objectmodel::Data<SReal> m_poisson_ratio;
    Eigen::Matrix3d m_C;//<<Elasticity Tensor
    SReal dt2;
};

#endif // PBDFEMCONSTRAINT_HPP

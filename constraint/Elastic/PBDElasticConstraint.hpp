#ifndef PBDELASTICCONSTRAINT_HPP
#define PBDELASTICCONSTRAINT_HPP

#include "../PBDBaseConstraint.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>

class PBDElasticConstraint : public PBDBaseConstraint
{
public:
    PBDElasticConstraint(sofa::simulation::Node* gnode = NULL)
        : PBDBaseConstraint(true),
          m_k(initData(&m_k,(SReal)1.0,"resilience","Resilience factor : bending, stretching etc...")){}
    virtual void bwdInit () override {}
protected:
    sofa::core::objectmodel::Data<SReal> m_k;
};

#endif // PBDELASTICCONSTRAINT_HPP

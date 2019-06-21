#ifndef PBDELASTICCONSTRAINT_HPP
#define PBDELASTICCONSTRAINT_HPP

#include "../PBDBaseConstraint.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>

class PBDElasticConstraint : public PBDBaseConstraint<sofa::defaulttype::Vec3Types>
{
public:
    PBDElasticConstraint(sofa::simulation::Node* gnode = NULL)
        : PBDBaseConstraint(),
          m_k(initData(&m_k,(SReal)1.0,"resilience","Resilience factor : bending, stretching etc...")){}
    /*
     * Init function of sofa. It's called after the first init of the tree.
     */
    virtual void bwdInit () override {}
protected:
    sofa::core::objectmodel::Data<SReal> m_k;
};

#endif // PBDELASTICCONSTRAINT_HPP

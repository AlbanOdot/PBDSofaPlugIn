#ifndef PBDELASTICCONSTRAINT_HPP
#define PBDELASTICCONSTRAINT_HPP

#include "../PBDBaseConstraint.hpp"
#include "../../InternalData/Data/PBDVertexMass.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>

class PBDElasticConstraint : public PBDConstraint<sofa::defaulttype::Vec3Types>
{
public:
    PBDElasticConstraint(sofa::simulation::Node* gnode = NULL)
        : PBDConstraint(),
          m_k(initData(&m_k,(SReal)1.0,"resilience","Resilience factor : bending, stretching etc...")){}

protected:
    PBDVertexMass<sofa::defaulttype::Vec3Types> m_mass;
    sofa::core::objectmodel::Data<SReal> m_k;
};

#endif // PBDELASTICCONSTRAINT_HPP

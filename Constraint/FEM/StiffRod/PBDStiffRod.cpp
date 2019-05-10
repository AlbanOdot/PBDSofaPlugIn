#include "PBDStiffRod.hpp"
#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>

int PBDStiffRodClass = sofa::core::RegisterObject("Constraint that correct beam.")
                       .add< PBDStiffRod >();

void PBDStiffRod::bwdInit()
{

}

void PBDStiffRod::solve(PBDObject &object, WriteCoord &p)
{

}

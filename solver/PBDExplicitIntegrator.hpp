#ifndef PBDEXPLICITEINTEGRATOR_HPP
#define PBDEXPLICITEINTEGRATOR_HPP

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/behavior/MultiVec.h>
#include <sofa/core/behavior/MultiMatrix.h>
#include "GaussSeidelSolver.hpp"

namespace sofa{
namespace core{
class PBDExplicitIntegrator : public virtual sofa::core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(PBDExplicitIntegrator, sofa::core::objectmodel::BaseObject);
    PBDExplicitIntegrator();

    virtual ~PBDExplicitIntegrator();

private:
    PBDExplicitIntegrator(const PBDExplicitIntegrator& n) ;
    PBDExplicitIntegrator& operator=(const PBDExplicitIntegrator& n) ;

public:
    // Main computation method.
    // Compute one step :  x = x + dt * dx
    void integrate(std::vector<SReal>& x, const std::vector<SReal>& dx, const SReal dt);

    //Compute one step in tmpposition :  tmpposition = position + dt * velocity
    void integrateTmp(std::vector<SReal>& tmpposition,const std::vector<SReal>& position, const std::vector<SReal>& velocities, const SReal dt);

    //Compute new velocity from external forces : velocity = velocity + dt * Forces_ext
    void integrateUniformExternalForces(std::vector<SReal>& velocities, const std::vector<SReal>& extForces, const SReal dt);

    //Call a solver (jk but will soon) to compute the constraints resolution
    // Solve only for the concerned points given by pointIdx (if void solves for all)
    void solveDistanceConstraint(std::vector<SReal>& position, const std::vector<SReal>& truth, const std::vector<uint>& pointIdx = std::vector<uint>());

    //Call a solver (jk but will soon) to compute the constraints resolution
    //Twice as fast as solving for all with a non empty pointIdx vector
    void solveDistanceConstraintAll(std::vector<SReal>& position, const std::vector<SReal>& truth, const std::vector<uint>& pointIdx);

    //Call a solver (jk but will soon) to compute the constraints resolution
    //Solve only for the concerned points given by pointIdx (if void solves for all)
    void solveFixedPointConstraint(std::vector<SReal>& position, const std::vector<SReal>& truth, const std::vector<uint>& pointIdx = std::vector<uint>());

    //Use the computed new position to compute set the enw velocity and position
    void PBDUpdate(const std::vector<SReal>& newPosition, std::vector<SReal>& velocity, std::vector<SReal>& position, const SReal dt);
protected:
    GaussSeidelSolver m_solver;
};
} //core namespace
} //sofa namespace
#endif //PBDEXPLICITINTEGRATOR_HPP

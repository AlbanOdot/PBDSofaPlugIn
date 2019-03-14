#ifndef PBDEXPLICITEINTEGRATOR_HPP
#define PBDEXPLICITEINTEGRATOR_HPP

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/behavior/MultiVec.h>
#include <sofa/core/behavior/MultiMatrix.h>
#include <sofa/defaulttype/VecTypes.h>

class PBDExplicitIntegrator : public virtual sofa::core::objectmodel::BaseObject
{
    typedef sofa::defaulttype::Vec3Types::Coord VecCoord;
    typedef sofa::defaulttype::Vec3Types::Deriv VecDeriv;
    typedef sofa::core::objectmodel::Data<VecCoord> Coordinates;
    typedef sofa::core::objectmodel::Data<VecDeriv> Derivatives;

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
    void integrate(Coordinates& x,
                   const Derivatives& dx,
                   const SReal dt);

    //Compute one step in tmpposition :  tmpposition = position + dt * velocity
    void integrateTmp(Coordinates& tmpposition,
                      const Coordinates& position,
                      const Derivatives& velocities,
                      SReal dt);

    //Compute new velocity from external forces : velocity = velocity + dt * Forces_ext
    void integrateUniformExternalForces(Derivatives& velocities,
                                        const Derivatives& extForces,
                                        const SReal dt);

    //Call a solver (jk but will soon) to compute the constraints resolution
    // Solve only for the concerned points given by pointIdx (if void solves for all)
    void solveDistanceConstraint(Coordinates& position,
                                 const Coordinates& truth,
                                 const std::vector<uint>& pointIdx = std::vector<uint>());

    //Call a solver (jk but will soon) to compute the constraints resolution
    //Twice as fast as solving for all with a non empty pointIdx vector
    void solveDistanceConstraintAll(Coordinates& position,
                                    const Coordinates& truth,
                                    const std::vector<uint>& pointIdx);

    //Call a solver (jk but will soon) to compute the constraints resolution
    //Solve only for the concerned points given by pointIdx (if void solves for all)
    void solveFixedPointConstraint(Coordinates& position,
                                   const Coordinates& truth,
                                   const std::vector<uint>& pointIdx = std::vector<uint>());

    //Use the computed new position to compute set the enw velocity and position
    void PBDUpdate(const Coordinates& newPosition,
                   Derivatives& velocity,
                   Coordinates& position,
                   const SReal dt);
protected:
    //GaussSeidelSolver m_solver;
};
#endif //PBDEXPLICITINTEGRATOR_HPP

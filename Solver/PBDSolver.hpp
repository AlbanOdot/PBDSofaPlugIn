#ifndef PBDSOLVER_H
#define PBDSOLVER_H
#include "config.h"

#include <sofa/core/behavior/OdeSolver.h>
#include <sofa/simulation/Node.h>
#include "../Constraint/PBDBaseConstraint.hpp"
#include <unordered_map>

/** Explicit time integrator using position difference
 */
using namespace sofa::defaulttype;
using namespace sofa::core::objectmodel;
using namespace sofa::component::container;
class PBDSolver
{
    typedef sofa::defaulttype::Vec3 Vec3;
    typedef sofa::helper::ReadAccessor  <sofa::core::objectmodel::Data<sofa::helper::vector<Vec3Types::Coord>>>       ReadCoord;
    typedef sofa::helper::WriteAccessor <sofa::core::objectmodel::Data<sofa::helper::vector<Vec3Types::Coord>>>       WriteCoord;
    typedef sofa::helper::WriteAccessor <sofa::core::objectmodel::Data<sofa::helper::vector<Vec3Types::Deriv>>>       WriteDeriv;
    typedef sofa::helper::ReadAccessor  <sofa::core::objectmodel::Data<sofa::helper::vector<RigidTypes::Coord>>>      ReadCoordR;
    typedef sofa::helper::WriteAccessor <sofa::core::objectmodel::Data<sofa::helper::vector<RigidTypes::Coord>>>      WriteCoordR;
    typedef sofa::helper::WriteAccessor <sofa::core::objectmodel::Data<sofa::helper::vector<RigidTypes::Deriv>>>      WriteDerivR;

public:
    PBDSolver(){}
    void integrate(SReal dt);
    void setupOrientations(SReal dt);
    void setupSolver(sofa::simulation::Node* node, int nbIter, SReal rayleigh);
    void solvePBDConstraints (const sofa::core::ExecParams* params);
    void solveSofaConstraints(const sofa::core::ExecParams* params, SReal dt, sofa::core::MultiVecCoordId freePosId, sofa::core::MultiVecDerivId freeVelId);
protected:
    static inline void computeVec3Integration(std::unordered_map<MechanicalObject<Vec3Types> * ,PBDObject<Vec3Types> *>& mechanicalObjects, SReal damping_times_inv_dt);
    static inline void computeRigidIntegration(sofa::simulation::Node* node, SReal damping_times_inv_dt);
private:

    SReal m_damping; ///< Rayleigh damping coefficient related to mass
    SReal m_damping_times_inv_dt;
    sofa::core::objectmodel::Data<bool> d_threadSafeVisitor;
    sofa::simulation::Node * m_node;
    std::vector<PBDBaseConstraint * > m_constraint;
    std::unordered_map<MechanicalObject<Vec3Types> * ,PBDObject<Vec3Types> *> m_v3m;
    std::unordered_map<MechanicalObject<RigidTypes> * ,PBDObject<RigidTypes> *> m_r3m;
    int m_nbIter;
};

#endif

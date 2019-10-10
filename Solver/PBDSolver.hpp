#ifndef PBDSOLVER_H
#define PBDSOLVER_H

#include "../Constraint/Other/PBDCollisionConstraint.hpp"
#include "../Constraint/PBDBaseConstraint.hpp"

#include <SofaTopologyMapping/Tetra2TriangleTopologicalMapping.h>
#include <sofa/core/CollisionModel.h>
#include <sofa/helper/map_ptr_stable_compare.h>
#include <sofa/core/behavior/OdeSolver.h>
#include <sofa/simulation/Node.h>
#include <unordered_map>
#include <sofa/core/objectmodel/BaseContext.h>
/**
 * Explicit time integrator using position difference
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
    void updateFreePositionsAndVelocities();
    void generateCollisions();
    bool generateCollisionConstraint(const std::pair< sofa::core::CollisionModel*, sofa::core::CollisionModel* >& collModels, sofa::core::collision::DetectionOutputVector* outputs);
    void clearCollisions() {for(auto& c : m_collision_constraints){delete c;}m_collision_constraints.clear ();}//delete and erase the vector content

protected:
    static inline void computeVec3Integration(std::unordered_map<MechanicalObject<Vec3Types> * ,PBDObject<Vec3Types> *>& mechanicalObjects, SReal damping_times_inv_dt);
    static inline void computeRigidIntegration(sofa::simulation::Node* node, SReal damping_times_inv_dt);
    static inline void collisionDetection(sofa::helper::vector<sofa::core::CollisionModel*>& cm);
    void findObjectInMaps(PBDObject<Rigid3Types>** pbdObj1r = nullptr,
                          PBDObject<Vec3Types>** pbdObj1v = nullptr,
                          BaseContext * cont = nullptr);
private:

    //Integration data
    SReal m_damping; ///< Rayleigh damping coefficient related to mass
    SReal m_damping_times_inv_dt;
    int m_nbIter;

    //Contexte related data
    sofa::simulation::Node * m_node;
    std::vector<PBDBaseConstraint * > m_constraints;

    //Mapping between sofa objects and PBD conception
    std::unordered_map<MechanicalObject<Vec3Types> * ,PBDObject<Vec3Types> *> m_v3m;
    std::unordered_map<MechanicalObject<RigidTypes> * ,PBDObject<RigidTypes> *> m_r3m;
    //Collision related data
    sofa::core::collision::NarrowPhaseDetection* m_narrowPhase;
    std::vector<PBDBaseCollisionConstraint *> m_collision_constraints;
};

#endif

#ifndef PBDCOLLISIONHANDLER_HPP
#define PBDCOLLISIONHANDLER_HPP
#include <SofaBaseCollision/ContactListener.h>
#include <SofaTopologyMapping/Tetra2TriangleTopologicalMapping.h>
#include <sofa/core/CollisionModel.h>
#include <sofa/simulation/Node.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/helper/map_ptr_stable_compare.h>
#include <SofaBaseTopology/EdgeSetTopologyContainer.h>
#include <SofaBaseTopology/PointSetTopologyContainer.h>
#include <SofaBaseTopology/TriangleSetTopologyContainer.h>

#include "../Solver/PBDSolver.hpp"

//This class only listen to the collision pipeline and tell the solver to create and solve the collision
class PBDCollisionHandler
{
public:
    void setPBDSolver(PBDSolver * s) { m_solver = s;}
    void handleCollision();
    void init();
protected:
    PBDCollisionHandler(  sofa::simulation::Node* n);
    virtual ~PBDCollisionHandler();

private:
    sofa::helper::vector<const sofa::helper::vector<sofa::core::collision::DetectionOutput>* > m_contactsVector;
    sofa::core::collision::NarrowPhaseDetection* m_narrowPhase;
    PBDSolver * m_solver;
    sofa::simulation::Node * m_node;
//    sofa::Data<sofa::component::topology::TriangleSetTopologyContainer *> m_trianglesContainer;
//    sofa::Data<sofa::component::topology::EdgeSetTopologyContainer *> m_edgesContainer;
//    sofa::Data<sofa::component::topology::PointSetTopologyContainer *> m_pointContainer;
};

#endif

#ifndef PBDCOLLISION_CONSTRAINT
#define PBDCOLLISION_CONSTRAINT
#include "../PBDBaseConstraint.hpp"
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/collision/DetectionOutput.h>
#include <SofaTopologyMapping/Tetra2TriangleTopologicalMapping.h>

class PBDBaseCollisionConstraint : public PBDBaseConstraint
{
public:
    PBDBaseCollisionConstraint(sofa::core::CollisionModel* colMod1 = nullptr,sofa::core::CollisionModel * colMod2 = nullptr,
                               sofa::core::collision::DetectionOutputVector * collData = nullptr,
                               uint nbIter =1) :
        PBDBaseConstraint(), m_colMod1(colMod1),m_colMod2(colMod2), m_detection(collData)
    {
        m_nbIter = nbIter;
        w1 = m_colMod1->isSimulated () ? 1.0 : 0.0;
        w2 = m_colMod2->isSimulated () ? 1.0 : 0.0;
    }
protected:
    SReal w1, w2;
    sofa::core::CollisionModel * m_colMod1;
    sofa::core::CollisionModel * m_colMod2;
    sofa::core::collision::DetectionOutputVector * m_detection;
};

template< class T, class U  >
class PBDCollisionConstraint : public PBDBaseCollisionConstraint
{
public:
    PBDCollisionConstraint(PBDObject<T>* obj1 = nullptr,PBDObject<U>* obj2 = nullptr,
                           sofa::core::CollisionModel* colMod1 = nullptr,sofa::core::CollisionModel * colMod2 = nullptr,
                           sofa::core::collision::DetectionOutputVector * collData = nullptr,
                           uint nbIter=1) :
        PBDBaseCollisionConstraint(colMod1,colMod2,collData,nbIter), m_obj1(obj1),m_obj2(obj2){}
    /*
     * Output : Solve the constraint adding in WriteCoord the computed displacement
     */
    virtual bool solve(sofa::simulation::Node * node);

private:
    void solvePP(const sofa::core::collision::DetectionOutput &);
    void solveEE(const sofa::core::collision::DetectionOutput &);
    void solveTT(const sofa::core::collision::DetectionOutput &);
    void solvePE(const sofa::core::collision::DetectionOutput &,bool swap = false);
    void solvePT(const sofa::core::collision::DetectionOutput &,bool swap = false);
    void solveET(const sofa::core::collision::DetectionOutput &,bool swap = false);
protected:
    PBDObject<T> * m_obj1;
    PBDObject<U> * m_obj2;

};
template class PBDCollisionConstraint<sofa::defaulttype::Vec3Types,sofa::defaulttype::Vec3Types>;
template class PBDCollisionConstraint<sofa::defaulttype::Vec3Types,sofa::defaulttype::Rigid3Types>;
template class PBDCollisionConstraint<sofa::defaulttype::Rigid3Types,sofa::defaulttype::Vec3Types>;
template class PBDCollisionConstraint<sofa::defaulttype::Rigid3Types,sofa::defaulttype::Rigid3Types>;
#endif

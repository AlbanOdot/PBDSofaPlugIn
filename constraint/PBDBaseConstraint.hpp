#ifndef PBDBASECONSTRAINT_HPP
#define PBDBASECONSTRAINT_HPP

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Data.h>
#include "../object/PBDObject.hpp"

class SOFA_CORE_API PBDBaseConstraint : public virtual sofa::core::objectmodel::BaseObject
{

public:
    SOFA_ABSTRACT_CLASS(PBDBaseConstraint, sofa::core::objectmodel::BaseObject);
    //SOFA_BASE_CAST_IMPLEMENTATION(PBDBaseConstraint)
    typedef sofa::defaulttype::Vec3Types::Coord       Coord;
    typedef sofa::helper::vector<Coord>               VecCoord;
    typedef sofa::core::objectmodel::Data<VecCoord>   Coordinates;
    typedef sofa::helper::ReadAccessor  <Coordinates> ReadCoord;
    typedef sofa::helper::WriteAccessor <Coordinates> WriteCoord;

    typedef sofa::defaulttype::Vec3Types::Deriv       Deriv;
    typedef sofa::helper::vector<Deriv>               VecDeriv;
    typedef sofa::core::objectmodel::Data<VecDeriv>   Derivatives;
    typedef sofa::helper::ReadAccessor  <Derivatives> ReadDeriv;
    typedef sofa::helper::WriteAccessor <Derivatives> WriteDeriv;

    typedef sofa::core::objectmodel::Data<sofa::helper::vector<uint>> IndexSet;
    typedef sofa::defaulttype::BaseMatrix Matrix;
protected:
    PBDBaseConstraint(bool optisolv = false)
        : m_indices(initData(&m_indices, sofa::helper::vector<uint>(), "indices", "ID of the vertices on wich this constraint is to apply")),
          m_nbIter(initData(&m_nbIter,(uint)1,"iter","Number of iteration for the solver")),
          m_hasOptimizedSolver(optisolv){}

    virtual ~PBDBaseConstraint() { }

public:

    inline bool hasOptiSolver() { return m_hasOptimizedSolver;}

    virtual Matrix * getConstraintMatrix() = 0;

    virtual void solve(PBDObject& object, WriteCoord& p) = 0;

public:
    bool m_hasOptimizedSolver;
    IndexSet m_indices; ///< Indices on wich to apply the constraint
    sofa::core::objectmodel::Data<unsigned int> m_nbIter;

};
#endif

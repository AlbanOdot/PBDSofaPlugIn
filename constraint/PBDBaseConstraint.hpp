#ifndef PBDANIMATIONLOOP_HPP
#define PBDANIMATIONLOOP_HPP

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Data.h>

class SOFA_CORE_API PBDBaseConstraint : public virtual sofa::core::objectmodel::BaseObject
{
    sofa::core::objectmodel::Data<std::vector<uint>> IndexSet;

public:
    SOFA_ABSTRACT_CLASS(PBDBaseConstraint, sofa::core::objectmodel::BaseObject);
    SOFA_BASE_CAST_IMPLEMENTATION(PBDBaseConstraint)

protected:
    PBDBaseConstraint(bool optisolv = false)
        : m_indices(initData(&m_indices, 0, "indices", "ID of the vertices on wich this constraint is to apply")),
          m_hasOptimizedSolver(optisolv){}

    virtual ~PBDBaseConstraint() { }

public:

    inline bool hasOptiSolver() { return m_hasOptimizedSolver;}

    /// Construct the Jacobian Matrix
    ///
    /// \param cId is the result constraint sparse matrix Id
    /// \param cIndex is the index of the next constraint equation: when building the constraint matrix, you have to use this index, and then update it
    /// \param cParams defines the state vectors to use for positions and velocities. Also defines the order of the constraint (POS, VEL, ACC)
    template<typename T>
    virtual void getConstraintMatrix(T& M) = 0;

    virtual void solve() = 0;

public:
    bool m_hasOptimizedSolver;
    IndexSet m_indices; ///< Indices on wich to apply the constraint

};
#endif

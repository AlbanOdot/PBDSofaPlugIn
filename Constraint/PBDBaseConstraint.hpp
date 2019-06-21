#ifndef PBDBASECONSTRAINT_HPP
#define PBDBASECONSTRAINT_HPP

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Data.h>
#include "../InternalData/PBDObject.hpp"

template < class T >
class SOFA_CORE_API PBDBaseConstraint : public virtual sofa::core::objectmodel::BaseObject
{

public:
    SOFA_ABSTRACT_CLASS(PBDBaseConstraint, sofa::core::objectmodel::BaseObject);
    //SOFA_BASE_CAST_IMPLEMENTATION(PBDBaseConstraint)
    typedef typename T::Coord       Coord;
    typedef sofa::helper::vector<Coord>               VecCoord;
    typedef sofa::core::objectmodel::Data<VecCoord>   Coordinates;
    typedef sofa::helper::ReadAccessor  <Coordinates> ReadCoord;
    typedef sofa::helper::WriteAccessor <Coordinates> WriteCoord;

    typedef typename T::Deriv       Deriv;
    typedef sofa::helper::vector<Deriv>               VecDeriv;
    typedef sofa::core::objectmodel::Data<VecDeriv>   Derivatives;
    typedef sofa::helper::ReadAccessor  <Derivatives> ReadDeriv;
    typedef sofa::helper::WriteAccessor <Derivatives> WriteDeriv;

    typedef sofa::core::objectmodel::Data<sofa::helper::vector<uint>> IndexSet;
    typedef sofa::defaulttype::BaseMatrix Matrix;
protected:
    PBDBaseConstraint()
        : m_indices(initData(&m_indices, sofa::helper::vector<uint>(), "indices", "ID of the vertices on wich this constraint is to apply")),
          m_nbIter(initData(&m_nbIter,(uint)1,"iter","Number of iteration for the solver")){}

    virtual ~PBDBaseConstraint() { }

public:
    /*
     * Inputs : PBDObject   -> Object on wich we will solve the constraint
     *          WriteCoord  -> Free positions on wich we apply the dispalcement
     *
     * Output : Solve the constraint adding in WriteCoord the computed displacement
     */
    virtual void solve(PBDObject<T>& object, WriteCoord& p) = 0;

    /*
     * Inputs : int -> Number of iterations
     *
     * Some constraints can have their own number of iterations. Do not use unless you know what you're doing.
     * This will most likely lead to unstabilities.
     */
    void setIterCount(int c) { m_nbIter.setValue (c);}

public:
    IndexSet m_indices; ///< Indices on wich to apply the constraint
    sofa::core::objectmodel::Data<unsigned int> m_nbIter;

};
#endif

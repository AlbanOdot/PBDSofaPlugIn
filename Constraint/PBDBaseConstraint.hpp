#ifndef PBDBASECONSTRAINT_HPP
#define PBDBASECONSTRAINT_HPP

#include <sofa/core/behavior/OdeSolver.h>
#include <sofa/core/objectmodel/Data.h>
#include "../InternalData/PBDObject.hpp"
#include <sofa/core/objectmodel/Link.h>


class PBDBaseConstraint : public virtual sofa::core::objectmodel::BaseObject
{

protected:
    typedef sofa::core::objectmodel::Data<sofa::helper::vector<uint>> IndexSet;
    typedef sofa::helper::vector<sofa::defaulttype::Vec3Types::Coord> VecCoord;
    typedef sofa::helper::vector<sofa::defaulttype::Rigid3Types::Coord> VecCoordR;
    typedef sofa::helper::ReadAccessor  <sofa::core::objectmodel::Data<sofa::helper::vector<sofa::defaulttype::Vec3Types::Coord>>>       ReadCoord;
    typedef sofa::helper::WriteAccessor <sofa::core::objectmodel::Data<sofa::helper::vector<sofa::defaulttype::Vec3Types::Coord>>>       WriteCoord;
    typedef sofa::helper::WriteAccessor <sofa::core::objectmodel::Data<sofa::helper::vector<sofa::defaulttype::Vec3Types::Deriv>>>       WriteDeriv;
    typedef sofa::helper::ReadAccessor  <sofa::core::objectmodel::Data<sofa::helper::vector<sofa::defaulttype::RigidTypes::Coord>>>      ReadCoordR;
    typedef sofa::helper::WriteAccessor <sofa::core::objectmodel::Data<sofa::helper::vector<sofa::defaulttype::RigidTypes::Coord>>>      WriteCoordR;
    typedef sofa::helper::WriteAccessor <sofa::core::objectmodel::Data<sofa::helper::vector<sofa::defaulttype::RigidTypes::Deriv>>>      WriteDerivR;

    PBDBaseConstraint(sofa::simulation::Node* gnode = nullptr):
        m_indices(initData(&m_indices, sofa::helper::vector<uint>(), "indices", "ID of the vertices on wich this constraint is to apply")),
        m_nbIter(initData(&m_nbIter,static_cast<uint>(1),"iter","Number of iteration for the solver"))
    {}

public:
    SOFA_ABSTRACT_CLASS(PBDBaseConstraint, sofa::core::objectmodel::BaseObject);
    /*
     * Output : Solve the constraint adding in WriteCoord the computed displacement
     */
    virtual bool solve(sofa::simulation::Node * node) = 0;
    /*
     * Inputs : int -> Number of iterations
     *
     * Some constraints can have their own number of iterations. Do not use unless you know what you're doing.
     * This will most likely lead to unstabilities.
     */
    void setIterCount(int c) { m_nbIter.setValue (static_cast<uint>(c));}
protected:
    IndexSet m_indices; ///< Indices on wich to apply the constraint
    sofa::core::objectmodel::Data<unsigned int> m_nbIter;
};



using namespace sofa::core::objectmodel;
template < class T >
class SOFA_CORE_API PBDConstraint : public virtual PBDBaseConstraint
{
    typedef sofa::core::objectmodel::Data<sofa::helper::vector<uint>> IndexSet;
    typedef typename sofa::core::objectmodel::BaseLink BaseLink;
public:
    SOFA_ABSTRACT_CLASS(PBDConstraint, PBDBaseConstraint);

protected:

    PBDConstraint(sofa::simulation::Node* gnode = nullptr): PBDBaseConstraint(gnode),
        m_mechanicalObject(initLink("attachedTo","Object on wich the constraint will apply")),
        m_topology(initLink("topology","Link to the topology relevant for this object"))
    {}
    virtual ~PBDConstraint() { }
public:
    inline sofa::component::container::MechanicalObject< T >*  mechanical() {return m_mechanicalObject.getValue();}
    inline sofa::core::topology::BaseMeshTopology * topology() {return m_topology.getValue ();}
    PBDObject<T> * getPBDObject() { return m_pbdObject;}
    void linkPBDObject(PBDObject<T>* obj) { m_pbdObject = obj;}
protected:
    /// Input Model, also called parent
    PBDObject<T>*    m_pbdObject;

    SingleLink< PBDConstraint,
                sofa::component::container::MechanicalObject< T >,
                BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> m_mechanicalObject;

    SingleLink< PBDConstraint,
                sofa::core::topology::BaseMeshTopology,
                BaseLink::FLAG_STRONGLINK|BaseLink::FLAG_STOREPATH> m_topology;
};

#endif

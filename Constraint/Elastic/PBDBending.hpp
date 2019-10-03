#ifndef PBDBENDING_HPP
#define PBDBENDING_HPP

#include "PBDElasticConstraint.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>

class PBDBending : public PBDElasticConstraint
{
public:
    typedef sofa::defaulttype::Vec3 Vec3;
    PBDBending(sofa::simulation::Node* gnode = nullptr) : PBDElasticConstraint(gnode),
        m_alpha_wann(initData(&m_alpha_wann,(SReal)0,"a1","Low frequency damping")),
        m_alpha_too(initData(&m_alpha_too,(SReal)1.0,"a2","High frequency damping"))
    {}
    /*
     * Output : Solve the constraint adding in WriteCoord the computed displacement
     */
    virtual bool solve(sofa::simulation::Node* node);
    /*
     * Init function of sofa. It's called after the first init of the tree.
     */
    virtual void bwdInit () override;
    /// Construction method called by ObjectFactory.
    template<class T>
    static typename T::SPtr create(T*, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg)
    {
        typename T::SPtr obj = sofa::core::objectmodel::New<T>();
        if (context) context->addObject(obj);
        if (arg) obj->parse(arg);
        return obj;
    }

protected:
    /*
     * Inputs : PBDObject   -> Object on wich we will solve the constraint
     *          uint        -> First vertex of the edge
     *          uint        -> Second vertex of the edge
     *          WriteCoord  -> Free positions on wich we apply the dispalcement
     *          ReadDeriv   -> Velocity of the object's vertices
     *
     * Output : Compute and apply the correction to the concerned vertices
     */
    bool correction(uint a, uint b, WriteCoord&p,WriteDeriv& v);
protected:
    Data<SReal> m_alpha_wann;
    Data<SReal> m_alpha_too;
    PBDVertexTopology<sofa::defaulttype::Vec3Types>  m_stretch_topology;
    PBDBendingTopology m_bending_topology;
    SReal m_K;
};

#endif // PBDBENDING_HPP

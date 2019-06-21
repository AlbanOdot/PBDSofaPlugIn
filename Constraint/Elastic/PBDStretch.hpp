#ifndef PBDSTRECH_HPP
#define PBDSTRECH_HPP

#include "PBDElasticConstraint.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>

class PBDStretch : public PBDElasticConstraint
{
public:
    PBDStretch(sofa::simulation::Node* gnode = NULL) : PBDElasticConstraint(){}
    /*
     * Inputs : PBDObject   -> Object on wich we will solve the constraint
     *          WriteCoord  -> Free positions on wich we apply the dispalcement
     *
     * Output : Solve the constraint adding in WriteCoord the computed displacement
     */
    virtual void solve(PBDObject<sofa::defaulttype::Vec3Types>& object, WriteCoord& p);

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
     * Inputs : uint                -> Index of the first point
     *          pair<uint,double>   -> Index and re'st distance of the neighbor
     *          WriteCoord          -> Free position (same as solve)
     */
    void correction (uint i, const std::pair<uint,double>& voisin, WriteCoord &p);
protected:
    double m_K;
};

#endif // PBDSTRECH_HPP

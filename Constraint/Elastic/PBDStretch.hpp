#ifndef PBDSTRECH_HPP
#define PBDSTRECH_HPP

#include "PBDElasticConstraint.hpp"
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/simulation/Node.h>

class PBDStretch : public PBDElasticConstraint
{
public:
    PBDStretch(sofa::simulation::Node* gnode = NULL) : PBDElasticConstraint()
    {

    }
    /*
     * Output : Solve the constraint adding in WriteCoord the computed displacement
     */
    virtual void solve(sofa::simulation::Node * node);

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

    /*
     * Inputs : uint                                            -> Index of the first point
     *          pair<uint,double>                               -> Index and re'st distance of the neighbor
     *          sofa::defaulttype::Vec3Types::VecCoord          -> Free position (same as solve)
     */
    void correction (uint i, const std::pair<uint,double>& voisin, WriteCoord& p, SReal w0);
protected:
    PBDVertexTopology<sofa::defaulttype::Vec3Types>  m_stretch_topology;
    double m_K;
};

#endif // PBDSTRECH_HPP

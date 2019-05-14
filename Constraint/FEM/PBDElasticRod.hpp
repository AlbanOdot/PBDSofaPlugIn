#ifndef PBDELASTICROD_HPP
#define PBDELASTICROD_HPP

#include "PBDFEMConstraint.hpp"

class PBDElasticRod : public PBDFEMConstraint
{
    typedef sofa::defaulttype::Vec3 vec3;
public:
    PBDElasticRod() : PBDFEMConstraint(),
        m_radius(initData(&m_radius,0.1,"radius","Diameter of the rod"))
    {}
    virtual void solve(PBDObject& object, WriteCoord& p) override;
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
private:
    vec3 m_stretchingAndShearingKs;
    vec3 m_bendingAndTwistingKs;
    sofa::core::objectmodel::Data<SReal> m_radius;
};



#endif // PBDElasticRod_HPP

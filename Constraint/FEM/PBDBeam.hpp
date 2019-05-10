#ifndef PBDBEAM_HPP
#define PBDBEAM_HPP

#include "PBDFEMConstraint.hpp"

class PBDBeam : public PBDFEMConstraint
{
    typedef sofa::defaulttype::Vec3 vec3;
    typedef Eigen::Quaternion<SReal> Quaternion;
    typedef Eigen::Matrix<double,6,6> Mat6;
    typedef Eigen::Matrix<double,6,1> Vec6;

public:
    PBDBeam() : PBDFEMConstraint(),
        m_radius(initData(&m_radius,0.1,"radius","Diameter of the rod")),
        m_Massq(initData(&m_Massq,1.0,"Qmass","Inertia tensor reduced to a mass"))
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
    sofa::core::objectmodel::Data<SReal> m_Massq;


};

#endif // PBDBEAM_HPP

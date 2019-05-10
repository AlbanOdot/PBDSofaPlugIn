#ifndef PBDSTIFFROD_HPP
#define PBDSTIFFROD_HPP

#include "../PBDFEMConstraint.hpp"

class PBDStiffRod : public PBDFEMConstraint
{
    typedef sofa::defaulttype::Vec3 vec3;
    typedef Eigen::Quaternion<SReal> Quaternion;
    typedef Eigen::Matrix<double,6,6> Mat6;
    typedef Eigen::Matrix<double,6,1> Vec6;

public:
    PBDStiffRod() : PBDFEMConstraint(),
        m_stretchingAndShearingKs(initData(&m_stretchingAndShearingKs,vec3(1.0,1.0,1.0),"stretch","Stretching and Shearing coefficients")),
        m_bendingAndTwistingKs(initData(&m_bendingAndTwistingKs,vec3(1.0,1.0,1.0),"bendingK","Bending and Twisting coefficients")),
        m_stiffnessCoefficientK(initData(&m_stiffnessCoefficientK,vec3(1e-8,1e-8,1e-8),"stiffnessK","Stiffness coefficients")),
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
    sofa::core::objectmodel::Data<vec3> m_stretchingAndShearingKs;
    sofa::core::objectmodel::Data<vec3> m_bendingAndTwistingKs;
    sofa::core::objectmodel::Data<vec3> m_stiffnessCoefficientK;
    sofa::core::objectmodel::Data<SReal> m_radius;

};

#endif // PBDStiffRod_HPP

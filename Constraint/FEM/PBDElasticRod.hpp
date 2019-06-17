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

    /*
     * Inputs : PBDObject   -> Object on wich we will solve the constraint
     *          WriteCoord  -> Free positions on wich we apply the dispalcement
     *
     * Output : Solve the constraint adding in WriteCoord the computed displacement
     */
    virtual void solve(PBDObject& object, WriteCoord& p) override;

    /*
     * Init function of sofa. It's called after the first init of the tree.
     */
    virtual void bwdInit () override;

    /*
     * Inputs : ElasticRodData          -> Data structure representing the rod informations
     *          vector<Quaternionr>     -> Free orientation
     *          PBDObject   -> Object on wich we will solve the constraint
     *          WriteCoord  -> Free positions on wich we apply the dispalcement
     *
     * Output : Compute and apply the correction
     */
    static void correction( ElasticRodData& eRod, std::vector<Quaternionr>& u, PBDObject& object, WriteCoord& p, const vec3& bending_twisting, const uint e);
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
    vec3 m_stretchingAndShearingKs;
    vec3 m_bendingAndTwistingKs;
    sofa::core::objectmodel::Data<SReal> m_radius;
};



#endif // PBDElasticRod_HPP

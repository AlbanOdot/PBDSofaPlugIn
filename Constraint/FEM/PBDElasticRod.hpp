#ifndef PBDELASTICROD_HPP
#define PBDELASTICROD_HPP

#include "PBDFEMConstraint.hpp"
#include "OrientedConstraint.hpp"
#include "../../InternalData/Data/FEM/PBDElasticRodData.hpp"

class PBDElasticRod : public OrientedConstraint
{
    typedef sofa::defaulttype::Vec3 vec3;
public:
    PBDElasticRod() : OrientedConstraint(),
        m_radius(initData(&m_radius,0.1,"radius","Diameter of the rod")){}

    /*
     * Inputs : PBDObject   -> Object on wich we will solve the constraint
     *          WriteCoord  -> Free positions on wich we apply the dispalcement
     *
     * Output : Solve the constraint adding in WriteCoord the computed displacement
     */
    virtual void solve(sofa::simulation::Node * node) override;

    /*
     * Init function of sofa. It's called after the first init of the tree.
     */
    virtual void bwdInit () override;

    /*
     * Inputs : ElasticRodData          -> Data structure representing the rod informations
     *          vector<Quaternion>     -> Free orientation
     *          PBDObject   -> Object on wich we will solve the constraint
     *          WriteCoord  -> Free positions on wich we apply the dispalcement
     *
     * Output : Compute and apply the correction
     */
    static inline void correction(sofa::defaulttype::RigidTypes::Coord& p0,sofa::defaulttype::RigidTypes::Coord& p1,
                                  const SReal& l,
                                  const SReal w0,const SReal w1,
                                  Quaternionr& q0,Quaternionr& q1,
                                  const SReal& wq0, const SReal& wq1,
                                  const Quaternionr& restDarboux,
                                  const Vec3& bending_twisting)
    {
        static const SReal eps = 1e-6;
        //  COMPUTE STRETCHING AND SHEARING
        vec3 d3(static_cast<SReal>(2.0) * (q0.x() * q0.z() + q0.w() * q0.y()),
                static_cast<SReal>(2.0) * (q0.y() * q0.z() - q0.w() * q0.x()),
                q0.w() * q0.w() - q0.x() * q0.x() - q0.y() * q0.y() + q0.z() * q0.z());	//third director d3 = q0 * e_3 * q0_conjugate

        vec3 gamma = (p1.getCenter ()- p0.getCenter ()) / l - d3;
        gamma     /= (w1 + w0) / l+ wq0 * static_cast<SReal>(4.0)*l + eps;

        p0.getCenter () += (wq0 * w0) * gamma;
        p1.getCenter () -= (wq1 * w1) * gamma;
        // Cs * q0 * e_3.conjugate (cheaper than quaternion product)
        //Quaternionr dq0 = Quaternionr(0.0, gamma.x(), gamma.y(), gamma.z()) * Quaternionr(q0.z(), -q0.y(), q0.x(), -q0.w());//Bending correction due to displacement in stretching
        //Equivalent to the one above without computing the zeros
        Quaternionr dq0 = { /*0.0*/           q0.y()*gamma[0] - q0.x()*gamma[1] + q0.w()*gamma[2],// Real part
                            q0.z()*gamma[0]   /*0.0*/         - q0.x()*gamma[2] - q0.w()*gamma[1],//i
                            q0.z()*gamma[1] - q0.y()*gamma[2]   /*0.0*/         + q0.w()*gamma[0],//j
                            q0.z()*gamma[2] + q0.y()*gamma[1] + q0.x()*gamma[0]   /*0.0*/};//k
        q0.coeffs() += (wq0 *static_cast<SReal>(2.0)* l) * dq0.coeffs ();
        q0.normalize ();

        // COMPUTE BENDING AND TWISTING
        Quaternionr omega   = q0.conjugate() * q1;   //darboux vector
        Quaternionr omega_plus;
        omega_plus.coeffs() = omega.coeffs() + restDarboux.coeffs(); //delta Omega with + Omega_0
        omega.coeffs()      = omega.coeffs() - restDarboux.coeffs(); //delta Omega with - Omega_0

        if (omega.squaredNorm() > omega_plus.squaredNorm())
            omega = omega_plus;

        for (uint i = 0; i < 3; i++)
            omega.coeffs()[i] *= bending_twisting[i] / (wq0 + wq1 + eps);
        omega.w() = 0.0;    //discrete Darboux vector does not have vanishing scalar part
        //Equivalent to "q0.coeffs() += wq0 * (q1 * omega).coeffs ()" without computing the zeros
        q0 = { /*0.0*/ q0.w() - wq0 * ( omega.x()*q1.x() + omega.y()*q1.y() + omega.z()*q1.z()),// Real part
               /*0.0*/ q0.x() + wq0 * ( omega.x()*q1.w() - omega.y()*q1.z() + omega.z()*q1.y()),//i
               /*0.0*/ q0.y() + wq0 * ( omega.x()*q1.z() + omega.y()*q1.w() - omega.z()*q1.x()),//j
               /*0.0*/ q0.z() + wq0 * (-omega.x()*q1.y() + omega.y()*q1.x() + omega.z()*q1.w())};
        q0.normalize ();
        //Equivalent to "q1.coeffs() -= wq1 * (q0 * omega).coeffs ()" without computing the zeros
        q1 = { /*0.0*/ q1.w() + wq1 * ( omega.x()*q0.x() + omega.y()*q0.y() + omega.z()*q0.z()),// Real part
               /*0.0*/ q1.x() - wq1 * ( omega.x()*q0.w() - omega.y()*q0.z() + omega.z()*q0.y()),//i
               /*0.0*/ q1.y() - wq1 * ( omega.x()*q0.z() + omega.y()*q0.w() - omega.z()*q0.x()),//j
               /*0.0*/ q1.z() - wq1 * (-omega.x()*q0.y() + omega.y()*q0.x() + omega.z()*q0.w())};
        q1.normalize ();
    }
    /// Construction method called by ObjectFactory.
    template<class T>
    static typename T::SPtr create(T*, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg)
    {
        typename T::SPtr obj = sofa::core::objectmodel::New<T>();
        if (context) context->addObject(obj);
        if (arg) obj->parse(arg);
        return obj;
    }

    virtual void draw(const sofa::core::visual::VisualParams* vparams) override;

protected:
    ElasticRodData m_elastic_rod;
    vec3 m_stretchingAndShearingKs;
    vec3 m_bendingAndTwistingKs;
    sofa::core::objectmodel::Data<SReal> m_radius;
};



#endif // PBDElasticRod_HPP

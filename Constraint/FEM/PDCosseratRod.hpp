#ifndef PDCOSSERATROD_HPP
#define PDCOSSERATROD_HPP

#include "OrientedConstraint.hpp"
#include "../../InternalData/Data/FEM/PDCosseratRodData.hpp"
/*
 *
 * This class implement the following paper
 * https://igl.ethz.ch/projects/cosserat-rods/CosseratRods-SCA2018.pdf
 *
 * Cosserat Rods with Projective Dynamics
 */
class PDCosseratRod : public OrientedConstraint
{
    typedef sofa::defaulttype::Vec3 vec3;
public:
    PDCosseratRod() : OrientedConstraint (),
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
    static inline void correction(const PDCosseratRodData& cRod,
                                  std::vector<Quaternionr>& u,
                                  const Quaternionr& restDarboux,
                                  const SReal w0, const SReal w1,
                                  WriteCoordR& p,
                                  const Eigen::Vector3d& bending_twisting,
                                  const uint e)
    {
        static const SReal eps = 1e-6;

        const uint a = cRod.beginIdx (e);
        const uint z = cRod.endIdx(e);

        //  COMPUTE STRETCHING AND SHEARING
        vec3 d3(static_cast<SReal>(2.0) * (u[a].x() * u[a].z() + u[a].w() * u[a].y()),
                static_cast<SReal>(2.0) * (u[a].y() * u[a].z() - u[a].w() * u[a].x()),
                u[a].w() * u[a].w() - u[a].x() * u[a].x() - u[a].y() * u[a].y() + u[a].z() * u[a].z());	//third director d3 = q0 * e_3 * q0_conjugate

        vec3 gamma = (p[z].getCenter ()- p[a].getCenter ()) / cRod.length(a) - d3;
        gamma     /= (w1 + w0) / cRod.length(a)+ cRod.wq(a) * static_cast<SReal>(4.0)*cRod.length(a) + eps;
        p[a].getCenter () += w0 * gamma;
        p[z].getCenter () -= w1 * gamma;
        // Cs * q * e_3.conjugate (cheaper than quaternion product)
        Quaternionr dq0 = Quaternionr(0.0, gamma.x(), gamma.y(), gamma.z()) * Quaternionr(u[a].z(), -u[a].y(), u[a].x(), -u[a].w());//Bending correction due to displacement in stretching
        u[a].coeffs() += cRod.wq(a) *(static_cast<SReal>(2.0)* cRod.length(a)) * dq0.coeffs ();

        // COMPUTE BENDING AND TWISTING
        Quaternionr omega   = u[a].conjugate() * u[z];   //darboux vector
        Quaternionr omega_plus;
        omega_plus.coeffs() = omega.coeffs() + restDarboux.coeffs(); //delta Omega with + Omega_0
        omega.coeffs()      = omega.coeffs() - restDarboux.coeffs(); //delta Omega with - Omega_0

        if (omega.squaredNorm() > omega_plus.squaredNorm())
            omega = omega_plus;

        for (uint i = 0; i < 3; i++)
            omega.coeffs()[i] *= bending_twisting[i] / (cRod.wq(a) + cRod.wq(z) + eps);
        omega.w() = 0.0;    //discrete Darboux vector does not have vanishing scalar part
        u[a].coeffs() += cRod.wq(a) * (u[z] * omega).coeffs ();
        u[z].coeffs() -= cRod.wq(z) * (u[a] * omega).coeffs ();
    }
    /*
     * Inputs : ElasticRodData          -> Data structure representing the rod informations
     *          vector<Quaternionr>     -> Free orientation
     *          vector<Vector3r>        -> Displacement
     *          vector<Quaternionr>     -> Rotation correction
     *          PBDObject               -> Object on wich we will solve the constraint
     *          WriteCoord              -> Free positions on wich we apply the dispalcement
     *          uint                    -> index of the segment
     *
     * Output : Compute and apply the correction
     */
    inline void solveLinearSystem(const PDCosseratRodData& cRod,
                                  std::vector<Quaternionr>& u,
                                  const SReal w0, const SReal w1,
                                  WriteCoordR& p,
                                  const uint e)
    {
        const uint a = cRod.beginIdx (e);
        const uint z = cRod.endIdx(e);

        //  COMPUTE STRETCHING AND SHEARING
        vec3 d3(static_cast<SReal>(2.0) * (u[a].x() * u[a].z() + u[a].w() * u[a].y()),
                static_cast<SReal>(2.0) * (u[a].y() * u[a].z() - u[a].w() * u[a].x()),
                u[a].w() * u[a].w() - u[a].x() * u[a].x() - u[a].y() * u[a].y() + u[a].z() * u[a].z());	//third director d3 = q0 * e_3 * q0_conjugate

        vec3 gamma = (p[z].getCenter ()- p[a].getCenter ()) / cRod.length(e) - d3;
        m_dx[a] += (w0 * cRod.ws(a) + 1e-6) * gamma;
        m_dx[z] -= (w1 * cRod.ws(z) + 1e-6) * gamma;
        // Cs * q * e_3.conjugate (cheaper than quaternion product)
        //The 0.5 is already computed inside wbt
        Quaternionr u0 = Quaternionr(0.0, gamma.x(), gamma.y(), gamma.z()) * Quaternionr(u[a].z(), -u[a].y(), u[a].x(), -u[a].w());
        m_dq[a].coeffs() = (cRod.wq(a) * static_cast<SReal>(2.0) * cRod.length(a) * cRod.wbt(a)) * u0.coeffs ();

        // COMPUTE BENDING AND TWISTING
        Quaternionr omega  = u[a].conjugate() * u[z];   //darboux vector
        omega.w() = 0.0;    //discrete Darboux vector does not have vanishing scalar part
        m_dq[a].coeffs () += (cRod.wq(a) * cRod.wbt(a)) * (m_dq[a] * omega).coeffs ();
        m_dq[z].coeffs () -= (cRod.wbt(z) * cRod.wq (z)) * (u[z] * omega).coeffs ();
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

private:
    std::vector<vec3> m_dx;
    std::vector<Quaternionr> m_dq;
    PDCosseratRodData m_cosserat_rod;
    vec3 m_stretchingAndShearingKs;
    vec3 m_bendingAndTwistingKs;
    sofa::core::objectmodel::Data<SReal> m_radius;
};



#endif // PDCosseratRod_HPP

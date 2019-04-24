#include "PBDBeam.hpp"
#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>
#include <sofa/helper/Quater.h>

int PBDBeamClass = sofa::core::RegisterObject("Constraint that correct beam.")
                   .add< PBDBeam >();

typedef sofa::defaulttype::Vec3 V3;

void PBDBeam::bwdInit ()
{

    auto node = dynamic_cast<sofa::simulation::Node*>(this->getContext());
    dt2 = node->getDt () * node->getDt ();
    SReal nu = m_poisson_ratio.getValue ();
    SReal E = m_young_modulus.getValue ();
    SReal mu = E / (2.0 * ( 1.0 + nu ));
    float d = m_radius.getValue ();
    float I = 0.25*M_PI*d*d*d*d;
    //This constant is here to make the unit of the scene in GPa -> 193 GPA = stainless steel
    m_bendingAndTwistingKs = vec3(I*E,I*E,I*(2.0*mu));
    m_bendingAndTwistingKs /= m_nbIter.getValue ();

}


void PBDBeam::solve(PBDObject &object, WriteCoord &p)
{
    static const SReal eps = 1e-6;
    const SReal invMass1 = object.invMass ();
    const SReal invMass0 = invMass1;
    const SReal invMassq0 = 1.0 / m_Massq.getValue ();
    const SReal invMassq1 = invMassq0;

    if(object.beam ().empty ())
        object.computeBeam();

    const auto& edges = object.sofaTopology ()->getEdges ();
    auto& beam = object.beam();
    for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
    {
        for(uint e = 0; e <= edges.size(); ++e)
        {
            BeamElement& be2 = beam[e+1];
            BeamElement& be = beam[e];
            vec3& p0 = p[be.current ()];
            vec3& p1 = p[be2.current()];
            Quaternion& q0 = be.m_q;
            Quaternion& q1 = be2.m_q;
            Quaternion& restDarbouxVector = be.restDarboux ();

            //  COMPUTE STRETCHING AND SHEARING
            vec3 d3(static_cast<SReal>(2.0) * (q0.x() * q0.z() + q0.w() * q0.y()),
                    static_cast<SReal>(2.0) * (q0.y() * q0.z() - q0.w() * q0.x()),
                    q0.w() * q0.w() - q0.x() * q0.x() - q0.y() * q0.y() + q0.z() * q0.z());	//third director d3 = q0 * e_3 * q0_conjugate

            vec3 gamma((p1 - p0) * be.invLength () - d3);
            gamma /= (invMass1 + invMass0) * be.invLength ()+ invMassq0 * static_cast<SReal>(4.0)*be.length() + eps;
            p0 += invMass0 * gamma;
            p1 -= invMass1 * gamma;

            //                               Cs                               *  q * e_3.conjugate (cheaper than quaternion product)
            Quaternion dq0 = Quaternion(0.0, gamma.x(), gamma.y(), gamma.z()) * Quaternion(q0.z(), -q0.y(), q0.x(), -q0.w());
            dq0.coeffs () *= static_cast<SReal>(2.0) *invMassq0 * be.length();
            q0.coeffs ()  += dq0.coeffs ();

            // COMPUTE BENDING AND TWISTING
            Quaternion omega = q0.conjugate() * q1;   //darboux vector
            Quaternion omega_plus;
            omega_plus.coeffs() = (omega.coeffs() + be.restDarboux ().coeffs()).eval();     //delta Omega with -Omega_0
            omega.coeffs() = omega.coeffs() - restDarbouxVector.coeffs();                 //delta Omega with + omega_0
            if (omega.squaredNorm() > omega_plus.squaredNorm())
                omega = omega_plus;

            for (int i = 0; i < 3; i++) omega.coeffs()[i] *= m_bendingAndTwistingKs[i] / (invMassq0 + invMassq1 + eps);
            omega.w() = 0.0;    //discrete Darboux vector does not have vanishing scalar part
            q0.coeffs () += invMassq0 * (q1 * omega).coeffs ();
            q0.normalize ();
            q1.coeffs () -= invMassq1 * (q0 * omega).coeffs ();
            q1.normalize ();
        }
    }

}

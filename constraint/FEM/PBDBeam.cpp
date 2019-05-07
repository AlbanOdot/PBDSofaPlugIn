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
    m_maxDisplacement = 0;

}


void PBDBeam::solve(PBDObject &object, WriteCoord &p)
{

    if(object.beam ().empty ())
    {
        object.computeBeam();
        object.applyFixedPoint(m_indices.getValue ());
    }
    auto& beam = object.beam();
    auto& u = object.freeOrientation ();
    static const SReal eps = 1e-6;
    const SReal invMass1 = object.invMass ();
    const SReal invMass0 = invMass1;
    uint  fixedIdx = 0;
    for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
    {
        for(uint e = 0; e < beam.size (); ++e )
        {

            //  COMPUTE STRETCHING AND SHEARING
            vec3 d3;
            d3[0] = static_cast<SReal>(2.0) * (u[e].x() * u[e].z() + u[e].w() * u[e].y());
            d3[1] = static_cast<SReal>(2.0) * (u[e].y() * u[e].z() - u[e].w() * u[e].x());
            d3[2] = u[e].w() * u[e].w() - u[e].x() * u[e].x() - u[e].y() * u[e].y() + u[e].z() * u[e].z();	//third director d3 = q0 * e_3 * q0_conjugate

            vec3 gamma = (p[beam[e].extremity(1)] - p[beam[e].extremity(0)]) / beam[e].length () - d3;
            gamma /= (invMass1 + invMass0) / beam[e].length ()+ beam[e].wq() * static_cast<SReal>(4.0)*beam[e].length() + eps;

            p[beam[e].extremity(0)] += invMass0 * gamma;
            p[beam[e].extremity(1)] -= invMass1 * gamma;

            //                               Cs                               *  q * e_3.conjugate (cheaper than quaternion product)
            Quaternion dq0 = Quaternion(0.0, gamma.x(), gamma.y(), gamma.z()) * Quaternion(u[e].z(), -u[e].y(), u[e].x(), -u[e].w());
            dq0.coeffs () *= static_cast<SReal>(2.0) * beam[e].wq() * beam[e].length();
            u[e].coeffs() += dq0.coeffs ();

            // COMPUTE BENDING AND TWISTING
            Quaternion omega = u[e].conjugate() * u[e+1];   //darboux vector
            Quaternion omega_plus;
            omega_plus.coeffs() = omega.coeffs() + beam[e].restDarboux ().coeffs();     //delta Omega with -Omega_0
            omega.coeffs() = omega.coeffs() - beam[e].restDarboux ().coeffs();          //delta Omega with + omega_0

            if (omega.squaredNorm() > omega_plus.squaredNorm())
                omega = omega_plus;

            for (uint i = 0; i < 3; i++)
                omega.coeffs()[i] *= m_bendingAndTwistingKs[i] / (beam[e].wq() + beam[e+1].wq() + eps);
            omega.w() = 0.0;    //discrete Darboux vector does not have vanishing scalar part
            u[e].coeffs() += beam[e].wq() * (u[e+1] * omega).coeffs ();
            u[e+1].coeffs() -=  beam[e+1].wq() * (u[e] * omega).coeffs ();
        }
        u[beam.size()-1] = u[beam.size()-2];
    }
//    SReal currentDisplacement = (object.rest ()[object.rest ().size () - 1] - p[object.rest ().size () - 1]).norm ();

//    m_maxDisplacement = currentDisplacement;
//    std::cout << "Displacement is now : "<<m_maxDisplacement<<std::endl;


}

inline void PBDBeam::solveCosserat(PBDObject& object, WriteCoord& p, uint e)
{
    auto& beam = object.beam();
    auto& u = object.freeOrientation ();



}
#include "PBDElasticRod.hpp"
#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>

int PBDElasticRodClass = sofa::core::RegisterObject("Constraint that correct elastic rod.")
                         .add< PBDElasticRod >();

typedef sofa::defaulttype::Vec3 V3;

void PBDElasticRod::bwdInit ()
{

    auto node = dynamic_cast<sofa::simulation::Node*>(this->getContext());
    dt2 = node->getDt () * node->getDt ();
    SReal nu = m_poisson_ratio.getValue ();
    SReal E = m_young_modulus.getValue ();
    SReal mu = E / (2.0 * ( 1.0 + nu ));
    SReal d = m_radius.getValue ();
    SReal I = 0.25*M_PI*d*d*d*d;
    //This constant is here to make the unit of the scene in GPa -> 193 GPA = stainless steel
    m_bendingAndTwistingKs = vec3(I*E,I*E,I*(2.0*mu));
    m_bendingAndTwistingKs /= m_nbIter.getValue ();

}


void PBDElasticRod::solve(PBDObject &object, WriteCoord &p)
{

    if(!object.hasDataType(PBDObject::ELASTICROD) || !object.hasDataType(PBDObject::ORIENTED))
    {
        if(!object.hasDataType(PBDObject::ORIENTED))
            object.computeOrientation ();
        object.computeElasticRod();
        object.elasticRod().applyFixedPoint(m_indices.getValue ());
        object.orientation ().setInertia ({{m_bendingAndTwistingKs[0],m_bendingAndTwistingKs[1], m_bendingAndTwistingKs[2]}});
    }
    auto& eRod = object.elasticRod ();
    auto& u    = object.orientation ().freeOrientation ();
    for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
    {
        for(uint e = 0; e < eRod.wq().size (); ++e )
        {
            correction(eRod,u,object,p,m_bendingAndTwistingKs,e);
        }
        u[u.size ()-1] = u[u.size()-2];
    }
}


void PBDElasticRod::correction( ElasticRodData& eRod, std::vector<Quaternionr>& u, PBDObject& object, WriteCoord& p,const vec3& bending_twisting, const uint e)
{
    static const SReal eps = 1e-6;

    uint a = eRod.beginIdx (e);
    uint z = eRod.endIdx(e);
    const SReal invMass1 = object.invMass(a);
    const SReal invMass0 = object.invMass(z);

    //  COMPUTE STRETCHING AND SHEARING
    vec3 d3;
    d3[0] = static_cast<SReal>(2.0) * (u[a].x() * u[a].z() + u[a].w() * u[a].y());
    d3[1] = static_cast<SReal>(2.0) * (u[a].y() * u[a].z() - u[a].w() * u[a].x());
    d3[2] = u[a].w() * u[a].w() - u[a].x() * u[a].x() - u[a].y() * u[a].y() + u[a].z() * u[a].z();	//third director d3 = q0 * e_3 * q0_conjugate

    vec3 gamma = (p[z] - p[a]) / eRod.length(e) - d3;
    gamma     /= (invMass1 + invMass0) / eRod.length(e)+ eRod.wq(e) * static_cast<SReal>(4.0)*eRod.length(e) + eps;
    p[a] += invMass0 * gamma;
    p[z] -= invMass1 * gamma;

    if( eRod.wq(a) > 0.0 )
    {
        // Cs * q * e_3.conjugate (cheaper than quaternion product)
        Quaternionr dq0 = Quaternionr(0.0, gamma.x(), gamma.y(), gamma.z()) * Quaternionr(u[a].z(), -u[a].y(), u[a].x(), -u[a].w());
        u[a].coeffs() += (static_cast<SReal>(2.0) * eRod.wq(e) * eRod.length(e)) * dq0.coeffs ();
    }

    // COMPUTE BENDING AND TWISTING
    Quaternionr omega    = u[a].conjugate() * u[z];   //darboux vector
    Quaternionr omega_plus;
    omega_plus.coeffs() = omega.coeffs() + object.orientation().restDarboux(a).coeffs(); //delta Omega with + Omega_0
    omega.coeffs()      = omega.coeffs() - object.orientation().restDarboux(a).coeffs(); //delta Omega with - Omega_0

    if (omega.squaredNorm() > omega_plus.squaredNorm())
        omega = omega_plus;

    for (uint i = 0; i < 3; i++)
        omega.coeffs()[i] *= bending_twisting[i] / (eRod.wq(a) + eRod.wq(z) + eps);
    omega.w() = 0.0;    //discrete Darboux vector does not have vanishing scalar part

    u[a].coeffs() +=  eRod.wq(a) * (u[z] * omega).coeffs ();
    u[z].coeffs() -=  eRod.wq(z) * (u[a] * omega).coeffs ();

}

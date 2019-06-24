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


void PBDElasticRod::solve(PBDObject<sofa::defaulttype::Rigid3Types> &object, WriteCoord &p)
{

    if(!object.hasDataType(PBDObject<sofa::defaulttype::Rigid3Types>::ELASTICROD) || !object.hasDataType(PBDObject<sofa::defaulttype::Rigid3Types>::ORIENTED))
    {
        if(!object.hasDataType(PBDObject<sofa::defaulttype::Rigid3Types>::ORIENTED))
            object.computeOrientation ();
        object.computeElasticRod();
        object.elasticRod().applyFixedPoint(m_indices.getValue ());
        object.orientation ().setInertia ({{m_bendingAndTwistingKs[0],m_bendingAndTwistingKs[1], m_bendingAndTwistingKs[2]}});
    }
    auto& eRod = object.elasticRod ();
    auto& u    = object.orientation().freeOrientation();
    for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
    {
        for(uint e = 0; e < eRod.wq().size (); ++e )
        {
            correction(eRod,u,object,p,m_bendingAndTwistingKs,e);
        }
        u[u.size ()-1] = u[u.size()-2];
    }
}

inline SReal squaredNorm( const Quaternion& q)
{
    return q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
}



void PBDElasticRod::correction( ElasticRodData& eRod, std::vector<Quaternion>& u, PBDObject<sofa::defaulttype::Rigid3Types>& object, WriteCoord& p,const vec3& bending_twisting, const uint e)
{
    static const SReal eps = 1e-6;

    uint a = eRod.beginIdx (e);
    uint z = eRod.endIdx(e);
    const SReal invMass1 = object.invMass(a);
    const SReal invMass0 = object.invMass(z);

    //  COMPUTE STRETCHING AND SHEARING
    // 3 = w, 0 = x, 1 = y, 2 = z
    vec3 d3;
    d3[0] = static_cast<SReal>(2.0) * (u[a][0] * u[a][2] + u[a][3] * u[a][1]);
    d3[1] = static_cast<SReal>(2.0) * (u[a][1] * u[a][2] - u[a][3] * u[a][0]);
    d3[2] = u[a][3] * u[a][3] - u[a][0] * u[a][0] - u[a][1] * u[a][1] + u[a][2] * u[a][2];	//third director d3 = q0 * e_3 * q0_conjugate

    vec3 gamma = (p[z].getCenter () - p[a].getCenter ()) / eRod.length(e) - d3;
    gamma     /= (invMass1 + invMass0) / eRod.length(e)+ eRod.wq(e) * static_cast<SReal>(4.0)*eRod.length(e) + eps;
    p[a].getCenter () += invMass0 * gamma;
    p[z].getCenter () -= invMass1 * gamma;

//    if( eRod.wq(a) > 0.0 )
//    {
//        // Cs * q * e_3.conjugate (cheaper than quaternion product)
//        Quaternion dq0 = Quaternion(0.0, gamma[0], gamma[1], gamma[2]) * Quaternion(u[a][2], -u[a][1], u[a][0], -u[a][3]);
//        u[a] += dq0 * (static_cast<SReal>(2.0) * eRod.wq(e) * eRod.length(e));
//    }

////    // COMPUTE BENDING AND TWISTING
////    Quaternion omega    = u[a].inverse ()* u[z];   //darboux vector
////    Quaternion omega_plus;
////    omega_plus = omega + object.orientation().restDarboux(a); //delta Omega with + Omega_0
////    omega = omega + (object.orientation().restDarboux(a) * -1); //delta Omega with - Omega_0

////    if (squaredNorm(omega) > squaredNorm (omega_plus))
////        omega = omega_plus;

////    for (uint i = 0; i < 3; i++)
////        omega[i] *= bending_twisting[i] / (eRod.wq(a) + eRod.wq(z) + eps);
////    omega[3] = 0.0;    //discrete Darboux vector does not have vanishing scalar part

////    u[a] += (u[z] * omega) * eRod.wq(a);
////    u[z] += (u[a] * omega) * eRod.wq(z) * -1.0;

}

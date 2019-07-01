#include "PDCosseratRod.hpp"
#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>

int PDCosseratRodClass = sofa::core::RegisterObject("Constraint that correct elastic rod.")
                         .add< PDCosseratRod >();

void PDCosseratRod::solve(PBDObject<sofa::defaulttype::RigidTypes> &object, WriteCoord &p)
{

    if(!object.hasDataType(PBDObject<sofa::defaulttype::RigidTypes>::COSSERATROD) || !object.hasDataType(PBDObject<sofa::defaulttype::RigidTypes>::ORIENTED))
    {
        if(!object.hasDataType(PBDObject<sofa::defaulttype::RigidTypes>::ORIENTED))
            object.computeOrientation ();

        object.computeCosseratRod ();
        object.cosseratRod().applyFixedPoint(m_indices.getValue ());
        object.cosseratRod().setupW(m_young_modulus.getValue (),m_poisson_ratio.getValue (), m_radius.getValue ());
        object.orientation ().setInertia ({{m_bendingAndTwistingKs[0],m_bendingAndTwistingKs[1], m_bendingAndTwistingKs[2]}});

        Quaternionr q; q.coeffs ().setZero ();
        vec3 x(0,0,0);
        for(uint e = 0; e < object.cosseratRod().wq().size (); ++e )
        {
            m_dq.emplace_back(q);
            m_dx.emplace_back(x);
        }
    }

    auto& cRod = object.cosseratRod ();
    auto& u    = object.orientation ().freeOrientation ();
    const uint nbElem = cRod.wq().size ();
    for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
    {
        //Compute first
        uint a = cRod.beginIdx (0);
        uint z = cRod.endIdx(0);
        m_dx[a].set(0,0,0);m_dx[z].set(0,0,0);
        m_dq[a].coeffs ().setZero ();m_dq[z].coeffs ().setZero ();
        //Ligne 9
        correction(cRod,u,object,p,m_bendingAndTwistingKs,0);
        //Ligne 10
        solveLinearSystem(cRod,u,m_dx,m_dq,object,p,0);

        //Compute the current one and integrate the previous one
        for(uint e = 1; e < nbElem; ++e )
        {
            uint a = cRod.beginIdx (e);
            uint z = cRod.endIdx(e);
            m_dx[a].set(0,0,0);m_dx[z].set(0,0,0);
            m_dq[a].coeffs ().setZero ();m_dq[z].coeffs ().setZero ();
            //Ligne 9
            correction(cRod,u,object,p,m_bendingAndTwistingKs,e);
            //Ligne 10
            solveLinearSystem(cRod,u,m_dx,m_dq,object,p,e);

            p[e].getCenter ()+= m_dx[e-1];
            u[e].coeffs () += m_dq[e-1].coeffs ();
        }

        //integrate the last one
        p[nbElem - 1 ].getCenter () += m_dx[nbElem - 1];
        u[nbElem - 1].coeffs () += m_dq[nbElem - 1].coeffs ();
    }
}


void PDCosseratRod::solveLinearSystem( PDCosseratRodData& cRod,
                                       std::vector<Quaternionr>& u,
                                       std::vector<vec3> dx,
                                       std::vector<Quaternionr> dq,
                                       PBDObject<sofa::defaulttype::RigidTypes>& object,
                                       WriteCoord& p,
                                       const uint e)
{

    const uint a = cRod.beginIdx (e);
    const uint z = cRod.endIdx(e);
    const SReal invMass1 = object.invMass(a);
    const SReal invMass0 = object.invMass(z);

    //  COMPUTE STRETCHING AND SHEARING
    vec3 d3(static_cast<SReal>(2.0) * (u[a].x() * u[a].z() + u[a].w() * u[a].y()),
            static_cast<SReal>(2.0) * (u[a].y() * u[a].z() - u[a].w() * u[a].x()),
            u[a].w() * u[a].w() - u[a].x() * u[a].x() - u[a].y() * u[a].y() + u[a].z() * u[a].z());	//third director d3 = q0 * e_3 * q0_conjugate

    vec3 gamma = (p[z].getCenter ()- p[a].getCenter ()) / cRod.length(e) - d3;
    dx[a] += (invMass0 * cRod.ws(a) + 1e-6) * gamma;
    dx[z] -= (invMass1 * cRod.ws(z) + 1e-6) * gamma;


    // COMPUTE BENDING AND TWISTING
    Quaternionr omega  = u[a].conjugate() * u[z];   //darboux vector
    omega.w() = 0.0;    //discrete Darboux vector does not have vanishing scalar part
    if( cRod.wq(a) > 0.0 )
    {
        // Cs * q * e_3.conjugate (cheaper than quaternion product)
        //The 0.5 is already computed inside wbt
        Quaternionr u0 = Quaternionr(0.0, gamma.x(), gamma.y(), gamma.z()) * Quaternionr(u[a].z(), -u[a].y(), u[a].x(), -u[a].w());
        dq[a].coeffs() = (static_cast<SReal>(2.0) * cRod.length(a) * cRod.wbt(a)) * u0.coeffs ();
        dq[a].coeffs () += cRod.wbt(a) * (dq[a] * omega).coeffs ();
    }
    dq[z].coeffs () -= (cRod.wbt(z) * cRod.wq (z)) * (u[z] * omega).coeffs ();
}

#include "PDCosseratRod.hpp"
#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>

int PDCosseratRodClass = sofa::core::RegisterObject("Constraint that correct elastic rod.")
                         .add< PDCosseratRod >();

typedef sofa::defaulttype::Vec3 V3;

void PDCosseratRod::bwdInit ()
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


void PDCosseratRod::solve(PBDObject &object, WriteCoord &p)
{

    if(!object.hasDataType(PBDObject::COSSERATROD) || !object.hasDataType(PBDObject::ORIENTED))
    {
        if(!object.hasDataType(PBDObject::ORIENTED))
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

    for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
    {
        uint a = cRod.beginIdx (0);
        uint z = cRod.endIdx(0);
        m_dx[a].set(0,0,0);m_dx[z].set(0,0,0);
        m_dq[a].coeffs ().setZero ();m_dq[z].coeffs ().setZero ();
        //Ligne 9
        correction(cRod,u,object,p,m_bendingAndTwistingKs,0);
        //Ligne 10
        solveLinearSystem(cRod,u,m_dx,m_dq,object,p,0);

        for(uint e = 1; e < cRod.wq().size (); ++e )
        {
            uint a = cRod.beginIdx (e);
            uint z = cRod.endIdx(e);
            m_dx[a].set(0,0,0);m_dx[z].set(0,0,0);
            m_dq[a].coeffs ().setZero ();m_dq[z].coeffs ().setZero ();
            //Ligne 9
            correction(cRod,u,object,p,m_bendingAndTwistingKs,e);
            //Ligne 10
            solveLinearSystem(cRod,u,m_dx,m_dq,object,p,e);

            p[e] += m_dx[e-1];
            u[e].coeffs () += m_dq[e-1].coeffs ();
        }

        u[u.size ()-1] = u[u.size()-2];
    }
}


void PDCosseratRod::solveLinearSystem( PDCosseratRodData& cRod,
                                       std::vector<Quaternionr>& u,
                                       std::vector<vec3> dx,
                                       std::vector<Quaternionr> dq,
                                       PBDObject& object,
                                       WriteCoord& p,
                                       const uint e)
{

    uint a = cRod.beginIdx (e);
    uint z = cRod.endIdx(e);
    const SReal invMass1 = object.invMass(a);
    const SReal invMass0 = object.invMass(z);

    //  COMPUTE STRETCHING AND SHEARING
    vec3 d3;
    d3[0] = static_cast<SReal>(2.0) * (u[a].x() * u[a].z() + u[a].w() * u[a].y());
    d3[1] = static_cast<SReal>(2.0) * (u[a].y() * u[a].z() - u[a].w() * u[a].x());
    d3[2] = u[a].w() * u[a].w() - u[a].x() * u[a].x() - u[a].y() * u[a].y() + u[a].z() * u[a].z();	//third director d3 = q0 * e_3 * q0_conjugate

    vec3 gamma = (p[z] - p[a]) / cRod.length(e) - d3;
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

#include "PBDStrainDynamic.hpp"
#include "../../Common/MathFunctions.hpp"
#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>

int PBDStrainDynamicClass = sofa::core::RegisterObject("Constraint that correct deformed Dynamic.")
                            .add< PBDStrainDynamic >();

using namespace sofa::defaulttype;

using namespace sofa::defaulttype;

void PBDStrainDynamic::bwdInit ()
{
    m_mass = PBDVertexMass<sofa::defaulttype::Vec3Types>(m_mechanicalObject.getValue (),m_topology.getValue ());
    m_basis = PBDTetrahedronBasis(m_mechanicalObject.getValue (),m_topology.getValue ());

    SReal nu = m_poisson_ratio.getValue ();
    SReal E = m_young_modulus.getValue ();
    m_lambda = (E*nu) / (( 1.0 + nu ) * ( 1.0 - 2.0 * nu ));
    m_mu = E / (2.0 * ( 1.0 + nu ));
    //C(i,i) defines x/y/z stretching while K[3] define the intensity of the shear motion
    m_C(0,0) = m_lambda + 2 * m_mu  ; m_C(0,1) = m_lambda         ; m_C(0,2) = m_lambda         ;
    m_C(1,0) = m_lambda           ; m_C(1,1) = m_lambda + 2 * m_mu; m_C(1,2) = m_lambda         ;
    m_C(2,0) = m_lambda           ; m_C(2,1) = m_lambda         ; m_C(2,2) = m_lambda + 2 * m_mu;
    dt2 = static_cast<sofa::simulation::Node*>(this->getContext())->getDt ();
    dt2 *= dt2;
    //Update shear so we don't have to compute it every time
    m_sshear[0] = m_shear.getValue()[0];//xy
    m_sshear[1] = m_shear.getValue()[1];//xz
    m_sshear[2] = m_shear.getValue()[2];//yz
    m_sstretch[0] = m_stretch.getValue()[0] * 0.5;//xy
    m_sstretch[1] = m_stretch.getValue()[1] * 0.5;//xz
    m_sstretch[2] = m_stretch.getValue()[2] * 0.5;//yz
}

void PBDStrainDynamic::solve(sofa::simulation::Node * node)
{
    static SReal eps = 1e-6;
    static SReal one_over_6 = 1.0/6.0;
    const auto& Dm_inv = m_basis.data ();
    const uint tetCount = m_topology.getValue()->getNbTetrahedra ();
    WriteCoord p = m_pbdObject->getFreePosition ();
    for(uint iter = 0; iter < m_nbIter.getValue(); ++iter)
    {
        for(uint tri = 0; tri < tetCount; ++tri)
        {
            const auto& t = m_topology.getValue()->getTetra (tri);
            Matrix3 Ds;
            Ds.x() = p[t[0]] - p[t[3]];
            Ds.y() = p[t[1]] - p[t[3]];
            Ds.z() = p[t[2]] - p[t[3]];

            const Matrix3& F(Ds * Dm_inv[tri].second);
            SReal volume = one_over_6 * dot(Ds.x(),Ds.y().cross(Ds.z()));
            Matrix3 S,piolaK;
            SReal energy;
            //Only goes in if the tetraheadron is inside-out or too flat
            if( volume * Dm_inv[tri].first < 0.0 )
            {
                computeGreenStrainAndPiolaInversion(F,S,piolaK,Dm_inv[tri].first,m_mu,m_lambda,energy,m_sstretch,m_sshear);
            }
            else
            {
                computeGreenStrainAndPiola(F,S,piolaK,Dm_inv[tri].first,m_mu,m_lambda,energy,m_sstretch,m_sshear);
            }

            Vec3 gradC[4];
            computeGradCGreen (Dm_inv[tri].first,Dm_inv[tri].second,piolaK,gradC);
            SReal sumGradSquared =
                    m_mass.w(t[0]) * gradC[0].norm2 ()
                    + m_mass.w(t[1]) * gradC[1].norm2 ()
                    + m_mass.w(t[2]) * gradC[2].norm2 ()
                    + m_mass.w(t[3]) * gradC[3].norm2 ();

            //Check if we reached an equilibrium
            if(sumGradSquared > eps){

                //                gradC[0][0] *= m_stretch.getValue()[0];gradC[0][1] *= m_sshear[0];            gradC[0][2] *= m_sshear[1];
                //                gradC[1][0] *= m_sshear[0];            gradC[1][1] *= m_stretch.getValue()[1];gradC[1][2] *= m_sshear[2];
                //                gradC[0][2] *= m_sshear[1];            gradC[1][2] *= m_sshear[2];            gradC[2][2] *= m_stretch.getValue()[2];

                SReal lagrangeMul = energy / sumGradSquared;
                p[t[0]] -= lagrangeMul * m_mass.w(t[0]) * gradC[0];
                p[t[1]] -= lagrangeMul * m_mass.w(t[1]) * gradC[1];
                p[t[2]] -= lagrangeMul * m_mass.w(t[2]) * gradC[2];
                p[t[3]] -= lagrangeMul * m_mass.w(t[3]) * gradC[3];
            }
        }
    }

}

void PBDStrainDynamic::computeGradCGreen(Real restVolume, const Matrix3 &invRestMat, const Matrix3 &piola, Vec3 *J)
{
    Matrix3 H(piola * invRestMat.transposed () * restVolume);
    J[0] = H.x();
    J[1] = H.y();
    J[2] = H.z();
    J[3] = -J[0] - J[1] - J[2];
}

void PBDStrainDynamic::computeGreenStrainAndPiola(const Matrix3 &F,Matrix3 &S,Matrix3 &P,SReal restVolume, SReal mu, SReal lambda,SReal& energy,Vec3& stretch, Vec3& shear)
{

    S(0, 0) = stretch[0]*(F(0, 0) * F(0, 0) + F(1, 0) * F(1, 0) + F(2, 0) * F(2, 0) - static_cast<Real>(1.0));		// xx
    S(1, 1) = stretch[1]*(F(0, 1) * F(0, 1) + F(1, 1) * F(1, 1) + F(2, 1) * F(2, 1) - static_cast<Real>(1.0));		// yy
    S(2, 2) = stretch[2]*(F(0, 2) * F(0, 2) + F(1, 2) * F(1, 2) + F(2, 2) * F(2, 2) - static_cast<Real>(1.0));		// zz
    S(0, 1) = shear[0]*(F(0, 0) * F(0, 1) + F(1, 0) * F(1, 1) + F(2, 0) * F(2, 1));			// xy
    S(0, 2) = shear[1]*(F(0, 0) * F(0, 2) + F(1, 0) * F(1, 2) + F(2, 0) * F(2, 2));			// xz
    S(1, 2) = shear[2]*(F(0, 1) * F(0, 2) + F(1, 1) * F(1, 2) + F(2, 1) * F(2, 2));			// yz
    S(1, 0) = S(0, 1);
    S(2, 0) = S(0, 2);
    S(2, 1) = S(1, 2);
    //compute Piola-Kirchhoff 1st stress tensor from St Venant Kirchhoff model
    // S = lambda*tr(E)*I + 2mu*E
    const Real trace = S(0, 0) + S(1, 1) + S(2, 2);
    const Real ltrace = lambda*trace;
    P = S * 2.0*mu;
    P(0, 0) += ltrace;
    P(1, 1) += ltrace;
    P(2, 2) += ltrace;
    P = F * P;

    //S symetric so we multiply by 2
    Real E2 = 2.0 * ( S(0,1) * S(0,1) + S(0,2) * S(0,2) + S(1,2) * S(1,2))
              + ( S(0,0) * S(0,0) + S(1,1) * S(1,1) + S(2,2) * S(2,2));
    //Isotropic material strain energy density function
    // psi = 0.5 * lambda * tr(E)^2 + mu*tr(E²)
    //Incompressible tr(E) = 0 -> psi = mu*tr(E^2)
    Real psi = mu*E2 + static_cast<Real>(0.5)*lambda * trace*trace;
    energy = restVolume * psi;
}

void PBDStrainDynamic::computeGreenStrainAndPiolaInversion(const Matrix3 &F, Matrix3 &S,Matrix3 &P,SReal restVolume, SReal mu, SReal lambda,SReal& energy,Vec3& stretch, Vec3& shear)
{
    //compute Piola-Kirchhoff 1st stress tensor from St Venant Kirchhoff model
    // S = lambda*tr(E)*I + 2mu*E
    Matrix3 U, VT;
    Vec3 hatF;
    MathFunctions::svdWithInversionHandling(F, hatF, U, VT);

    // Clamp small singular values
    static Real minXVal = 0.577;

    for (unsigned char j = 0; j < 3; j++)
    {
        if (hatF[j] < minXVal)
            hatF[j] = minXVal;
    }

    // S for hatF
    const Vec3 SHatF( stretch[0]*(hatF[0]*hatF[0] - 1.0), stretch[0]*(hatF[1]*hatF[1] - 1.0),  stretch[0]*(hatF[2]*hatF[2] - 1.0));

    const Real trace = SHatF[0] + SHatF[1] + SHatF[2];
    const Real ltrace = lambda*trace;
    Vec3 PVec = SHatF * 2.0*mu;
    PVec[0] += ltrace;
    PVec[1] += ltrace;
    PVec[2] += ltrace;
    PVec[0] = hatF[0] * PVec[0];
    PVec[1] = hatF[1] * PVec[1];
    PVec[2] = hatF[2] * PVec[2];

    Matrix3 epsVT, sigVT;
    epsVT.x () = SHatF[0] * VT.x ();
    epsVT.y () = SHatF[1] * VT.y ();
    epsVT.z () = SHatF[2] * VT.z ();
    sigVT.x () = PVec[0] * VT.x ();
    sigVT.y () = PVec[1] * VT.y ();
    sigVT.z () = PVec[2] * VT.z ();

    S = U*epsVT;
    S(1, 0) = S(0, 1) = shear[0] * S(0, 1);
    S(2, 0) = S(0, 2) = shear[1] * S(0, 2);
    S(2, 1) = S(1, 2) = shear[2] * S(1, 2);

    P = U*sigVT;

    //S symetric so we multiply by 2
    Real E2 = 2.0 * ( S(0,1) * S(0,1) + S(0,2) * S(0,2) + S(1,2) * S(1,2))
              + ( S(0,0) * S(0,0) + S(1,1) * S(1,1) + S(2,2) * S(2,2));
    //Isotropic material strain energy density function
    // psi = 0.5 * lambda * tr(E)^2 + mu*tr(E²)
    //Incompressible tr(E) = 0 -> psi = mu*tr(E^2)
    Real psi = mu*E2 + static_cast<Real>(0.5)*lambda * trace*trace;
    energy = restVolume * psi;
}

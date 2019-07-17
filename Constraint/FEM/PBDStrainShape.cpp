#include "PBDStrainShape.hpp"

#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>
#include "../../Common/MathFunctions.hpp"

int PBDStrainShapeClass = sofa::core::RegisterObject("Constraint that correct deformed shape.")
                          .add< PBDStrainShape >();

using namespace sofa::defaulttype;

void PBDStrainShape::bwdInit ()
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
}

void PBDStrainShape::solve(sofa::simulation::Node * node)
{
    static SReal eps = 1e-6;
    static SReal one_over_6 = 1.0/6.0;
    const auto& Dm_inv = m_basis.data ();
    const uint tetCount = m_topology.getValue()->getNbTetrahedra ();
    WriteCoord p = m_pbdObject->getFreePosition ();
    for(uint iter = 0; iter < m_nbIter.getValue(); ++iter)
    {
        for(uint i = 0; i < tetCount; ++i)
        {
            const auto& t = m_topology.getValue()->getTetra (i);
            //First compute Ds to get F
            Matrix3 Ds;
            Ds.x() = p[t[0]] - p[t[3]];
            Ds.y() = p[t[1]] - p[t[3]];
            Ds.z() = p[t[2]] - p[t[3]];

            const Matrix3& F = Ds*Dm_inv[i].second;
            SReal volume = one_over_6 * dot(Ds.x(),Ds.y().cross(Ds.z()));
            Matrix3 S,epsilon;
            SReal energy;
            //Only goes in if the tetraheadron is inside-out or too flat
            if( volume / Dm_inv[i].first < 0.08 )
            {
                computeGreenStrainAndPiolaStressInversion(F,Dm_inv[i].first,m_mu,m_lambda,epsilon,S,energy);
            }
            else
            {
                computeGreenStrainAndPiolaStress(F,Dm_inv[i].first,m_mu,m_lambda,epsilon,S,energy);
            }

            Vec3 gradC[4];
            computeGradCGreen (Dm_inv[i].first,Dm_inv[i].second,S,gradC);

            SReal sumGradSquared =   m_mass.w(t[0]) * gradC[0].norm2 ()
                    + m_mass.w(t[1]) * gradC[1].norm2 ()
                    + m_mass.w(t[2]) * gradC[2].norm2 ()
                    + m_mass.w(t[3]) * gradC[3].norm2 ();

            //Check if we reached an equilibrium
            if(sumGradSquared > eps){
                SReal lagrangeMul = energy / sumGradSquared;
                p[t[0]] -= lagrangeMul * m_mass.w(t[0]) * gradC[0];
                p[t[1]] -= lagrangeMul * m_mass.w(t[1]) * gradC[1];
                p[t[2]] -= lagrangeMul * m_mass.w(t[2]) * gradC[2];
                p[t[3]] -= lagrangeMul * m_mass.w(t[3]) * gradC[3];
            }
        }

    }
}


void PBDStrainShape::computeGradCGreen(Real restVolume, const Matrix3 &invRestMat, const Matrix3 &sigma, Vec3 *J)
{
    Matrix3 H(sigma * invRestMat.transposed () * restVolume);
    J[0] = H.x();
    J[1] = H.y();
    J[2] = H.z();
    J[3] = -J[0] - J[1] - J[2];
}

void PBDStrainShape::computeGreenStrainAndPiolaStress(const Matrix3 &F, const Real restVolume,
                                                      const Real mu, const Real lambda,
                                                      Matrix3 &epsilon, Matrix3 &sigma, Real &energy)
{

    epsilon(0, 0) = static_cast<Real>(0.5)*(F(0, 0) * F(0, 0) + F(1, 0) * F(1, 0) + F(2, 0) * F(2, 0) - static_cast<Real>(1.0));		// xx
    epsilon(1, 1) = static_cast<Real>(0.5)*(F(0, 1) * F(0, 1) + F(1, 1) * F(1, 1) + F(2, 1) * F(2, 1) - static_cast<Real>(1.0));		// yy
    epsilon(2, 2) = static_cast<Real>(0.5)*(F(0, 2) * F(0, 2) + F(1, 2) * F(1, 2) + F(2, 2) * F(2, 2) - static_cast<Real>(1.0));		// zz
    epsilon(0, 1) = static_cast<Real>(0.5)*(F(0, 0) * F(0, 1) + F(1, 0) * F(1, 1) + F(2, 0) * F(2, 1));			// xy
    epsilon(0, 2) = static_cast<Real>(0.5)*(F(0, 0) * F(0, 2) + F(1, 0) * F(1, 2) + F(2, 0) * F(2, 2));			// xz
    epsilon(1, 2) = static_cast<Real>(0.5)*(F(0, 1) * F(0, 2) + F(1, 1) * F(1, 2) + F(2, 1) * F(2, 2));			// yz
    epsilon(1, 0) = epsilon(0, 1);
    epsilon(2, 0) = epsilon(0, 2);
    epsilon(2, 1) = epsilon(1, 2);
    //compute Piola-Kirchhoff 1st stress tensor from St Venant Kirchhoff model
    // S = lambda*tr(E)*I + 2mu*E
    const Real trace = epsilon(0, 0) + epsilon(1, 1) + epsilon(2, 2);
    const Real ltrace = lambda*trace;
    sigma = epsilon * 2.0*mu;
    sigma(0, 0) += ltrace;
    sigma(1, 1) += ltrace;
    sigma(2, 2) += ltrace;
    sigma = F * sigma;

    //epsilon symetric so we multiply by 2
    Real E2 = 2.0 * ( epsilon(0,1) * epsilon(0,1) + epsilon(0,2) * epsilon(0,2) + epsilon(1,2) * epsilon(1,2))
              + ( epsilon(0,0) * epsilon(0,0) + epsilon(1,1) * epsilon(1,1) + epsilon(2,2) * epsilon(2,2));
    //Isotropic material strain energy density function
    // psi = 0.5 * lambda * tr(E)^2 + mu*tr(E²)
    //Incompressible tr(E) = 0 -> psi = mu*tr(E^2)
    Real psi = mu*E2 + static_cast<Real>(0.5)*lambda * trace*trace;
    energy = restVolume * psi;
}

void PBDStrainShape::computeGreenStrainAndPiolaStressInversion(const Matrix3 &F, const Real restVolume,
                                                               const Real mu, const Real lambda,
                                                               Matrix3 &epsilon, Matrix3 &sigma, Real &energy)
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

    // epsilon for hatF
    const Vec3 epsilonHatF(	0.5*(hatF[0]*hatF[0] - 1.0),  0.5*(hatF[1]*hatF[1] - 1.0), 0.5*(hatF[2]*hatF[2] - 1.0));

    const Real trace = epsilonHatF[0] + epsilonHatF[1] + epsilonHatF[2];
    const Real ltrace = lambda*trace;
    Vec3 sigmaVec = epsilonHatF * 2.0*mu;
    sigmaVec[0] += ltrace;
    sigmaVec[1] += ltrace;
    sigmaVec[2] += ltrace;
    sigmaVec[0] = hatF[0] * sigmaVec[0];
    sigmaVec[1] = hatF[1] * sigmaVec[1];
    sigmaVec[2] = hatF[2] * sigmaVec[2];

    Matrix3 epsVT, sigVT;
    epsVT.x () = epsilonHatF[0] * VT.x ();
    epsVT.y () = epsilonHatF[1] * VT.y ();
    epsVT.z () = epsilonHatF[2] * VT.z ();
    sigVT.x () = sigmaVec[0] * VT.x ();
    sigVT.y () = sigmaVec[1] * VT.y ();
    sigVT.z () = sigmaVec[2] * VT.z ();

    epsilon = U*epsVT;
    sigma = U*sigVT;

    //epsilon symetric so we multiply by 2
    Real E2 = 2.0 * ( epsilon(0,1) * epsilon(0,1) + epsilon(0,2) * epsilon(0,2) + epsilon(1,2) * epsilon(1,2))
              + ( epsilon(0,0) * epsilon(0,0) + epsilon(1,1) * epsilon(1,1) + epsilon(2,2) * epsilon(2,2));
    //Isotropic material strain energy density function
    // psi = 0.5 * lambda * tr(E)^2 + mu*tr(E²)
    //Incompressible tr(E) = 0 -> psi = mu*tr(E^2)
    Real psi = mu*E2 + static_cast<Real>(0.5)*lambda * trace*trace;
    energy = restVolume * psi;
}

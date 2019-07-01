#include "PBDStrainShape.hpp"

#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>
#include "../../Common/MathFunctions.hpp"

int PBDStrainShapeClass = sofa::core::RegisterObject("Constraint that correct deformed shape.")
                          .add< PBDStrainShape >();

using namespace sofa::defaulttype;

void PBDStrainShape::bwdInit ()
{
    SReal p = m_poisson_ratio.getValue ();
    SReal E = m_young_modulus.getValue ();
    m_lambda = (E*p) / (( 1.0 + p ) * ( 1.0 - 2.0 * p ));
    m_mu = E / (2.0 * ( 1.0 + p ));
    //C(i,i) defines x/y/z stretching while K[3] define the intensity of the shear motion
    m_C(0,0) = m_lambda + 2 * m_mu  ; m_C(0,1) = m_lambda         ; m_C(0,2) = m_lambda         ;
    m_C(1,0) = m_lambda           ; m_C(1,1) = m_lambda + 2 * m_mu; m_C(1,2) = m_lambda         ;
    m_C(2,0) = m_lambda           ; m_C(2,1) = m_lambda         ; m_C(2,2) = m_lambda + 2 * m_mu;
    dt2 = static_cast<sofa::simulation::Node*>(this->getContext())->getDt ();
    dt2 *= dt2;
}

void PBDStrainShape::solve(PBDObject<sofa::defaulttype::Vec3Types> &object, WriteCoord &x)
{
    static SReal eps = 1e-6;
    if(!object.hasDataType(PBDObject<sofa::defaulttype::Vec3Types>::TETRAHEDRON))
    {
        object.computeTetrahedraBasis ();
    }
    const auto& Dm_inv = object.tetrahedraBases ().data ();
    const uint tetCount = object.sofaTopology ()->getNbTetrahedra ();

    for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
    {
        for(uint i = 0; i < tetCount; ++i)
        {
            const auto& t = object.sofaTopology ()->getTetra (i);
            //First compute Ds to get F
            const auto& r1 = x[t[0]] - x[t[3]];
            const auto& r2 = x[t[1]] - x[t[3]];
            const auto& r3 = x[t[2]] - x[t[3]];
            Matrix3 Ds;
            Ds.x () = r1;
            Ds.y () = r2;
            Ds.z () = r3;

            const Matrix3& F = Ds*Dm_inv[i].second;

            //From F compute Green Strain Tensor { 0.5 (F^t * F - I) }
            Matrix3 epsilon = 0.5 * ( F.transposed() * F); //What happend if F is orthogonal ? :-> if F is orthogonal then Ds is a pure rotation of Dm hence no strain is applied
            epsilon(0,0) -= 0.5;epsilon(1,1) -= 0.5;epsilon(2,2) -= 0.5;

            //compute Piola-Kirchhoff stress tensor
            const Matrix3& P = F*m_C*epsilon;

            //Then the gradient and divergence (which are actually displacement since m_C is premultiplied by dt^2)
            const Matrix3& displacement = dt2 * Dm_inv[i].first * (P * Dm_inv[i].second.transposed());

            x[t[0]] += object.invMass(t[0]) * displacement.x ();
            x[t[1]] += object.invMass(t[1]) * displacement.y ();
            x[t[2]] += object.invMass(t[2]) * displacement.z ();
            x[t[3]] -= object.invMass(t[3]) * displacement.x () + displacement.y () + displacement.z ();
        }
    }
}


void PBDStrainShape::computeGradCGreen(Real restVolume, const Matrix3 &invRestMat, const Matrix3 &sigma, Vec3 *J)
{
    Matrix3 H(sigma * invRestMat.transposed () * restVolume);

    J[0][0] = H(0, 0);
    J[1][0] = H(0, 1);
    J[2][0] = H(0, 2);

    J[0][1] = H(1, 0);
    J[1][1] = H(1, 1);
    J[2][1] = H(1, 2);

    J[0][2] = H(2, 0);
    J[1][2] = H(2, 1);
    J[2][2] = H(2, 2);

    J[3] = -J[0] - J[1] - J[2];
}

void PBDStrainShape::computeGreenStrainAndPiolaStress(const Matrix3 &F, const Real restVolume,
                                                      const Real mu, const Real lambda,
                                                      Matrix3 &epsilon, Matrix3 &sigma, Real &energy)
{
    // Determine \partial x/\partial m_i

    // epsilon = 1/2 F^T F - I
    epsilon(0, 0) = static_cast<Real>(0.5)*(F(0, 0) * F(0, 0) + F(1, 0) * F(1, 0) + F(2, 0) * F(2, 0) - static_cast<Real>(1.0));		// xx
    epsilon(1, 1) = static_cast<Real>(0.5)*(F(0, 1) * F(0, 1) + F(1, 1) * F(1, 1) + F(2, 1) * F(2, 1) - static_cast<Real>(1.0));		// yy
    epsilon(2, 2) = static_cast<Real>(0.5)*(F(0, 2) * F(0, 2) + F(1, 2) * F(1, 2) + F(2, 2) * F(2, 2) - static_cast<Real>(1.0));		// zz
    epsilon(0, 1) = static_cast<Real>(0.5)*(F(0, 0) * F(0, 1) + F(1, 0) * F(1, 1) + F(2, 0) * F(2, 1));			// xy
    epsilon(0, 2) = static_cast<Real>(0.5)*(F(0, 0) * F(0, 2) + F(1, 0) * F(1, 2) + F(2, 0) * F(2, 2));			// xz
    epsilon(1, 2) = static_cast<Real>(0.5)*(F(0, 1) * F(0, 2) + F(1, 1) * F(1, 2) + F(2, 1) * F(2, 2));			// yz
    epsilon(1, 0) = epsilon(0, 1);
    epsilon(2, 0) = epsilon(0, 2);
    epsilon(2, 1) = epsilon(1, 2);

    // P(F) = F(2 mu E + lambda tr(E)I) => E = green strain
    const Real trace = epsilon(0, 0) + epsilon(1, 1) + epsilon(2, 2);
    const Real ltrace = lambda*trace;
    sigma = epsilon * 2.0*mu;
    sigma(0, 0) += ltrace;
    sigma(1, 1) += ltrace;
    sigma(2, 2) += ltrace;
    sigma = F * sigma;

    Real psi = 0.0;
    for (unsigned char j = 0; j < 3; j++)
        for (unsigned char k = 0; k < 3; k++)
            psi += epsilon(j, k) * epsilon(j, k);
    psi = mu*psi + static_cast<Real>(0.5)*lambda * trace*trace;
    energy = restVolume * psi;
}

void PBDStrainShape::computeGreenStrainAndPiolaStressInversion(const Matrix3 &F, const Real restVolume,
                                                               const Real mu, const Real lambda,
                                                               Matrix3 &epsilon, Matrix3 &sigma, Real &energy)
{
    Matrix3 U, VT;
    Vec3 hatF;
    MathFunctions::svdWithInversionHandling(F, hatF, U, VT);

//    // Clamp small singular values
//    const Real minXVal = static_cast<Real>(0.577);

//    for (unsigned char j = 0; j < 3; j++)
//    {
//        if (hatF[j] < minXVal)
//            hatF[j] = minXVal;
//    }

    // epsilon for hatF
    Vec3 epsilonHatF(	static_cast<Real>(0.5)*(hatF[0]*hatF[0] - static_cast<Real>(1.0)),
            static_cast<Real>(0.5)*(hatF[1]*hatF[1] - static_cast<Real>(1.0)),
            static_cast<Real>(0.5)*(hatF[2]*hatF[2] - static_cast<Real>(1.0)));

    const Real trace = epsilonHatF[0] + epsilonHatF[1] + epsilonHatF[2];
    const Real ltrace = lambda*trace;
    Vec3 sigmaVec = epsilonHatF * 2.0*mu;
    sigmaVec[0] += ltrace;
    sigmaVec[1] += ltrace;
    sigmaVec[2] += ltrace;
    sigmaVec[0] = hatF[0] * sigmaVec[0];
    sigmaVec[1] = hatF[1] * sigmaVec[1];
    sigmaVec[2] = hatF[2] * sigmaVec[2];

    Matrix3 sigmaDiag, epsDiag;

    sigmaDiag.x() = Vec3(sigmaVec[0], 0.0, 0.0);
    sigmaDiag.y() = Vec3(0.0, sigmaVec[1], 0.0);
    sigmaDiag.z() = Vec3(0.0, 0.0, sigmaVec[2]);

    epsDiag.x() = Vec3(epsilonHatF[0], 0.0, 0.0);
    epsDiag.y() = Vec3(0.0, epsilonHatF[1], 0.0);
    epsDiag.z() = Vec3(0.0, 0.0, epsilonHatF[2]);

    epsilon = U*epsDiag*VT;
    sigma = U*sigmaDiag*VT;
    Real psi = 0.0;
    for (unsigned char j = 0; j < 3; j++)
        for (unsigned char k = 0; k < 3; k++)
            psi += epsilon(j, k) * epsilon(j, k);
    psi = mu*psi + static_cast<Real>(0.5)*lambda * trace*trace;
    energy = restVolume * psi;
}

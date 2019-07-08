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
    static SReal eps = 1e-8;
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
            Matrix3 Ds;
            Ds.x() = x[t[0]] - x[t[3]];
            Ds.y() = x[t[1]] - x[t[3]];
            Ds.z() = x[t[2]] - x[t[3]];

            const Matrix3& F = Ds*Dm_inv[i].second;
            SReal volume = -0.166666666666666666 * dot(Ds.x(),Ds.y().cross(Ds.z()));
            Matrix3 S,epsilon;
            SReal energy;
            //Only goes in if the tetraheadron is inside-out
            if( volume / Dm_inv[i].first < 0.2 )
            {
                computeGreenStrainAndPiolaStressInversion(F,Dm_inv[i].first,m_mu,m_lambda,epsilon,S,energy);
            }
            else
            {
                computeGreenStrainAndPiolaStress(F,Dm_inv[i].first,m_mu,m_lambda,epsilon,S,energy);
            }

            Vec3 gradC[4];
            computeGradCGreen (Dm_inv[i].first,Dm_inv[i].second,S,gradC);

            SReal sumGradSquared =   object.invMass(t[0]) * gradC[0].norm2 ()
                    + object.invMass(t[1]) * gradC[1].norm2 ()
                    + object.invMass(t[2]) * gradC[2].norm2 ()
                    + object.invMass(t[3]) * gradC[3].norm2 ();

            //Check if we reached an equilibrium
            if(sumGradSquared > eps){
                SReal lagrangeMul = energy / sumGradSquared;
                x[t[0]] -= lagrangeMul * object.invMass(t[0]) * gradC[0];
                x[t[1]] -= lagrangeMul * object.invMass(t[1]) * gradC[1];
                x[t[2]] -= lagrangeMul * object.invMass(t[2]) * gradC[2];
                x[t[3]] -= lagrangeMul * object.invMass(t[3]) * gradC[3];
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

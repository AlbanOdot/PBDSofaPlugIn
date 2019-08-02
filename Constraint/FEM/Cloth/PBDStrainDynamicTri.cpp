#include "PBDStrainDynamicTri.hpp"

#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>
#include "../../../Common/MathFunctions.hpp"

int PBDStrainDynamicTriClass = sofa::core::RegisterObject("Constraint that correct deformed Dynamic.")
                               .add< PBDStrainDynamicTri >();

using namespace sofa::defaulttype;

void PBDStrainDynamicTri::bwdInit ()
{
    m_mass = PBDVertexMass<sofa::defaulttype::Vec3Types>(m_mechanicalObject.getValue (),m_topology.getValue ());
    m_basis = PBDTriDynamicBasis(m_mechanicalObject.getValue (),m_topology.getValue ());

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
    m_stretch_and_shear[0] = m_sshear.getValue ()[0];
    m_stretch_and_shear[1] = m_sshear.getValue ()[1];
    m_stretch_and_shear[2] = m_sshear.getValue ()[2];
}

void PBDStrainDynamicTri::solve(sofa::simulation::Node * node)
{
    static SReal eps = 1e-6;
    const auto& Dm_inv = m_basis.data ();
    const uint triCount = m_topology.getValue()->getNbTriangles ();
    WriteCoord p = m_pbdObject->getFreePosition ();
    for(uint iter = 0; iter < m_nbIter.getValue(); ++iter)
    {
        for(uint tri = 0; tri < triCount; ++tri)
        {

            const auto& t = m_topology.getValue()->getTriangle(tri);
            // Determine \partial x/\partial m_i
            sofa::defaulttype::Mat<3,2,SReal> F;
            const Vec3 p13 = p[t[0]] - p[t[2]];
            const Vec3 p23 = p[t[1]] - p[t[2]];
            F(0,0) = p13[0] * Dm_inv[tri].second(0,0) + p23[0] * Dm_inv[tri].second(1,0);
            F(0,1) = p13[0] * Dm_inv[tri].second(0,1) + p23[0] * Dm_inv[tri].second(1,1);
            F(1,0) = p13[1] * Dm_inv[tri].second(0,0) + p23[1] * Dm_inv[tri].second(1,0);
            F(1,1) = p13[1] * Dm_inv[tri].second(0,1) + p23[1] * Dm_inv[tri].second(1,1);
            F(2,0) = p13[2] * Dm_inv[tri].second(0,0) + p23[2] * Dm_inv[tri].second(1,0);
            F(2,1) = p13[2] * Dm_inv[tri].second(0,1) + p23[2] * Dm_inv[tri].second(1,1);

            // epsilon = F^T * F - I
            //Here epsilon == S because if you do the math we have the same computations
            Matrix2 S;
            S(0,0) = (F(0,0) * F(0,0) + F(1,0) * F(1,0) + F(2,0) * F(2,0) - static_cast<Real>(1.0)) * m_stretch_and_shear[0];		// xx
            S(1,1) = (F(0,1) * F(0,1) + F(1,1) * F(1,1) + F(2,1) * F(2,1) - static_cast<Real>(1.0)) * m_stretch_and_shear[1];		// yy
            S(0,1) = (F(0,0) * F(0,1) + F(1,0) * F(1,1) + F(2,0) * F(2,1)) * m_stretch_and_shear[2];			// xy
            S(1,0) = S(0,1);

            SReal psi = 0.5* (S(0,0) * S(0,0) + S(0,1) * S(0,1) + S(1,0) * S(1,0) + S(1,1) * S(1,1)) ;
            SReal energy =  Dm_inv[tri].first*psi;

            // compute gradient -> F * S == The first Piola-Kirchhoff stress tensor
            sofa::defaulttype::Mat<3,2,SReal> H = Dm_inv[tri].first * F * S * Dm_inv[tri].second.transposed ();
            Vec3 gradC3 = -H.col(0) - H.col(1);

            Real sum_normGradC = m_mass.w(t[0]) * H.col(0).norm2() + m_mass.w(t[1]) *H.col(1).norm2() + m_mass.w(t[2]) * gradC3.norm2();

            // exit early if required
            if (sum_normGradC > eps)
            {
                // compute scaling factor
                const Real s = energy / sum_normGradC;

                // update positions
                p[t[0]] -= (s*m_mass.w(t[0])) * H.col(0);
                p[t[1]] -= (s*m_mass.w(t[1])) * H.col(1);
                p[t[2]] -= (s*m_mass.w(t[2])) * gradC3;

            }
        }

    }
}


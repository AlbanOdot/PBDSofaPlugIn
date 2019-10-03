#include "PBDFEMTri.hpp"

#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>
#include "../../Common/MathFunctions.hpp"
#include <sofa/core/visual/VisualParams.h>
#include <sofa/defaulttype/VecTypes.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/core/behavior/ForceField.inl>
#include <SofaBaseTopology/TopologyData.inl>
#include <omp.h>

int PBDFEMTriClass = sofa::core::RegisterObject("Constraint that correct deformed shape.")
                     .add< PBDFEMTri >();

using namespace sofa::defaulttype;
using namespace	sofa::component::topology;
using namespace sofa;
using namespace core::topology;


void PBDFEMTri::bwdInit ()
{
    m_mass = PBDVertexMass<sofa::defaulttype::Vec3Types>(m_mechanicalObject.getValue (),m_topology.getValue ());
    m_basis = PBDTriBasis(m_mechanicalObject.getValue (),m_topology.getValue ());

    SReal nu = m_poisson_ratio.getValue ();
    SReal E = m_young_modulus.getValue ();
    m_lambda = (E*nu) / (( 1.0 + nu ) * ( 1.0 - 2.0 * nu ));
    m_mu = E / (2.0 * ( 1.0 + nu ));
}

bool PBDFEMTri::solve(sofa::simulation::Node * node)
{
    static SReal eps = 1e-6;
    static SReal one_over_6 = 1.0/6.0;
    const auto& Dm_inv = m_basis.data ();
    const uint triCount = m_topology.getValue()->getNbTriangles ();
    WriteCoord p = m_pbdObject->getFreePosition ();
    bool modification = false;
    for(uint tri = 0; tri < triCount; ++tri)
    {
        const auto& t = m_topology.getValue()->getTriangle (tri);
        const Matrix2& invRestMat = Dm_inv[tri].second;
        const SReal area = Dm_inv[tri].first;
        Vec3& p0 = p[t[0]];
        Vec3& p1 = p[t[1]];
        Vec3& p2 = p[t[2]];
        // Determine \partial x/\partial m_i
        Mat<3,2,SReal> F;
        const Vec3& p13 = p0 - p2;
        const Vec3& p23 = p1 - p2;
        F(0,0) = p13[0] * invRestMat(0,0) + p23[0] * invRestMat(1,0);
        F(0,1) = p13[0] * invRestMat(0,1) + p23[0] * invRestMat(1,1);
        F(1,0) = p13[1] * invRestMat(0,0) + p23[1] * invRestMat(1,0);
        F(1,1) = p13[1] * invRestMat(0,1) + p23[1] * invRestMat(1,1);
        F(2,0) = p13[2] * invRestMat(0,0) + p23[2] * invRestMat(1,0);
        F(2,1) = p13[2] * invRestMat(0,1) + p23[2] * invRestMat(1,1);

        // epsilon = 0.5(F^T * F - I)
        Matrix2 epsilon;
        epsilon(0,0) = static_cast<Real>(0.5)*(F(0,0) * F(0,0) + F(1,0) * F(1,0) + F(2,0) * F(2,0) - static_cast<Real>(1.0));		// xx
        epsilon(1,1) = static_cast<Real>(0.5)*(F(0,1) * F(0,1) + F(1,1) * F(1,1) + F(2,1) * F(2,1) - static_cast<Real>(1.0));		// yy
        epsilon(0,1) = static_cast<Real>(0.5)*(F(0,0) * F(0,1) + F(1,0) * F(1,1) + F(2,0) * F(2,1));			// xy
        epsilon(1,0) = epsilon(0,1);

        // P(F) = det(F) * C*E * F^-T => E = green strain
        Matrix2 sigma;
        const Real trace = epsilon(0, 0) + epsilon(1, 1);
        const Real ltrace = m_lambda*trace;
        sigma = epsilon * 2.0*m_mu;
        sigma(0, 0) += ltrace;
        sigma(1, 1) += ltrace;
        const Mat<3,2,SReal>& piolaKirchhoffStress = F * sigma;

        Real psi = 0.0;
        for (unsigned char j = 0; j < 2; j++)
            for (unsigned char k = 0; k < 2; k++)
                psi += epsilon(j,k) * sigma(j,k);
        psi = static_cast<Real>(0.5)*psi;
        Real energy = area*psi;

        // compute gradient
        const Mat<3,2,SReal>& H = area * piolaKirchhoffStress * invRestMat.transposed ();
        Vec3 gradC[3];
        for (unsigned char j = 0; j < 3; ++j)
        {
            gradC[0][j] = H(j,0);
            gradC[1][j] = H(j,1);
        }
        gradC[2] = -gradC[0] - gradC[1];


        Real sum_normGradC = m_mass.w(t[0]) * gradC[0].norm2();
        sum_normGradC += m_mass.w(t[1]) * gradC[1].norm2();
        sum_normGradC += m_mass.w(t[2]) * gradC[2].norm2();

        // exit early if required
        if (fabs(sum_normGradC) > eps)
        {
            // compute scaling factor
            const Real s = energy / sum_normGradC;
            // update positions
            p0 -= (s*m_mass.w(t[0])) * gradC[0];
            p1 -= (s*m_mass.w(t[1])) * gradC[1];
            p2 -= (s*m_mass.w(t[2])) * gradC[2];

            modification = true;
        }
        modification |= false;
    }
    return modification;
}

#include "PBDStrainTri.hpp"
#include "../../Common/MathFunctions.hpp"
#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/defaulttype/VecTypes.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/core/behavior/ForceField.inl>
#include <SofaBaseTopology/TopologyData.inl>

int PBDStrainTriClass = sofa::core::RegisterObject("Constraint that correct deformed Dynamic.")
                        .add< PBDStrainTri >();

using namespace sofa::defaulttype;
using namespace	sofa::component::topology;
using namespace sofa;
using namespace core::topology;
void PBDStrainTri::bwdInit ()
{
    m_mass = PBDVertexMass<sofa::defaulttype::Vec3Types>(m_mechanicalObject.getValue (),m_topology.getValue ());
    m_basis = PBDTriDynamicBasis(m_mechanicalObject.getValue (),m_topology.getValue ());

    SReal nu = m_poisson_ratio.getValue ();
    SReal E = m_young_modulus.getValue ();
    m_lambda = (E*nu) / (( 1.0 + nu ) * ( 1.0 - 2.0 * nu ));
    m_mu = E / (2.0 * ( 1.0 + nu ));
    dt2 = static_cast<sofa::simulation::Node*>(this->getContext())->getDt ();
    dt2 *= dt2;
    //Update shear so we don't have to compute it every time
    m_sshear[0] = m_shear.getValue()[0];//xx
    m_sshear[1] = m_shear.getValue()[1];//yy
    m_sshear[2] = m_shear.getValue()[2];//xy
}

bool PBDStrainTri::solve(sofa::simulation::Node * node)
{
    static SReal eps = 1e-6;
    static bool normalizeShear = false;
    static bool normalizeStretch = false;
    const auto& Dm_inv = m_basis.data ();
    const uint triCount = m_topology.getValue()->getNbTriangles ();
    WriteCoord p = m_pbdObject->getFreePosition ();
    bool modification = false;

    for(uint tri = 0; tri < triCount; ++tri)
    {
        const auto& t = m_topology.getValue()->getTriangle (tri);
        Vec3& p0 = p[t[0]];
        Vec3& p1 = p[t[1]];
        Vec3& p2 = p[t[2]];

        const Matrix2& invRestMat = Dm_inv[tri].second;
        Vec3 c[2];
        c[0] = Vec3(invRestMat(0, 0), invRestMat(1, 0), 0.0);
        c[1] = Vec3(invRestMat(0, 1), invRestMat(1, 1), 0.0);

        Vec3 r[3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j <= i; j++) {

                r[0] = Vec3(p1[0] - p0[0], p2[0] - p0[0], 0.0);
                r[1] = Vec3(p1[1] - p0[1], p2[1] - p0[1], 0.0);
                r[2] = Vec3(p1[2] - p0[2], p2[2] - p0[2], 0.0);


                Real Sij = 0.0;
                for (int k = 0; k < 3; k++)
                    Sij += dot(r[k],c[i]) * dot(r[k],c[j]);

                Vec3 d[3];
                d[0] = Vec3(0.0, 0.0, 0.0);

                for (int k = 0; k < 2; k++) {
                    d[k+1]  = Vec3(dot(r[0],c[i]), dot(r[1],c[i]), dot(r[2],c[i])) * invRestMat(k, j);
                    d[k+1] += Vec3(dot(r[0],c[j]), dot(r[1],c[j]), dot(r[2],c[j])) * invRestMat(k, i);
                    d[0] -= d[k+1];
                }

                if (i != j && normalizeShear) {
                    Real fi2 = 0.0;
                    Real fj2 = 0.0;
                    for (int k = 0; k < 3; k++) {
                        fi2 += dot(r[k],c[i]) * dot(r[k],c[i]);
                        fj2 += dot(r[k],c[j]) * dot(r[k],c[j]);
                    }
                    Real fi = sqrt(fi2);
                    Real fj = sqrt(fj2);

                    d[0] = Vec3(0.0, 0.0, 0.0);
                    Real s = Sij / (fi2*fi*fj2*fj);
                    for (int k = 0; k < 2; k++) {
                        d[k+1] /= fi * fj;
                        d[k+1] -= fj*fj * Vec3(dot(r[0],c[i]), dot(r[1],c[i]), dot(r[2],c[i])) * invRestMat(k, i) * s;
                        d[k+1] -= fi*fi * Vec3(dot(r[0],c[j]), dot(r[1],c[j]), dot(r[2],c[j])) * invRestMat(k, j) * s;
                        d[0] -= d[k+1];
                    }
                    Sij = Sij / (fi * fj);
                }

                Real lambda =
                        m_mass.w(t[0]) * d[0].norm2() +
                        m_mass.w(t[1]) * d[1].norm2() +
                        m_mass.w(t[2]) * d[2].norm2();

                if (lambda <= eps)
                    continue;

                if (i == 0 && j == 0) {
                    if (normalizeStretch) {
                        Real s = sqrt(Sij);
                        lambda = static_cast<Real>(2.0) * s * (s - static_cast<Real>(1.0)) / lambda * m_sshear[0];
                    }
                    else {
                        lambda = (Sij - static_cast<Real>(1.0)) / lambda * m_sshear[0];
                    }
                }
                else if (i == 1 && j == 1) {
                    if (normalizeStretch) {
                        Real s = sqrt(Sij);
                        lambda = static_cast<Real>(2.0) * s * (s - static_cast<Real>(1.0)) / lambda * m_sshear[1];
                    }
                    else {
                        lambda = (Sij - static_cast<Real>(1.0)) / lambda * m_sshear[1];
                    }
                }
                else {
                    lambda = Sij / lambda * m_sshear[2];
                }

                p0 -= lambda * m_mass.w(t[0]) * d[0];
                p1 -= lambda * m_mass.w(t[1]) * d[1];
                p2 -= lambda * m_mass.w(t[2]) * d[2];
                modification = true;
            }
            modification |= false;
        }
    }
    return  modification;
}


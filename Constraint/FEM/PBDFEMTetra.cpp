#include "PBDFEMTetra.hpp"

#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>
#include "../../Common/MathFunctions.hpp"
#include <sofa/core/visual/VisualParams.h>
#include <sofa/defaulttype/VecTypes.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/core/behavior/ForceField.inl>
#include <SofaBaseTopology/TopologyData.inl>

int PBDFEMTetraClass = sofa::core::RegisterObject("Constraint that correct deformed shape.")
                          .add< PBDFEMTetra >();

using namespace sofa::defaulttype;
using namespace	sofa::component::topology;
using namespace sofa;
using namespace core::topology;


void PBDFEMTetra::bwdInit ()
{
    m_mass = PBDVertexMass<sofa::defaulttype::Vec3Types>(m_mechanicalObject.getValue (),m_topology.getValue ());
    m_basis = FEMTetraData(m_mechanicalObject.getValue (),m_topology.getValue ());

    SReal nu = m_poisson_ratio.getValue ();
    SReal E = m_young_modulus.getValue ();
    m_lambda = (E*nu) / (( 1.0 + nu ) * ( 1.0 - 2.0 * nu ));
    m_mu = E / (2.0 * ( 1.0 + nu ));
    dt2 = static_cast<sofa::simulation::Node*>(this->getContext())->getDt ();
    dt2 *= dt2;
}

void PBDFEMTetra::solve(sofa::simulation::Node * node)
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
            // Determine \partial x/\partial m_i
            const Matrix3& invRestMat = Dm_inv[i].second;
            Matrix3 F;
            const Vec3 p14 = p[t[0]] - p[t[3]];
            const Vec3 p24 = p[t[1]] - p[t[3]];
            const Vec3 p34 = p[t[2]] - p[t[3]];
            F(0, 0) = p14[0]*invRestMat(0, 0) + p24[0]*invRestMat(1, 0) + p34[0]*invRestMat(2, 0);
            F(0, 1) = p14[0]*invRestMat(0, 1) + p24[0]*invRestMat(1, 1) + p34[0]*invRestMat(2, 1);
            F(0, 2) = p14[0]*invRestMat(0, 2) + p24[0]*invRestMat(1, 2) + p34[0]*invRestMat(2, 2);

            F(1, 0) = p14[1]*invRestMat(0, 0) + p24[1]*invRestMat(1, 0) + p34[1]*invRestMat(2, 0);
            F(1, 1) = p14[1]*invRestMat(0, 1) + p24[1]*invRestMat(1, 1) + p34[1]*invRestMat(2, 1);
            F(1, 2) = p14[1]*invRestMat(0, 2) + p24[1]*invRestMat(1, 2) + p34[1]*invRestMat(2, 2);

            F(2, 0) = p14[2]*invRestMat(0, 0) + p24[2]*invRestMat(1, 0) + p34[2]*invRestMat(2, 0);
            F(2, 1) = p14[2]*invRestMat(0, 1) + p24[2]*invRestMat(1, 1) + p34[2]*invRestMat(2, 1);
            F(2, 2) = p14[2]*invRestMat(0, 2) + p24[2]*invRestMat(1, 2) + p34[2]*invRestMat(2, 2);
            SReal volume = one_over_6 * dot(p14,p24.cross(p34));
            Matrix3 piolaKirchhoff,epsilon;
            SReal energy;
            //Only goes in if the tetraheadron is inside-out or too flat (i.e 20% volume remaining)
            if( volume / Dm_inv[i].first < 0.2 )
            {
                computeGreenStrainAndPiolaStressInversion(F,Dm_inv[i].first,m_mu,m_lambda,epsilon,piolaKirchhoff,energy);
            }
            else
            {
                computeGreenStrainAndPiolaStress(F,Dm_inv[i].first,m_mu,m_lambda,epsilon,piolaKirchhoff,energy);
            }
            Vec3 gradC[4];
            computeGradCGreen (Dm_inv[i].first,Dm_inv[i].second,piolaKirchhoff,gradC);
            SReal sumGradSquared =
                    m_mass.w(t[0]) * gradC[0].norm2 ()
                    + m_mass.w(t[1]) * gradC[1].norm2 ()
                    + m_mass.w(t[2]) * gradC[2].norm2 ()
                    + m_mass.w(t[3]) * gradC[3].norm2 ();

            //Check if we reached an equilibrium
            if(sumGradSquared > eps){
                SReal lagrangeMul = energy / sumGradSquared;
                p[t[0]] -= (lagrangeMul * m_mass.w(t[0])) * gradC[0];
                p[t[1]] -= (lagrangeMul * m_mass.w(t[1])) * gradC[1];
                p[t[2]] -= (lagrangeMul * m_mass.w(t[2])) * gradC[2];
                p[t[3]] -= (lagrangeMul * m_mass.w(t[3])) * gradC[3];
            }

            if(m_volumeConservation.getValue())
            {
                Matrix3 Ds;
                Ds.x() = p[t[0]] - p[t[3]];
                Ds.y() = p[t[1]] - p[t[3]];
                Ds.z() = p[t[2]] - p[t[3]];

                Vec3 p2cp3 = Ds.y().cross(Ds.z());
                //            actual volume - rest volume
                SReal C = dot(Ds.x(),p2cp3) - 6.0*Dm_inv[i].first;
                std::cout << " C : "<<C<<std::endl;
                //gradient
                gradC[1] = p2cp3;
                gradC[2] = Ds.z().cross(Ds.x());
                gradC[3] = Ds.x().cross(Ds.y());
                gradC[0] = - gradC[1] - gradC[2] - gradC[3];
                SReal sumGradSquared =
                        m_mass.w(t[0]) * gradC[0].norm2 ()
                        + m_mass.w(t[1]) * gradC[1].norm2 ()
                        + m_mass.w(t[2]) * gradC[2].norm2 ()
                        + m_mass.w(t[3]) * gradC[3].norm2 ();
                if(sumGradSquared > eps)
                {
                    SReal lambda = C/sumGradSquared;
                    p[t[0]] -= (lambda * m_mass.w(t[0])) * gradC[0];
                    p[t[1]] -= (lambda * m_mass.w(t[1])) * gradC[1];
                    p[t[2]] -= (lambda * m_mass.w(t[2])) * gradC[2];
                    p[t[3]] -= (lambda * m_mass.w(t[3])) * gradC[3];
                }
            }

        }
    }
}


void PBDFEMTetra::computeGradCGreen(Real restVolume, const Matrix3 &invRestMat, const Matrix3 &sigma, Vec3 *J)
{
    Matrix3 H;
    Matrix3 T;
    T = invRestMat.transposed ();
    H = sigma * T * restVolume;

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

void PBDFEMTetra::computeGreenStrainAndPiolaStress(const Matrix3 &F, const Real restVolume,
                                                      const Real mu, const Real lambda,
                                                      Matrix3 &epsilon, Matrix3 &sigma, Real &energy)
{

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
    //compute Piola-Kirchhoff 1st stress tensor from St Venant Kirchhoff model
    // S = lambda*tr(E)*I + 2mu*E
    const Real trace = epsilon(0, 0) + epsilon(1, 1) + epsilon(2, 2);
    const Real ltrace = lambda*trace;
    sigma = epsilon * 2.0*mu;
    sigma(0, 0) += ltrace;
    sigma(1, 1) += ltrace;
    sigma(2, 2) += ltrace;
    sigma = F * sigma;

    //Isotropic material strain energy density function
    // psi = 0.5 * lambda * tr(E)^2 + mu*tr(E²)
    //Incompressible tr(E) = 0 -> psi = mu*tr(E^2)
    Real psi = 0.0;
    for (unsigned char j = 0; j < 3; j++)
        for (unsigned char k = 0; k < 3; k++)
            psi += epsilon(j, k) * epsilon(j, k);
    psi = mu*psi + static_cast<Real>(0.5)*lambda * trace*trace;
    energy = restVolume * psi;
}

void PBDFEMTetra::computeGreenStrainAndPiolaStressInversion(const Matrix3 &F, const Real restVolume,
                                                               const Real mu, const Real lambda,
                                                               Matrix3 &epsilon, Matrix3 &sigma, Real &energy)
{
    //compute Piola-Kirchhoff 1st stress tensor from St Venant Kirchhoff model
    // S = lambda*tr(E)*I + 2mu*E
    Matrix3 U, VT;
    Vec3 hatF;
    MathFunctions::svdWithInversionHandling(F, hatF, U, VT);

    // Clamp small singular values
    const Real minXVal = static_cast<Real>(0.577);

    for (unsigned char j = 0; j < 3; j++)
    {
        if (hatF[j] < minXVal)
            hatF[j] = minXVal;
    }

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
    energy = restVolume*psi;

}

void PBDFEMTetra::draw(const sofa::core::visual::VisualParams *vparams)
{

    //	unsigned int i;
    if (!vparams->displayFlags().getShowForceFields()) return;
    if (!m_mechanicalObject.getValue ()) return;

    vparams->drawTool()->saveLastState();

    const VecCoord& x = m_mechanicalObject.getValue ()->read(core::ConstVecCoordId::position())->getValue();

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,true);


    std::vector< Vector3 > points[4];
    for(Topology::TetrahedronID i = 0 ; i<m_topology->getNbTetrahedra();++i)
    {
        const Tetrahedron t=m_topology->getTetrahedron(i);

        uint a = t[0];
        uint b = t[1];
        uint c = t[2];
        uint d = t[3];
        Vec3Types::Coord center = (x[a]+x[b]+x[c]+x[d])*0.125;
        Vec3Types::Coord pa = (x[a]+center)*(Real)0.666667;
        Vec3Types::Coord pb = (x[b]+center)*(Real)0.666667;
        Vec3Types::Coord pc = (x[c]+center)*(Real)0.666667;
        Vec3Types::Coord pd = (x[d]+center)*(Real)0.666667;

        points[0].push_back(pa);
        points[0].push_back(pb);
        points[0].push_back(pc);

        points[1].push_back(pb);
        points[1].push_back(pc);
        points[1].push_back(pd);

        points[2].push_back(pc);
        points[2].push_back(pd);
        points[2].push_back(pa);

        points[3].push_back(pd);
        points[3].push_back(pa);
        points[3].push_back(pb);
    }

    Vec<4,float> color1;
    Vec<4,float> color2;
    Vec<4,float> color3;
    Vec<4,float> color4;

    color1 = Vec<4,float>(0.047, 0.372, 0.074,1.0);
    color2 = Vec<4,float>(0.047, 0.372, 0.074,1.0);
    color3 = Vec<4,float>(0.047, 0.372, 0.074,1.0);
    color4 = Vec<4,float>(0.8,0.8,0.8,1.0);

    vparams->drawTool()->drawTriangles(points[0], color1);
    vparams->drawTool()->drawTriangles(points[1], color2);
    vparams->drawTool()->drawTriangles(points[2], color3);
    vparams->drawTool()->drawTriangles(points[3], color4);

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,false);

    vparams->drawTool()->restoreLastState();
}

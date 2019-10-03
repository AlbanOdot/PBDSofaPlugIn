#include "PBDStrainEnergy.hpp"
#include "../../Common/MathFunctions.hpp"
#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/defaulttype/VecTypes.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/core/behavior/ForceField.inl>
#include <SofaBaseTopology/TopologyData.inl>

int PBDStrainEnergyClass = sofa::core::RegisterObject("Constraint that correct deformed Dynamic.")
                            .add< PBDStrainEnergy >();

using namespace sofa::defaulttype;
using namespace	sofa::component::topology;
using namespace sofa;
using namespace core::topology;
void PBDStrainEnergy::bwdInit ()
{
    m_mass = PBDVertexMass<sofa::defaulttype::Vec3Types>(m_mechanicalObject.getValue (),m_topology.getValue ());
    m_basis = StrainEnergyData(m_mechanicalObject.getValue (),m_topology.getValue ());

    SReal nu = m_poisson_ratio.getValue ();
    SReal E = m_young_modulus.getValue ();
    m_lambda = (E*nu) / (( 1.0 + nu ) * ( 1.0 - 2.0 * nu ));
    m_mu = E / (2.0 * ( 1.0 + nu ));
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

bool PBDStrainEnergy::solve(sofa::simulation::Node * node)
{
    static SReal eps = 1e-6;
    static bool normalizeShear = false;
    static bool normalizeStretch = false;
    const auto& Dm_inv = m_basis.data ();
    const uint tetCount = m_topology.getValue()->getNbTetrahedra ();
    WriteCoord p = m_pbdObject->getFreePosition ();
    bool modification = false;

        for(uint tri = 0; tri < tetCount; ++tri)
        {
            const auto& t = m_topology.getValue()->getTetra (tri);
            Vec3& p0 = p[t[0]];
            Vec3& p1 = p[t[1]];
            Vec3& p2 = p[t[2]];
            Vec3& p3 = p[t[3]];

            const Matrix3& invRestMat = Dm_inv[tri].second;
            Vec3 c[3];
            c[0] = invRestMat.col(0);
            c[1] = invRestMat.col(1);
            c[2] = invRestMat.col(2);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j <= i; j++) {

                    Matrix3 ptin;
                    ptin.x() = p1 - p0;
                    ptin.y() = p2 - p0;
                    ptin.z() = p3 - p0;
                    const Matrix3& P = ptin.transposed ();

                    Vec3 fi = P * c[i];
                    Vec3 fj = P * c[j];

                    Real Sij = dot(fi,fj);

                    Real wi,wj,s1,s3;
                    if (normalizeShear && i != j) {
                        wi = fi.norm();
                        wj = fj.norm();
                        s1 = static_cast<Real>(1.0) / (wi*wj);
                        s3 = s1 * s1 * s1;
                    }

                    Vec3 d[4];
                    d[0] = Vec3(0.0, 0.0, 0.0);

                    for (int k = 0; k < 3; k++) {
                        d[k+1] = fj * invRestMat(k,i) + fi * invRestMat(k,j);

                        if (normalizeShear && i != j) {
                            d[k+1] = s1 * d[k] - Sij*s3 * (wj*wj * fi*invRestMat(k,i) + wi*wi * fj*invRestMat(k,j));
                        }

                        d[0] -= d[k+1];
                    }

                    if (normalizeShear && i != j)
                        Sij *= s1;

                    Real lambda =
                            m_mass.w(t[0]) * d[0].norm2() +
                            m_mass.w(t[1]) * d[1].norm2() +
                            m_mass.w(t[2]) * d[2].norm2() +
                            m_mass.w(t[3]) * d[3].norm2();

                    if (lambda >= eps){		// foo: threshold should be scale dependent

                        if (i == j) {	// diagonal, stretch
                            if (normalizeStretch)  {
                                Real s = sqrt(Sij);
                                lambda = static_cast<Real>(2.0) * s * (s - static_cast<Real>(1.0)) / lambda;// * stretchStiffness[i];
                            }
                            else {
                                lambda = (Sij - static_cast<Real>(1.0)) / lambda  * m_sstretch[i];
                            }
                        }
                        else {		// off diagonal, shear
                            lambda = Sij / lambda * m_sshear[i + j - 1];
                        }
                        p0 -= lambda * m_mass.w(t[0]) * d[0];
                        p1 -= lambda * m_mass.w(t[1]) * d[1];
                        p2 -= lambda * m_mass.w(t[2]) * d[2];
                        p3 -= lambda * m_mass.w(t[3]) * d[3];
                        modification = true;
                    }
                    modification |= false;
                }
            }
        }
        return  modification;
}


void PBDStrainEnergy::draw(const sofa::core::visual::VisualParams *vparams)
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

    color1 = Vec<4,float>(0.607, 0.011, 0.098,1.0);
    color2 = Vec<4,float>(0.607, 0.011, 0.098,1.0);
    color3 = Vec<4,float>(0.607, 0.011, 0.098,1.0);
    color4 = Vec<4,float>(0.8,0.8,0.8,1.0);

    vparams->drawTool()->drawTriangles(points[0], color1);
    vparams->drawTool()->drawTriangles(points[1], color2);
    vparams->drawTool()->drawTriangles(points[2], color3);
    vparams->drawTool()->drawTriangles(points[3], color4);

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,false);

    vparams->drawTool()->restoreLastState();
}

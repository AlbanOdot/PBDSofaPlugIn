#include "PBDElasticRod.hpp"
#include <sofa/core/ObjectFactory.h>
#include <unsupported/Eigen/MatrixFunctions>

#include <sofa/core/visual/VisualParams.h>
#include <sofa/defaulttype/VecTypes.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/core/behavior/ForceField.inl>
#include <SofaBaseTopology/TopologyData.inl>

int PBDElasticRodClass = sofa::core::RegisterObject("Constraint that correct elastic rod.")
                         .add< PBDElasticRod >();

typedef sofa::defaulttype::Vec3 V3;
using namespace sofa::defaulttype;
using namespace sofa;
using namespace core::topology;
using namespace	sofa::component::topology;
void PBDElasticRod::bwdInit ()
{
    m_mass = PBDVertexMass<sofa::defaulttype::RigidTypes>(m_mechanicalObject.getValue (), m_topology.getValue ());
    m_orientation = PBDOrientation(m_mechanicalObject.getValue(), m_topology.getValue());
    m_elastic_rod = PBDElasticRodData(m_mechanicalObject.getValue (), m_topology.getValue());

    auto node = dynamic_cast<sofa::simulation::Node*>(this->getContext());
    dt2 = node->getDt () * node->getDt ();
    SReal nu = m_poisson_ratio.getValue ();
    SReal E = m_young_modulus.getValue ();
    SReal mu = E / (2.0 * ( 1.0 + nu ));
    SReal d = m_radius.getValue () * 1e-1;
    SReal I = 0.25*M_PI*d*d*d*d;
    m_bendingAndTwistingKs = vec3(I*E,I*E,I*2.0*mu);
    m_bendingAndTwistingKs /= m_nbIter.getValue ();
    m_elastic_rod.applyFixedPoint(m_indices.getValue ());
    m_orientation.setInertia ({{m_bendingAndTwistingKs[0],m_bendingAndTwistingKs[1], m_bendingAndTwistingKs[2]}});
}


bool PBDElasticRod::solve(sofa::simulation::Node* node)
{

    WriteCoordR p = m_pbdObject->getFreePosition ();
    const auto& edges = m_topology.getValue ()->getEdges ();
    uint size = edges.size();
    uint step = (size % 2) == 0 ? 2 : 1;
    bool modification = false;

    for(uint e = 0 ; e < size ; e += step )
    {
        uint nextQuat = e < size - 1  ? e+1 : e;
        const sofa::core::topology::Topology::Edge& edge = edges[e];
        if(m_elastic_rod.color (e) == PBDBeamElement::RED)
            modification |= correction(p[edge[0]],p[edge[1]],
                    m_elastic_rod.length (e),
                    m_mass.w(edge[0]),m_mass.w (edge[1]),//3 firsts row are for the stretch
                    m_orientation.freeOrientation (e),m_orientation.freeOrientation (nextQuat),//The 4 next are for the bend
                    m_elastic_rod.wq (e),m_elastic_rod.wq (nextQuat),
                    m_orientation.restDarboux (e),
                    m_bendingAndTwistingKs
                    );
        uint opposite = size - 1 - e;
        nextQuat = opposite < size - 1  ? opposite+1 : opposite;
        const sofa::core::topology::Topology::Edge& op = edges[opposite];
        if(m_elastic_rod.color (opposite) == PBDBeamElement::BLACK)
            modification |= correction(p[op[0]],p[op[1]],
                    m_elastic_rod.length (opposite),
                    m_mass.w(op[0]),m_mass.w (op[1]),//3 firsts row are for the stretch
                    m_orientation.freeOrientation (opposite),m_orientation.freeOrientation (nextQuat),//The 4 next are for the bend
                    m_elastic_rod.wq (opposite),m_elastic_rod.wq (nextQuat),
                    m_orientation.restDarboux (opposite),
                    m_bendingAndTwistingKs
                    );

    }
    return modification;

}

void PBDElasticRod::draw(const sofa::core::visual::VisualParams* vparams)
{
    //	unsigned int i;
    vparams->drawTool()->saveLastState();

    const auto& x = m_mechanicalObject.getValue ()->read(sofa::core::ConstVecCoordId::position())->getValue();

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,true);

    const auto& edges = m_topology.getValue ()->getEdges ();
    uint size = edges.size();
    for(uint e = 0 ; e < size ; ++e )
    {
        //This will compute the bending deformation.
        uint nextQuat = e < size - 1  ? e+1 : e;
        Quaternionr omega   = m_orientation.freeOrientation (e).conjugate() * m_orientation.freeOrientation (nextQuat);   //darboux vector
        omega.coeffs()      = omega.coeffs() - m_orientation.restDarboux(e).coeffs(); //delta Omega with - Omega_0
        float coef = std::min(omega.squaredNorm ()* 1e5,1.0);
        float B = coef < 0.50 ? 1.0-2.0*coef : 0.0;
        float R = coef > 0.50 ? 2.0*coef : 0.0;
        float G = coef >= 0.25 && coef <= 0.75 ? (coef <= 0.5 ? 4.0*(coef-0.25)  : 1.0-4.0*(coef-0.5)): 0.0;
        Vec<4,float> color(R,G,B,1.0f);
        vparams->drawTool()->drawCylinder (x[edges[e][0]].getCenter (),x[edges[e][1]].getCenter (),0.025,color);
    }
    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,false);

    vparams->drawTool()->restoreLastState();
}

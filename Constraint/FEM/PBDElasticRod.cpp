#include "PBDElasticRod.hpp"
#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>

int PBDElasticRodClass = sofa::core::RegisterObject("Constraint that correct elastic rod.")
                         .add< PBDElasticRod >();

typedef sofa::defaulttype::Vec3 V3;

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
    SReal d = m_radius.getValue () * 1e-2;
    SReal I = 0.25*M_PI*d*d*d*d;
    m_bendingAndTwistingKs = vec3(I*E,I*E,I*2.0*mu);
    m_bendingAndTwistingKs /= m_nbIter.getValue ();
    m_elastic_rod.applyFixedPoint(m_indices.getValue ());
    m_orientation.setInertia ({{m_bendingAndTwistingKs[0],m_bendingAndTwistingKs[1], m_bendingAndTwistingKs[2]}});
}


void PBDElasticRod::solve(sofa::simulation::Node* node)
{

    WriteCoordR p = m_pbdObject->getFreePosition ();
    const auto& edges = m_topology.getValue ()->getEdges ();
    uint size = edges.size();
    for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
    {
        for(uint e = 0 ; e < size ; ++e )
        {
            const sofa::core::topology::Topology::Edge& edge = edges[e];
            if(m_elastic_rod.color (e) == PBDBeamElement::RED)
                correction(m_elastic_rod,
                           m_orientation.freeOrientation(),
                           m_orientation.restDarboux(e),
                           m_mass.w (edge[0]),
                           m_mass.w (edge[1]),
                           p,
                           m_orientation.inertia (e),
                           edge);
            uint opposite = size - 1 - e;
            const sofa::core::topology::Topology::Edge& op = edges[opposite];
            if(m_elastic_rod.color (opposite) == PBDBeamElement::BLACK)
                correction(m_elastic_rod,
                           m_orientation.freeOrientation(),
                           m_orientation.restDarboux(opposite),
                           m_mass.w (op[0]),
                           m_mass.w (op[1]),
                           p,
                           m_orientation.inertia (opposite),
                           op);

        }
    }
}



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
    SReal d = m_radius.getValue ();
    SReal I = 0.25*M_PI*d*d*d*d;
    m_bendingAndTwistingKs = vec3(I*E,I*E,I*2.0*mu);
    m_bendingAndTwistingKs /= m_nbIter.getValue ();
    m_elastic_rod.applyFixedPoint(m_indices.getValue ());
    m_orientation.setInertia ({{m_bendingAndTwistingKs[0],m_bendingAndTwistingKs[1], m_bendingAndTwistingKs[2]}});
}


void PBDElasticRod::solve(sofa::simulation::Node* node)
{

    WriteCoordR p = m_pbdObject->getFreePosition ();
    for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
    {
        for(uint e = 0 ; e < m_elastic_rod.wq().size (); ++e )
        {
            if(m_elastic_rod.color (e) == PBDBeamElement::RED)
                correction(m_elastic_rod,
                           m_orientation.freeOrientation(),
                           m_orientation.restDarboux(),
                           m_mass.w (m_elastic_rod.beginIdx (e)),
                           m_mass.w (m_elastic_rod.endIdx (e)),
                           p,
                           m_orientation.inertia (e),
                           e);
        }
        for(int e = static_cast<int>(m_elastic_rod.wq().size () - 1) ; e >= 0; --e )
        {
            if(m_elastic_rod.color (e) == PBDBeamElement::BLACK)
                correction(m_elastic_rod,
                           m_orientation.freeOrientation(),
                           m_orientation.restDarboux(),
                           m_mass.w (m_elastic_rod.beginIdx (e)),
                           m_mass.w (m_elastic_rod.endIdx (e)),
                           p,
                           m_orientation.inertia (e),
                           e);
        }
    }
}



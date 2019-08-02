#include "PDCosseratRod.hpp"
#include <sofa/core/ObjectFactory.h>
#include <Eigen/MatrixFunctions>

int PDCosseratRodClass = sofa::core::RegisterObject("Constraint that correct cosserat rod.")
                         .add< PDCosseratRod >();

void PDCosseratRod::bwdInit ()
{

    m_mass = PBDVertexMass<sofa::defaulttype::RigidTypes>(m_mechanicalObject.getValue (), m_topology.getValue ());
    m_orientation = PBDOrientation(m_mechanicalObject.getValue(), m_topology.getValue());
    m_cosserat_rod = PDCosseratRodData(m_mechanicalObject.getValue (), m_topology.getValue());

    auto node = dynamic_cast<sofa::simulation::Node*>(this->getContext());
    dt2 = node->getDt () * node->getDt ();
    SReal nu = m_poisson_ratio.getValue ();
    SReal E = m_young_modulus.getValue ();
    SReal mu = E / (2.0 * ( 1.0 + nu ));
    SReal d = m_radius.getValue ();
    SReal I = 0.25*M_PI*d*d*d*d;
    m_bendingAndTwistingKs = vec3(I*E,I*E,I*2.0*mu);
    m_bendingAndTwistingKs /= m_nbIter.getValue ();

    m_cosserat_rod.applyFixedPoint(m_indices.getValue ());
    m_cosserat_rod.setupW(m_young_modulus.getValue (),m_poisson_ratio.getValue (), m_radius.getValue ());
    m_orientation.setInertia ({{m_bendingAndTwistingKs[0],m_bendingAndTwistingKs[1], m_bendingAndTwistingKs[2]}});

    Quaternionr q; q.coeffs ().setZero ();
    for(uint e = 0; e < m_cosserat_rod.wq().size (); ++e )
    {
        m_dq.emplace_back(q);
    }

}

void PDCosseratRod::solve(sofa::simulation::Node * node)
{

    WriteCoordR p = m_pbdObject->getFreePosition ();
    uint size = m_cosserat_rod.wq().size ();
    for(uint iter = 0; iter < m_nbIter.getValue (); ++iter)
    {
        for(auto& dq : m_dq){dq.coeffs ().setZero ();}//RESET des valeurs du solver global
        for(uint e = 0 ; e < size; ++e )
        {
            if(m_cosserat_rod.color (e) == PBDBeamElement::RED)
            {
                correction(m_cosserat_rod,
                           m_orientation.freeOrientation(),
                           m_orientation.restDarboux(m_cosserat_rod.beginIdx (e)),
                           m_mass.w (m_cosserat_rod.beginIdx (e)),
                           m_mass.w (m_cosserat_rod.endIdx (e)),
                           p,
                           m_orientation.inertia (e),
                           e);
                //Ligne 10
                solveLinearSystem(m_cosserat_rod,
                                  m_orientation.freeOrientation(),
                                  m_orientation.restDarboux(m_cosserat_rod.beginIdx (e)),
                                  m_mass.w (m_cosserat_rod.beginIdx (e)),
                                  m_mass.w (m_cosserat_rod.endIdx (e)),
                                  p,
                                  m_orientation.inertia (e),
                                  e);
                std::cout << "ROUGE : "<<e<<std::endl;
            }
            uint opposite = size - 1 - e;
            if(m_cosserat_rod.color (opposite) == PBDBeamElement::BLACK)
            {
                correction(m_cosserat_rod,
                           m_orientation.freeOrientation(),
                           m_orientation.restDarboux(m_cosserat_rod.beginIdx (opposite)),
                           m_mass.w (m_cosserat_rod.beginIdx (opposite)),
                           m_mass.w (m_cosserat_rod.endIdx (opposite)),
                           p,
                           m_orientation.inertia (opposite),
                           opposite);
                //Ligne 10
                solveLinearSystem(m_cosserat_rod,
                                  m_orientation.freeOrientation(),
                                  m_orientation.restDarboux(m_cosserat_rod.beginIdx (opposite)),
                                  m_mass.w (m_cosserat_rod.beginIdx (opposite)),
                                  m_mass.w (m_cosserat_rod.endIdx (opposite)),
                                  p,
                                  m_orientation.inertia (opposite),
                                  opposite);
                std::cout << "NOIR : "<<opposite<<std::endl;
            }
        }
        std::cout<<std::endl;
        for(uint e = 0 ; e < size; ++e )
        {
            uint a = m_cosserat_rod.beginIdx (e);
            m_orientation.freeOrientation(a).coeffs () += m_dq[a].coeffs ();
             m_orientation.freeOrientation(a).normalize ();
        }
    }
}
//        for(int e = static_cast<int>(m_cosserat_rod.wq().size () - 1) ; e >= 0; --e )
//        {
//            if(m_cosserat_rod.color (e) == PBDBeamElement::BLACK)
//            {
//                correction(m_cosserat_rod,
//                           m_orientation.freeOrientation(),
//                           m_orientation.restDarboux(m_cosserat_rod.beginIdx (e)),
//                           m_mass.w (m_cosserat_rod.beginIdx (e)),
//                           m_mass.w (m_cosserat_rod.endIdx (e)),
//                           p,
//                           m_orientation.inertia (e),
//                           e);
//                //Ligne 10
//                solveLinearSystem(m_cosserat_rod,
//                                  m_orientation.freeOrientation(),
//                                  m_orientation.restDarboux(m_cosserat_rod.beginIdx (e)),
//                                  m_mass.w (m_cosserat_rod.beginIdx (e)),
//                                  m_mass.w (m_cosserat_rod.endIdx (e)),
//                                  p,
//                                  m_orientation.inertia (e),
//                                  e);
//            }
//            //We can optimize this by putting the last loop into this one since all of the reds arer already computed
//        }
